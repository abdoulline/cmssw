#include "ZdcHitReconstructor_Run3.h"
#include "DataFormats/Common/interface/EDCollection.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "CalibFormats/HcalObjects/interface/HcalCoderDb.h"
#include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"

#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalSeverityLevelComputer.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalSeverityLevelComputerRcd.h"
#include "Geometry/CaloTopology/interface/HcalTopology.h"

#include <iostream>

#include "CalibFormats/HcalObjects/interface/HcalCalibrationWidths.h"
#include "CondFormats/HcalObjects/interface/HcalPedestal.h"

namespace zdchelper {
  void setZDCSaturation(ZDCRecHit rh, QIE10DataFrame& digi, int maxValue) {
    unsigned int digisize = digi.size();
    for (unsigned int i = 0; i < digisize; i++) {
      if (digi[i].adc() >= maxValue) {
        rh.setFlagField(1, HcalCaloFlagLabels::ADCSaturationBit);
        break;
      }
    }
  }

}  // namespace zdchelper

/*  Zdc Hit reconstructor allows for CaloRecHits with status words */
ZdcHitReconstructor_Run3::ZdcHitReconstructor_Run3(edm::ParameterSet const& conf)

    : reco_(conf.getParameter<bool>("correctForTimeslew"),
            conf.getParameter<bool>("correctForPhaseContainment"),
            conf.getParameter<double>("correctionPhaseNS"),
            conf.getParameter<int>("recoMethod")),
      saturationFlagSetter_(nullptr),
      det_(DetId::Hcal),
      correctTiming_(conf.getParameter<bool>("correctTiming")),
      setNoiseFlags_(conf.getParameter<bool>("setNoiseFlags")),
      setHSCPFlags_(conf.getParameter<bool>("setHSCPFlags")),
      setSaturationFlags_(conf.getParameter<bool>("setSaturationFlags")),
      setTimingTrustFlags_(conf.getParameter<bool>("setTimingTrustFlags")),
      dropZSmarkedPassed_(conf.getParameter<bool>("dropZSmarkedPassed")),
      ignoreRPD_(conf.getParameter<bool>("ignoreRPD")) {
  tok_input_QIE10 = consumes<QIE10DigiCollection>(conf.getParameter<edm::InputTag>("digiLabelQIE10ZDC"));

  std::string subd = conf.getParameter<std::string>("Subdetector");

  if (setSaturationFlags_) {
    const edm::ParameterSet& pssat = conf.getParameter<edm::ParameterSet>("saturationParameters");
    maxADCvalue_ = pssat.getParameter<int>("maxADCvalue");
  }
  if (!strcasecmp(subd.c_str(), "ZDC")) {
    det_ = DetId::Calo;
    subdet_ = HcalZDCDetId::SubdetectorId;
    produces<ZDCRecHitCollection>();
  } else if (!strcasecmp(subd.c_str(), "CALIB")) {
    subdet_ = HcalOther;
    subdetOther_ = HcalCalibration;
    produces<HcalCalibRecHitCollection>();
  } else {
    std::cout << "ZdcHitReconstructor_Run3 is not associated with a specific subdetector!" << std::endl;
  }

  // ES tokens
  htopoToken_ = esConsumes<HcalTopology, HcalRecNumberingRecord, edm::Transition::BeginRun>();
  paramsToken_ = esConsumes<HcalLongRecoParams, HcalLongRecoParamsRcd, edm::Transition::BeginRun>();
  conditionsToken_ = esConsumes<HcalDbService, HcalDbRecord>();
  qualToken_ = esConsumes<HcalChannelQuality, HcalChannelQualityRcd>(edm::ESInputTag("", "withTopo"));
  sevToken_ = esConsumes<HcalSeverityLevelComputer, HcalSeverityLevelComputerRcd>();
}

ZdcHitReconstructor_Run3::~ZdcHitReconstructor_Run3() { delete saturationFlagSetter_; }

void ZdcHitReconstructor_Run3::beginRun(edm::Run const& r, edm::EventSetup const& es) {
  const HcalTopology& htopo = es.getData(htopoToken_);
  const HcalLongRecoParams& p = es.getData(paramsToken_);
  longRecoParams_ = std::make_unique<HcalLongRecoParams>(p);
  longRecoParams_->setTopo(&htopo);
}

void ZdcHitReconstructor_Run3::endRun(edm::Run const& r, edm::EventSetup const& es) {}

void ZdcHitReconstructor_Run3::produce(edm::Event& e, const edm::EventSetup& eventSetup) {
  // get conditions
  const HcalDbService* conditions = &eventSetup.getData(conditionsToken_);
  const HcalChannelQuality* myqual = &eventSetup.getData(qualToken_);
  const HcalSeverityLevelComputer* mySeverity = &eventSetup.getData(sevToken_);

  // define vectors to pass noiseTS and signalTS
  std::vector<unsigned int> mySignalTS;
  std::vector<unsigned int> myNoiseTS;

  if (det_ == DetId::Calo && subdet_ == HcalZDCDetId::SubdetectorId) {
    edm::Handle<QIE10DigiCollection> digi;
    e.getByToken(tok_input_QIE10, digi);

    // create empty output
    auto rec = std::make_unique<ZDCRecHitCollection>();
    rec->reserve(digi->size());

    // testing QEI10 conditions
    for (auto it = digi->begin(); it != digi->end(); it++) {
      QIE10DataFrame QIE10_i = static_cast<QIE10DataFrame>(*it);
      HcalZDCDetId cell = QIE10_i.id();
      bool isRPD = cell.section() == 4;
      if (isRPD && ignoreRPD_)
        continue;
      if (cell.section() == 1 && cell.channel() > 5)
        continue;  // ignore extra EM channels

      DetId detcell = (DetId)cell;

      // check on cells to be ignored and dropped: (rof,20.Feb.09)
      const HcalChannelStatus* mydigistatus = myqual->getValues(detcell.rawId());
      if (mySeverity->dropChannel(mydigistatus->getValue()))
        continue;
      if (dropZSmarkedPassed_)
        if (QIE10_i.zsMarkAndPass())
          continue;

      const HcalCalibrations& calibrations = conditions->getHcalCalibrations(cell);
      const HcalQIECoder* channelCoder = conditions->getHcalCoder(cell);
      const HcalQIEShape* shape = conditions->getHcalShape(channelCoder);
      HcalCoderDb coder(*channelCoder, *shape);

      const HcalLongRecoParam* myParams = longRecoParams_->getValues(detcell);
      mySignalTS.clear();
      myNoiseTS.clear();
      mySignalTS = myParams->signalTS();
      myNoiseTS = myParams->noiseTS();

      // pass the effective pedestals to rec hit since both ped value and width used in subtraction of pedestals
      const HcalPedestal* effPeds = conditions->getEffectivePedestal(cell);
      rec->push_back(reco_.reconstruct(QIE10_i, myNoiseTS, mySignalTS, coder, calibrations, *effPeds, isRPD));

      // saturationFlagSetter_ doesn't work with QIE10
      // created new function zdchelper::setZDCSaturation to work with QIE10
      (rec->back()).setFlags(0);
      if (setSaturationFlags_)
        zdchelper::setZDCSaturation(rec->back(), QIE10_i, maxADCvalue_);
    }
    // return result
    e.put(std::move(rec));
  }  // else if (det_==DetId::Calo...)

}  // void HcalHitReconstructor::produce(...)

//define this as a plug-in
DEFINE_FWK_MODULE(ZdcHitReconstructor_Run3);