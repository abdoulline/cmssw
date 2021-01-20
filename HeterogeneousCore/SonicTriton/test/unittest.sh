#!/bin/bash

LOCALTOP=$1

# the test is not possible if:
# 1. avx instructions not supported (needed for singularity on CPU)
# 2. singularity not found or not usable
# so just return true in those cases

if grep -q avx /proc/cpuinfo; then
	echo "has avx"
else
	echo "missing avx"
	exit 0
fi

if type singularity >& /dev/null; then
	echo "has singularity"
else
	echo "missing singularity"
	exit 0
fi

cmsRun ${LOCALTOP}/src/HeterogeneousCore/SonicTriton/test/tritonTest_cfg.py modules=TritonGraphProducer,TritonGraphFilter,TritonGraphAnalyzer maxEvents=1 unittest=1 verbose=1
CMSEXIT=$?

STOP_COUNTER=0
while ! LOGFILE="$(ls -rt ${LOCALTOP}/log_triton_server_instance*.log 2>/dev/null | tail -n 1)" && [ "$STOP_COUNTER" -lt 5 ]; do
	STOP_COUNTER=$((STOP_COUNTER+1))
	sleep 5
done

if [ -n "$LOGFILE" ]; then
	echo -e '\n=====\nContents of '$LOGFILE':\n=====\n'
	cat "$LOGFILE"
fi

exit $CMSEXIT
