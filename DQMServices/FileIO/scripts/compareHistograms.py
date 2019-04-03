#!/bin/env python

from __future__ import print_function
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import os
import sys
import argparse
import root_numpy
import numpy as np
from blacklist import get_blacklist

def create_dif(base_file_path, pr_file_path, pr_number, test_number, cmssw_version, output_dir_path):
   base_file = ROOT.TFile(base_file_path, 'read')
   pr_file = ROOT.TFile(pr_file_path, 'read')

   if base_file.IsOpen():
      print('Baseline file successfully opened', file=sys.stderr)
   else:
      print('Unable to open base file', file=sys.stderr)
      return

   if pr_file.IsOpen():
      print('PR file successfully opened', file=sys.stderr)
   else:
      print('Unable to open PR file', file=sys.stderr)
      return

   run_nr = get_run_nr(pr_file_path)
      
   # Get list of paths (lists of directories)
   base_flat_dict = flatten_file(base_file, run_nr)
   pr_flat_dict = flatten_file(pr_file, run_nr)

   # Paths that appear in both baseline and PR data. (Intersection)
   shared_paths = list(set(pr_flat_dict).intersection(set(base_flat_dict)))

   # Paths that appear only in PR data. (Except)
   only_pr_paths = list(set(pr_flat_dict).difference(set(base_flat_dict)))

   # Paths that appear only in baseline data. (Except)
   only_base_paths = list(set(base_flat_dict).difference(set(pr_flat_dict)))

   # Histograms pointed to by these paths will be written to baseline output
   paths_to_save_in_base = []

   # Histograms pointed to by these paths will be written to pr output
   paths_to_save_in_pr = []

   # Make comparison
   compare(shared_paths, pr_flat_dict, base_flat_dict, paths_to_save_in_pr, paths_to_save_in_base)
   
   # Collect paths that have to be written to baseline output file
   for path in only_base_paths:
      item = base_flat_dict[path]

      if item == None:
         continue

      paths_to_save_in_base.append(path)
   
   # Collect paths that have to be written to PR output file
   for path in only_pr_paths:
      item = pr_flat_dict[path]

      if item == None:
         continue

      paths_to_save_in_pr.append(path)

   base_output_filename = get_output_filename(pr_file_path, pr_number, test_number, cmssw_version, False)
   pr_output_filename = get_output_filename(pr_file_path, pr_number, test_number, cmssw_version, True)

   # Write baseline output
   save_paths(base_flat_dict, paths_to_save_in_base, os.path.join(output_dir_path, 'base', base_output_filename))
   
   # Write PR output
   save_paths(pr_flat_dict, paths_to_save_in_pr, os.path.join(output_dir_path, 'pr', pr_output_filename))

   ROOT.gROOT.GetListOfFiles().Remove(base_file)
   ROOT.gROOT.GetListOfFiles().Remove(pr_file)

   pr_file.Close()
   base_file.Close()

   # Info about changed, added and removed elements
   nr_of_changed_elements = len(set(paths_to_save_in_base).intersection(set(paths_to_save_in_pr)))
   nr_of_removed_elements = len(paths_to_save_in_base) - nr_of_changed_elements
   nr_of_added_elements = len(paths_to_save_in_pr) - nr_of_changed_elements

   print('Base output file. PR output file. Changed elements, removed elements, added elements:')
   print(base_output_filename)
   print(pr_output_filename)
   print('%s %s %s' % (nr_of_changed_elements, nr_of_removed_elements, nr_of_added_elements))

def compare(shared_paths, pr_flat_dict, base_flat_dict, paths_to_save_in_pr, paths_to_save_in_base):
   # Collect paths that have to be written to both output files
   for path in shared_paths:
      pr_item = pr_flat_dict[path]
      base_item = base_flat_dict[path]

      if pr_item == None or base_item == None:
         continue

      if pr_item.InheritsFrom('TH1') and base_item.InheritsFrom('TH1'):
         # Compare bin by bin
         pr_array = root_numpy.hist2array(pr_item)
         base_array = root_numpy.hist2array(base_item)
         if pr_array.shape != base_array.shape or not np.allclose(pr_array, base_array, equal_nan=True):
            paths_to_save_in_pr.append(path)
            paths_to_save_in_base.append(path)
            continue
      else:
         # Compare non histograms
         if pr_item != base_item:
            paths_to_save_in_pr.append(path)
            paths_to_save_in_base.append(path)

def flatten_file(file, run_nr):
   result = {} 
   for key in file.GetListOfKeys():
      try:
         traverse_till_end(key.ReadObj(), [], result, run_nr)
      except:
         pass
   
   return result

def traverse_till_end(node, dirs_list, result, run_nr):
   new_dir_list = dirs_list + [get_node_name(node)]
   if hasattr(node, 'GetListOfKeys'): 
      for key in node.GetListOfKeys():
         traverse_till_end(key.ReadObj(), new_dir_list, result, run_nr)
   else:
      path = tuple(new_dir_list)
      if path not in get_blacklist(run_nr):
         result[path] = node

def get_node_name(node):
   if node.InheritsFrom('TObjString'):
      # Strip out just the name from a tag (<name>value</name>)
      return node.GetName().split('>')[0][1:]
   else:
      return node.GetName()

def save_paths(flat_dict, paths, result_file_path):
   if len(paths) == 0:
      print('No differences were observed - output will not be written', file=sys.stderr)
      return

   # Make sure output dir exists
   result_dir = os.path.dirname(result_file_path)
   if not os.path.exists(result_dir):
      os.makedirs(result_dir)
   
   result_file = ROOT.TFile(result_file_path, 'recreate')

   if not result_file.IsOpen():
      print('Unable to open %s output file' % result_file_path, file=sys.stderr)
      return

   for path in paths:
      save_to_file(flat_dict, path, result_file)

   result_file.Close()
   print('Output written to %s file' % result_file_path, file=sys.stderr)

# Saves file from flat_dict in the same dir of currently open file for writing
def save_to_file(flat_dict, path, output_file):
   histogram = flat_dict[path]

   current = output_file

   # Last item is filename. No need to create dir for it
   for directory in path[:-1]:
      current = create_dir(current, directory)
      current.cd()

   histogram.Write()

# Create dir in root file if it doesn't exist
def create_dir(parent_dir, name):
   dir = parent_dir.Get(name)
   if not dir:
      dir = parent_dir.mkdir(name)
   return dir

def get_output_filename(input_file_path, pr_number, test_number, cmssw_version, isPr):
   # Samples of correct output file format:
   # DQM_V0001_R000320822__wf136_892_pr__CMSSW_10_4_0_pre3-PR25518-1234__DQMIO.root
   # When run number is 1 we have to use RelVal naming pattern:
   # DQM_V0002_R000000001__RelVal_wf136_892_pr__CMSSW_10_4_0_pre3-PR25518-1234__DQMIO.root

   input_file_name = os.path.basename(input_file_path)

   run = input_file_name.split('_')[2]
   workflow = os.path.basename(os.path.dirname(input_file_path)).split('_')[0].replace('.', '_')
   if not workflow:
      workflow = 'Unknown'

   relval_prefix = ''
   if run == 'R000000001':
      relval_prefix = 'RelVal_'

   baseOrPr = 'base'
   if isPr:
      baseOrPr = 'pr'

   return 'DQM_V0001_%s__%swf%s_%s__%s-PR%s-%s__DQMIO.root' % (run, relval_prefix, workflow, baseOrPr, cmssw_version, pr_number, test_number)

def get_run_nr(file_path):
   return os.path.basename(file_path).split('_')[2].lstrip('R').lstrip('0')

if __name__ == '__main__':
   parser = argparse.ArgumentParser(description="This tool compares DQM monitor elements found in base-file with the ones found in pr-file."
      "Comparison is done bin by bin and output is written to a root file containing only the changes.")
   parser.add_argument('-b', '--base-file', help='Baseline IB DQM root file', required=True)
   parser.add_argument('-p', '--pr-file', help='PR DQM root file', required=True)
   parser.add_argument('-n', '--pr-number', help='PR number under test', default='00001')
   parser.add_argument('-t', '--test-number', help='Unique test number to distinguish different comparisons of the same PR.', default='1')
   parser.add_argument('-r', '--release-format', help='Release format in this format: CMSSW_10_5_X_2019-02-17-0000', default=os.environ['CMSSW_VERSION'])
   parser.add_argument('-o', '--output-dir', help='Comparison root files output directory', default='dqmHistoComparisonOutput')
   args = parser.parse_args()

   cmssw_version = '_'.join(args.release_format.split('_')[:4])
   
   create_dif(args.base_file, args.pr_file, args.pr_number, args.test_number, cmssw_version, args.output_dir)
