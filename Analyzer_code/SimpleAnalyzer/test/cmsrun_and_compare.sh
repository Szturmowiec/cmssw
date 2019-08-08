#!/bin/bash

if [ "$#" -lt 7 ]; then
    echo "Wrong number of parameters - in order to run the script, you need to pass the following arguments in this order:"
    printf "Path to the CMSSW directory\nPath to the initial calibration file\nPath to the desired geometry file\nvalidOOT parameter for the cmsRun\n"
    printf "Debug parametr for the cmsRun\nPath to the python script for comparing the outputs\nPath to the original output root file name\n"
    exit 2
fi

path=$(realpath $6)
original_root=$(realpath $7)
cd $1$"/src"
if [ $? -ne 0 ]; then
    echo "Wrong CMSSW path or the CMSSW used doesn't contain the Analyzer_code/SimpleAnalyzer package."
    exit 2
fi

eval `scramv1 runtime -sh`
if [ $? -ne 0 ]; then
    echo "Couldn't set up the CMSSW environment correctly."
    exit 2
fi

echo "Compiling..."
scram b -j 8
if [ $? -ne 0 ]; then
    echo "Compilation failed."
    exit 2
fi
echo "Compilation completed successfully."

echo "Running timing analysis..."
cmsRun Analyzer_code/SimpleAnalyzer/test/test.py calibFile=$2 geometryFile=$3 validOOT=$4 debug=$5 outputFile="out.root"
  if [ $? -ne 0 ]; then
      echo "Error while executing the cmsRun command."
      exit 2
  fi
echo "Timing analysis done."

echo "Comparsion script output:"
echo "export PYTHONPATH=/usr/lib64/root" >> ~/.bashrc
python $path out.root $original_root
