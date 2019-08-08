#!/bin/bash

if [ "$#" -lt 8 ]; then
    echo "Wrong number of parameters - in order to run the script, you need to pass the following arguments in this order:"
    printf "Number of iterations\nPath to the CMSSW directory\nPath to the initial calibration file\nPath to the desired geometry file\nvalidOOT parameter for the cmsRun\n"
    printf "Debug parametr for the cmsRun\nPath to the JSON generator\nPath to the python script for generating the resolution plots\n"
    exit 2
fi

json_producer=$(realpath $7)
plot=$(realpath $8)
src_dir=$(realpath $2$"/src")
cal=$(realpath $3)
cd $2$"/src/Analyzer_code/SimpleAnalyzer/test"
if [ $? -ne 0 ]; then
    echo "Wrong CMSSW path or the CMSSW used doesn't contain the Analyzer_code/SimpleAnalyzer package."
    exit 2
fi

echo "Removing the old JSON files..."
rm -f *.json
cd ../../..
echo "JSON files from the previous run have been successfully deleted."

eval `scramv1 runtime -sh`
if [ $? -ne 0 ]; then
    echo "Couldn't set up the CMSSW environment correctly."
    exit 2
fi

for i in $(seq 0 $1)
do
  echo "Running timing analysis..."
  if [ $i -eq 0 ]
  then
    cmsRun Analyzer_code/SimpleAnalyzer/test/test.py calibFile=$cal geometryFile=$4 validOOT=$5 debug=$6 outputFile="out.root"
    if [ $? -ne 0 ]; then
        echo "Error while executing the cmsRun command."
        exit 2
    fi
  else
    x=$(($i-1))
    f="Calibration_L"$x$".json"
    cmsRun Analyzer_code/SimpleAnalyzer/test/test.py calibFile="Analyzer_code/SimpleAnalyzer/test/"$f geometryFile=$4 validOOT=$5 debug=$6 outputFile="out.root"
    if [ $? -ne 0 ]; then
        echo "Error while executing the cmsRun command."
        exit 2
    fi
  fi

  echo "Timing analysis done."
  let l=$i
  if [ $l -gt 2 ]
  then
    l=2
  fi

  echo "Generating JSON file from the analysis ouput..."
  if [ $? -ne 0 ]; then
      echo "Can't open the JSON generator - wrong path specified."
      exit 2
  fi

  root -b -q $json_producer$'++("'$src_dir$'/out.root",'$l',"'$src_dir$'/Analyzer_code/SimpleAnalyzer/test/Calibration_L'$i$'.json",'0')'
  if [ $? -ne 0 ]; then
      echo "Error while generating JSON file."
      exit 2
  fi
  cd $src_dir
  echo "JSON file successfully generated."
done

echo "Drawing resolution plots..."
output=$(date "+%F_%T")
python3 $plot Analyzer_code/SimpleAnalyzer/test $output
if [ $? -ne 0 ]; then
    echo "Error while drawing plots."
    exit 2
fi
echo "Plots successfully created."
