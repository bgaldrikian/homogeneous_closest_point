#!/bin/bash
FNAME="hcp_tests_auto.txt"
FPATH=$(pwd)/$FNAME
pushd tests > /dev/null
echo Writing output to $FNAME
echo HCP tests automatic build and unit test execution > $FPATH
date >> $FPATH
echo >> $FPATH
echo OSTYPE=$OSTYPE > $FPATH
echo >> $FPATH
echo Running lscpu: >> $FPATH
lscpu >> $FPATH
echo >> $FPATH
echo Building tests... | tee -a $FPATH
echo >> $FPATH
make >> $FPATH
echo >> $FPATH
echo Running ./bin/unit_tests.exe... | tee -a $FPATH
echo >> $FPATH
./bin/unit_tests.exe >> $FPATH
echo >> $FPATH
echo Completed $(date) >> $FPATH
echo Done.
popd > /dev/null
