#!/bin/bash
FNAME="hcp_examples_auto.txt"
FPATH=$(pwd)/$FNAME
pushd examples > /dev/null
echo Writing output to $FNAME
echo HCP examples automatic build and run > $FPATH
date >> $FPATH
echo >> $FPATH
echo OSTYPE=$OSTYPE > $FPATH
echo >> $FPATH
echo Running lscpu: >> $FPATH
lscpu >> $FPATH
echo >> $FPATH
echo Building examples... | tee -a $FPATH
echo >> $FPATH
make >> $FPATH
for i in {1..8}
do
	echo >> $FPATH
	echo Running ./bin/ex$i.exe: >> $FPATH
	echo >> $FPATH
	./bin/ex$i.exe >> $FPATH
done
echo >> $FPATH
echo Completed $(date) >> $FPATH
echo Done.
popd > /dev/null
