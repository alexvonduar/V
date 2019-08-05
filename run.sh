#!/bin/bash

#echo $#
if [ $# -gt 0 ]
then
	#echo $1
	for file in $1/*.pgm
	do
		./build/test_v $file
	done
fi

