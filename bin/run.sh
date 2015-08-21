#!/bin/bash
# first argument: executable name
# second argument: input file name
# third argument: inputfile id number
nice nohup ./"$1" "$2" > out"$3".txt &
