#!/bin/sh

#cd Dropbox/qBio/HIV_c #move to the correct directory

echo What should we name the output pdf? #ask what to name the output

read fn #read the filename

g++ -o hiv_c HIV.cpp  -lm #compile the HIV code

./hiv_c #run the compiled code

#run the python plotting script
if [ -z "$fn" ]
then
	echo you did not enter name, thus we get testah
	python HIVplot.py testah
	open testah.pdf
else
	python HIVplot.py $fn
	open $fn.pdf
fi
