#!/bin/bash

#wget https://github.com/dealii/dealii/releases/download/v8.4.1/dealii-8.4.1.dmg
#hdiutil attach dealii-8.4.1.dmg
#./deal.II
#sudo apt-get install libdeal.ii-dev libdeal.ii-doc g++ cmake 

PRG=$PWD/programs

if [ ! -d $PRG ] 
then
    echo "create folder `$PRG`"
    mkdir $PRG
fi

# astyle
if [ ! -d $PRG/astyle ]
then
    echo "Downloading and installing astyle."
    mkdir $PRG/astyle
    wget http://downloads.sourceforge.net/project/astyle/astyle/astyle%202.04/astyle_2.04_linux.tar.gz  > /dev/null
    tar xfz astyle_2.04_linux.tar.gz -C $PRG > /dev/null
    cd $PRG/astyle/build/gcc
    make -j4 > /dev/null
fi

