#!/bin/bash

#wget https://github.com/dealii/dealii/releases/download/v8.4.1/dealii-8.4.1.dmg
#hdiutil attach dealii-8.4.1.dmg
#./deal.II
#sudo apt-get install libdeal.ii-dev libdeal.ii-doc g++ cmake 


git clone https://github.com/dealii/dealii
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=. ../dealii -DCMAKE_BUILD_TYPE=Release
make -j2 install
#make test
