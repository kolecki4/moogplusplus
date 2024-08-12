#! /bin/bash
ARGC=$#

pip install astropy numpy astroquery scipy

g++ RunAbundanceOnGoodLines.cpp -o RunAbundanceOnGoodLines -lgsl -lgslcblas -O2 -Wall -std=c++20

gfortran linemake/linemake.f -o linemake/linemake -O2 -ffixed-line-length-none

echo "Downloading mooglists"
curl "https://drive.usercontent.google.com/download?id=1zZEePXb1mW17dlddEtWwYrrRPksNKfek&confirm=xxx" -o mooglists.tar.gz
echo "Downloading isochrones"
curl "https://drive.usercontent.google.com/download?id=1FS8aIk9e-7hbHdNQTs0e2u-8vTNzM0T9&confirm=xxx" -o isochrones.tar.gz
echo "Downloading PHOENIX data"
curl "https://drive.usercontent.google.com/download?id=17Mpy5yd1zTaYQgtBdA_H2ae8P17CWLvM&confirm=xxx" -o PHOENIX.tar.gz


echo "Extracting files from tar.gz archives"
tar -xf mooglists.tar.gz
tar -xf isochrones.tar.gz
tar -xf PHOENIX.tar.gz


