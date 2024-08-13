#! /bin/bash
ARGC=$#

echo "Compiling C++ code..."
g++ RunAbundanceOnGoodLines.cpp -o RunAbundanceOnGoodLines -lgsl -lgslcblas -O2 -Wall -std=c++20


if [ $? -ne 0 ];
then
echo "Something went wrong when compiling the C++ code"
echo "Check that you have g++ installed and that you've pointed the compiler to the correct directory where libgsl is located, using the compiler flag -L[/path/to/libgsl]"
exit 1
fi;



pip show numpy 
if [ $? -eq 1 ];
then
echo "Installing numpy..."
pip install numpy
fi;

pip show astropy
if [ $? -eq 1 ];
then
echo "Installing astropy..."
pip install astropy
fi;

pip show scipy
if [ $? -eq 1 ];
then
echo "Installing scipy..."
pip install scipy
fi;

pip show astroquery
if [ $? -eq 1 ];
then
echo "Installing astroquery..."
pip install astroquery
fi;

gfortran linemake/linemake.f -o linemake/linemake -O2 -ffixed-line-length-none
if [ $? -ne 0 ];
then
echo "Something went wrong when compiling Linemake"
echo "Check that you have gfortran installed"
exit 1
fi;


cd moog
make
cd ..


echo "Downloading mooglists"
curl "https://drive.usercontent.google.com/download?id=1zZEePXb1mW17dlddEtWwYrrRPksNKfek&confirm=xxx" -o mooglists.tar.gz
echo "Downloading isochrones"
curl "https://drive.usercontent.google.com/download?id=1FS8aIk9e-7hbHdNQTs0e2u-8vTNzM0T9&confirm=xxx" -o isochrones.tar.gz
echo "Downloading PHOENIX data"
curl "https://drive.usercontent.google.com/download?id=17Mpy5yd1zTaYQgtBdA_H2ae8P17CWLvM&confirm=xxx" -o PHOENIX.tar.gz


echo "Extracting mooglists..."
tar -xf mooglists.tar.gz
echo "Extracting isochrones..."
tar -xf isochrones.tar.gz
echo "Extracting PHOENIX data..."
tar -xf PHOENIX.tar.gz
echo "DONE\n"
echo "Try executing the following command:"
echo ""
echo "  ./runStar.sh \"tau Ceti\" example-tauCet/ -0.41"
exit 0
