HydroBlocks
==========

HydroBlocks relies on a number python 2.7 libraries. To make this straightforward use conda (http://conda.pydata.org/miniconda.html) to install all the packages. Here are the steps to install the model.


1. conda install gdal
2. conda install scikit-learn
3. conda install netCDF4
4. conda install matplotlib
5. git clone https://github.com/chaneyn/HydroBlocks.git
6. cd HydroBlocks
7. python setup.py
8. cd ..
9. git clone https://github.com/chaneyn/geospatialtools.git
10. cd geospatialtools
11. python setup.py install
12. cd ..

If you want to enable the water management module
1.  conda install -c snorfalorpagus lpsolve=5.5.2.3
2.  conda install -c conda-forge glpk
3. conda install -c pywr pywr=0.2


Note: If you have the intel mkl library on your machine, we recommend that you set the mkl_flag variable in setup.py to True.

To run the model on a test dataset, do the following:

1. wget https://www.dropbox.com/s/tw4z4rf9ol6p24x/HB_sample.tar.gz?dl=0
2. tar -xvzf HB_sample.tar.gz
3. cd HB_sample
4. python ../HydroBlocks/Preprocessing/Driver.py metadata.json
5. python ../HydroBlocks/HydroBlocks/Driver.py metadata.json


