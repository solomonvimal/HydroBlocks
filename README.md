HydroBlocks
==========

HydroBlocks relies on a number python libraries. To make this straightforward use conda (http://conda.pydata.org/miniconda.html) to install all the packages. Here are the steps to install the model.


1. conda install gdal
2. conda install scikit-learn
3. conda install netCDF4
4. conda install matplotlib
5. git clone https://github.com/chaneyn/HydroBlocks.git
6. cd HydroBlocks
7. python setup.py

If you want to enable the water management module
8.  conda install -c snorfalorpagus lpsolve=5.5.2.3
9.  conda install -c conda-forge glpk
10. conda install -c pywr pywr=0.2

Note: If you have the intel mkl library on your machine, we recommend that you set the mkl_flag variable in setup.py to True.

To run the model on a test dataset, do the following:

1. wget https://www.dropbox.com/s/384ug2i70j033zh/LittleWashita.tar.gz?dl=0
2. tar -xvzf LittleWashita.tar.gz
3. cd LittleWashita
4. python ../HydroBlocks/Preprocessing/Driver.py metadata_littlewashita.json
5. python ../HydroBlocks/HydroBloks/Driver.py metadata_littlewashita.json


