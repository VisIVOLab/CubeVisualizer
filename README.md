# CubeVisualizer


To generate makefile:
```
cd build
cmake -DCFITSIO_INC_PATH=<CFITSIO INCLUDE PATH> -DCFITSIO_LIB_PATH=<CFITSIO LIB PATH> -DBOOST_INC_PATH=<BOOST PATH> .. 
```

To generate Xcode project:
```
cd build
cmake -G Xcode  -DCFITSIO_INC_PATH=<CFITSIO INCLUDE PATH> -DCFITSIO_LIB_PATH=<CFITSIO LIB PATH> -DBOOST_INC_PATH=<BOOST PATH> .. 
```
