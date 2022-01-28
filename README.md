# CubeVisualizer


To generate makefile:
```
$ cd build
$ cmake \
    -DVTK_DIR=<VTK DIR> \
    -DCFITSIO_INC_PATH=<CFITSIO INCLUDE PATH> \
    -DCFITSIO_LIB_PATH=<CFITSIO LIB PATH> \
    -DBOOST_INC_PATH=<BOOST PATH> \
    -DOPENMP_LIB_PATH=<OpenMP LIB PATH> \
    -DOPENMP_INC_PATH=<OpenMP INCLUDE PATH> \
    .. 
```

To generate Xcode project add the following argument:

`-G Xcode`

On an Apple Silicon Mac add the following argument to build for x86_64 arch:

`-DCMAKE_OSX_ARCHITECTURES=x86_64`
