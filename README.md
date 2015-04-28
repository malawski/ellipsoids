# ellipsoids
Dense packing of ellipsoids

## Compiling

GCC:

```
g++  ellipsoids-openmp.cpp -fopenmp -o ellipsoids-openmp
```

Intel Compiler:

```
module add intel
icc ellipsoids-openmp.cpp -fopenmp -o ellipsoids-openmp
```




## Execution of single program

```
./ellipsoids-openmp ell.dat el.log 1
```
