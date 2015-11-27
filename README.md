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

## Using HyperFLow
### Generate workflow description

```
cd ellipsoids/js
npm install
node generator_template.js /path/to/tests /path/to/ellipsoids/ellipsoids-openmp 3 /path/to/ellipsoids/mean-pack-ell > wf.json
```

### Execute workflow locally

```
hflow run wf.json 
```

### Use amqpExecutor

```
node generator_template.js /path/to/tests /path/to/ellipsoids/ellipsoids-openmp 3 /path/to/ellipsoids/mean-pack-ell amqpCommand > amqpwf.json
```
