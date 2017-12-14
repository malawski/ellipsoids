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

## Publications

When using this code, please cite:

1. Bargieł, M., Szczygłowski, Ł., Trzcionkowski, R., Malawski, M., Towards large-scale parallel simulated packings of ellipsoids with OpenMP and HyperFlow, CGW workshop'15 : October 26–28, 2015 Kraków, Poland. Proceedings: eds. Marian Bubak, Michał Turała, Kazimierz Wiatr, Kraków, Academic Computer Centre CYFRONET AGH, 2015. ISBN: 978-83-61433-14-9, 117–118.

2. Bargieł, M., Geometrical properties of simulated packings of ellipsoids, CGW workshop'14 : October 27–29, 2014, Krakow, Poland. Proceedings: eds. Marian Bubak, Michał Turała, Kazimierz Wiatr, Kraków: Academic Computer Centre CYFRONET AGH, 2014. ISBN: 978-83-61433-10-1, 123–124.
