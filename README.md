# Real Space Renormalization Group

## Quick Start

```bash
# install
$ git clone https://github.com/refraction-ray/realspace-RG.git
$ cd realspace-RG
$ make 
# usage
$ bin/rsrg -i utils/input.txt -o utils/output.txt
$ mpirun -np 3 bin/rsrgmpi -i utils/input.txt -o utils/output.txt
```

## Install and requirements

This project requires C++11 support from compiler, LAPACK support and mpi support.  The external library [Armadillo](http://arma.sourceforge.net/) is also included for matrix calculations. 

To install and configure LAPACK and mpi, please refer to their detailed manuals (there might be various implementations of them, and they are in general well configured for supercomputers).

For Armadillo library, it is not necessary to install it as a shared library and call it by `-larmadillo` when compile. Instead, including the path of Armadillo header file on compiling is enough. Namely `gcc -I <arma header file path> -DARMA_DONT_USE_WRAPPER`, and some more flags on shared library may need depending on the system. On macOS, you may need `-framework Accelerate`, on tianhe-2 SC, one should include `-mkl=parallel  -openmp`, on some linux, one should use `-lblas -llapack`, etc.

The default makefile is designed for macOS with Aramdillo installed as shared library. For other systems, some modifications on makefile should be applied based on the discussion above.

## System support

Both non-parallel version and parallel version of the program are tested on MacBook Pro with macOS 10.11 and Tianhe-2 supercomputer in Guangzhou.

With at most some small modifications on makefile, the project should work on all *nix system. Windows is not supported for now.

## Physics picture

## Usage

## Development

## References