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

This project requires C++11 support from compiler, LAPACK and mpi support are also required.  The external library [Armadillo](http://arma.sourceforge.net/) is also included for matrix calculations. 

To install and configure LAPACK and mpi, please refer to their detailed manuals (there might be various implementations of them, and they are in general well configured for supercomputers).

For Armadillo library, it is not necessary to install it as a shared library and call it by `-larmadillo` when compile. Instead, including the path of Armadillo header file on compiling is enough. Namely `gcc -I <arma header file path> -DARMA_DONT_USE_WRAPPER`, and some more flags on shared library may need depending on the system. On macOS, you may need `-framework Accelerate`, on tianhe-2 SC, one should include `-mkl=parallel  -openmp`, on some linux, one should use `-lblas -llapack`, etc.

The default makefile is designed for macOS with Aramdillo installed as shared library. For other systems, some modifications on makefile should be applied based on the discussion above.

## System support

Both non-parallel version and parallel version of the program are tested on MacBook Pro with macOS 10.11 and Tianhe-2 supercomputer in Guangzhou.

With at most some small modifications on makefile, the project should work on all *nix system. Windows is not supported for now.

## Physics picture

The original idea of this algorithm is from these two papers: 1) A. C. Potter, R. Vasseur, and S. A. Parameswaran, Phys. Rev. X 5, 031033 (2015) 2)  P. T. Dumitrescu, R. Vasseur, and A. C. Potter, Phys. Rev. Lett. 119, 110604 (2017). They proposed a new approach of RG in real space to address the universality of  MBL-thermal criticality. By adopted their method and made some crucial updates, we generalized the novel RSRG approach to more general settings and specifically studied the MBL transition induced by quasi-periodic potentials. This method can reveal the resonace structure and entanglement entropy of systems near MBL critical points and one can extract critical exponents from the results. This project is the backend of the paper **arXiv: 1805.05958**. If you utilize this project in your study, please cite this paper.

## Usage

One can only run the program as a demo on his or her desktop or laptop. To get some more serious results, a supercomputer is necessary for this program. And the mpi version of the program is specifically designed for large supercomputers, i.e. `bin/rsrgmpi` compiled from `src/main_mpi.cpp`.

## Development

## How to cite

S.-X. Zhang and H. Yao, *Universal properties of many-body localization transitions in quasiperiodic systems*, arXiv: 1805.05958 (2018).

## References

*We list arXiv version of the papers mentioned above such that everyone can access them.*

* [Universal properties of many-body localization transitions in quasiperiodic systems](https://arxiv.org/abs/1805.05958)

* [Universal properties of many-body delocalization transitions](https://arxiv.org/abs/1501.03501)

* [Scaling Theory of Entanglement at the Many-Body Localization Transition](https://arxiv.org/abs/1701.04827)