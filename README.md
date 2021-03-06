# Real Space Renormalization Group

## Quick Start

```bash
# install
$ git clone https://github.com/refraction-ray/realspace-RG.git
$ cd realspace-RG
$ make 
# usage
$ bin/rsrg -i utils/input.txt -o utils/output.txt
$ mpirun -np 3 bin/rsrgmpi -i utils/input.txt -o utils/output.txt -r -d ./
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

The original idea of this algorithm is from these two papers: 1) A. C. Potter, R. Vasseur, and S. A. Parameswaran, Phys. Rev. X **5**, 031033 (2015) 2)  P. T. Dumitrescu, R. Vasseur, and A. C. Potter, Phys. Rev. Lett. **119**, 110604 (2017). They proposed a new approach of RG in real space to address the universality of  MBL-thermal criticality. By adopted their method and made some crucial updates, we generalized the novel RSRG approach to more general settings and specifically studied the MBL transition induced by quasi-periodic potentials. This method can reveal the resonace structure and entanglement entropy of systems near MBL critical points and one can extract critical exponents from the results. This project is the backend of the paper Phys. Rev. Lett. **121**, 206601 (2018). If you utilize this project in your study, please cite this paper.

## Usage

One can only run the program as a demo on his or her desktop or laptop. To get some more serious results, a supercomputer is necessary for this program. And the mpi version of the program is specifically designed for large supercomputers, i.e. `bin/rsrgmpi` compiled from `src/main_mpi.cpp`.

basic options: `-i <input/file/path>`, `-o <output/file/path>`, for file io pathes. If no path specified, the default path is `input.txt` and `output.txt`, if there is file with the same name, they will be override.

 `-q` is set by default, which is suitable for quasiperiodic system (and system with both random and quasiperiodic potentials as well as simple random systems).

`-d` is for further detail enumeration of observables instead of simple average information. `-d <path>` can generate multiple txt files with detail values of observables for each disorder configuration. Such data can be further analysed on its distributions.

The Hamiltonian used in current version has NNN hopping and both type of potentials (random vs. quasiperiodic). The measures input used in current version is default `get_measures()` which includes three quantities:  average-length ratio, max-length ratio and entanglement entropy density.

The format of the input.txt, every line is for one process, the format is `NN-hopping qp-potential-amplitude qp-potential-wavevector qp-potential-sample-window NNN-hopping random-onsite-potential system-size interaction repeat-times` which is separated by space. Note the random potential term can be negative -W. This means the potential window is from -W/2 to W/2 instead of 0 to W for positive values. The format of output file is similar, with several terms more append on each line. They are `mean-average-length-ratio std-average-length-ratio mean-max-length-ratio std-max-length-ratio mean-entanglement-entropy std-entanglement-entropy`. The format of dist file, arranged three data a group are values of mean-length, max-length and entanglement-entropy for each disorder configuration.

### New format on pseudorandom models (experimental)

The model is a two band model with random potential but share the similar universal behavior with QP models in terms of MBL criticality, see arXiv:1906.00971 for details.

The format of input file `lower-band-width`, `upper-band-width`,`middle-gap-width`,`no-transfer-probability`,`system-size`,`interaction`,`repeat-times`. To explore this model, use the command option `-p` instead of `-q` (by default) or `-r`.

### New format on model with thermal bath inclusion (experimental)

The format is `NN-hopping qp-potential-amplitude qp-potential-wavevector qp-potential-sample-window NNN-hopping random-onsite-potential thermal-inclusion-size system-size interaction repeat-times` . The model is the same as the quasiperiodic model except the fact that the first `thermal-inclusion-size` of sites share the same chemical potential and thus consist of a thermal bath. Use `-t` flag  for this model. Note: the short options `q,p,t,g` are exclusive and only one can be used for the computations.

### New format on model with general form of on-site potential and NN/NNN hoppings (experimental)

The format is still one line for each Hamiltonian paramters with the form `hm_param, size, interaction, repeat-times`. In terms of hm_param, to describe more general systems, we divide it into three parts which control on-site potential, NN hopping and NNN hopping, respectively. The form of each one is like `c1,c2,c3,c4,c5,c6,c7,c8,c9…`, representing the form $$c_1+c_2*U[0,1]+c_3\cos(c_4*i+c_5)+c_6\cos(c_7*i+c_8)$$, the sum of cosine ended when 0 is in the position of $$c_{3n}$$. And the sequence begins to describe the next term. For example, `-w/2,w,0,0,0,t,k,2pi,0,0,0`, represents the system with on-site potential $$[-W/2,W/2]$$, also with quasiperiodic NN hopping by $$t\cos(ki+\phi)$$, where $$\phi$$ is a fixed random number distributed in $$[0,2\pi]$$ for each Hamiltonian.

To read such data form and do the real space RG calculation, you need `-g` flag in cli.

## Workflow

To use it on supercomputers with slurm. One can first write a bash script as `run.sh`

```bash
#!/bin/bash
srun -n 80 -N 4 bin/rsrgmpi -i utils/input.txt -o utils/output.txt -q -d data/
```

where 80 should be the number of lines in inputfile **plus 1**, and 4 is the number of total nodes you want to utilize. To submit the task, just use `sbatch -N 4 run.sh` .

If one want to utilize the massive data with `-d` option for each disorder configuration, one should utilize the python3 script in `utils/` to merge files and generate simple data report for further regression or visulization analysis. One may need further hack on the python script to make it suitable for certain type of filename.

For more scripts in mathematica to generate input files and extract content from output files, together with some basic critical behavior analysis, see `utils/`.

## Development

This project is desinged with cosideration on extensibility and scalability. So it is very ease to add new types of Hamiltonians and new types of observables. Hack `src/hamiltonian.cpp` for the former one and hack `src/model.cpp` for the latter one.

## How to cite

S.-X. Zhang and H. Yao, Phys. Rev. Lett. **121**, 206601 (2018).

## References

*We list arXiv version of the papers mentioned above so that everyone can access them.*

* [Universal properties of many-body localization transitions in quasiperiodic systems](https://arxiv.org/abs/1805.05958)
* [Universal properties of many-body delocalization transitions](https://arxiv.org/abs/1501.03501)
* [Scaling Theory of Entanglement at the Many-Body Localization Transition](https://arxiv.org/abs/1701.04827)
* [Strong and Weak Many-Body Localizations](https://arxiv.org/abs/1906.00971)