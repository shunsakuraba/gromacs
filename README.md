# Implementation of the ZMM on GROMACS (2016+ version)

This is a patched version of GROMACS that implements[[1]](#1) the Zero-multipole summation methods.[[2]](#2)

## Branches

* zmm-2016/zmm-2018/zmm-2019
    - This branch. Implementations are based on GROMACS 2016/2018/2019.
* zq-nonzeroalpha-50
    - Zero-quadrupole implementation on top of GROMACS 5.0.
* zo-50
    - Zero-octupole implementation on top of GROMACS 5.0.

## Installation

A typical installation procedure is as follows:

````sh
git clone --branch zmm-2016 https://github.com/shunsakuraba/gromacs.git
cd gromacs
# from now on, consult standard GROMACS installation procedure.
# for example:
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$HOME/opt/gromacs-2016-zmm ..
make -j8 && make install
````

## Usage

Add following terms to the mdp and run `gmx grompp`  as an usual GROMACS run.

````
Coulombtype = Zero-multipole
zmm-degree = 2 ; 0 (Wolf method), 1 (up to dipole cancellation), 2 (quadupole), 3 (octupole) are supported
zmm-alpha = 1.0 ; unit is nm^-1.
````

One of my recommended combination of parameters is:
````
zmm-degree = 2
zmm-alpha = 0
rcoulomb = 1.2
rvdw = 1.2
````
Which corresponds to "ZMM 1.2 2 0.0" in Wang et al.,[[3]](#3) Table 2; with these parameters, Wang et al. showed that the excess chemical potential, constant volume molar heat capacity, constant pressure heat capacity, isothermal compressibility, thermal expansion coefficient, diffusion constant, and viscocity of TIP3P water molecules are all within 95% confidence interval of the smooth particle mesh Ewald method. I recommend this parameter combination because it is fast AND accurate.

## Limitations

* The cases with `zmm-alpha` ≠ 0 are NOT supported with GPU. If you want to run `zmm-alpha` ≠ 0 with GPU-equipped machine, specify `-nb cpu` to the `mdrun` arguments to disable the GPU.

* ZMM with dipole cancellation and ZMM with `zmm-alpha` ≠ 0 are slightly slower than the published version. This is a trade-off we paid to simplify the code change. For the implementation we measured please consult branches for GROMACS 5.0.

* These branches are frequently rebased for merging upstream gromacs updates. Furture bug-fixes in zmm branches may also be merged and rebased into existing commits, which may cause a problem during `git fetch` or `git pull`. It is advised to use `-f` option in such a case.

<a id="1">[1]</a> Shun Sakuraba, Ikuo Fukuda. Performance evaluation of the zero-multipole summation method in modern molecular dynamics software. Journal of Computational Chemistry 39 (20), 1551-1560 (2018). [Link](https://doi.org/10.1002/jcc.25228), [arXiv preprint](https://arxiv.org/pdf/1704.07071).

<a id="2">[2]</a> Ikuo Fukuda. Zero-multipole summation method for efficiently estimating electrostatic interactions in molecular system. Journal of Chemical Physics, 139, 174107 (2013). [Link](https://doi.org/10.1063/1.4827055). Ikuo Fukuda, Kamiya Narutoshi, and Haruki Nakamura. The zero-multipole summation method for estimating electrostatic interactions in molecular dynamics: Analysis of the accuracy and application to liquid systems. Journal of Chemical Physics, 140, 194307 (2014). [Link](https://doi.org/10.1063/1.4875693).

<a id="3">[3]</a> Han Wang, Haruki Nakamura, Ikuo Fukuda. A Critical Appraisal of the Zero-Multipole Method: Structural, Thermodynamic, Dielectric, and Dynamical Properties of a Water System. Journal of Chemical Physics, 144, 114503 (2016). [Link](https://doi.org/10.1063/1.4943956).
