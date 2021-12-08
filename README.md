fesom2-kernels

# Intro

This repository contains kernels extracted from https://github.com/FESOM/fesom2 (233d72b62585c490e27c1212437d0ba283754169) to ease:
- to compare different optimisation
- to compare different OpenACC porting strategies

The authors acknowledge the funding from the ESiWACE2 project. The ESiWACE2 project has received funding from the European Union’s Horizon 2020 research and innovation programme under grant agreement No 823988.

# Usage

This repository contains binary input files. Those files have not been commited to the history but are tracked using git-annex.

You can use DataLad, a tool wrapping Git and git-annex. Install it via your pacakge manager or using pip
````sh
pip install datalad
````

The steps are:
1. clone this repository
2. enable the github LFS
3. get the input files

````sh
datalad clone https://github.com/ESiWACE-S1/fesom2-kernels
cd fesom2-kernels
git annex enableremote esiwace-s1-lfs
datalad get intputs/
````

Happy benchmarking
