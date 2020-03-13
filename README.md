# PoF

Calculates the probability of failure (i.e. the complement of structural robustness) for metabolic networks. 


### Prerequisites

GCC 4.8.1 or newer 


### Installing

Clone and build by

```
git clone https://github.com/julibeg/PoF.git
cd PoF
make
```

### Example
To test the installation, go into `/test_files` and run the analysis on compressed MCSs of the *E. coli* model *i*JO1366 
via 
```
PoFcalc -m iJO1366.mcs.comp.binary -c iJO1366.num_comp_rxns -r 2583 -d 5
```
The command line arguments are explained under `PoFcalc -h`.

## License

The program relies on some Boost libraries (https://www.boost.org/) as well as Luigi Pertoldi's progressbar 
(https://github.com/gipert/progressbar). 
The corresponding licenses can be found in the respective subdirectories in `/include`. The rest of the code in `/src` is
licensed under GPL v3 (see LICENSE file in root directory).
