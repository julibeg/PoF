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

## License

The program relies on some Boost libraries (https://www.boost.org/) as well as Luigi Pertoldi's progressbar 
(https://github.com/gipert/progressbar). 
The corresponding licenses can be found in the respective subdirectories in `/include`. The rest of the code in `/src` is
licensed under GPL v3 (see LICENSE file in root directory).
