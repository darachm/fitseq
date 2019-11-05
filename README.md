[![MIT License][license-shield]][license-url]
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python3][python-shield]][https://www.python.org/]

## FitSeq

### 1. What is FitSeq?

FitSeq is a Python-based fitness estimation tool for pooled amplicon sequencing studies. FitSeq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


### 2. How to install FitSeq?
* Tested on MacBook Pro (3.1 GHz Intel Core i5), with Python 3.7.4.
* Clone this repository by running `git clone https://github.com/FangfeiLi05/FitSeq.git` in terminal.
* `cd` to the root directory of the project (the folder containing README.md)
* Install dependencies by running `pip install -r requirements.txt` in terminal.

### 3. How to use FitSeq?

#### 3.1. Evolution Simulations
`evo_simulator.py` models competative pooled growth of a population of genotypes with different fitnesses. This simulation may include many sources of noise, including growth noise, noise from cell transfers, DNA extraction, PCR, and sequencing.

##### OPTIONS
* `--input` or `-i`: a .csv file, with
  + 1st column: fitness of each genotype, [x1, x2, ...]
  + 2nd column: initial cell number of each genotype at generation 0, [n1, n2, ...]
* `--t_seq` or `-t`: sequenced time-points (`format: 0 t1 t2 ...`)
* `--read_num_average_seq` or `-r`: average number of reads per genotype per sequencing time-point (`format: 0 r1 r2 ...`)
* `--noise_option` or `-n`: five types of possible noise (cell growth, bottleneck transfer, DNA extraction, PCR, sequencing (`default: growth bottleneck_transfer DNA_extraction PCR sequencing`)
* `--dna_copies` or `-d`: average copy number of genome DNA per genotype as template in PCR (`default: 500`)
* `--pcr_cycles` or `-p`: number of cycles in PCR (`default: 25`) 
* `--output_filename` or `-o`: prefix of output .csv files (`default: output`)

##### OUTPUTS
* `output_filename_EvoSimulation_Read_Number.csv`: read number per genotype per sequencing time-point
* `output_filename_EvoSimulation_Mean_Fitness.csv`: mean fitness per sequencing time-point
* `output_filename_EvoSimulation_Input_Log.csv`: a record of all inputs

##### EXAMPLE
```
python evo_simulator.py --help

python evo_simulator.py -i Input.csv -t 0 3 6 9 12 -r 50 50 50 50 50 -o result
```      

#### 3.2. Fitness Estimation
`fitseq.py` estimates the fitness of each genotype from read-count time-series data.

##### OPTIONS
* `--input` or `-i`: a .csv file, with each column being the read number per genotype at each sequencing time-point
* `--t_seq` or `-t`: sequenced time-points (`format: 0 t1 t2 ...`)
* `--max_iter_num` or `-m`: maximum number of iterations in the optimization (`default: 10`)
* `--kappa` or `-k`: a noise parameter that characterizes the total noise introduced by growth, cell transfer, DNA extraction, PCR, and sequencing. (To measure kappa empirically, see the reference: Levy et al. Nature 2015 519, 181-186.) (`default: 2.5`)
* `--regression_num` or `-g`: number of points used in initial linear regression (`default: 2`)
* `--output_filename` or `-o`: prefix of output .csv files (`default: output`)

##### OUTPUTS
* `output_filename_FitSeq_Result.csv`: a .csv file, with
  + 1st column: estimated fitness of each genotype, [x1, x2, ...]
  + 2nd column: log likelihood value of each genotype, [f1, f2, ...]
  + 3rd column: estimated mean fitness per sequencing time-point, [x_mean(0), x_mean(t1), ...]
  + 4th column: estimated reads number per genotype per sequencing time-point

##### EXAMPLE
```
python fitseq.py --help

python fitseq.py -i result_EvoSimulation_Read_Number.csv -t 0 3 6 9 12 -o result
```  



