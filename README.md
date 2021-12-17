[![Python 3.8](https://img.shields.io/badge/python-3.8-green.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

# fitseq

This is a utility for analyzing counts data of lineages to best estimate
individual and population-average lineage fitness.

<!--
-->

# Versions

- PyFitSeq is a Python-based fitness estimation tool for pooled amplicon 
    sequencing studies. 
    The conceptual/math work and first implementation is described in the paper
    [Unbiased Fitness Estimation of Pooled Barcode or Amplicon Sequencing Studies](https://doi.org/10.1016/j.cels.2018.09.004),
    and this code is [available here](https://github.com/sashaflevy/Fit-Seq).
- This was rewritten in python, [available here](https://github.com/FangfeiLi05/PyFitSeq)
    and is a python-translated version of the MATLAB tool 
    FitSeq above.
- This repo is a fork of that python version to fix some bugs and tweak the 
    speed, flexibility, and interface.
    **Wrightian fitness does not yet work in this version**. Sorry.

If you use this software, please reference: [F. Li, et al. Unbiased Fitness Estimation of Pooled Barcode or Amplicon Sequencing Studies. Cell Systems, 7: 521-525 (2018)](https://doi.org/10.1016/j.cels.2018.09.004)

# Installation

## With pip

You can install from this git repo directly as:

    python3 -m pip install git+https://github.com/darachm/fitseq.git

or from PyPi with:

    python3 -m pip install fitseq

Install the latest development branch with something like:

    python3 -m pip install git+https://github.com/darachm/fitseq.git@dev

Test intstallationg with:

    fitseq.py -h

## Or don't install, use a container

Use the `Dockerfile` in this repo like so:

    TODO gotta write it

On a multi-user HPC? Want to get an achive-ready monolithic stable container?
[Singularity](https://sylabs.io/guides/3.8/user-guide/quick_start.html#quick-installation-steps)
is a container system for scientific multi-user HPC computing and archiving.
You can build your own container from the Singularity file in this repo using
a command like:

    singularity build fitseq.sif Singularity.fitseq

or, download from the github container registry like so

    TODO gotta figure this out

Built on a Thinkpad t480, so AVX instruction set.


# Usage 

The `fitseq.py` script functions to estimate fitnesses of lineages in a pool.
There is also a script `evo_simulator.py` that can perform simulations of 
competitive pooled growth of lineages, in order to generate synthetic data for
benchmarking.

## `fitseq.py` - estimate fitnesses from counts data

This tool expects a comma-separated table (CSV) of your best estimate of
lineage counts of the lineage, with one column per timepoint. Each lineage is
a row, and outputs are in the same order as the input.

For an example using data distributed in this repo, try:

    python3 pyfitseq.py \
        --input testing/data/ppiseq_test_counts_1000.csv \
        --processes 8 \
        --t-seq 0 1 2 3 4 \
        --min-iter 10 \
        --max-iter-num 100 \
        --min-step 0.001 \
        --output-mean-fitness test_output_means.csv \
        -o test_output.csv

This reads an input file at 
`testing/data/ppiseq_test_counts_1000.csv`, and uses 8 processes.
It assumes each sample is 1 "generation" of growth.
It does a mandatory 10 iterations of burn-in to stabilize the estimates, then
proceeds until the sum of negative log likelihood doesn't improve by at least
0.1% over the previous step, at 100 iterations max.
Then it writes the mean fitness values to that CSV, 
and the rest to `test_output.csv`.

### Options, from `pyfitseq.py -h`

#### input

* `-i` INPUT, `--input` INPUT The path to a header`-less` CSV file, where each
  column contains the count of each lineage (each row is a lineage) at that
  sample/timepoint. REQUIRED
* `--t-seq` [T_SEQ [T_SEQ ...]], `-t` [T_SEQ [T_SEQ ...]] The estimated
  "generations" of growth elapse at each sampled timepoint. This is useful for
  scaling the fitness or using unevenly spaced timepoints. REQUIRED

#### output

* `-o` OUTPUT, `--output` OUTPUT The path (default STDOUT) from which to output
  the fitnesses and errors and likelihoods and estimated reads. CSV format.
  (default: STDOUT, so that you can pipe it into other tools)
* `--output-mean-fitness` OUTPUT_MEAN_FITNESS, `-om` OUTPUT_MEAN_FITNESS The
  path (default None) to which to write the mean fitnessescalculated per
  sample. 

#### parallelism

* `-p` PROCESSES, `--processes` PROCESSES Number of processes to launch with
  multiprocessing
* `--max-chunk-size` MAX_CHUNK_SIZE The max chunksize for parallelism,
  automatically set to a roughly even split of lineages per chunk. 
  Tune if you want to.

#### optimization stopping control

* `--min-iter` MIN_ITER   Force FitSeq to run at least this many iterations in
  the optimization (default: 10)
* `--max-iter-num` MAX_ITER_NUM, `-m` MAX_ITER_NUM Maximum number of iterations
  in the optimization (of optimizing population average fitness) (default: 100)
* `--minimum-step-size` MINIMUM_STEP_SIZE, `--min-step` MINIMUM_STEP_SIZE Set a
  minimum fracitonal step size for improvement, if below this then the
  optimization iterations terminate. (default: 0.0001)

#### tuning details
* `--fitness-type` {m,w}, `-f` {m,w} 
  SORRY but **Wrightian fitness does not yet work in this version**, 
  so just don't set the `--fitness_type`, or set to `m`. Sorry.
  Maybe we'll re-implement Wrightian fitness (w). Maybe.
* `-k` KAPPA, `--kappa` KAPPA a noise parameter that characterizes the total
  noise introduced. For estimation, see doi:10.1038/nature14279 (default: 2.5)
* `--gtol` GTOL The gradient tolerance parameter for the BFGS opitmization, 
  default (from SciPy) is 1e-5
* `-g` REGRESSION_NUM, `--regression-num` REGRESSION_NUM number of points used
  in the initial linear`-regression-based` fitness estimate (default: 2)


<!--
A walk-through of an old version is included as a jupyter notebook in a previous version of the software [here](https://github.com/FangfeiLi05/PyFitSeq/blob/master/PyFitSeq_Walk_Through.ipynb).
-->


## Evolution Simulator

`evo_simulator.py` simulates competitive pooled growth of lineages.
This simulation includes sampling noise from growth, 
cell transfers (bottlenecks), DNA extraction, PCR, and sequencing.
For example:

    python evo_simulator.py -i input_EvoSimulation.csv \
        -t 0 3 6 9 12 -r 50 50 50 50 50 \
        -o output

    python evo_simulator.py -i input_EvoSimulation.csv \
        -t 0 2 4 6 8 -r 75 75 75 75 50 \
        -n DNA_extraction PCR sequencing -d 300 -p 27 -f w \
        -o output

### Options

* `--input` or `-i`: a .csv file, with
  + 1st column of .csv: fitness of each genotype, [x1, x2, ...]
  + 2nd column .csv: initial cell number of each genotype at generation 0, 
    [n1, n2, ...]
* `--t_seq` or `-t`: time-points evaluated in number of generations 
    (format: 0 t1 t2 ...)
* `--read_num_average_seq` or `-r`: average number of reads per genotype 
    for each time-point (format: 0 r1 r2 ...)
* `--noise_option` or `-n`: which types of noise to include in the simulation, 
    default is all sources of noise 
    (default: growth bottleneck_transfer DNA_extraction PCR sequencing)
* `--dna_copies` or `-d`: average genome copy number per genotype used as 
    template in PCR (default: 500)
* `--pcr_cycles` or `-p`: number of cycles of PCR (default: 25) 
* `--fitness_type` or `-f`: type of fitness: 
    Wrightian fitness (w), or Malthusian fitness (m)' (default: m)
* `--output_filename` or `-o`: prefix of output .csv files (default: output),
    this tool *AUTOMATICALLY* generates files named, for a `-o` option of 
    `output`:
    + `output_filename_EvoSimulation_Read_Number.csv`: 
        read number per genotype for each time-point
    + `output_filename_EvoSimulation_Mean_Fitness.csv`: 
        mean fitness for each time-point
    + `output_filename_EvoSimulation_Input_Log.csv`: 
        a record of all inputs

See `python evo_simulator.py --help` for a reminder...



