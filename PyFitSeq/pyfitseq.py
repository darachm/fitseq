#!/usr/bin/env python3

import numpy as np
import pandas as pd
import math
import argparse
import itertools
import sys
from scipy.stats import linregress
from scipy.optimize import minimize
from scipy.optimize import Bounds
from tqdm import tqdm
from scipy.misc import derivative
from multiprocessing import Pool
import itertools

x0_global = None
read_num_measure_global = None
kappa_global = None
read_depth_seq_global = None
t_seq_global = None
seq_num_global = None
sum_term_global = None
fitness_type_global = None



def estimate_parameters(x,processes,total_reads,max_chunk_size):
    """Estimate parameters?
    This copied over from the old old old PyFitSeq - dunno if still relevant
    but it's missing in this version !!!
    
    A SUB-FUNCTION CALLED BY MAIN FUNCTION main() TO CALCULATE THE LOG
    LIKELIHOOD VALUE OF EACH GENOTYPE GIVEN ITS FITNESS, THE ESTIMATED READ 
    NUMBER PER GENOTYPE PER SEQUENCING TIME-POINT, AND THE ESTIMATED MEAN 
    FITNESS PER SEQUENCING TIME-POINT
    
    INPUTS ( NOT ANY more apparently....)
        * x: fitness of each genotype, [x1, x2, ...] 
        * read_num_seq: read number per genotype at each sequencing time-point 
        * t_seq: sequenced time-points in number of generations, 
            [0, t1, t2, ...] 
        * kappa: a noise parameter that characterizes the total noise introduced
            by growth, cell transfer, DNA extraction, PCR, and sequencing 
            (To measure kappa empirically, see the reference: 
                [S. F. Levy, et al. Quantitative Evolutionary Dynamics Using 
                High-resolution Lineage Tracking. Nature, 519: 181â186 (2015)].
            ) .  (default: 2.5) 
        * fitness_type: type of fitness: Wrightian fitness (w), or 
            Malthusian fitness (m)' (default: m)
    
    OUTPUTS
        * estimate_parameters_output: log likelihood value of each genotype, 
            estimated reads number per genotype per sequencing time-point,
            estimated mean fitness per sequencing time-point, 
            [x_mean(t0),x_mean(t1),...]
    """


    global read_num_measure_global
    global read_num_measure_original
    global read_depth_seq_global
    global t_seq_global
    global kappa_global
    global fitness_type_global
    global seq_num_global
    global fitness_type_global

    read_num_theory = 1e-1*np.ones(read_num_measure_global.shape, dtype=float)
    read_num_theory[:,0] = read_num_measure_global[:,0]  

    x_mean = np.zeros(seq_num_global, dtype=float)
    sum_term = np.zeros(seq_num_global, dtype=float)
    
    if fitness_type_global == 'm':  
        for k in range(1, seq_num_global):

            freq_of_lineage = (
                    read_num_measure_original[:, k] / 
                        np.sum(read_num_measure_original[:, k])
                    )
            x_mean[k] = np.average(x, weights=freq_of_lineage)

            sum_term[k] = (
                    (t_seq_global[k]-t_seq_global[k-1]) * 
                        (x_mean[k]+x_mean[k-1]) / 2
                    )

            tempt = (
                    read_num_measure_original[:, k-1] * 
                        np.exp(
                            (t_seq_global[k]-t_seq_global[k-1]) * 
                                x - sum_term[k]
                            )
                    )

            read_num_theory[:,k] = ( tempt /
                    read_depth_seq_global[k-1]*read_depth_seq_global[k]
                    )

    elif fitness_type_global == 'w':
        for k in range(1, seq_num_global):

            freq_of_lineage = (
                    read_num_measure_global[:, k] / 
                        np.sum(read_num_measure_global[:, k])
                    )
            x_mean[k] = np.maximum( np.average(x, weights=freq_of_lineage) , 0)

            if x_mean[k] != x_mean[k-1]:
                sum_term[k] = ((x_mean[k]+1)*np.log(x_mean[k]+1) - (x_mean[k-1]+1)*np.log(x_mean[k-1]+1) 
                               - (x_mean[k]-x_mean[k-1])) * (t_seq_global[k]-t_seq_global[k-1])/(x_mean[k]-x_mean[k-1])
            else:
                sum_term[k] = (t_seq_global[k] - t_seq_global[k-1]) * np.log(1 + x_mean[k-1])
                
            tempt = (
                read_num_measure_global[:,k-1] * 
                    np.exp( (t_seq_global[k]-t_seq_global[k-1]) * 
                        np.log(1+x) - sum_term[k]
                        )
                )
            read_num_theory[:,k] = tempt/read_depth_seq_global[k-1]*read_depth_seq_global[k]
    
            #x_mean[k] = np.maximum(np.dot(x, read_num_theory[:, k]) / np.sum(read_num_theory[:, k]),0)
            if x_mean[k] != x_mean[k-1]:
                sum_term[k] = ((x_mean[k]+1)*np.log(x_mean[k]+1) - (x_mean[k-1]+1)*np.log(x_mean[k-1]+1) 
                               - (x_mean[k]-x_mean[k-1])) * (t_seq_global[k]-t_seq_global[k-1])/(x_mean[k]-x_mean[k-1])
            else:
                sum_term[k] = (t_seq_global[k] - t_seq_global[k-1]) * np.log(1 + x_mean[k-1])

    if processes > 1:
        pool_obj = Pool(processes)
        other_result = pool_obj.starmap(
                calculate_likelihood_of_fitness_vector, 
                tqdm(
                    [ (x0_global[i],read_num_measure_global[i,:],kappa_global,total_reads,sum_term) 
                        for i in range(read_num_measure_global.shape[0]) ]
                    ) ,
                chunksize=np.minimum(
                        max_chunk_size,
                        int(len(x)/processes)+1
                        )
                )
    else:
        other_result = list(itertools.starmap(
                calculate_likelihood_of_fitness_vector, 
                tqdm(
                    [ (x0_global[i],read_num_measure_global[i,:],kappa_global,total_reads,sum_term) 
                        for i in range(read_num_measure_global.shape[0]) ]
                    ) ))

    parameter_output = {'Likelihood_Log': other_result,
                        'Estimated_Read_Number': read_num_theory,
                        'Estimated_Mean_Fitness': x_mean, 
                        'Sum_Term': sum_term}

    return parameter_output
 
    

##################################################        
def predict_counts(fitness,observations,total_reads,sum_term):
    """predict expected counts?
    """
    global t_seq_global
    global seq_num_global
    global fitness_type_global

    number_of_timepoints = len(observations)   

    read_num_lineage_theory = 1e-1 * np.ones(number_of_timepoints, dtype=float)
    read_num_lineage_theory[0] = observations[0]

    if fitness_type_global == 'm':
        for k in range(1, number_of_timepoints):
            tempt = (
                    observations[k-1] * 
                        np.exp(
                            (t_seq_global[k]-t_seq_global[k-1]) * 
                            fitness - sum_term[k]
                            )
                    )
# wait a sec, so this is predicting from the observed previous timepoint at every step????? that seems odd,maybe wrong
            read_num_lineage_theory[k] = (
                    tempt / total_reads[k-1] * 
                        total_reads[k]
                    )
    
    elif fitness_type_global == 'w':
        for k in range(1, number_of_timepoints):  
            tempt = observations[k-1] * np.exp((t_seq_global[k]-t_seq_global[k-1])*np.log(1+fitness) 
                                                                  - sum_term[k])
            read_num_lineage_theory[k] = tempt/total_reads[k-1]*total_reads[k]
    
    return read_num_lineage_theory



def calculate_likelihood_of_fitness_vector(fitness,observations,kappa,
        total_reads,sum_term):
    """given a fitness value, calculate the likelihood of that

    Arguments:
    fitness -- fitness to calc likelihood for
    observations -- the counts to calc likelihood for
    kappa -- that kappa parameter for noise
    """

    # generate expected counts


    expected_counts = predict_counts(fitness,observations,
        total_reads,sum_term)

    number_of_timepoints = len(observations)

    likelihood_log_seq_lineage = np.zeros(number_of_timepoints, dtype=float)
    
    read_threshold = 20
    read_threshold_2 = 10

    positions_to_consider = np.where(observations[:-1] >= read_threshold)[0]

    likelihood_log_seq_lineage[positions_to_consider + 1] = (
            0.25 * np.log(expected_counts[positions_to_consider + 1])
                - 0.5 * np.log(4 * np.pi * kappa)
                - 0.75 * np.log(observations[positions_to_consider + 1])
                - ( np.sqrt(observations[positions_to_consider + 1]) - 
                        np.sqrt(expected_counts[positions_to_consider + 1])
                    ) ** 2 / kappa
            )

    pos = np.where(observations[:-1] < read_threshold)[0]

    pos_p1 = np.where(
            observations[pos + 1] >= read_threshold_2
            )[0]
    pos_p2 = np.where(
            observations[pos + 1] < read_threshold_2
            )[0]
    pos2 = pos[pos_p1]
    pos3 = pos[pos_p2]

    likelihood_log_seq_lineage[pos2 + 1] = (
            np.multiply(
                observations[pos2 + 1],
                np.log(expected_counts[pos2 + 1])
                ) - 
                expected_counts[pos2 + 1] - 
                np.multiply(
                    observations[pos2 + 1], 
                    np.log(observations[pos2 + 1])
                    ) + 
                observations[pos2 + 1] - 
                0.5 * np.log(2 * np.pi * 
                observations[pos2 + 1])
            )

    factorial_tempt = [
            float(math.factorial(i)) for i in 
                observations[pos3 + 1].astype(int)
            ]

    likelihood_log_seq_lineage[pos3 + 1] = (
            np.multiply(
                observations[pos3 + 1],
                np.log(expected_counts[pos3 + 1])
                ) - 
                expected_counts[pos3 + 1] - 
                np.log(factorial_tempt)
            )

    likelihood_log_lineage = np.sum(likelihood_log_seq_lineage)

    return -likelihood_log_lineage


##################################################
def fun_x_est_lineage(i,tolerance):
    global x0_global
    global read_num_measure_global   
    global kappa_global
    global read_depth_seq_global
    global t_seq_global
    global seq_num_global
    global sum_term_global
    global fitness_type_global
    
    # x0_global is the currently worked on fitnesses
    optimization_result = minimize(
            fun=calculate_likelihood_of_fitness_vector, 
            x0=x0_global[i],
            args=(
                read_num_measure_global[i,:] ,
                kappa_global,
                read_depth_seq_global,
                sum_term_global
                ),
            method='BFGS',
            options={'gtol':tolerance}
            )

    return optimization_result['x'][0]

   
    
##################################################
def main():
    """
    """
    global x0_global
    global read_num_measure_global
    global read_num_measure_original
    global kappa_global
    global read_depth_seq_global
    global t_seq_global
    global seq_num_global
    global sum_term_global
    global fitness_type_global
    
    parser = argparse.ArgumentParser(description='Estimate fitness of each genotype in a competitive pooled growth experiment', 
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-p', '--processes', type=int, default=1,
        help='Number of processes to launch with multiprocessing')

    parser.add_argument('--max-chunk-size', type=int, default=None,
        help=('The max chunksize for parallelism, automatically set to '
            'a roughly even split of lineages per chunk. Tune if you want to.')
        )
    
    parser.add_argument('-i', '--input', type=str, required=True, 
        help=('The path to a header-less CSV file, where each column '
            'contains the count of each lineage (each row is a lineage) '
            'at that sample/timepoint.')
        )
    
    parser.add_argument('--t-seq', '-t', nargs='*', required=True, type=float,
        help=('The estimated "generations" of growth elapse at each sampled '
            'timepoint. This is useful for scaling the fitness or using '
            'unevenly spaced timepoints')
        )


    parser.add_argument('-o', '--output', type=str, default=sys.stdout, 
        help='The path (default STDOUT) from which to output the fitnesses '
            'and errors and likelihoods and estimated reads. CSV format.')

    parser.add_argument('--output-mean-fitness','-om', type=str, 
        default=None,
        help='The path (default None) to which to write the mean fitnesses'
            'calculated per sample.')


    parser.add_argument('--min-iter', type=int, default=10,
        help='Force FitSeq to run at least this many iterations in the '
            'optimization')

    parser.add_argument('--max-iter-num', '-m', type=int, default=100,
        help=('Maximum number of iterations in the optimization '
            '(of optimizing population average fitness)')
        )

    parser.add_argument('--minimum-step-size', '--min-step',
        type=float, default=0.0001,
        help=('Set a minimum fracitonal step size for improvement, if below '
            'this then the optimization iterations terminate.')
        )
    

    parser.add_argument('--fitness-type', '-f', type=str, default='m', 
        choices = ['m', 'w'],
        help=('SORRY no choice, only Malthusian fitness (m) works. '
            'But in later verions, '
            'maybe we\'ll re-implement Wrightian fitness (w).')
        )

    parser.add_argument('-k', '--kappa', type=float, default=2.5,
        help=('a noise parameter that characterizes the total '
             'noise introduced. For estimateion, see doi:10.1038/nature14279')
        )

    parser.add_argument('--gtol', type=float, default=1e-5,
        help='The gradient tolerance parameter for the BFGS opitmization, '
            'default (from SciPy) is 1e-5')
    
    parser.add_argument('-g', '--regression-num', type=int, default=2,
        help='number of points used in the initial '
             'linear-regression-based fitness estimate')
    

    args = parser.parse_args()
    read_num_measure_global = np.array(pd.read_csv(args.input, header=None), dtype=float)
    t_seq_global = np.array(args.t_seq, dtype=float)
    max_iter_num = args.max_iter_num
    min_iter = args.min_iter
    kappa_global = args.kappa
    regression_num = args.regression_num
    fitness_type_global = args.fitness_type
    minimum_step_size = args.minimum_step_size

    lineages_num, seq_num_global = read_num_measure_global.shape

    max_chunk_size = args.max_chunk_size
    if max_chunk_size is None:
        max_chunk_size = int(lineages_num/args.processes)+1
    else:
        max_chunk_size = int(np.minimum(max_chunk_size,lineages_num))
    

    if fitness_type_global == 'w':
        exit("Wrightian fitness does not yet work in this version")
        print('Estimating Wrightian fitness for %d lineages...' %lineages_num,file=sys.stderr)
    elif fitness_type_global == 'm':
        print('Estimating Malthusian fitness for %d lineages...' %lineages_num,file=sys.stderr)  

    read_num_measure_original = read_num_measure_global
    read_num_measure_global[read_num_measure_global < 1] = 0.1
        # This is where the minimum read is set to 0.1, so that later
        # log values do not error out
    read_depth_seq_global = np.sum(read_num_measure_original, axis=0)
 
    read_freq_seq = read_num_measure_global / read_depth_seq_global
    if fitness_type_global == 'm':
        if regression_num == 2:
            x0_tempt = np.true_divide(read_freq_seq[:, 1] - read_freq_seq[:, 0], t_seq_global[1] - t_seq_global[0])
        else:
            x0_tempt = [regression_output.slope for i in range(lineages_num) for regression_output in
                        [linregress(t_seq[0:regression_num], np.log(read_freq_seq[i, 0:regression_num]))]]
        x0 = x0_tempt #- np.dot(read_freq_seq[:, 0], x0_tempt)  # normalization

    elif fitness_type_global == 'w':
        if regression_num == 2:
            x0_tempt = np.power(np.true_divide(read_freq_seq[:, 1], read_freq_seq[:, 0]), 1 
                                / (t_seq_global[1] - t_seq_global[0])) - 1
        else:
            x0_tempt = np.exp([regression_output.slope for i in range(lineages_num) for regression_output in 
                               [linregress(t_seq_global[0:regression_num], np.log(read_freq_seq[i, 0:regression_num]))]]) - 1
        x0 = (1 + x0_tempt) / (1 + np.dot(read_freq_seq[:, 0], x0_tempt)) - 1  # normalization

    x0_global = x0


    print(r'-- Estimating initial guesses of global parameters ',file=sys.stderr)
    parameter_output = estimate_parameters(x0_global,args.processes,
            read_depth_seq_global,
            max_chunk_size
            )
    x_mean_global = parameter_output['Estimated_Mean_Fitness']
    sum_term_global = parameter_output['Sum_Term']
    likelihood_log = parameter_output['Likelihood_Log']

    likelihood_log_sum_iter = [np.sum(likelihood_log)]

    for k_iter in range(max_iter_num):   

        if fitness_type_global == 'w':
            x0_global[x0_global <= -1] = -1 + 1e-7
         
        print(r'-- Optimizing fitness for every lineage with global parms',file=sys.stderr)
        if args.processes > 1:
            with Pool(args.processes) as pool_obj:
                x0_global = np.array(
                        pool_obj.starmap(   
                            fun_x_est_lineage, 
                            tqdm([ (i,args.gtol) 
                                    for i in range(lineages_num) ]),
                            chunksize=np.minimum(
                                    max_chunk_size,
                                    int(len(x0_global)/args.processes)+1
                                )
                            )
                        )
        else:
            x0_global = np.array(
                    list(
                        itertools.starmap(fun_x_est_lineage, 
                            tqdm([ (i,args.gtol) 
                                    for i in range(lineages_num) ])
                            )
                        )
                    )

        print(r'-- Re-estimating global parms',file=sys.stderr)
        parameter_output = estimate_parameters(x0_global,args.processes,
                read_depth_seq_global,
                max_chunk_size
                )
        x_mean_global = parameter_output['Estimated_Mean_Fitness']
        sum_term_global = parameter_output['Sum_Term']
        likelihood_log = parameter_output['Likelihood_Log']

        print(r'-- Average fitnesses ', x_mean_global,file=sys.stderr)

        likelihood_log_sum_iter.append(np.sum(likelihood_log))
        print(r'-- log likelihood after iteration %i: %.4f' 
            %(k_iter+1, likelihood_log_sum_iter[-1]) ,
            file=sys.stderr)

        if (    k_iter >= min_iter and 
                (likelihood_log_sum_iter[-2] / likelihood_log_sum_iter[-1]) - 1 <= minimum_step_size
                ):
            break

    print(r'-- Calculating second derivatives around final fitness estimates',file=sys.stderr)
    # estimation error
    if args.processes > 1:
        with Pool(args.processes) as pool_obj:
            second_derivative = pool_obj.starmap(
                        derivative, 
                        tqdm( [ ( calculate_likelihood_of_fitness_vector,
                                    x0_global[i], 1e-6, 2,
                                    ( read_num_measure_global[i,:],
                                        kappa_global, 
                                        read_depth_seq_global,
                                        sum_term_global )
                                    ) for i in range(lineages_num) ] ) ,
                        chunksize=np.minimum(
                                max_chunk_size,
                                int(lineages_num/args.processes)+1
                                )
                            )
    else:
        second_derivative = list(itertools.starmap(
                        derivative, 
                        tqdm( [ ( calculate_likelihood_of_fitness_vector,
                                    x0_global[i], 1e-6, 2,
                                    ( read_num_measure_global[i,:],
                                        kappa_global, 
                                        read_depth_seq_global,
                                        sum_term_global ) 
                                    ) for i in range(lineages_num) 
                                ] 
                            ) 
                        )
                    )

    estimation_error = np.array(
            [ 1/np.sqrt(i) 
                if type(i) is np.double and i > 0 and np.sqrt(i) is not None
                else np.nan 
                for i in second_derivative
                ]
            )

    print(r'-- Writing outputs',file=sys.stderr)

    read_num_theory = parameter_output['Estimated_Read_Number']
    if fitness_type_global == 'm':
        x_opt = x0_global #- np.dot(read_num_theory[:, 0], x0_global) / np.sum(read_num_theory[:, 0]) # normalization
    elif fitness_type_global == 'w':
        x_opt = (1 + x0_global) / (1 + np.dot(read_num_theory[:, 0], x0_global)) - 1  # normalization

    fitseq_output = {'Estimated_Fitness': x_opt,
                     'Estimation_Error': estimation_error,
                     'Likelihood_Log': likelihood_log}
    for k in range(seq_num_global):
        fitseq_output['Estimated_Read_Number_t%d' % k] = read_num_theory[:, k].astype(float)

    pd.DataFrame(fitseq_output).to_csv(args.output,index=False)

    pd.DataFrame(
        {'Samples':list(range(seq_num_global)),
            'Estimate_Mean_Fitness':x_mean_global}
        ).to_csv(args.output_mean_fitness,index=False)
    
    print('Finished!',file=sys.stderr)

if __name__ == "__main__":
    main()
        
     

                
                
                

