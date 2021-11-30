#!/usr/bin/env python3

import numpy as np
import pandas as pd
import math
import argparse
import itertools
import csv
from scipy.stats import linregress
from scipy.optimize import minimize
from scipy.optimize import Bounds
from tqdm import tqdm
from scipy.misc import derivative
from multiprocess import Pool, Process


x0_global = None
read_num_measure_global = None
kappa_global = None
read_num_lineage_measure_global = None
read_depth_seq_global = None
t_seq_global = None
seq_num_global = None
sum_term_global = None
fitness_type_global = None



##################################################
def estimate_parameters(x):
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

            x_mean[k] = (
                    np.dot(x, read_num_measure_global[:, k]) / 
                        read_depth_seq_global[k]
                    )

            sum_term[k] = (
                    (t_seq_global[k]-t_seq_global[k-1]) * 
                        (x_mean[k]+x_mean[k-1])/2
                    )
            
            tempt = (
                    read_num_measure_global[:, k-1] * 
                        np.exp(
                            (t_seq_global[k]-t_seq_global[k-1])*x - 
                                sum_term[k]
                            )
                    )

            read_num_theory[:,k] = (
                    tempt / read_depth_seq_global[k-1] * 
                        read_depth_seq_global[k]     
                    )

            x_mean[k] = (
                    np.dot(x, read_num_theory[:, k]) / 
                        np.sum(read_num_theory[:, k])
                    )

            sum_term[k] = (
                    (t_seq_global[k]-t_seq_global[k-1]) * 
                        (x_mean[k]+x_mean[k-1]) / 2
                    )
            
    elif fitness_type_global == 'w':
        for k in range(1, seq_num_global):
            x_mean[k] = np.maximum(np.dot(x, read_num_measure_global[:, k]) / read_depth_seq_global[k], 0)
            if x_mean[k] != x_mean[k-1]:
                sum_term[k] = ((x_mean[k]+1)*np.log(x_mean[k]+1) - (x_mean[k-1]+1)*np.log(x_mean[k-1]+1) 
                               - (x_mean[k]-x_mean[k-1])) * (t_seq_global[k]-t_seq_global[k-1])/(x_mean[k]-x_mean[k-1])
            else:
                sum_term[k] = (t_seq_global[k] - t_seq_global[k-1]) * np.log(1 + x_mean[k-1])
                
            tempt = read_num_measure_global[:,k-1] * np.exp((t_seq_global[k]-t_seq_global[k-1])*np.log(1+x) - sum_term[k])
            read_num_theory[:,k] = tempt/read_depth_seq_global[k-1]*read_depth_seq_global[k]
    
            x_mean[k] = np.maximum(np.dot(x, read_num_theory[:, k]) / np.sum(read_num_theory[:, k]),0)
            if x_mean[k] != x_mean[k-1]:
                sum_term[k] = ((x_mean[k]+1)*np.log(x_mean[k]+1) - (x_mean[k-1]+1)*np.log(x_mean[k-1]+1) 
                               - (x_mean[k]-x_mean[k-1])) * (t_seq_global[k]-t_seq_global[k-1])/(x_mean[k]-x_mean[k-1])
            else:
                sum_term[k] = (t_seq_global[k] - t_seq_global[k-1]) * np.log(1 + x_mean[k-1])

    print(x_mean)
                
    likelihood_log_seq = np.zeros(read_num_measure_global.shape, dtype=float)
    
    pos1_r, pos1_c = np.where(read_num_measure_global[:, :-1] >= 20)
    likelihood_log_seq[pos1_r, pos1_c + 1] = (0.25 * np.log(read_num_theory[pos1_r, pos1_c + 1])
                                              - 0.5 * np.log(4 * np.pi * kappa_global)
                                              - 0.75 * np.log(read_num_measure_global[pos1_r, pos1_c + 1])
                                              - (np.sqrt(read_num_measure_global[pos1_r, pos1_c + 1])
                                                 - np.sqrt(read_num_theory[pos1_r, pos1_c + 1])) ** 2 / kappa_global)

    pos_r, pos_c = np.where(read_num_measure_global[:, :-1] < 20)
    pos_p1 = np.where(read_num_measure_global[pos_r, pos_c + 1] >= 10)[0]
    pos_p2 = np.where(read_num_measure_global[pos_r, pos_c + 1] < 10)[0]
    pos2_r = pos_r[pos_p1]
    pos2_c = pos_c[pos_p1]
    pos3_r = pos_r[pos_p2]
    pos3_c = pos_c[pos_p2]

    likelihood_log_seq[pos2_r, pos2_c + 1] = (np.multiply(read_num_measure_global[pos2_r, pos2_c + 1],
                                                          np.log(read_num_theory[pos2_r, pos2_c + 1]))
                                              - read_num_theory[pos2_r, pos2_c + 1]
                                              - np.multiply(read_num_measure_global[pos2_r, pos2_c + 1],
                                                            np.log(read_num_measure_global[pos2_r, pos2_c + 1]))
                                              + read_num_measure_global[pos2_r, pos2_c + 1]
                                              - 0.5 * np.log(2 * np.pi * read_num_measure_global[pos2_r, pos2_c + 1]))

    factorial_tempt = [float(math.factorial(i)) for i in read_num_measure_global[pos3_r, pos3_c + 1].astype(int)]
    likelihood_log_seq[pos3_r, pos3_c + 1] = (np.multiply(read_num_measure_global[pos3_r, pos3_c + 1],
                                                          np.log(read_num_theory[pos3_r, pos3_c + 1]))
                                              - read_num_theory[pos3_r, pos3_c + 1] - np.log(factorial_tempt))

    likelihood_log = np.sum(likelihood_log_seq, axis=1)
    
    parameter_output = {'Likelihood_Log': likelihood_log,
                        'Estimated_Read_Number': read_num_theory,
                        'Estimated_Mean_Fitness': x_mean, 
                        'Sum_Term': sum_term}

    return parameter_output
 
    

##################################################        
def fun_read_num_lineage_theory(x):
    global read_num_lineage_measure_global
    global read_depth_seq_global
    global t_seq_global
    global seq_num_global
    global sum_term_global #seq_num_global, sum_term_global[0]=0
    global fitness_type_global
    
    read_num_lineage_theory = 1e-1*np.ones(seq_num_global, dtype=float)
    read_num_lineage_theory[0] = read_num_lineage_measure_global[0]
        
    if fitness_type_global == 'm':
        for k in range(1, seq_num_global):
            tempt = read_num_lineage_measure_global[k-1] * np.exp((t_seq_global[k]-t_seq_global[k-1])*x - sum_term_global[k])
            read_num_lineage_theory[k] = tempt/read_depth_seq_global[k-1]*read_depth_seq_global[k]
    
    elif fitness_type_global == 'w':
        for k in range(1, seq_num_global):  
            tempt = read_num_lineage_measure_global[k-1] * np.exp((t_seq_global[k]-t_seq_global[k-1])*np.log(1+x) 
                                                                  - sum_term_global[k])
            read_num_lineage_theory[k] = tempt/read_depth_seq_global[k-1]*read_depth_seq_global[k]
    
    return read_num_lineage_theory



##################################################
def fun_likelihood_lineage_opt(x):
    global kappa_global
    global read_num_lineage_measure_global
    global read_depth_seq_global
    global t_seq_global
    global seq_num_global
    global sum_term_global
    global fitness_type_global
    
    read_num_lineage_theory = fun_read_num_lineage_theory(x)
    
    likelihood_log_seq_lineage = np.zeros(seq_num_global, dtype=float)
    
    read_threshold = 1
    read_threshold_2 = 1
    pos1 = np.where(read_num_lineage_measure_global[:-1] >= read_threshold)[0]
    likelihood_log_seq_lineage[pos1 + 1] = (
            0.25 * np.log(read_num_lineage_theory[pos1 + 1])
                - 0.5 * np.log(4 * np.pi * kappa_global)
                - 0.75 * np.log(read_num_lineage_measure_global[pos1 + 1])
                - ( np.sqrt(read_num_lineage_measure_global[pos1 + 1]) - 
                        np.sqrt(read_num_lineage_theory[pos1 + 1])
                    ) ** 2 / kappa_global
            )

    pos = np.where(read_num_lineage_measure_global[:-1] < read_threshold)[0]

    pos_p1 = np.where(
            read_num_lineage_measure_global[pos + 1] >= read_threshold_2
            )[0]
    pos_p2 = np.where(
            read_num_lineage_measure_global[pos + 1] < read_threshold_2
            )[0]
    pos2 = pos[pos_p1]
    pos3 = pos[pos_p2]

    likelihood_log_seq_lineage[pos2 + 1] = (
            np.multiply(
                read_num_lineage_measure_global[pos2 + 1],
                np.log(read_num_lineage_theory[pos2 + 1])
                ) - 
                read_num_lineage_theory[pos2 + 1] - 
                np.multiply(
                    read_num_lineage_measure_global[pos2 + 1], 
                    np.log(read_num_lineage_measure_global[pos2 + 1])
                    ) + 
                read_num_lineage_measure_global[pos2 + 1] - 
                0.5 * np.log(2 * np.pi * 
                read_num_lineage_measure_global[pos2 + 1])
            )
    
    factorial_tempt = [
            float(math.factorial(i)) for i in 
                read_num_lineage_measure_global[pos3 + 1].astype(int)
            ]

    likelihood_log_seq_lineage[pos3 + 1] = (
            np.multiply(
                read_num_lineage_measure_global[pos3 + 1],
                np.log(read_num_lineage_theory[pos3 + 1])
                ) - 
                read_num_lineage_theory[pos3 + 1] - 
                np.log(factorial_tempt)
            )

    likelihood_log_lineage = np.sum(likelihood_log_seq_lineage)
    
    return -likelihood_log_lineage



##################################################
def fun_x_est_lineage(i):
    global x0_global
    global read_num_measure_global   
    global kappa_global
    global read_num_lineage_measure_global
    global read_depth_seq_global
    global t_seq_global
    global seq_num_global
    global sum_term_global
    global fitness_type_global
    
    read_num_lineage_measure_global = read_num_measure_global[i,:]
    
    #opt_output_lineage = minimize(fun_likelihood_lineage_opt, x0_global[i], method='BFGS',
    #                              options={'disp': False, 'maxiter': 500})  # how to add boundry
    
    opt_output_lineage = minimize(fun_likelihood_lineage_opt, x0_global[i], method='Nelder-Mead', 
                                  options={'ftol': 1e-8, 'disp': False, 'maxiter': 500})

    return opt_output_lineage['x'][0]

   
    
##################################################
def main():
    """
    """
    global x0_global
    global read_num_measure_global
    global kappa_global
    global read_num_lineage_measure_global
    global read_depth_seq_global
    global t_seq_global
    global seq_num_global
    global sum_term_global
    global fitness_type_global
    
    parser = argparse.ArgumentParser(description='Estimate fitness of each genotype in a competitive pooled growth experiment', 
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-p', '--processes', type=int, default=1,
        help='Number of processes to launch with multiprocessing')
    
    parser.add_argument('-i', '--input', type=str, required=True, 
        help='a .csv file: with each column being the read number per '
            'genotype at each sequenced time-point')
    
    parser.add_argument('-t', '--t_seq', nargs='*', type=float,
        help='sequenced time-points in number of generations')
    
    parser.add_argument('-m', '--max_iter_num', type=int, default=10,
                        help='maximum number of iterations in the optimization')

    parser.add_argument('--min-iter', type=int, default=0,
        help='Force FitSeq to run at least this many iterations in the '
            'optimization')
    parser.add_argument('--min-step','--minimum-step-size', 
        type=float, default=0.0001,
        help='Set a minimum fracitonal step size for improvement, if below '
            'this then the optimization iterations terminate.')
    
    parser.add_argument('-k', '--kappa', type=float, default=2.5,
        help='a noise parameter that characterizes the total '
             'noise introduced by growth, cell transfer, DNA '
             'extraction, PCR, and sequencing (To measure kappa '
             'empirically, see the reference: [S. F. Levy et al. '
             'Quantitative Evolutionary Dynamics Using '
             'High-resolution Lineage Tracking. Nature, 519: '
             '181–186 (2015)].)')
    
    parser.add_argument('-g', '--regression_num', type=int, default=2,
        help='number of points used in the initial '
             'linear-regression-based fitness estimate')
    
    parser.add_argument('-f', '--fitness_type', type=str, default='m', 
        choices = ['m', 'w'],
        help='type of fitness: Wrightian fitness (w), or '
             'Malthusian fitness (m)')
    
    parser.add_argument('-o', '--output_filename', type=str, default='output', 
                        help='prefix of output .csv files')

    
    args = parser.parse_args()
    read_num_measure_global = np.array(pd.read_csv(args.input, header=None), dtype=float)
    t_seq_global = np.array(args.t_seq, dtype=float)
    max_iter_num = args.max_iter_num
    min_iter = args.min_iter
    kappa_global = args.kappa
    regression_num = args.regression_num
    fitness_type_global = args.fitness_type
    output_filename = args.output_filename
    minimum_step_size = args.min_step
    
    lineages_num, seq_num_global = read_num_measure_global.shape

    if fitness_type_global == 'w':
        print('Estimating Wrightian fitness for %d lineages...' %lineages_num)
    elif fitness_type_global == 'm':
        print('Estimating Malthusian fitness for %d lineages...' %lineages_num)  


    ##################################################
    read_num_measure_global[read_num_measure_global < 1] = 0.1
        # This is where the minimum read is set to 0.1, so that later
        # log values do not error out
    read_depth_seq_global = np.sum(read_num_measure_global, axis=0)
 
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

    
    ##################################################
    likelihood_log_sum_iter = []
    for k_iter in range(max_iter_num):   
        if k_iter == 0:
            x0_global = x0
        else:
            x0_global = opt_result
         
        if fitness_type_global == 'w':
            x0_global[x0_global <= -1] = -1 + 1e-7
    
        parameter_output = estimate_parameters(x0_global)
        x_mean_global = parameter_output['Estimated_Mean_Fitness']
        sum_term_global = parameter_output['Sum_Term']
        likelihood_log = parameter_output['Likelihood_Log']
   
        likelihood_log_sum_iter.append(np.sum(likelihood_log))
        print(r'-- log likelihood after iteration %i: %.4f' %(k_iter+1, likelihood_log_sum_iter[-1]))

        if (    k_iter>=1 and 
                k_iter >= min_iter and 
                (likelihood_log_sum_iter[-2] / likelihood_log_sum_iter[-1]) - 1 <= minimum_step_size
                ):
            break

        if args.processes > 1:
            pool_obj = Pool(args.processes)
            opt_result = pool_obj.map(fun_x_est_lineage, tqdm(range(lineages_num)))
            opt_result = np.array(opt_result)
            print(np.subtract(x0_global,opt_result))
            print(np.mean(np.subtract(x0_global,opt_result)))
            print(np.median(np.subtract(x0_global,opt_result)))
        else:
            opt_result = list(map(fun_x_est_lineage, tqdm(range(lineages_num))))
            opt_result = np.array(opt_result)
            print(np.subtract(x0_global,opt_result))
            print(np.mean(np.subtract(x0_global,opt_result)))
            print(np.median(np.subtract(x0_global,opt_result)))

    ################################################## estimation error
    second_derivative = np.zeros(lineages_num, dtype=float)
    for i in range(lineages_num):
        read_num_lineage_measure_global = read_num_measure_global[i,:]
        second_derivative[i] = derivative(fun_likelihood_lineage_opt, x0_global[i], dx=1e-6, n=2)
    estimation_error = 1/np.sqrt(second_derivative)

    read_num_theory = parameter_output['Estimated_Read_Number']
    if fitness_type_global == 'm':
        x_opt = x0_global #- np.dot(read_num_theory[:, 0], x0_global) / np.sum(read_num_theory[:, 0]) # normalization
    elif fitness_type_global == 'w':
        x_opt = (1 + x0_global) / (1 + np.dot(read_num_theory[:, 0], x0_global)) - 1  # normalization

    fitseq_output = {'Estimated_Fitness': x_opt,
                     'Estimation_Error': estimation_error,
                     'Likelihood_Log': likelihood_log,
                     'Estimated_Mean_Fitness': x_mean_global}

    for k in range(seq_num_global):
        fitseq_output['Estimated_Read_Number_t%d' % k] = read_num_theory[:, k].astype(float)

    tempt = list(itertools.zip_longest(*list(fitseq_output.values())))
    with open(output_filename + '_FitSeq.csv', 'w') as f:
        w = csv.writer(f)
        w.writerow(fitseq_output.keys())
        w.writerows(tempt)
    
    print('Finished!')

if __name__ == "__main__":
    main()
        
     

                
                
                

