import pandas as pd
from niapy.task import Task, OptimizationType

import analysis
import clean_files as cf
import find_triades as ft
import encoder
import os
import configparser

import util
import problem
import algorithm
from task import TaskModified

config = configparser.ConfigParser()
config.read(os.path.join(os.pardir, 'config.ini'))
config = config['default']

transpath = config['transformed_location']
output = config['output_location']
encoded_directory = config['encoded_location']

ga_output = config['ga_output_location']
ga_most_common = config['ga_most_common']
ga_enzyme_common = config['ga_enzyme_common']

output_analysis = config['analysis_output_location']
best_occurrences = config['best_occurrences']
similarity = config['similarity']
fitness = config['fitness']
plots = config['plots']

HEADER = ['Nuc', 'Acid', 'Base', 'D1', 'D2', 'fitness']

if __name__ == "__main__":
    # check if files have been transformed
    if not os.path.isdir(transpath) or not os.listdir(transpath):
        cf.clean_files()

    # find triads and make csv files
    if not os.path.isdir(output) or not os.listdir(output):
        ft.find_triads()

    triads_protein = util.store_triads_protein(output)  # triads by protein

    # found triads encoding to files
    if not os.path.isdir(encoded_directory) or not os.listdir(encoded_directory):
        encoder.old_encode(triads=triads_protein)

    # load all encoded triads to dataframe and calculate occurrence counts
    triads_df = util.read_triads_df(encoded_directory)
    triads_count = triads_df.groupby(list(triads_df.columns)).size().reset_index(name='Count')

    triads_dict = util.read_triads_dict(encoded_directory)
    triads_dict_count = {}
    for protein in triads_dict.keys():
        triads_dict_count[protein] = triads_dict[protein].groupby(
            list(triads_dict[protein].columns)).size().reset_index(name='Count')

    # - - - - - - - - - -

    # run genetic algorithm
    util.create_folder(ga_output)
    util.create_folder(fitness)
    util.create_folder(plots)
    file_most_common = open(os.path.join(ga_output, ga_most_common), 'w+')

    most_common_df = pd.DataFrame()
    for i in range(3):
        task = TaskModified(problem=problem.MostCommonPattern(dimension=5, triads_count=triads_count, method='old'),
                            max_evals=1000, optimization_type=OptimizationType.MAXIMIZATION, enable_logging=True)

        algo = algorithm.GeneticAlgorithmModified(type='most_common', iteration=str(i).zfill(2), population_size=100,
                                                  crossover=algorithm.single_point_crossover,
                                                  mutation=algorithm.old_mutation,
                                                  crossover_rate=0.9, mutation_rate=0.01,
                                                  initialization_function=algorithm.population_init_mixed,
                                                  individual_type=problem.TriadIndividual)

        best = algo.run(task=task)
        print(best)
        task.convergence_data(algo_type=algo.type, iteration=str(i).zfill(2), output_directory=fitness, x_axis="evals")

        most_common_df = pd.concat([most_common_df, util.get_iteration_info(population=algo.population_list,
                                                                            header=HEADER)])
    most_common_df.to_csv(file_most_common, header=HEADER, index=False)
    file_most_common.close()

    #####

    file_enzyme_common = open(os.path.join(ga_output, ga_enzyme_common), 'w+')
    enzyme_common_df = pd.DataFrame()

    for i in range(3):
        task = TaskModified(problem=problem.EnzymeCommonPattern(dimension=5, triads_count=triads_dict_count,
                                                                triads_count_dict=triads_dict_count, method='old'),
                            max_evals=250, optimization_type=OptimizationType.MAXIMIZATION, enable_logging=True)

        algo = algorithm.GeneticAlgorithmModified(type='enzyme_common', iteration=str(i).zfill(2), population_size=100,
                                                  crossover=algorithm.single_point_crossover,
                                                  mutation=algorithm.old_mutation,
                                                  crossover_rate=0.9, mutation_rate=0.01,
                                                  initialization_function=algorithm.population_init_mixed,
                                                  individual_type=problem.TriadIndividual)

        algo.run(task=task)
        task.convergence_data(algo_type=algo.type, iteration=str(i).zfill(2), output_directory=fitness, x_axis="evals")

        enzyme_common_df = pd.concat([enzyme_common_df, util.get_iteration_info(population=algo.population_list,
                                                                                header=HEADER)])
    enzyme_common_df.to_csv(file_enzyme_common, header=HEADER, index=False)
    file_enzyme_common.close()

    # results analysis
    util.create_folder(output_analysis)
    util.create_folder(best_occurrences)
    util.create_folder(similarity)

    analysis.store_triad_count(output, output_analysis, similar=True)
    analysis.store_best_individual_occurrences(ga_output, HEADER, best_occurrences)
    analysis.store_similarity_best(ga_output, similarity)  # similarity between best individuals
