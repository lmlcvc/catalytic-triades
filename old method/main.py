import itertools

import pandas as pd
from niapy.algorithms.basic import GeneticAlgorithm
from niapy.algorithms.basic.ga import multi_point_crossover, uniform_mutation, two_point_crossover
from niapy.task import Task, OptimizationType

import clean_files as cf
import find_triades as ft
import encoder
import os
import configparser

import util
import problem
import algorithm

config = configparser.ConfigParser()
config.read(os.path.join(os.pardir, 'config.ini'))
config = config['default']

transpath = config['transformed_location']
output = config['output_location']
encoded_directory = config['encoded_location']

if __name__ == "__main__":
    # check if files have been transformed
    if not os.path.isdir(transpath) or not os.listdir(transpath):
        cf.clean_files()

    # find triades and make csv files
    if not os.path.isdir(output) or not os.listdir(output):
        ft.find_triads()

    triads_protein = util.store_triads_protein(output)  # triads by protein

    # found triads encoding to files
    if not os.path.isdir(encoded_directory) or not os.listdir(encoded_directory):
        encoder.old_encode(triads=triads_protein)

    # load all encoded triads to dataframe and calculate occurrence counts
    triads_df = util.read_triads_df(encoded_directory)
    triads_count = triads_df.groupby(list(triads_df.columns)).size().reset_index(name='Count')

    for i in range(5):
        task = Task(problem=problem.MostCommonPattern(dimension=5, triads_count=triads_count), max_evals=10000,
                    optimization_type=OptimizationType.MAXIMIZATION, enable_logging=True)
        algo = algorithm.GeneticAlgorithmModified(population_size=100, crossover=algorithm.single_point_crossover,
                                                  mutation=algorithm.old_mutation,
                                                  crossover_rate=0.9, mutation_rate=0.01,
                                                  initialization_function=problem.population_init,
                                                  individual_type=problem.TriadIndividual)
        best = algo.run(task=task)
        print('%s -> %s' % (best[0], best[1]))
