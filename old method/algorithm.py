import configparser
import itertools
import math
import operator
import os
import random

import numpy as np
import pandas as pd
from niapy.algorithms import Individual, default_individual_init
from niapy.algorithms.basic import GeneticAlgorithm
from niapy.algorithms.basic.ga import tournament_selection, uniform_crossover, uniform_mutation
from niapy.util import objects_to_array

import util
from problem import TriadIndividual

config = configparser.ConfigParser()
config.read(os.path.join(os.pardir, 'config.ini'))
config = config['default']

output_analysis = config['analysis_output_location']
final_population = config['final_population']
encoded_directory = config['encoded_location']  # TODO: adapt to other methods if necessary

HEADER = ['Nuc', 'Acid', 'Base', 'D1', 'D2', 'fitness']
HEADER_NO_FITNESS = ['NUC', 'ACID', 'BASE', 'Dist_Nuc_Acid', 'Dist_Acid_Base']  # TODO: change depending on columns nr


# TODO: split code into multiple methods
def population_init_mixed(task, population_size, rng, all_distances=False, angles=False, distance_categories=20,
                          angle_categories=20, **kwargs):
    triads_df = util.read_triads_df(encoded_directory)
    n_triads = triads_df.shape[0]

    if n_triads > 0.2 * population_size:
        n_triads = math.floor(0.2 * population_size)

    n_permutations = population_size - n_triads

    # generate all possible permutations
    options = [[0, 1], [0, 1], [0, 1], list(range(0, distance_categories)),
               list(range(0, distance_categories))]  # atom types and 2 distances

    if all_distances:  # append list for third distance
        options.append(range(0, distance_categories))

    if angles:  # append lists for angles
        options.append(range(0, angle_categories))
        options.append(range(0, angle_categories))
        options.append(range(0, angle_categories))

    permutations = [list(a) for a in itertools.product(*options)]
    permutations_df = pd.DataFrame(permutations, columns=HEADER_NO_FITNESS)

    pop = [triads_df.sample(n=n_triads).to_numpy(), permutations_df.sample(n=n_permutations).to_numpy()]
    pop = np.concatenate(pop, axis=0)

    pop = [TriadIndividual(i) for i in pop]
    fpop = np.asarray([x.f for x in pop])

    return pop, fpop


def population_init_random(task, population_size, rng, all_distances=False, angles=False, distance_categories=20,
                           angle_categories=20, **kwargs):
    # generate all possible permutations
    options = [[0, 1], [0, 1], [0, 1], list(range(0, distance_categories)),
               list(range(0, distance_categories))]  # atom types and 2 distances

    if all_distances:  # append list for third distance
        options.append(range(0, distance_categories))

    if angles:  # append lists for angles
        options.append(range(0, angle_categories))
        options.append(range(0, angle_categories))
        options.append(range(0, angle_categories))

    permutations = [list(a) for a in itertools.product(*options)]
    permutations_df = pd.DataFrame(permutations, columns=HEADER)

    pop = permutations_df.sample(n=population_size).to_numpy()  # pick random n=population_size elements

    fpop = np.apply_along_axis(task.eval, 1, pop)

    pop = list(zip(pop, fpop))
    pop = [TriadIndividual(i[0], i[1]) for i in pop]

    fpop = np.asarray([x.f for x in pop])
    return pop, fpop


def single_point_crossover(pop, ic, _cr, rng, task, new_pop, algorithm):
    r"""Single point crossover method.
    Args:
        pop (list): Current population.
        ic (int): Index of current individual.
        _cr (float): Crossover probability.
        rng (numpy.random.Generator): Random generator.
        task
        new_pop (list): Population newly created through crossover
        algorithm
    Returns:
        numpy.ndarray: New genotype or old individual if no crossover.
    """

    if rng.random() < _cr:
        io = ic
        while io == ic:
            io = rng.integers(len(pop))

        x = pop[ic]
        y = pop[io]

        crossover_index = random.randint(0, len(pop[ic].x) - 2)  # -2 to ensure at least one gene is changed

        # crossover
        tmp = pop[ic].copy()
        x[0:crossover_index] = y[0:crossover_index]
        y[0: crossover_index] = tmp[0:crossover_index]

        # mutation
        x.x = algorithm.mutation(pop, x, algorithm.mutation_rate, task, rng)
        y.x = algorithm.mutation(pop, x, algorithm.mutation_rate, task, rng)

        # new fitness
        x.evaluate(task, rng)
        y.evaluate(task, rng)

        # append offspring to new population list
        new_pop.append(x)
        new_pop.append(y)

    else:
        pass


def old_mutation(pop, individual, mr, task, rng, distance_categories=20, angle_categories=20):
    r"""Custom mutation method for 'old method' implementation.
    Args:
        pop (numpy.ndarray[Individual]): Current population.
        individual (TriadIndividual): Current individual.
        mr (float): Mutation probability.
        task (Task): Optimization task.
        rng (numpy.random.Generator): Random generator.
        distance_categories (int)
        angle_categories (int)
    Returns:
        numpy.ndarray: Mutated individual if rng.random < mr, else pop[ic].
    """
    if rng.random() < mr:
        new_gene = -1
        gene_index = random.randrange(0, task.dimension)
        # individual = pop[ic]

        if 0 <= gene_index < 2:  # NUC and ACID
            new_gene = 1 - individual[gene_index]  # swap 1 with 0 and vice versa
        elif gene_index == 2:  # BASE
            current_gene = individual[gene_index]
            other_genes = list(range(0, current_gene)) + list(range(current_gene + 1, 2))
            new_gene = random.choice(other_genes)
        elif 3 <= gene_index < 5:  # distances
            current_gene = individual[gene_index]
            other_genes = list(range(0, current_gene)) + list(range(current_gene + 1, distance_categories))
            new_gene = random.choice(other_genes)
        elif gene_index >= 5:  # angles
            current_gene = individual[gene_index]
            other_genes = list(range(0, current_gene)) + list(range(current_gene + 1, angle_categories))
            new_gene = random.choice(other_genes)

        if new_gene <= -1:
            raise ValueError("Can't mutate on index -1")
        else:
            individual[gene_index] = new_gene

    else:
        pass

    return individual.x


class GeneticAlgorithmModified(GeneticAlgorithm):

    def __init__(self, type, iteration, population_size=100, tournament_size=2, mutation_rate=0.25, crossover_rate=0.25,
                 selection=tournament_selection, crossover=uniform_crossover, mutation=uniform_mutation, *args,
                 **kwargs):
        """Initialize GeneticAlgorithm.
        Args:
            type (String): Algorithm name.
            iteration (str): Ordinal number of algorithm iteration in main.
            population_size (Optional[int]): Population size.
            tournament_size (Optional[int]): Tournament selection.
            mutation_rate (Optional[int]): Mutation rate.
            crossover_rate (Optional[float]): Crossover rate.
            selection (Optional[Callable[[numpy.ndarray[Individual], int, int, Individual, numpy.random.Generator], Individual]]): Selection operator.
            crossover (Optional[Callable[[numpy.ndarray[Individual], int, float, numpy.random.Generator], Individual]]): Crossover operator.
            mutation (Optional[Callable[[numpy.ndarray[Individual], int, float, Task, numpy.random.Generator], Individual]]): Mutation operator.
        See Also:
            * :func:`niapy.algorithms.Algorithm.set_parameters`
            * selection:
                * :func:`niapy.algorithms.basic.tournament_selection`
                * :func:`niapy.algorithms.basic.roulette_selection`
            * Crossover:
                * :func:`niapy.algorithms.basic.uniform_crossover`
                * :func:`niapy.algorithms.basic.two_point_crossover`
                * :func:`niapy.algorithms.basic.multi_point_crossover`
                * :func:`niapy.algorithms.basic.crossover_uros`
            * Mutations:
                * :func:`niapy.algorithms.basic.uniform_mutation`
                * :func:`niapy.algorithms.basic.creep_mutation`
                * :func:`niapy.algorithms.basic.mutation_uros`
        """
        super().__init__(population_size,
                         initialization_function=kwargs.pop('initialization_function', default_individual_init),
                         *args, **kwargs)
        self.type = type
        self.iteration = iteration
        self.tournament_size = tournament_size
        self.mutation_rate = mutation_rate
        self.crossover_rate = crossover_rate
        self.selection = selection
        self.crossover = crossover
        self.mutation = mutation

        self.population_list = []
        self.individual_type = TriadIndividual

    def run_task(self, task):
        r"""Start the optimization.
        Args:
            task (Task): Task with bounds and objective function for optimization.
        Returns:
            Tuple[numpy.ndarray, float]:
                1. Best individuals components found in optimization process.
                2. Best fitness value found in optimization process.
        See Also:
            * :func:`niapy.algorithms.Algorithm.iteration_generator`
        """
        algo, xb, fxb = self.iteration_generator(task), None, 0.0
        while not task.stopping_condition():
            xb, fxb = next(algo)
            task.next_iter()
        return xb, fxb

    def run_iteration(self, task, population, population_fitness, best_x, best_fitness, **params):
        r"""Core function of GeneticAlgorithm algorithm.
        Args:
            task (Task): Optimization task.
            population (numpy.ndarray): Current population.
            population_fitness (numpy.ndarray): Current populations fitness/function values.
            best_x (numpy.ndarray): Global best individual.
            best_fitness (float): Global best individuals function/fitness value.
            **params (Dict[str, Any]): Additional arguments.
        Returns:
            Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray, float, Dict[str, Any]]:
                1. New population.
                2. New populations function/fitness values.
                3. New global best solution
                4. New global best solution's fitness/objective value
                5. Additional arguments.
        """

        new_pop = []
        for i in range(self.population_size):
            ind_tmp = self.selection(population, i, self.tournament_size, best_x, self.rng)
            ind = TriadIndividual(ind_tmp.x)

            self.crossover(population, i, self.crossover_rate, self.rng, task=task,
                           new_pop=new_pop, algorithm=self)

        double_list = [population, new_pop]
        population_double = [item for sublist in double_list for item in sublist]

        best_x, best_fitness = self.get_best(ind, ind.f, best_x, best_fitness)

        population_reduced = sorted(population_double, key=operator.attrgetter('f'))[-self.population_size:]

        population_parameter_array = []
        for i in range(len(population_reduced)):
            individual = [parameter for parameter in population_reduced[i]]
            individual.append(population_reduced[i].f)

            population_parameter_array.append(individual)

        self.population_list = population_reduced
        # print(len(self.population_list))

        util.get_iteration_info(population=self.population_list,
                                header=HEADER, task=task)
        print(self)

        return population_reduced, np.asarray([i.f for i in population_reduced]), best_x, best_fitness, {}
