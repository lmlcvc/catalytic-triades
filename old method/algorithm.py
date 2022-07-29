import random

import numpy as np
from niapy.algorithms import Individual, default_individual_init
from niapy.algorithms.basic import GeneticAlgorithm
from niapy.algorithms.basic.ga import tournament_selection, uniform_crossover, uniform_mutation


def single_point_crossover(pop, ic, _cr, rng):
    r"""Single point crossover method.
    Args:
        pop (numpy.ndarray[Individual]): Current population.
        ic (int): Index of current individual.
        _cr (float): Crossover probability.
        rng (numpy.random.Generator): Random generator.
    Returns:
        numpy.ndarray: New genotype or old individual if no crossover.
    """

    if rng.random() < _cr:
        io = ic
        while io == ic:
            io = rng.integers(len(pop))

        x = pop[ic].x

        crossover_index = random.randint(0, len(pop[ic].x) - 2)  # -2 to ensure at least one gene is changed

        tmp = x.copy()
        x[0:crossover_index] = pop[io].x[0:crossover_index]
        pop[io].x[0:crossover_index] = tmp[0:crossover_index]

        return np.asarray(x)
    else:
        return pop[ic].x


def old_mutation(pop, ic, mr, task, rng, distance_categories=20, angle_categories=20):
    r"""Custom mutation method for 'old method' implementation.
    Args:
        pop (numpy.ndarray[Individual]): Current population.
        ic (int): Index of current individual.
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
        individual = pop[ic]

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

        return individual.x
    else:
        return pop[ic].x


class GeneticAlgorithmModified(GeneticAlgorithm):

    def __init__(self, population_size=25, tournament_size=2, mutation_rate=0.25, crossover_rate=0.25,
                 selection=tournament_selection, crossover=uniform_crossover, mutation=uniform_mutation, *args,
                 **kwargs):
        """Initialize GeneticAlgorithm.
        Args:
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
                         individual_type=kwargs.pop('individual_type', Individual),
                         initialization_function=kwargs.pop('initialization_function', default_individual_init),
                         *args, **kwargs)
        self.tournament_size = tournament_size
        self.mutation_rate = mutation_rate
        self.crossover_rate = crossover_rate
        self.selection = selection
        self.crossover = crossover
        self.mutation = mutation

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
                4. New global best solutions fitness/objective value
                5. Additional arguments.
        """
        print('---------- ITERATION ----------')

        new_pop = np.empty(self.population_size, dtype=object)
        for i in range(self.population_size):

            ind_tmp = self.selection(population, i, self.tournament_size, best_x, self.rng)
            ind = self.individual_type(ind_tmp.x, ind_tmp.f)
            ind.x = self.crossover(population, i, self.crossover_rate, self.rng)
            ind.x = self.mutation(population, i, self.mutation_rate, task, self.rng)
            ind.evaluate(task, rng=self.rng)
            new_pop[i] = ind
            if new_pop[i].f < best_fitness:
                best_x, best_fitness = self.get_best(new_pop[i], new_pop[i].f, best_x, best_fitness)

        return new_pop, np.asarray([i.f for i in new_pop]), best_x, best_fitness, {}