import configparser
import itertools
import math
import os

import numpy as np
import pandas as pd
from niapy.algorithms import Individual
from niapy.problems.problem import Problem

__all__ = ['MostCommonPattern']

from numpy import shape

import util

HEADER = ['NUC', 'ACID', 'BASE', 'Dist_Nuc_Acid', 'Dist_Acid_Base']  # TODO: change depending on columns nr

config = configparser.ConfigParser()
config.read(os.path.join(os.pardir, 'config.ini'))
config = config['default']
encoded_directory = config['encoded_location']  # TODO: adapt to other methods if necessary


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
    permutations_df = pd.DataFrame(permutations, columns=HEADER)

    pop = [triads_df.sample(n=n_triads).to_numpy(), permutations_df.sample(n=n_permutations).to_numpy()]
    pop = np.concatenate(pop, axis=0)

    fpop = np.apply_along_axis(task.eval, 1, pop)

    pop = list(zip(pop, fpop))
    pop = [TriadIndividual(i[0], i[1]) for i in pop]

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


class TriadIndividual(Individual):
    def __new__(cls, x, f):
        return super(TriadIndividual, cls).__new__(cls)

    def __init__(self, x, f, **kwargs):
        super().__init__(x, **kwargs)
        self.x = x
        self.f = f

    def fitness(self):
        return self.f


class MostCommonPattern(Problem):

    def __init__(self, lower=0, upper=1000, dimension=5, triads_count=None, method='old', *args, **kwargs):
        r"""Initialize MostCommonPattern problem.
        Args:
            dimension (Optional[int]): Dimension of the problem.
            lower (Optional[Union[float, Iterable[float]]]): Lower bounds of the problem.
            upper (Optional[Union[float, Iterable[float]]]): Upper bounds of the problem.
        See Also:
            :func:`niapy.problems.Problem.__init__`
        """

        super().__init__(dimension, lower, upper, *args, **kwargs)
        self.triads_count = triads_count.astype(int)
        self.method = method

    def _evaluate(self, x):

        x_df = pd.DataFrame(data=x).T.astype(int)
        x_df.columns = HEADER

        merged_df = self.triads_count.merge(x_df, on=x_df.columns.tolist(), how='inner')

        if merged_df.empty:
            return 0
        else:
            return merged_df['Count'].item()


class EnzymeCommonPattern(Problem):

    def __init__(self, lower=0, upper=1000, dimension=5, triads_count=None, triads_count_dict=None, method='old', *args,
                 **kwargs):
        r"""Initialize EnzymeCommonPattern problem.
        Args:
            dimension (Optional[int]): Dimension of the problem.
            lower (Optional[Union[float, Iterable[float]]]): Lower bounds of the problem.
            upper (Optional[Union[float, Iterable[float]]]): Upper bounds of the problem.
        See Also:
            :func:`niapy.problems.Problem.__init__`
        """

        super().__init__(dimension, lower, upper, *args, **kwargs)
        self.triads_count = triads_count
        self.triads_count_dict = triads_count_dict
        self.method = method

    def _evaluate(self, x):

        x_df = pd.DataFrame(data=x).T.astype(int)
        x_df.columns = HEADER

        merged_df = pd.DataFrame()

        for protein in self.triads_count_dict.keys():
            protein = protein.replace(".csv", "")
            print(type(self.triads_count[protein]))
            print(self.triads_count[protein])
            print(type(self.triads_count_dict[protein]))
            print(self.triads_count_dict[protein])
            merged_df = pd.concat(
                [self.triads_count_dict[protein].merge(x_df, on=x_df.columns.tolist(), how='inner'), merged_df])

        if merged_df.empty:
            return 0
        else:
            return merged_df.shape[0]
