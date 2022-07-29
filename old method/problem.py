import itertools

import numpy as np
import pandas as pd
from niapy.algorithms import Individual
from niapy.problems.problem import Problem

__all__ = ['MostCommonPattern']

from numpy import shape

HEADER = ['NUC', 'ACID', 'BASE', 'Dist_Nuc_Acid', 'Dist_Acid_Base']  # TODO: change depending on columns nr


def population_init(task, population_size, rng, all_distances=False, angles=False, distance_categories=20,
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

    print(permutations_df)

    # pick random n=population_size elements
    print(permutations_df.values)
    print(shape(permutations_df.values))

    pop = permutations_df.sample(n=population_size).to_numpy()

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
        self.__name__ = x
        self.f = f


class MostCommonPattern(Problem):

    def __init__(self, lower=0, upper=1000, dimension=5, triads_count=None, *args, **kwargs):
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

    def _evaluate(self, x):

        x_df = pd.DataFrame(data=x).T.astype(int)
        x_df.columns = HEADER

        merged_df = self.triads_count.merge(x_df, on=x_df.columns.tolist(), how='inner')

        if merged_df.empty:
            return 0
        else:
            print('Merged, Count = ' + str(merged_df['Count'].item()))
            return merged_df['Count'].item()
