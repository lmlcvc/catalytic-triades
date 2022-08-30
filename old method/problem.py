import configparser
import itertools
import math
import os

import numpy as np
import pandas as pd
from niapy.algorithms import Individual
from niapy.problems.problem import Problem

__all__ = ['MostCommonPattern']

from niapy.util import objects_to_array

from numpy import shape
from numpy.random import default_rng

import util


HEADER = ['NUC', 'ACID', 'BASE', 'Dist_Nuc_Acid', 'Dist_Acid_Base']


class TriadIndividual(Individual):
    def __init__(self, x=None, task=None, e=True, rng=None, **kwargs):
        r"""Initialize new individual.
        Parameters:
            task (Optional[Task]): Optimization task.
            rand (Optional[numpy.random.Generator]): Random generator.
            x (Optional[numpy.ndarray]): Individuals components.
            e (Optional[bool]): True to evaluate the individual on initialization. Default value is True.
        """
        super().__init__(x=None, task=None, e=False, **kwargs)
        self.f = 0.0
        if x is not None:
            self.x = x if isinstance(x, np.ndarray) else np.asarray(x)
        if e and task is not None:
            self.evaluate(task, rng)



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
            merged_df = pd.concat(
                [self.triads_count_dict[protein].merge(x_df, on=x_df.columns.tolist(), how='inner'), merged_df])

        if merged_df.empty:
            return 0
        else:
            return merged_df.shape[0]
