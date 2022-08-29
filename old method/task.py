import os

import numpy as np
from matplotlib import pyplot as plt, ticker
from niapy.problems import Problem
from niapy.task import Task, OptimizationType, logger
from niapy.util import limit
from niapy.util.factory import get_problem


class TaskModified(Task):
    r"""Class representing an optimization task.
    Date:
        2019
    Author:
        Klemen Berkovič and others
    Attributes:
        problem (Problem): Optimization problem.
        dimension (int): Dimension of the problem.
        lower (numpy.ndarray): Lower bounds of the problem.
        upper (numpy.ndarray): Upper bounds of the problem.
        range (numpy.ndarray): Search range between upper and lower limits.
        optimization_type (OptimizationType): Optimization type to use.
        iters (int): Number of algorithm iterations/generations.
        evals (int): Number of function evaluations.
        max_iters (int): Maximum number of algorithm iterations/generations.
        max_evals (int): Maximum number of function evaluations.
        cutoff_value (float): Reference function/fitness values to reach in optimization.
        x_f (float): Best found individual function/fitness value.
    """

    def __init__(self, problem=None, dimension=None, lower=None, upper=None,
                 optimization_type=OptimizationType.MINIMIZATION, repair_function=limit, max_evals=np.inf,
                 max_iters=np.inf, cutoff_value=None, enable_logging=False):
        super().__init__(problem, dimension, lower, upper, optimization_type, repair_function, max_evals, max_iters,
                         cutoff_value, enable_logging)

    def plot_convergence(self, algo_type, iteration, output_directory, x_axis='iters', title='Convergence Graph'):
        """Plot a simple convergence graph.
        Args:
            algo_type (str)
            iteration (str)
            output_directory (str)
            x_axis (Literal['iters', 'evals']): Quantity to be displayed on the x-axis. Either 'iters' or 'evals'.
            title (str): Title of the graph.
        """
        x, fitness = self.convergence_data(x_axis)
        fitness *= -1
        _, ax = plt.subplots()
        ax.plot(x, fitness)
        ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
        if x_axis == 'iters':
            plt.xlabel('Iterations')
        else:
            plt.xlabel('Fitness Evaluations')
        plt.ylabel('Fitness')
        plt.title(title)
        plt.savefig(os.path.join(output_directory, algo_type + iteration + ".png"))
