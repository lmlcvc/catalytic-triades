import configparser
import os

import config as config
import numpy as np
from matplotlib import pyplot as plt, ticker
from niapy.problems import Problem
from niapy.task import Task, OptimizationType, logger
from niapy.util import limit
from niapy.util.factory import get_problem

import util


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
        r"""Initialize task class for optimization.
                Args:
                    problem (Union[str, Problem]): Optimization problem.
                    dimension (Optional[int]): Dimension of the problem. Will be ignored if problem is instance of the `Problem` class.
                    lower (Optional[Union[float, Iterable[float]]]): Lower bounds of the problem. Will be ignored if problem is instance of the `Problem` class.
                    upper (Optional[Union[float, Iterable[float]]]): Upper bounds of the problem. Will be ignored if problem is instance of the `Problem` class.
                    optimization_type (Optional[OptimizationType]): Set the type of optimization. Default is minimization.
                    repair_function (Optional[Callable[[numpy.ndarray, numpy.ndarray, numpy.ndarray, Dict[str, Any]], numpy.ndarray]]): Function for repairing individuals components to desired limits.
                    max_evals (Optional[int]): Number of function evaluations.
                    max_iters (Optional[int]): Number of generations or iterations.
                    cutoff_value (Optional[float]): Reference value of function/fitness function.
                    enable_logging (Optional[bool]): Enable/disable logging of improvements.
                """
        super().__init__(problem, dimension, lower, upper, optimization_type, repair_function, max_evals, max_iters,
                         cutoff_value, enable_logging)
        if isinstance(problem, str):
            params = dict(dimension=dimension, lower=lower, upper=upper)
            params = {key: val for key, val in params.items() if val is not None}
            self.problem = get_problem(problem, **params)
        elif isinstance(problem, Problem):
            self.problem = problem
            if dimension is not None or lower is not None or upper is not None:
                logger.warning(
                    'An instance of the Problem class was passed in, `dimension`, `lower` and `upper` parameters will be ignored.')
        else:
            raise TypeError('Unsupported type for problem: {}'.format(type(problem)))

        self.optimization_type = optimization_type
        self.dimension = self.problem.dimension
        self.lower = self.problem.lower
        self.upper = self.problem.upper
        self.range = self.upper - self.lower
        self.repair_function = repair_function

        self.iters = 0
        self.evals = 0

        self.cutoff_value = -np.inf * optimization_type.value if cutoff_value is None else cutoff_value
        self.enable_logging = enable_logging
        self.x_f = np.inf
        self.max_evals = max_evals
        self.max_iters = max_iters
        self.n_evals = []
        self.fitness_evals = []  # fitness improvements at self.n_evals evaluations
        self.fitness_iters = []  # best fitness at each iteration

    def stopping_condition(self):
        r"""Check if optimization task should stop.
        Returns:
            bool: `True` if number of function evaluations or number of algorithm iterations/generations or reference values is reach else `False`.
        """
        return (self.evals >= self.max_evals) or (self.iters >= self.max_iters)

    def eval(self, x):
        r"""Evaluate the solution A.
        Args:
            x (numpy.ndarray): Solution to evaluate.
        Returns:
            float: Fitness/function values of solution.
        """

        self.evals += 1
        x_f = self.problem.evaluate(x) * self.optimization_type.value

        if x_f < self.x_f * self.optimization_type.value:
            self.x_f = x_f * self.optimization_type.value
            self.n_evals.append(self.evals)
            self.fitness_evals.append(x_f)
            if self.enable_logging:
                logger.info('evals:%d => %s' % (self.evals, self.x_f))
        return x_f

    def convergence_data(self, algo_type, iteration, output_directory, x_axis='iters'):
        r"""Get values of x and y-axis for plotting covariance graph.
        Args:
            algo_type (str)
            iteration (str)
            output_directory (str)
            x_axis (Literal['iters', 'evals']): Quantity to be displayed on the x-axis. Either 'iters' or 'evals'.
        Returns:
            Tuple[np.ndarray, np.ndarray]:
                1. array  of function evaluations.
                2. array of fitness values.
        """
        if x_axis == 'iters':
            return np.arange(self.iters), np.array(self.fitness_iters)
        else:  # x_axis == 'evals'
            util.store_fitness_trends(fitness=self.fitness_evals,
                                      filepath=os.path.join(output_directory, algo_type + iteration + ".csv"))

            r1, r2 = [], []
            for i, v in enumerate(self.n_evals):
                r1.append(v)
                r2.append(self.fitness_evals[i])
                if i >= len(self.n_evals) - 1:
                    break
                diff = self.n_evals[i + 1] - v
                if diff <= 1:
                    continue
                for j in range(diff - 1):
                    r1.append(v + j + 1)
                    r2.append(self.fitness_evals[i])
            return np.array(r1), np.array(r2)

    def plot_convergence(self, algo_type, iteration, output_directory, x_axis='iters', title='Convergence Graph'):
        """Plot a simple convergence graph.
        Args:
            algo_type (str)
            iteration (str)
            output_directory (str)
            x_axis (Literal['iters', 'evals']): Quantity to be displayed on the x-axis. Either 'iters' or 'evals'.
            title (str): Title of the graph.
        """
        config = configparser.ConfigParser()
        config.read(os.path.join(os.pardir, 'config.ini'))
        config = config['default']

        x, fitness = self.convergence_data(algo_type, iteration, config["fitness"], x_axis="evals")
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
