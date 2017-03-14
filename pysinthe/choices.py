import itertools
import bisect
import random


def weighted_choice(population, weights=None, *, cum_weights=None):
    """
    Adapted from Python 3.6 source code.

    Return a k sized list of population elements chosen with replacement.
    If the relative weights or cumulative weights are not specified,
    the selections are made with equal probability.
    """
    if cum_weights is None:
        if weights is None:
            _int = int
            total = len(population)
            return population[_int(random.random() * total)]
        cum_weights = list(itertools.accumulate(weights))
    elif weights is not None:
        raise TypeError('Cannot specify both weights and cumulative weights')
    if len(cum_weights) != len(population):
        raise ValueError('The number of weights does not match the population')
    total = cum_weights[-1]
    return population[bisect.bisect(cum_weights, random.random() * total)]
