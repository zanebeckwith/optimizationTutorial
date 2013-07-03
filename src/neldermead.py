"""
This file implements a function to find the location of the minimum
for a given well-behaved non-linear function in arbitrary dimensions
using the Nelder-Mead (downhill simplex) method.

Originally written for the UNC-CH CAP-REU program, 
summer of 2013.

Author: Zane Beckwith
Adapted from code provided by Dr. Fabian Heitsch.

"""

import numpy as np
import sys

def max_distance(points):
  """
  Takes a list of 
  Numpy arrays and finds
  the largest (Euclidean) distance between them

  function: max_distance
  input :
  """
  from math import sqrt
  distances = []
  for (index1, p1,) in enumerate(points):
    for p2 in points[index1:]:
      distances.append(sqrt(sum((p1 - p2) ** 2)))

  return max(distances)

def minimize(func, initial_simplex, max_iterations, tolerance):
    """
    Finds the location of the minimum of given function
    using the Nelder-Mead method.
  
    Arguments:
    initial_simplex   -- list with n+1 members; each member
                        is a Numpy array of n values, 
                        where n is the dimension of func()
    func              -- the objective function. Must accept length-n Numpy array
    max_iterations    -- maximum number of iterations used
    tolerance         -- maximum distance between points in final simplex
    Return value :
    xmin              -- Location of the minimum as Numpy array
    """

    alpha = 1.0
    gamma = 2.0
    rho = 0.5
    sigma = 0.5

    values_and_points = [ [func(p), p] for p in initial_simplex ]

    iterations = 0

    while max_distance([ p[1] for p in values_and_points ]) > tolerance:
      if iterations == max_iterations:
        print 'minimize: reached maximum number of iterations'
        break
      iterations += 1
      values_and_points.sort()
      lowest_values_and_points = values_and_points[:-1]
      highest_value_and_point = values_and_points[-1]
      center_of_gravity = sum( \
                             [np.multiply(p[0], p[1]) for p in lowest_values_and_points] \
                             ) / sum([p[0] for p in lowest_values_and_points])
      reflected_point = center_of_gravity + \
          alpha * (center_of_gravity - highest_value_and_point[1])
      reflected_value = func(reflected_point)
      if values_and_points[0][0] <= reflected_value <= values_and_points[-1][0]:
        values_and_points[-1] = [reflected_value, reflected_point]
        continue
      if reflected_value < values_and_points[0][0]:
        expanded_point = center_of_gravity + \
            gamma * (center_of_gravity - highest_value_and_point[1])
        expanded_value = func(expanded_point)
        if expanded_value < reflected_value:
          values_and_points[-1] = [expanded_value, expanded_point]
        else:
          values_and_points[-1] = [reflected_value, reflected_point]
        continue
      contracted_point = values_and_points[-1][1] + \
          rho * (center_of_gravity - highest_value_and_point[1])
      contracted_value = func(contracted_point)
      if contracted_value < values_and_points[-1][0]:
        values_and_points[-1] = [contracted_value, contracted_point]
        continue
      for p in values_and_points:
        p[1] = values_and_points[0][1] + \
            np.multiply(sigma, p[1] - values_and_points[0][1])
        p[0] = func(p[1])

    return sorted(values_and_points)[0][1]


if __name__ == '__main__':

    def parabola(params):
        return (params[0] - 0.5) ** 2 + (params[1] + 0.3) ** 2

    known_min = np.array([0.5, -0.3])
    assert parabola(known_min) == 0.0
    i_simp = [np.array([-2.0, -2.0]), np.array([0.0, 2.0]), np.array([2.0, -2.0])]
    max_iters = 500
    tol = 1e-05
    min_pos = minimize(parabola, i_simp, max_iters, tol)
    assert max_distance([min_pos, known_min]) < tol

