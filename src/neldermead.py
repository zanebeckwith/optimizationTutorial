#!/usr/bin/env python
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


def max_distance(points):
    """
    Takes a list of
    Numpy arrays and finds
    the largest (Euclidean) distance between them

    Argument:
    points  --  list of Numpy arrays
    Return value:
            --  float: maximum distance
                between points
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
    func              -- objective function. Must accept length-n Numpy array
    max_iterations    -- int: maximum number of iterations used
    tolerance         -- float: max. distance between points in final simplex
    Return value :
    xmin              -- Numpy array: location of the minimum
    """

    # Set constants controlling:
    alpha = 1.0     # reflection
    gamma = 2.0     # expansion
    rho = 0.5       # contraction
    sigma = 0.5     # reduction

    # Find function values for initial simplex points,
    # and combine these together with the points.
    #
    # The result is a list of
    # [function value, point]
    values_and_points = [[func(p), p] for p in initial_simplex]

    # Initialize iteration counter
    iterations = 0

    while max_distance([p[1] for p in values_and_points]) > tolerance:
        # Increment iteration counter.
        # If we've already done to many, quit loop.
        if iterations == max_iterations:
            print 'minimize: reached maximum number of iterations'
            break
        iterations += 1

        # Sort function values and points
        # by function value.
        #
        # Separate the highest-value point
        # from the others.
        values_and_points.sort()
        lowest_values_and_points = values_and_points[:-1]
        highest_value_and_point = values_and_points[-1]

        # Find the center of gravity of the lowest-value points.
        #
        # The center of gravity is the weighted average
        # (weighted by function value) of the points.
        #
        # Movements of the simplex's points will be relative to
        # this point.
        center_of_gravity = sum(
            [np.multiply(p[0], p[1]) for p in lowest_values_and_points]
            ) / sum([p[0] for p in lowest_values_and_points])

        # Reflect simplex.
        #
        # Reflect highest-value point across center of gravity
        # of other points.
        #
        # If this new point is better than the old highest-value point,
        # but is still not the absolute lowest-value point,
        # replace the old highest-value point with the new point
        # and restart the loop.
        reflected_point = center_of_gravity + \
            alpha * (center_of_gravity - highest_value_and_point[1])
        reflected_value = func(reflected_point)
        if (values_and_points[0][0] <=
                reflected_value <=
                values_and_points[-1][0]):
            values_and_points[-1] = [reflected_value, reflected_point]
            continue

        # Expand simplex.
        #
        # If the reflected point found above is the lowest-value
        # point, move it even more in the same direction (gamma > alpha)
        # to see if it will keep getting better. If it does,
        # keep this new farther point and restart the loop;
        # if not, just keep the original reflected point
        # and restart the loop.
        #
        # 'Keeping' a point means replacing the highest-value point.
        if reflected_value < values_and_points[0][0]:
            expanded_point = center_of_gravity + \
                gamma * (center_of_gravity - highest_value_and_point[1])
            expanded_value = func(expanded_point)
            if expanded_value < reflected_value:
                values_and_points[-1] = [expanded_value, expanded_point]
            else:
                values_and_points[-1] = [reflected_value, reflected_point]
            continue

        # Contract the simplex.
        #
        # If we're still in the loop, the reflected point
        # must be higher-value than any of the points other
        # than the highest-value point (and it may be higher than that);
        # thus, back-track from the reflection direction.
        # As usual, if that's a good choice, keep it and restart the loop.
        contracted_point = values_and_points[-1][1] + \
            rho * (center_of_gravity - highest_value_and_point[1])
        contracted_value = func(contracted_point)
        if contracted_value < values_and_points[-1][0]:
            values_and_points[-1] = [contracted_value, contracted_point]
            continue

        # Reduce the simplex.
        #
        # It's generally rare that we get this far in the loop.
        # If we do, none of the movements have been successful,
        # so we just contract ALL points except the lowest-value
        # one around that lowest-value point. We hit the end of
        # the loop, and start over again.
        for p in values_and_points:
            p[1] = values_and_points[0][1] + \
                np.multiply(sigma, p[1] - values_and_points[0][1])
            p[0] = func(p[1])

    # Once the loop is finished, either because we found a
    # good simplex or because we ran out of iterations,
    # find the lowest-value point in the final simplex
    # and return it.
    return sorted(values_and_points)[0][1]


if __name__ == '__main__':

    # Define a simple function with a well-known minimum,
    # and make sure we find it.

    def parabola(params):
        return (params[0] - 0.5) ** 2 + (params[1] + 0.3) ** 2

    known_min = np.array([0.5, -0.3])
    assert parabola(known_min) == 0.0

    i_simp = [np.array([-2.0, -2.0]),
              np.array([0.0, 2.0]),
              np.array([2.0, -2.0])
              ]
    max_iters = 500
    tol = 1e-05
    min_pos = minimize(parabola, i_simp, max_iters, tol)
    assert max_distance([min_pos, known_min]) < tol
