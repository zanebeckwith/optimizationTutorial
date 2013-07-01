#!/usr/bin/env python

"""
This file defines a python function 
that is intended to be used as part of
a tutorial on mathematical optimization.

Author: Zane Beckwith

In the toy example given in the tutorial,
we imagine that this function actually represents
simulation software for nuclear cross-section
calculations. A driver program ('find_best_params.py')
uses this function along with a Nelder-Mead
optimizer method (included in this directory)
to find the 'parameter values' 
('energy' and 'helicity')
that give the largest cross-section.
"""

from math import sqrt

def simulate(energy, helicity):
	return sqrt(\
				(energy**3 - 3*energy*(helicity**2) - 1)**2 + \
				(3*helicity*(energy**2) - helicity**3)**2 \
			) - \
			energy + 10
