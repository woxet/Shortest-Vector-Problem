#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 2022

@author: Jordan
"""

from fpylll import IntegerMatrix, LLL, GSO, SVP, FPLLL
import numpy as np
import random
import math
import time

FPLLL.set_random_seed(time.time())
random.seed(time.time())

# Initiation of the matrix
# In : n the lattice length
# Out: B a basis
def init_matrix(n):
	B = IntegerMatrix(n, n)
	B.randomize("uniform", bits=8)
	#print("Basis\n",B, end="\n\n")
	return B

# Reduction of the base by using Gram-Schmidt and LLL
# In : B a basis
# Out: L a reduced basis from B
def reduce_base(B):
	M = GSO.Mat(B) 			# Using Gram-Schmidt
	_ = M.update_gso()

	L = LLL.reduction(B)	# Using LLL-reduction
	return L


global n
n = 10
B = init_matrix(n)
B = reduce_base(B)

S = SVP.shortest_vector(B)

for i in range(n):
	print(type(S[i]))


