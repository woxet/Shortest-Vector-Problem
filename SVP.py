#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 2022

@author: Jordan
"""

from fpylll import IntegerMatrix, LLL, GSO, SVP, FPLLL
import numpy as np
import random
from time import time

global n # Lattice length
FPLLL.set_random_seed(time())
random.seed(time())

# Initiation of the matrix
# In : n the lattice length
# Out : B a basis
def init_matrix(n):
	B = IntegerMatrix(n, n)
	B.randomize("uniform", bits=8)
	print("Basis\n",B, end="\n\n")
	return B

# Reduction of the base by using Gram-Schmidt and LLL
# In : B a basis
# Out : L a reduced basis from B
def reduce_base(B):
	M = GSO.Mat(B) 			# Using Gram-Schmidt
	_ = M.update_gso()

	L = LLL.reduction(B)	# Using LLL-reduction
	print("Reduced Basis\n",L, end="\n\n")
	print("Is reduced ?",LLL.is_reduced(L))
	return L

#######################################
##############  CODE  #################
#######################################

n = 2
B = init_matrix(n)
L = reduce_base(B)

#S = SVP.shortest_vector(L)
#print(S)




