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

# Depth-first search of the lattice [-5;+5]
# In : B a reduced basis
# Out : vMin a list of vectors with a minimal norm
def depth_first_search(B):
	vMin = []
	vMax = IntegerMatrix(1,n)
	# Choose randomly a vector in the base and define it as the largest
	r = random.randint(0,n-1)
	print("Random :", r)
	vMax = np.copy(B[r])
	print("vMax =", vMax)
	#print(vMax[0].norm()**2)

	# Randomly create a vector in the lattice
	alea = IntegerMatrix(1, n)
	for k in range(3):
		r = random.randint(0,n-1)
		alea[0].addmul(B[k])
	print("alea =", alea)
	#print(type(alea))

	neighborList = neighbor(B, alea) # We build the neighbors 
	print(type(neighborList[0]))

	return vMin

def neighbor(B, v):
	nList = [] # List that contains the neighbors of the vector in parameter
	for i in range(n):
		for k in range(-5,5):
			neigh = v[0]
			#neigh = np.copy(v) # A neighbor 
			neigh.addmul(B[i],k)
			nList.append(neigh)
	return nList

#######################################
##############  CODE  #################
#######################################

n = 10
B = init_matrix(n)
L = reduce_base(B)
vMin = depth_first_search(L)

#S = SVP.shortest_vector(L)
#print(S)




