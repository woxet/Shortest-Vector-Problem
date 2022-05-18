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

import csv,math
from pathlib import Path

global n # Lattice length
FPLLL.set_random_seed(time())
random.seed(time())

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
	#print("Reduced Basis\n",L, end="\n\n")
	#print("Is reduced ?",LLL.is_reduced(L))
	return L

# Depth-first search of the lattice [-5;+5]
# In : B a reduced basis
# Out: vMin a list of vectors with a minimal norm
def depth_first_search(B):
	vMin = []
	borne = 0
	for k in range(n): # Search for the largest vector of the basis
		if borne < B[k].norm():
			borne = B[k].norm()
		#borne += B[k].norm() 
	#borne = borne/n
	
	neigh = []
	v = IntegerMatrix(1, n)
	x = []
	for i in range(n):
		x.append(0)

	for k in range(3):
		r = random.randint(0,n-1)
		x[r] += 1
		v[0].addmul(B[r])

	neigh.append(v[0])
	vMin = neighbor(B, x, v, neigh, vMin, borne)
	print("vMin len =", len(vMin))
	return vMin # We build the neighbors 

def isIn(x,List):
	for i in range(len(List)):
		elem = List[i]
		if np.array_equal(x[0], elem):
			#print("x:",x, "\n", "elem:", elem,"\n\n")
			return True
	return False

# Generate neighbors of the vector
# In : B a reduced basis
# Out: neigh a list of candidates neighbors with a minimal norm | vMin a list of good vectors
def neighbor(B, x, v, neighPrev, vMin, borne):
	for i in range(n):
		if x[i] < -5 or x[i] > 5:
			#print(x)
			return vMin
		x2 = x.copy()
		w = IntegerMatrix(1, n)
		#print(np.dot(v[0],B[i]))
		if np.dot(v[0],B[i]) > 0:
			w[0].addmul(v[0])
			w[0].addmul(B[i], x=-1)
			x2[i] -= 1
			#print(w[0])
		else:
			w[0].addmul(v[0])
			w[0].addmul(B[i])
			x2[i] += 1
		
		if w[0].norm() != 0 and not isIn(w,neighPrev):
			if w[0].norm() <= v[0].norm():
				#print(w[0])
				neighPrev.append(w[0])
				if len(neighPrev) > 100000:
					neighPrev.pop(0)
				if w[0].norm() <= borne:
					print(w[0])
					vMin.append(w[0])
				#print("on descend\n\n")
				vMin = neighbor(B, x2, w, neighPrev, vMin, borne)

	return vMin

####################################
############### CODE ###############
####################################

n = 40
essais = 1
somme = 0
vMin = []
B = init_matrix(n)
L = reduce_base(B)
#print(L)

vMin = depth_first_search(L)

print("vMin len =", len(vMin))

sv = vMin[0]
for k in range(1,len(vMin)):
	if sv.norm() > vMin[k].norm():
		sv = vMin[k]

print("sv :", sv)

S = SVP.shortest_vector(L)
print(S)

 


