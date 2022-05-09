#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 2022

@author: Jordan
"""

from fpylll import IntegerMatrix, LLL, GSO, SVP, FPLLL
import numpy as np
import random
import time

import csv,math
from pathlib import Path

global n # Lattice length
FPLLL.set_random_seed(48623)
random.seed(53265)

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
	
	neigh = []
	v = IntegerMatrix(1, n)
	x = []
	for i in range(n):
		x.append(0)

	for k in range(3):
		r = random.randint(0,n-1)
		x[r] += 1
		v[0].addmul(B[r])

	neigh, vMin = neighbor(B, x, v, neigh, borne)
	print("vMin len =", len(vMin))
	return vMin # We build the neighbors 

def isIn(x,List):
	for i in range(len(List)):
		elem = List[i]
		if np.array_equal(x[0], elem):
			return True
		"""
		for j in range(n):
			if j == n-1 and x[0,j] == elem[j]:
				#print("x:",x, "\n", "elem:", elem,"\n\n")
				return True 
			if x[0,j] != elem[j]:
				break
		"""
	return False

# Generate neighbors of the vector
# In : B a reduced basis
# Out: neigh a list of candidates neighbors with a minimal norm | vMin a list of good vectors
def neighbor(B, x, v, neighPrev, borne):
	vMin = []
	vMin2 = []
	neigh = list(neighPrev)
	neigh2 = []

	neigh.append(v[0])

	for i in range(n):
		if x[i] < -3 or x[i] > 3:
			#print(x)
			return neigh2, vMin2
		x2 = list(x)
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
		
		if w[0].norm() != 0 and not isIn(w,neigh):
			if w[0].norm() <= v[0].norm():
				#print(w[0])
				neigh.append(w[0])
				if w[0].norm() <= borne:
					print(w[0])
					vMin.append(w[0])
				#print("on descend\n\n")
				neigh2, vMin2 = neighbor(B, x2, w, neigh, borne)
		#print(x)
	#print(len(neigh))

	return neigh+neigh2, vMin+vMin2

#####################################
############### STATS ###############
#####################################

def stats(somme,essais):
	moy = somme/essais
	var = (1/essais)*(somme**2)-(moy**2)
	ec = math.sqrt(var)
	print(moy)
	print(ec)
	with open('stats.csv','a') as fichier:
		writer = csv.writer(fichier)
		elem = [moy,ec]
		writer.writerow(elem)

####################################
############### CODE ###############
####################################

n = 15
essais = 1
somme = 0

"""
fileName = r"stats.csv"
fileObj = Path(fileName)
if not fileObj.is_file():
	with open('stats.csv','w') as fichier:
		writer = csv.writer(fichier)
		elem = ['Moyenne','Ecart-Type']
		writer.writerow(elem)
"""
for x in range(essais):
	B = init_matrix(n)
	L = reduce_base(B)
	vMin = depth_first_search(L)
	#somme += len(vMin)

print("vMin len =", len(vMin))
for i in range(len(vMin)):
	print(vMin[i],end="\n\n")
#stats(somme,essais)

S = SVP.shortest_vector(L)
print(S)

 


