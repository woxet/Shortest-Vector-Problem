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

global n # Lattice length
global timer_start
global S

FPLLL.set_random_seed(time.time())
random.seed(time.time())

# Initiation of the matrix
# In : n the lattice length
# Out: B a basis
def init_matrix(n):
	B = IntegerMatrix(n, n)
	B.randomize("uniform", bits=8)

	# We sort the base
	mini = 0
	tmp = IntegerMatrix(1,n)
	for x in range(n):
		for y in range(x+1,n):
			if B[y].norm() < B[mini].norm():
				mini = y
		for z in range(n):
			tmp[0,z] = B[x,z]
			B[x,z] = B[mini,z]
			B[mini,z] = tmp[0,z]
	return B

def moyenne(M):
	somme = 0
	for i in range(n):
		somme += M[i].norm()
	moy = somme/n
	return moy

def norm(v):
	somme = 0
	for i in range(n):
		somme += v[i]**2
	return math.sqrt(somme)

###################################
########## PERMUTATIONS ###########
###################################

# Swap 2 rows of a matrix
# In : M a matrix, i and j the indexes of the lines to swap
# Out : None
def swap_rows(M,i,j):
	tmp = IntegerMatrix(1,n)
	for k in range(n):
		tmp[0,k] = int(M[int(i),k])
	for k in range(n):
		M[i,k] = int(M[j,k])
	for k in range(n):
		M[j,k] = tmp[0,k]

# Swap 2 columns of a matrix
# In : M a matrix, i and j the indexes of the columns to swap
# Out : None
def swap_cols(M,i,j):
	tmp = IntegerMatrix(n,1)
	for k in range(n):
		tmp[k,0] = M[k,i]
	for k in range(n):
		M[k,i] = M[k,j]
	for k in range(n):
		M[k,j] = tmp[k,0]

# Transform a Column-Style Matrix into a Row-Style Matrix
# In : C the matrix you want to transform
# Out : L the matrix converted
def matrix_cols_into_rows(C):
	L = IntegerMatrix(C.ncols,C.nrows)
	for i in range(C.ncols):
		for j in range(C.nrows):
			L[i,j] = C[j,i]
	return L

# Transform a Row-Style Matrix into a Column-Style Matrix
# In : L the matrix you want to transform
# Out : C the matrix converted
def matrix_rows_into_cols(L):
	C = IntegerMatrix(L.ncols,L.nrows)
	for i in range(L.ncols):
		for j in range(L.nrows):
			C[i,j] = L[j,i]
	return C

# Returns the first non-zero pivot
# In : A a matrix, x an index
# Out : x or i the index of the row swapped
def pivot(A,x):
	p = A[x,x]
	col = 0
	while p == 0:
		for i in range(x,A.nrows):
			if A[i,x] != 0:
				swap_rows(A,x,i)
				return i
		col += 1
		swap_cols(A,x,x+col)
		p = A[x,x]
	else:
		return x

###################################
############### SVP ###############
###################################

# Reduction of the base by using Gram-Schmidt and LLL
# In : B a basis
# Out: L a reduced basis from B
def reduce_base(B):
	M = GSO.Mat(B) 			# Using Gram-Schmidt
	_ = M.update_gso()

	L = LLL.reduction(B)	# Using LLL-reduction
	return L

# Search an element in a list
# In : x the element to found, List the list where search in
# Out : boolean
def isIn(x,List):
	for i in range(len(List)):
		elem = List[i]
		if np.array_equal(x[0], elem):
			return True
	return False

# Depth-first search of the lattice [-5;+5]
# In : B a reduced basis
# Out: vMin a list of vectors with a minimal norm
def depth_first_search(B):
	vMin = []
	vMinX = []
	borne = 0
	for k in range(n): # Search for the largest vector of the basis
		if borne < B[k].norm():
			borne = B[k].norm()
	
	neigh = []
	v = IntegerMatrix(1, n)
	x = IntegerMatrix(1, n)

	for k in range(3):
		r = random.randint(0,n-1)
		x[0,r] += 1
		v[0].addmul(B[r])

	neigh.append(v[0])
	vMin, vMinX = neighbor(B, x, v, neigh, vMin, vMinX, borne)
	return vMin, vMinX # We build the neighbors 

# Generate neighbors of the vector
# In : B a reduced basis
# Out: neigh a list of candidates neighbors with a minimal norm | vMin a list of good vectors
def neighbor(B, x, v, neighPrev, vMin, vMinX, borne):
	for i in range(n):
		if (x[0,i] < -5 or x[0,i] > 5 or len(vMin) > 1.25*n):
			return vMin, vMinX
		x2 = IntegerMatrix(x)
		w = IntegerMatrix(1,n)
		if np.dot(v[0],B[i]) > 0:
			w[0].addmul(v[0])
			w[0].addmul(B[i], x=-1)
			x2[0,i] -= 1
		else:
			w[0].addmul(v[0])
			w[0].addmul(B[i])
			x2[0,i] += 1
		
		if w[0].norm() != 0 and not isIn(w,neighPrev):
			if w[0].norm() <= v[0].norm():
				neighPrev.append(w[0])
				if len(neighPrev) > 1000000:
					neighPrev.pop(0)
				if w[0].norm() <= borne:
					#print(w[0])
					vMin.append(w[0])
					vMinX.append(x2[0])
					if len(vMin) == 1.25*n:
						print("J'ai trouvé",len(vMin),"\"bon(s)\" vecteurs en",time.time()-timer_start)
				vMin, vMinX = neighbor(B, x2, w, neighPrev, vMin, vMinX, borne)

	return vMin, vMinX

####################################
############# HERMITE ##############
####################################

# Copy the lane i of A in B
# In : A and B a matrix, i an index 
# Out : None
def line_copy(A,B,i):
	for x in range(n):
		B[i,x] = A[i,x]

# Sort the list by the norm of his vectors
# In : L a list
# Out : None
def sortByNorm(L,Lx):
	mini = 0
	for x in range(len(L)):
		for y in range(x,len(L)):
			if L[y].norm() < L[mini].norm():
				mini = y
		tmp = L[x]
		tmpx = Lx[x]
		L[x] = L[mini]
		Lx[x] = Lx[mini]
		L[mini] = tmp
		Lx[mini] = tmpx

# We search if our vectors are linearly dependent by using the Hermite normal form
# In : M a matrix
# Out : N a matrix of vectors linearly independent
def HNF(M,Mx):
	N = IntegerMatrix(n,n)
	M2 = IntegerMatrix(M)
	for i in range(n):
		p = pivot(M2,i)
		swap_rows(M,p,i)
		line_copy(M,N,i)
		for j in range(i+1,n):
			d = math.gcd(abs(M2[i,i]),abs(M2[j,i]))
			a = int(M2[i,i]/d)
			b = int(M2[j,i]/d)
			for k in range(i,M2.nrows):
				M2[j,k] = a*M2[j,k]-b*M2[i,k]
	return N

####################################
############### MAIN ###############
####################################

n = 50

print("Création de la base..")
B = init_matrix(n)
print("Réduction de la base..")
B = reduce_base(B)
print("Base réduite.\nDébut du parcours du graphe..")

for i in range(2):
	# We test if the vector found is the same as the one proposed by BKZ
	S = SVP.shortest_vector(B)
	print(S)
	print("norme sv BKZ:", norm(S))

	timer_start = time.time()
	vMin, vMinX = depth_first_search(B)
	print("vMin len =", len(vMin))	
	sortByNorm(vMin, vMinX)
	print("sv :", vMin[0])
	print("norme sv:", vMin[0].norm())

	# Creation of the matrix m x n formed from vMin + the vectors of the base
	print("Création de la matrice m x n..")
	A = IntegerMatrix((len(vMin)+n), n) 	
	Ax = IntegerMatrix((len(vMin)+n), n)
	for x in range(len(vMin)):
		elem = vMin[x]
		elemx = vMinX[x]
		for y in range(n):
			A[x,y] = elem[y]
			Ax[x,y] = elemx[y]
	for x in range(n):
		Ax[len(vMin)+x,x] = 1
		for y in range(n):
			A[len(vMin)+x,y] = B[x,y]

	print("Début HNF..")
	A = matrix_cols_into_rows(A)
	Ax = matrix_cols_into_rows(Ax)

	A = HNF(A,Ax)

	A = matrix_rows_into_cols(A)
	Ax = matrix_rows_into_cols(Ax)
	print("Fin HNF.\nRéducton de la nouvelle base..")

	A = reduce_base(A)
	print("Base réduite.")

	print("Moyenne des normes des vecteurs de la base initiale :",moyenne(B))
	print("Moyenne des normes des vecteurs de la nouvelle base :",moyenne(A))