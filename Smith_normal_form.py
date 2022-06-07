#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 16 2022

@author: Jordan
"""

from fpylll import IntegerMatrix, FPLLL, GSO
import numpy as np
import random
import math

from time import time
FPLLL.set_random_seed(time())

###########################################
############## PERMUTATION ################
###########################################
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

###########################################
###########################################

# Initiation of the matrix
# In : n the lattice length
# Out: B a basis
def init_matrix(n,m):
	B = IntegerMatrix(n,m)
	for i in range(n):
		l = IntegerMatrix(1,m)
		for j in range(m):
			r = random.randint(-5,5)
			l[0,j] = r
		B[i].addmul(l[0])
	return B

# Returns the first non-zero pivot
# In : A a matrix, x an index
# Out : None
def pivot(A,x):
	p = A[x,x]
	for i in range(x,A.nrows):
		if A[i,x] != 0:
			swap_rows(A,x,i)
			break

def improving_pivot(A):
	for p in range(n):
		pivot(A,p)
		for j in range(p+1,n):
			d = math.gcd(abs(A[p,p]), abs(A[j,p]))
			a = int(A[p,p]/d)
			b = int(A[j,p]/d)
			for k in range(j,n):
				A[j,p] = a*A[j,p]-b*A[p,p]
	print(A,"\n")
	"""for p in range(n):
		for i in range(n-1,p,-1):
			d = math.gcd(abs(A[p,p]), abs(A[p,i]))
			a = int(A[p,p]/d)
			b = int(A[p,i]/d)
			A[p,i] = a*A[p,i]-b*A[p,p]

	print(A,"\n")"""

global k,n
k = 15
n = 10

A = init_matrix(k,n)
print(A,"\n")
#pivot(A)
improving_pivot(A)
d = IntegerMatrix(1,n)
for x in range(n):
	d[0,x] = A[x,x]
print(d)