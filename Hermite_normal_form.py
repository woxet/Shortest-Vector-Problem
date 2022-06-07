#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 31 2022

@author: Jordan
"""

from fpylll import IntegerMatrix, FPLLL, GSO
import numpy as np
import random
import math
import time

FPLLL.set_random_seed(time.time())

global n,m

# Initiation of the matrix
# In : n the lattice length
# Out: B a basis
def init_matrix(n,m):
	B = IntegerMatrix(n,m)
	for i in range(n):
		l = IntegerMatrix(1,m)
		for j in range(m):
			r = random.randint(-1,1)
			l[0,j] = r
		B[i].addmul(l[0])
	return B

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

# Returns the first non-zero pivot
# In : A a matrix, x an index
# Out : x or i the index swapped
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

def HNF(M):
	for i in range(n):
		p = pivot(M,i)
		for j in range(i+1,n):
			d = math.gcd(abs(M[i,i]),abs(M[j,i]))
			a = int(M[i,i]/d)
			b = int(M[j,i]/d)
			for k in range(i,M.nrows):
				M[j,k] = a*M[j,k]-b*M[i,k]
	return M

n = 10
m = 20

M = init_matrix(n,m)
print(M,"\n")
M = HNF(M)

print(M)