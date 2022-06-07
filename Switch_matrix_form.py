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
			r = random.randint(-5,5)
			l[0,j] = r
		B[i].addmul(l[0])
	#print("Basis\n",B, end="\n\n")
	return B

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

n = 10
m = 20

C = init_matrix(n,m)
print(C,"\n")
L = matrix_cols_into_rows(C)
print(L,"\n")
C = matrix_rows_into_cols(L)
print(C,"\n")