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
import time

FPLLL.set_random_seed(time.time())

global n,m,dec

# Initiation of the matrix
# In : n the lattice length
# Out: B a basis
def init_matrix(m,n):
	B = IntegerMatrix(m,n)
	for i in range(m):
		l = IntegerMatrix(1,n)
		for j in range(n):
			r = random.randint(-5,5)
			l[0,j] = r
		B[i].addmul(l[0])
	#print("Basis\n",B, end="\n\n")
	return B

# Swap 2 rows of a matrix
# In : M a matrix, i and j the indexes of the lines to swap
# Out : None
def swap_rows(M,i,j):
	tmp = IntegerMatrix(1,n)
	for k in range(n):
		tmp[0,k] = M[i,k]
	for k in range(n):
		M[i,k] = int(M[j,k])
	for k in range(n):
		M[j,k] = tmp[0,k]

def pivot(A,x):
	if A[x+dec,x] == 0:
		for i in range(x+dec,n):
			if A[i,x] != 0:
				swap_rows(A,x,i)
				break

def HNF(M):
	for i in range(m):
		pivot(M,i)
		p = M[i+dec,i]
		for j in range(1,n):
			#print(M,"\n")
			
			d = math.gcd(abs(p),abs(M[i+dec,j]))
			#print(d,abs(p),abs(M[i+dec,j]))
			a = int(p/d)
			b = int(M[i+dec,j]/d)
			#print(M,"\n")
			for h in range(j,n):
				print(a,M[i+dec,i],b,M[i+dec,h])
				M[i+dec,h] = a*M[i+dec,i]-b*M[i+dec,h]
			#print(M,"\n\n")
	return M

n = 10
m = 20
dec = 0

M = init_matrix(m,n)
print(M,"\n")
M = HNF(M)

print(M)


