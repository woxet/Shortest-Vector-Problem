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
def swap_rows(M,i,j):
	for l in range(n):
		tmp = M[i,l]
		M[i,l] = M[j,l]
		M[j,l] = tmp

def swap_columns(M,i,j):
	for l in range(k):
		tmp = M[l,i]
		M[l,i] = M[l,j]
		M[l,j] = tmp

def add_to_row(M,p,d):
	a = int(M[p,p]/d)
	b = int(M[p+1,p+1]/d)
	for x in range(p+1,n):
		print(M[p+1,x])
		M[p+1,x] = a*M[p+1,x]-b*M[p,x]
		print(M[p+1,x],"\n")

def add_to_column(M,p,d):
	a = int(M[p,p]/d)
	b = int(M[p+1,p]/d)
	for x in range(p+1,k):
		M[x,p] = a*M[x,p]-b*M[p,p]

def change_sign_row(M,x):
	for y in range(n):
		M[x][y] = - M[x][y]

def change_sign_column(M,x):
	for y in range(k):
		M[y][x] = - M[y][x]

###########################################
###########################################

def init_matrix(k,n):
	B = IntegerMatrix(k, n)
	for i in range(k):
		l = IntegerMatrix(1,n)
		for j in range(n):
			r = random.randint(-5,5)
			l[0,j] = r
		B[i].addmul(l[0])
	#print("Basis\n",B, end="\n\n")
	return B

# On prend la plus petite entrée de A non-zero
def pivot(A,x):
	if A[x,x] == 0:
		for i in range(n):
			if A[i,0] != 0:
				swap_rows(A,x,i)
				break

def improving_pivot(A):
	for p in range(n):
		pivot(A,p)
		for i in range(p+1,k):
			d = math.gcd(abs(A[p,p]), abs(A[i,p]))
			a = int(A[p,p]/d)
			b = int(A[i,p]/d)
			A[i,p] = a*A[i,p]-b*A[p,p]

	print(A,"\n")

global k,n
k = 15
n = 10

A = init_matrix(k,n)
print(A,"\n")
#pivot(A)
improving_pivot(A)