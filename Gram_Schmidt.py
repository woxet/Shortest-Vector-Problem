#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 10 2022

@author: Jordan
"""

from fpylll import IntegerMatrix, FPLLL, GSO
import numpy as np

from time import time
FPLLL.set_random_seed(time())

def init_matrix(n):
	B = IntegerMatrix(n, n)
	B.randomize("uniform", bits=8)
	#print("Basis\n",B, end="\n\n")
	return B

def proj(u,v):
	p = IntegerMatrix(1,n)
	coeff = np.dot(u,v)/(v.norm()**2)
	#print(coeff)
	p[0].addmul(v, x = int(coeff))
	return p[0]

def gram_schmidt(B):
	Bgs = IntegerMatrix(n,n)
	Bgs[0].addmul(B[0])
	for i in range(n):
		sum_proj = IntegerMatrix(1,n)
		for j in range(1,i):
			sum_proj[0].addmul(proj(B[i], Bgs[j]))
		Bgs[i].addmul(B[i])
		Bgs[i].addmul(sum_proj[0], x=-1)
	return Bgs

global n
n = 80

B = init_matrix(n)
GS = gram_schmidt(B)

#print(B, end="\n\n")
#print(GS, end="\n\n")
