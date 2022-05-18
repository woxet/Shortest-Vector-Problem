#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 12 2022

@author: Jordan
"""

from fpylll import IntegerMatrix, FPLLL, GSO
import numpy as np
import random

from time import time
FPLLL.set_random_seed(time())

def vector_sum(u,v, coeff):
	w = IntegerMatrix(1,n)
	w[0].addmul(u)
	for i in range(n):
		w[0,i] = int(v[i]*coeff)
	#print("w :",w)
	return w[0]

def init_matrix(n):
	B = IntegerMatrix(n, n)
	for i in range(n):
		l = IntegerMatrix(1,n)
		for j in range(n):
			r = random.randint(0,5)
			l[0,j] = r
		B[i].addmul(l[0])
	#print("Basis\n",B, end="\n\n")
	return B

def gram_schmidt(B):
	Bgs = IntegerMatrix(n,n)
	Bgs[0].addmul(B[0])

	for i in range(1,n):
		prod_norm2 = 1
		u = IntegerMatrix(1,n)
		for k in range(i):
			prod_norm2 = prod_norm2*B[i].norm()**2
			u[0].addmul(B[k], x=int(prod_norm2))
		Bgs[i].addmul(u[0])

		v = IntegerMatrix(1,n)
		for j in range(1,i):
			prod_norm2 = 1
			for k in range(1,j):
				prod_norm2 = prod_norm2*Bgs[k].norm()**2
			coeff = prod_norm2*np.dot(B[i],Bgs[j])
			print(prod_norm2,np.dot(B[i],Bgs[j]))
			v[0].addmul(Bgs[j], x=int(coeff))
		Bgs[i].addmul(v[0],x=-1)

	return Bgs

global n
n = 10

B = init_matrix(n)
print(B)
print(gram_schmidt(B))