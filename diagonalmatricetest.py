# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 13:58:25 2021

@author: janousu
"""
import numpy as np

# Penta-diagonal matrix setup

w0 = np.random.randint(10, size=(16,12))
w = np.ravel(w0)
e = np.random.randint(10, size=(16,12))
e = np.ravel(e)
no = np.random.randint(10, size=(16,12))
no = np.ravel(no)
s = np.random.randint(10, size=(16,12))
s = np.ravel(s)

cols = w0.shape[1]
rows = w0.shape[0]
a = np.zeros((cols*rows, cols*rows))

implic = 1
alfa = 2


def diagmatrix(a=a, cols=cols, rows=rows, implic=implic, w=w, e=e, no=no, s=s, alfa=alfa):
    i,j = np.indices(a.shape)
    n = cols*rows
    a[i==j] = (w + e + no + s)   # Diagonal
    a[i==j+1] = w[1:]  # West element
    a[i==j-1] = e[:-1]  # East element
    a[i==j+cols] = no[cols:]  # North element
    a[i==j-cols] = s[:n-cols]  # South element
    return a
    

test = diagmatrix(a=a, cols=cols,implic=implic, w=w, e=e, no=no, s=s, alfa=alfa)


Htmp1 = np.linalg.multi_dot([np.linalg.inv(a),w])


inv_a = np.linalg.inv(a)


from numpy.linalg import multi_dot

# Prepare some data
A = np.random.randint(10, size=(16,12))
B = np.random.randint(10, size=(16,12))
C = np.random.randint(10, size=(16,12))
D = np.random.randint(10, size=(16,12))
# the actual dot multiplication
multi_dot([A, B, C, D])


a = np.array([[1,2],[3,4]]) 
b = np.array([[11,12],[13,14]]) 
c = np.array([[21,22],[23,24]]) 
temp = np.dot(a,b)
np.dot(temp,c)

multi_dot([a,b,c])
