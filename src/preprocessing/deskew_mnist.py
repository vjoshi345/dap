# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 18:10:58 2017

@author: Varun Joshi
Downloaded from: https://fsix.github.io/mnist/Deskewing.html
"""

import numpy as np
from scipy.ndimage import interpolation
from numpy import genfromtxt
import matplotlib.pyplot as plt

def moments(image):
    c0,c1 = np.mgrid[:image.shape[0],:image.shape[1]] # A trick in numPy to create a mesh grid
    totalImage = np.sum(image) #sum of pixels
    m0 = np.sum(c0*image)/totalImage #mu_x
    m1 = np.sum(c1*image)/totalImage #mu_y
    m00 = np.sum((c0-m0)**2*image)/totalImage #var(x)
    m11 = np.sum((c1-m1)**2*image)/totalImage #var(y)
    m01 = np.sum((c0-m0)*(c1-m1)*image)/totalImage #covariance(x,y)
    mu_vector = np.array([m0,m1]) # Notice that these are \mu_x, \mu_y respectively
    covariance_matrix = np.array([[m00,m01],[m01,m11]]) # Do you see a similarity between the covariance matrix
    return mu_vector, covariance_matrix

def deskew(image):
    c,v = moments(image)
    alpha = v[0,1]/v[0,0]
    affine = np.array([[1,0],[alpha,1]])
    ocenter = np.array(image.shape)/2.0
    offset = c-np.dot(affine,ocenter)
    return interpolation.affine_transform(image,affine,offset=offset)


data = genfromtxt('C:\CMU\CMU-Spring-2016\DAP\working-directory\dap\data\mnist2.csv', delimiter=',')
print type(data)
d, n = data.shape[0], data.shape[1]
print d, n

#print np.array_equal(data[:, 150], data[:, 150].reshape((28,28), order='F').flatten(order='F'))

deskewed = np.zeros(shape=(d,n))
for i in range(n):
    deskewed[:, i] = deskew(data[:, i].reshape((28,28), order='F')).flatten(order='F')

#plt.imshow(deskewed[:, 5].reshape((28,28), order='F'), cmap='gray')
np.savetxt(fname='mnist2-deskewed.csv', X=deskewed, delimiter=',', fmt='%6.6f')

for i in range(50):
    plt.subplot(1, 2, 1)
    plt.imshow(data[:, i].reshape((28,28), order='F'), cmap='gray')
    
    plt.subplot(1, 2, 2)
    plt.imshow(deskew(data[:, i].reshape((28,28), order='F')), cmap='gray')
    
    name = 'Deskewed-mnist2\image' + str(i) + '.png'
    plt.savefig(name)













