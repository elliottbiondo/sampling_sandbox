import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from numpy.random import random

fig = plt.figure()
ax = plt.axes(projection='3d')

#z = np.linspace(0, 1, 100)
#x = z * np.sin(20 * z)
#y = z * np.cos(20 * z)

a=[]
with open("unstr.out", 'r') as f:
   for line in f.readlines():
       a.append([float(x) for x in line.split()])


p = np.array(a)
ax.scatter(p[:,0],p[:,1],p[:,2], c=p[:,4])

plt.draw()

plt.show()
