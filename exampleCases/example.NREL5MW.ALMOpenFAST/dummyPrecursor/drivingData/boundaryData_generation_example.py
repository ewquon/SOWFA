#!/usr/bin/env python
import numpy as np
import datatools.SOWFA.constant.boundaryData as bd

# TODO: setup 1D times array
times = [0,99999.9999999]

# TODO: setup 1D x,y,z arrays
Ly,Lz = 1000., 500.
ds = 10.
x = np.array([0])
y = np.arange(0,Ly+ds,ds)
z = np.arange(0,Lz+ds,ds)
patch = bd.CartesianPatch(x,y,z,dpath='boundaryData',name='west')
patch.write_points()

# TODO: generate U,T arrays with shape==(time, height[, direction])
T = 300.0 * np.ones((len(times),len(z)))
U = np.zeros((len(times),len(z),3))
U[:,:,0] = 8.0  # uniform westerly flow
patch.write_profiles(times,z,U=U,T=T)

