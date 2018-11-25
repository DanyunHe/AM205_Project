# -*- coding: utf-8 -*-
"""
Created on Sat Nov 24 20:41:01 2018

@author: jiayi
"""
import numpy as np
import matplotlib.pyplot as plt

#Set up
dx = 1
dt = 0.005
H = 25
L = 150

c = dx/dt
r = 0.6

n = int(H/dx) + 1
m = int(L/dx) + 1
plane = np.zeros((n,m))
plane[0, :] = 1
plane[-1, :] = 1


e = [np.array([0, 0]), np.array([1, 0]), np.array([0, 1]), np.array([-1, 0]), \
     np.array([0, -1]), np.array([1, 1]), np.array([-1, 1]), np.array([-1, -1]), \
     np.array([1, -1])]

def density(f0, f1, f2, f3, f4, f5, f6, f7, f8, s):
    dens = f0[:, :, s] + f1[:, :, s] + f2[:, :, s] + f3[:, :, s] + f4[:, :, s] \
           + f5[:, :, s] + f6[:, :, s] + f7[:, :, s] + f8[:, :, s]
    return dens

def velocity(f0, f1, f2, f3, f4, f5, f6, f7, f8, dens, s):
    vel = np.zeros((n, m, 2))
    for i in range(0, n):
        for j in range(0, m):
            velo = 1/dens[i, j] * c * (f0[i, j, s] * e[0] + f1[i, j, s] * e[1] \
                   + f2[i, j, s] * e[2] + f3[i, j, s] * e[3] + f4[i, j, s] * e[4] \
                   + f5[i, j, s] * e[5] + f6[i, j, s] * e[6] + f7[i, j, s] * e[7] \
                   + f8[i, j, s] * e[8])
            vel[i, j, 0] = velo[0]
            vel[i, j, 1] = velo[1]
    return vel


def feq(dens, vel):
    w = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36]
    f_eq = np.zeros((n, m, 9))
    for i in range(0, 9):
        for j in range(0, n):
            for k in range(0, m):
                s = w[i] * (3 * (e[i] @ vel[j, k, :])/c + \
                    9/2 * (e[i] @ vel[j, k, :])**2/(c**2) - \
                    3/2 * (vel[j, k, :] @ vel[j, k, :])/(c**2))
                f_eq[j, k, i] = w[i] * dens[j, k] + dens[j, k] * s
    return f_eq



'''
'''
def solve_t(t, p0, p1):
    it = int(t/dt)
    if it < t/dt:
        it+=1
        
    f0 = np.zeros((n, m ,3))
    f0[:, :, :] = 1/9
    f1 = np.zeros((n, m ,3))
    f1[:, :, :] = 1/9
    f2 = np.zeros((n, m ,3))
    f2[:, :, :] = 1/9
    f3 = np.zeros((n, m ,3))
    f3[:, :, :] = 1/9
    f4 = np.zeros((n, m ,3))
    f4[:, :, :] = 1/9
    f5 = np.zeros((n, m ,3))
    f5[:, :, :] = 1/9
    f6 = np.zeros((n, m ,3))
    f6[:, :, :] = 1/9
    f7 = np.zeros((n, m ,3))
    f7[:, :, :] = 1/9
    f8 = np.zeros((n, m ,3))
    f8[:, :, :] = 1/9
    
    dens = density(f0, f1, f2, f3, f4, f5, f6, f7, f8, 0)
    vel = velocity(f0, f1, f2, f3, f4, f5, f6, f7, f8, dens, 0)
    f_eq = feq(dens, vel)
    
    
    for i in range(0, it):
        for j in range(0, n):
            for k in range(0, m):
                
                #upper wall BC:
                if j == 0 and k != 0 and k != m-1:
                    f0[j, k, 1] = f0[j, k, 0]
                    f1[j, k, 1] = f1[j, k-1, 0]
                    f2[j, k, 1] = f2[j+1, k, 0]
                    f3[j, k, 1] = f3[j, k+1, 0]
                    f4[j, k, 1] = f2[j, k, 0]  #bounce back
                    f5[j, k, 1] = f5[j+1, k-1, 0]
                    f6[j, k, 1] = f6[j+1, k+1, 0]
                    f7[j, k, 1] = f5[j, k, 0]   #bounce back
                    f8[j, k, 1] = f6[j, k, 0]   #bounce back
                    
                #lower wall BC:
                elif j == n-1 and k != 0 and k != m-1:
                    f0[j, k, 1] = f0[j, k, 0]
                    f1[j, k, 1] = f1[j, k-1, 0]
                    f2[j, k, 1] = f4[j, k, 0]   #bounce back
                    f3[j, k, 1] = f3[j, k+1, 0]
                    f4[j, k, 1] = f4[j-1, k, 0]
                    f5[j, k, 1] = f7[j, k, 0]   #bounce back
                    f6[j, k, 1] = f8[j, k, 0]   #bounce back
                    f7[j, k, 1] = f7[j-1, k+1, 0]
                    f8[j, k, 1] = f8[j-1, k-1, 0]
                
                
                #inlet BC: p0, vy=0
                elif k == 0 and j != 0 and j != n-1:
                    f0[j, k, 1] = f0[j, k, 0]
                    f2[j, k, 1] = f2[j+1, k, 0]
                    f6[j, k, 1] = f6[j+1, k+1, 0]
                    f3[j, k, 1] = f3[j, k+1, 0]
                    f7[j, k, 1] = f7[j-1, k+1, 0]
                    f4[j, k, 1] = f4[j-1, k, 0]
                    vx = c * (1- 1/p0 * (f0[j, k, 1] + f2[j, k, 1] + \
                         f4[j, k, 1] + 2 * (f3[j, k, 1] + f6[j, k, 1] + \
                           f7[j, k, 1])))
                    v = np.array([vx, 0])
                    s1 = 1/9 * (3 * (e[1] @ v) /c + 9/2 * (e[1] @ v)**2 / (c**2)\
                                - 3/2 * (v @ v)/(c**2))
                    s3 = 1/9 * (3 * (e[3] @ v) /c + 9/2 * (e[3] @ v)**2 / (c**2)\
                                - 3/2 * (v @ v)/(c**2))
                    f_eq1 =  1/9 * p0 + p0 * s1
                    f_eq3 =  1/9 * p0 + p0 * s3
                    
                    f1[j, k, 1] = f3[j, k, 1] + f_eq1 - f_eq3
                    f5[j, k, 1] = 1/2 * (p0 * vx/c - f1[j, k, 1] - f2[j, k, 1] + \
                                   f3[j, k, 1] + f4[j, k, 1] + 2*f7[j, k, 1])
                    f8[j, k, 1] = 1/2 * (p0 * vx/c - f1[j, k, 1] + f2[j, k, 1] + \
                                   f3[j, k, 1] - f4[j, k, 1] + 2*f6[j, k, 1])
                
                
                #outlet BC: p1, vy=0
                elif k == m-1 and j != 0 and j != n-1:
                    f0[j, k, 1] = f0[j, k, 0]
                    f2[j, k, 1] = f2[j+1, k, 0]
                    f5[j, k, 1] = f5[j+1, k-1, 0]
                    f1[j, k, 1] = f1[j, k-1, 0]
                    f8[j, k, 1] = f8[j-1, k-1, 0]
                    f4[j, k, 1] = f4[j-1, k, 0]
                    vx = c * (-1 +  1/p1 * (f0[j, k, 1] + f2[j, k, 1] + \
                         f4[j, k, 1] + 2 * (f1[j, k, 1] + f5[j, k, 1] + \
                           f8[j, k, 1])))
                    v = np.array([vx, 0])
                    s1 = 1/9 * (3 * (e[1] @ v) /c + 9/2 * (e[1] @ v)**2 / (c**2)\
                                - 3/2 * (v @ v)/(c**2))
                    s3 = 1/9 * (3 * (e[3] @ v) /c + 9/2 * (e[3] @ v)**2 / (c**2)\
                                - 3/2 * (v @ v)/(c**2))
                    f_eq1 =  1/9 * p1 + p1 * s1
                    f_eq3 =  1/9 * p1 + p1 * s3
                    
                    f3[j, k, 1] = f1[j, k, 1] + f_eq3 - f_eq1
                    f6[j, k, 1] = 1/2 * (-p1*vx/c + f1[j, k, 1] - f2[j, k, 1] - \
                                   f3[j, k, 1] + f4[j, k, 1] + 2*f8[j, k, 1])
                    f7[j, k, 1] = 1/2 * (-p1 * vx/c + f1[j, k, 1] + f2[j, k, 1] - \
                                   f3[j, k, 1] - f4[j, k, 1] + 2*f5[j, k, 1])
                 
                    
                #corner BC
                elif j == 0 and k == 0:
                    f0[j, k, 1] = f0[j, k, 0]
                    f2[j, k, 1] = f2[j+1, k, 0]
                    f6[j, k, 1] = f6[j+1, k+1, 0]
                    f3[j, k, 1] = f3[j, k+1, 0]
                    
                    f1[j, k, 1] = f3[j, k, 1]
                    f8[j, k, 1] = f6[j, k, 1]
                    f4[j, k, 1] = f2[j, k, 1]
                    
                    f5[j, k, 1] = 1/2 * (p0 - f0[j, k, 1] - f1[j, k, 1] - \
                                  f2[j, k, 1] - f3[j, k, 1] - f4[j, k, 1] - \
                                  f6[j, k, 1] - f8[j, k, 1])
                    f7[j, k, 1] = f5[j, k, 1]
                    
                elif j == 0 and k == m-1:
                    f0[j, k, 1] = f0[j, k, 0]
                    f2[j, k, 1] = f2[j+1, k, 0]
                    f5[j, k, 1] = f5[j+1, k-1, 0]
                    f1[j, k, 1] = f1[j, k-1, 0]
                    
                    f3[j, k, 1] = f1[j, k, 1]
                    f7[j, k, 1] = f5[j, k, 1]
                    f4[j, k, 1] = f2[j, k, 1]
                    
                    f6[j, k, 1] = 1/2 * (p1 - f0[j, k, 1] - f1[j, k, 1] - \
                                  f2[j, k, 1] - f3[j, k, 1] - f4[j, k, 1] - \
                                  f5[j, k, 1] - f7[j, k, 1])
                    f8[j, k, 1] = f6[j, k, 1]
                
                elif j == n-1 and k == 0:
                    f0[j, k, 1] = f0[j, k, 0]
                    f3[j, k, 1] = f3[j, k+1, 0]
                    f7[j, k, 1] = f7[j-1, k+1, 0]
                    f4[j, k, 1] = f4[j-1, k, 0]
                    
                    f1[j, k, 1] = f3[j, k, 1]
                    f5[j, k, 1] = f7[j, k, 1]
                    f2[j, k, 1] = f4[j, k, 1]
                    
                    f6[j, k, 1] = 1/2 * (p0 - f0[j, k, 1] - f1[j, k, 1] - \
                                  f2[j, k, 1] - f3[j, k, 1] - f4[j, k, 1] - \
                                  f5[j, k, 1] - f7[j, k, 1])
                    f8[j, k, 1] = f6[j, k, 1]
                    
                    
                elif j == n-1 and k == m-1:
                    f0[j, k, 1] = f0[j, k, 0]
                    f1[j, k, 1] = f1[j, k-1, 0]
                    f8[j, k, 1] = f8[j-1, k-1, 0]
                    f4[j, k, 1] = f4[j-1, k, 0]
                    
                    f3[j, k, 1] = f1[j, k, 1]
                    f6[j, k, 1] = f8[j, k, 1]
                    f2[j, k, 1] = f4[j, k, 1]
                    
                    f5[j, k, 1] = 1/2 * (p1 - f0[j, k, 1] - f1[j, k, 1] - \
                                  f2[j, k, 1] - f3[j, k, 1] - f4[j, k, 1] - \
                                  f6[j, k, 1] - f8[j, k, 1])
                    f7[j, k, 1] = f5[j, k, 1]
                    
                else: 
                    f0[j, k, 1] = f0[j, k, 0]
                    f1[j, k, 1] = f1[j, k-1, 0]
                    f2[j, k, 1] = f2[j+1, k, 0]
                    f3[j, k, 1] = f3[j, k+1, 0]
                    f4[j, k, 1] = f4[j-1, k, 0]
                    f5[j, k, 1] = f5[j+1, k-1, 0]
                    f6[j, k, 1] = f6[j+1, k+1, 0]
                    f7[j, k, 1] = f7[j-1, k+1, 0]
                    f8[j, k, 1] = f8[j-1, k-1, 0]
        dens = density(f0, f1, f2, f3, f4, f5, f6, f7, f8, 1)
        vel = velocity(f0, f1, f2, f3, f4, f5, f6, f7, f8, dens, 1)
        f_eq = feq(dens, vel)
        
        f0[:, :, 2] = f0[:, :, 1] - 1/r * (f0[:, :, 1] - f_eq[:, :, 0])
        f1[:, :, 2] = f1[:, :, 1] - 1/r * (f1[:, :, 1] - f_eq[:, :, 1])
        f2[:, :, 2] = f2[:, :, 1] - 1/r * (f2[:, :, 1] - f_eq[:, :, 2])
        f3[:, :, 2] = f3[:, :, 1] - 1/r * (f3[:, :, 1] - f_eq[:, :, 3])
        f4[:, :, 2] = f4[:, :, 1] - 1/r * (f4[:, :, 1] - f_eq[:, :, 4])
        f5[:, :, 2] = f5[:, :, 1] - 1/r * (f5[:, :, 1] - f_eq[:, :, 5])
        f6[:, :, 2] = f6[:, :, 1] - 1/r * (f6[:, :, 1] - f_eq[:, :, 6])
        f7[:, :, 2] = f7[:, :, 1] - 1/r * (f7[:, :, 1] - f_eq[:, :, 7])
        f8[:, :, 2] = f8[:, :, 1] - 1/r * (f8[:, :, 1] - f_eq[:, :, 8])
        
        f0[:, :, 0] = f0[:, :, 2]
        f1[:, :, 0] = f1[:, :, 2]
        f2[:, :, 0] = f2[:, :, 2]
        f3[:, :, 0] = f3[:, :, 2]
        f4[:, :, 0] = f4[:, :, 2]
        f5[:, :, 0] = f5[:, :, 2]
        f6[:, :, 0] = f6[:, :, 2]
        f7[:, :, 0] = f7[:, :, 2]
        f8[:, :, 0] = f8[:, :, 2]
        
        
    dens = density(f0, f1, f2, f3, f4, f5, f6, f7, f8, 2)
    vel = velocity(f0, f1, f2, f3, f4, f5, f6, f7, f8, dens, 2)
    result = {"density": dens, "velocity": vel}
    return result
  
   '''test'''
        
p0 = 1.0019
p1 = 0.9981
test = solve_t(6, p0, p1)    
     

      
test_dens = test["density"]
test_vel = test["velocity"]
test_vel_abs = np.zeros((n, m))
for j in range(0, n):
    for k in range(0, m):
        test_vel_abs[j, k] = np.sqrt(test_vel[j, k, 0] ** 2)


plt.figure(figsize = (15, 8))
plt.imshow(test_dens)
plt.colorbar()
plt.show()    


plt.figure(figsize = (15, 8))
plt.imshow(test_vel_abs)
plt.colorbar()
plt.show() 


exact = [(p1 - p0)/(2*((2*r-1)/6*c*dx)*test_dens[0, 0]*L)*y*(y-H) for y in range(0, n)]
plt.plot(test_vel[:, 0, 0], np.arange(0, n, 1), "--")
#plt.plot(exact, np.arange(0, n, 1))
 
'''test2'''
test2 = solve_t(1, p0, p1)    
     

      
test2_dens = test2["density"]
test2_vel = test2["velocity"]
test2_vel_abs = np.zeros((n, m))
for j in range(0, n):
    for k in range(0, m):
        test2_vel_abs[j, k] = np.sqrt(test2_vel[j, k, 0] ** 2)


plt.figure(figsize = (15, 8))
plt.imshow(test2_dens)
plt.colorbar()
plt.show()    


plt.figure(figsize = (15, 8))
plt.imshow(test2_vel_abs)
plt.colorbar()
plt.show() 

plt.plot(test2_vel[:, 0, 0], np.arange(0, n, 1), "--")

'''test end'''


    
    
    
    
    
                
                
                
 








