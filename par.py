"""
Par.py

Written by Jacob N. McLane

Simple code for grabing projected seperation and position
angle from a system.
"""

from math import *
import numpy as np
xp=632  #Target Position in x
yp=510  #Target Position in y
xcen=511.5    #center of motion for body in x
ycen=511.5     #center of motion for body in y
scale=0.009942 #plate scale

X=xp-xcen  #Distace in x between body and center of mass
Y=yp-ycen  #Distace in y between body and center of mass

r=np.hypot(X,Y) #Distace between body and center of mass

PS=r*scale #Projected Seperation

#Calculate position angle
    
if X <= 0:
        PA=abs(degrees(atan2(X,Y)))
  
else :
        PA =360-degrees(atan2(X,Y))

#Report Results
print('Position Angle is '+str(PA)+' degrees.')
print('Projected Seperation is '+str(PS)+' as.')
