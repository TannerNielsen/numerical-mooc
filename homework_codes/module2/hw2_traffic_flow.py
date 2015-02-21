#!/usr/bin/python

# MAE 6286 Practical Numerical Methods with Python
# Homework 2
#
# Solution to traffic flow equations-
# Velocity: v = vmax(1 - rho/rhomax)
# Traffic Flux: F = v*rho
# Non-linear convection of car density: d_rho/dt + (d_F/d_rho)*(d_rho/dx) = 0 or by the chain rule d_rho/dt + d_F/dx = 0
#
# Tanner Nielsen 2/15/2015

import numpy as np
import matplotlib.pyplot as plt
import pylab
import sys
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

def vel(vmax,rhomax,rho):

  v = vmax*(1.0-rho/rhomax)

  return v

def main():

  #vmax = 80.0  # max velocity (km/hr)
  vmax = 136.0  # max velocity (km/hr)
  L = 11.0  # length of road (km)
  rhomax = 250  # max density of cars (cars/km)
  nx = 51
  dt = 0.001
  dx = L/(nx-1)
  print 'dx: ', dx, 'dt: ', dt
  nt = 101

  x = np.linspace(0,L,nx)
  #rho = np.ones(nx)*10
  rho = np.ones(nx)*20
  rho[10:20] = 50

  # initial velocity at t=0
  print 'min velocity (m/s) at t=0: ', min(vel(vmax,rhomax,rho))*(1000.0/3600.0)

  plt.plot(x, rho, color='#003366', ls='--', lw=3)
  plt.ylim(0,75);

  for n in range(1,nt):  
    rhon = rho.copy() 
    for i in range(1,nx):
        #rho[0] = 10.0
        rho[0] = 20.0 
        rho[i] = rhon[i] - dt/dx*((vmax*rhon[i]*(1.0-rhon[i]/rhomax)) - (vmax*rhon[i-1]*(1.0-rhon[i-1]/rhomax)))
#    # central different in space
#    for i in range(1,nx-1):
#        #rho[0] = 10.0
#        rho[0] = 20.0
#        rho[i] = rhon[i] - dt/(2.0*dx)*((vmax*rhon[i+1]*(1.0-rhon[i+1]/rhomax)) - (vmax*rhon[i-1]*(1.0-rhon[i-1]/rhomax)))
#        rho[nx-1] = rhon[nx-1] - dt/dx*((vmax*rhon[nx-1]*(1.0-rhon[nx-1]/rhomax)) - (vmax*rhon[nx-2]*(1.0-rhon[nx-2]/rhomax)))

    print 'Time: ', dt*n*60.0
    if n*dt*60 == 3.0:
      print 'Avg. velocity (m/s) at t=3 minutes: ', np.average(vel(vmax,rhomax,rho))*(1000.0/3600.0)
      print 'Min. velocity (m/s) at t=3 minutes: ', min(vel(vmax,rhomax,rho))*(1000.0/3600.0)
    if n*dt*60 == 6.0:
      print 'Min. velocity (m/s) at t=6 minutes: ', min(vel(vmax,rhomax,rho))*(1000.0/3600.0) 
    plt.plot(x, rho, color='red', ls='-', lw=3)
    plt.ylim(0,75);
    plt.xlabel(r'Length of Road (km)', fontsize=18)
    plt.ylabel(r'Density (cars/km)', fontsize=18)
    plt.title(r'Traffic Flow', fontsize=18)
    plt.text(6, 60, 'Time: '+str(dt*n*60.0)+' min.', fontsize=15)
    pylab.savefig('./plots/forwardTime_backwardSpace_'+str(n)+'.png')
    #pylab.savefig('./plots/forwardTime_centralSpace_'+str(n)+'.png')
    plt.close();

    plt.figure(2)
    plt.plot(x, vel(vmax,rhomax,rho)*(1000.0/3600.0), color='red', ls='-', lw=3)
    plt.ylim(0,75);
    plt.xlabel(r'Length of Road (km)', fontsize=18)
    plt.ylabel(r'Velocity (m/s)', fontsize=18)
    plt.title(r'Traffic Flow', fontsize=18)
    plt.text(6, 60, 'Time: '+str(dt*n*60.0)+' min.', fontsize=15)
    pylab.savefig('./plots/velocity/forwardTime_backwardSpace_'+str(n)+'.png')
    #pylab.savefig('./plots/forwardTime_centralSpace_'+str(n)+'.png')
    plt.close(2);

  print 'Done.'


if __name__ == '__main__':
  main()
