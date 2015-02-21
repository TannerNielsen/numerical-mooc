#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import pylab


m_s = 50.0  # [kg] weight of the rocket shell
g = 9.81  # [m/s^2] gravitational constant
rho = 1.091  # [kg/m^3] average air density
r = 0.5  # [m] radius of the rocket
v_e = 325.0  # [m/s] exhaust speed
C_D = 0.15  # drag coefficient
m_p0 = 100.0 # initial weight of rocket propellant
A = np.pi*r**2  # calculate the cross sectional area of the rocket
h = 0.0  # [m] initial height
v = 0.0  # [m/s] initial velocity

T = 40.0

def m_p(t):
  # returns the propellant mass as a function of time
  if (t <= 5.0):
    m_p = m_p0 - 20.0*t
  else:
    m_p = 0.0
  return m_p

def m_p_dot(t):
  if (t <= 5.0):
    m_p_dot = 20.0
  else:
    m_p_dot = 0.0

  return m_p_dot


def f(u,t):

  h = u[0]
  v = u[1]
  return np.array([v,-g + m_p_dot(t)*v_e/(m_s+m_p(t)) - 0.5*rho*v*abs(v)*A*C_D/(m_s+m_p(t))])


def euler_step(u,t,dt):

  return u + dt*f(u,t)
  

def main():

#  t = 0.0  # time

  dt = 0.1 # timestep
  N = int(T/dt)

  t = np.linspace(0.0,T,N)

  u = np.empty((N,2))
  u[0] = np.array([h,v])

  test = []

  # initialize variables
  max_v = 0.0 # max velocity
  t_max_v = 0.0  # time at max velocity
  h_max_v = 0.0  # height at max velocity
  max_h = 0.0  # max height
  t_max_h = 0.0  # time at max height
#  t_impact   # time at impact to ground
#  v_impact   # velocity at impact



  for n in range(N-1):
    u[n+1] = euler_step(u[n],t[n],dt)
    if u[n+1,1] > max_v:
      max_v = u[n+1,1]
      t_max_v = t[n+1]
      h_max_v = u[n+1,0]
    if u[n+1,0] > max_h:
      max_h = u[n+1,0]
      t_max_h = t[n+1]
    if u[n+1,0] <= 0.0:
      t_impact = t[n+1]
      v_impact = u[n+1,1]
    #print u

  for i in range(N):
    print u[i,0],u[i,1],t[i]
    
  print 'max velocity = ', max_v, max(u[:,1])
  print 't at max velocity = ', t_max_v
  print 'height at max velocity = ', h_max_v
  print ' '
  print 'max height = ', max(u[:,0]), max_h
  print 't at max height = ', t_max_h
  print ' '
  print 'time of impact = ', t_impact
  print 'velocity at impact = ', v_impact


  plt.plot(t,u[:,0])
  plt.plot(t_max_h,max_h,'ro-')
  plt.annotate('max height = %.2f' % max_h, xy=(t_max_h+0.5, max_h), xytext=(t_max_h+10, max_h-5),
            arrowprops=dict(facecolor='black', shrink=0.1),
            )
  plt.xlabel(r'time (sec)', fontsize=18)
  plt.ylabel(r'height (m)', fontsize=18)
  plt.title('Rocket Height', fontsize=18)
  #pylab.show()

  plt.figure(2)
  plt.plot(t,u[:,1])
  plt.plot(t_max_v,max_v, 'ro-')
  plt.xlabel(r'time (sec)', fontsize=18)
  plt.ylabel(r'velocity (m/s)', fontsize=18)
  plt.title('Rocket Velocity', fontsize=18)
  pylab.show()

if __name__ == '__main__':
  main()
