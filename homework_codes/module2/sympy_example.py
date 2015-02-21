#!/usr/bin/python

import numpy
import sympy
#%matplotlib inline
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

from sympy import init_printing
init_printing()

from sympy.utilities.lambdify import lambdify

def main():

#  x = sympy.symbols('x')

#  print x

  

  x, nu, t = sympy.symbols('x nu t')
  phi = sympy.exp(-(x-4*t)**2/(4*nu*(t+1))) + \
  sympy.exp(-(x-4*t-2*numpy.pi)**2/(4*nu*(t+1)))
  print phi 
  print 

  
  phiprime = phi.diff(x)
  print phiprime
  print

  u = -2*nu*(phiprime/phi)+4
  print(u)

  ufunc = lambdify((t, x, nu), u)
  print("The value of u at t=1, x=4, nu=3 is {}.".format(ufunc(1,4,3)))

  print
  y = sympy.symbols('y')
  z = ((sympy.cos(y)**2)*(sympy.sin(y)**3))/((4*y**5)*sympy.exp(y))
  print(z)
  print

  zprime = z.diff(y)
  zfunc = lambdify((y), zprime)
  print("The value of zprime at y=2.2 is {}.".format(zfunc(2.2)))


if __name__ == '__main__':
  main()
