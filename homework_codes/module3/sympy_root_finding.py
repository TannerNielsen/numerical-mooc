#!/usr/bin/python

import sympy
sympy.init_printing()

def main():

  u_max, u_star, rho_max, rho_star, A, B, rho = sympy.symbols('u_max u_star rho_max rho_star A B rho')

  eq1 = sympy.Eq( 0, u_max*rho_max*(1 - A*rho_max-B*rho_max**2) )
  eq2 = sympy.Eq( 0, u_max*(1 - 2*A*rho_star-3*B*rho_star**2) )
  eq3 = sympy.Eq( u_star, u_max*(1 - A*rho_star - B*rho_star**2) )
  eq4 = sympy.Eq(eq2.lhs - 3*eq3.lhs, eq2.rhs - 3*eq3.rhs)
  eq4.simplify()
  
  rho_sol = sympy.solve(eq4,rho_star)[0]
  B_sol = sympy.solve(eq1,B)[0]
  quadA = eq2.subs([(rho_star, rho_sol), (B,B_sol)])
  quadA.simplify()

  A_sol = sympy.solve(quadA, A)
  #print A_sol[0]
  #print A_sol[1]

  aval = A_sol[0].evalf(subs={u_star: 1.5, u_max:2.0, rho_max:15.0} )
  bval = B_sol.evalf(subs={rho_max:15.0, A:aval} )
  print aval
  #print A_sol[1].evalf(subs={u_star: 0.7, u_max:1.0, rho_max:10.0} )
  print bval

  eq5 = sympy.Eq( 0, u_max*(1 - 2*aval*rho-3*bval*rho**2) )
  rho_min = sympy.solve(eq5, rho)
  print 'rho_min= ', rho_min



if __name__ == '__main__':
  main()
