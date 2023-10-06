#!/usr/bin/env python3

import numpy as np
import scipy.optimize as opt

__all__ = ["make_integer"]


def _residualsfrominteger(alpha, vec_in):
  """
  Return vec_in*alpha-round(vec_in*alpha)
  """
  
  aux1=vec_in*alpha

  aux2=np.empty(len(aux1), dtype=np.float)
  np.rint(aux1, aux2)

  return aux2-aux1


def make_integer(vec_in):
  """
  Given a numpy vector "vec_in", compute alpha such that all components of 
  alpha*vec_in are on average integer and return 
  alpha, round(alpha*vec)
  """

  #to have a first guess
 
  length=100
  x=np.linspace(0.2, 2.0, length) 
  y=np.empty(length)

  for i in range(length):
    y[i]=np.linalg.norm(_residualsfrominteger(x[i], vec_in))

  #to get the index corresponding to min(y)
  mask = (y<min(y)*(1.0+1.0e-15))
  aux=x[mask]    
  param=aux[0]

  ris = opt.leastsq( _residualsfrominteger, param, args=(vec_in))
  return ris[0][0], np.rint(ris[0][0]*vec_in)


#***************************
# unit testing

if __name__=="__main__":
 
  print("**********************")
  print("UNIT TESTING")
  print()

  scale_in=np.random.uniform(0.8, 1.2)
  print("Result has to be "+str(scale_in))
  print("")

  length=100
  vec0=np.random.random_integers(-6, 6, length)

  vec1=vec0/scale_in
 
  scale_out, vec_int = make_integer(vec1)

  print("Resut="+str(scale_out))
