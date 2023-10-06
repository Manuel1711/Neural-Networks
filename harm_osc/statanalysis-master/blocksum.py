#!/usr/bin/env python3

import numpy as np

__all__ = ["blocksum"]


def blocksum(vec_in, block):
  """Block sums of a vector.

  Given a numpy vector "vec_in", it return the vector of block sums using
  blocksize "block" for blocking. 
  The returned vector has length int(len(vec_in)/block)
  """

  end =  block * int(len(vec_in)/block)
  ris=np.sum(np.reshape(vec_in[:end],(-1, block)), axis=1)
  return ris



#***************************
# unit testing

if __name__=="__main__":
  
  print("**********************")
  print("UNIT TESTING")
  print()

  vec=np.array([1, 2, 3, 4, 5], dtype=np.float)

  print("Original vector (vec)")
  print(vec)
  print()

  print("blocksum(vec, 1)")
  print(blocksum(vec, 1))
  print()

  print("blocksum(vec, 2)")
  print(blocksum(vec, 2))
  print()

  print("blocksum(vec, 3)")
  print(blocksum(vec, 3))
  print("**********************")

