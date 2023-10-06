#!/usr/bin/env python3

import numpy as np
import blocksum as bs
import progressbar as pb

__all__ = ["bootstrap_for_primary", "bootstrap_for_secondary"]

#***************************
#library functions

def bootstrap_for_primary(func, vec_in, block, samples):
  """Bootstrap for primary observables.

  Given a numpy vector "vec_in", compute 
  <func(vec_in)> 
  using blocksize "block" for blocking
  and "samples" resamplings.
  """

  if not isinstance(block, int):
    print("ERROR: blocksize has to be an integer!")
    sys.exit(1)

  if block<1:
    print("ERROR: blocksize has to be positive!")
    sys.exit(1)

  numblocks=int(len(vec_in)/block)
  end =  block * numblocks

  # cut vec_in to have a number of columns multiple of "block" and apply "func" 
  data=func(vec_in[:end])  

  block_sum_data=bs.blocksum(data, block)/float(block)

  # generate bootstrap samples
  aux=len(block_sum_data)
  bootsample=np.random.choice(block_sum_data,(samples, aux) )
  
  # sum up the samples
  risboot=np.sum(bootsample, axis=1)/len(block_sum_data)

  ris=np.mean(risboot)
  err=np.std(risboot, ddof=1)
 
  return ris, err


def bootstrap_for_secondary(func2, block, samples, show_progressbar, *args):
  """Bootstrap for secondary observables.
  
  Every element of *arg is a list of two element of the form
  args[i]=[func_i, vec_i]
  and the final result is 
  func2(<func_0(vec_0)>, ..,<func_n(vec_n)>) 
  with blocksize "block" for blocking
  and "samples" resampling
  show_progressbar: if =1 show the progressbar
  """

  if not isinstance(block, int):
    print("ERROR: blocksize has to be an integer!")
    sys.exit(1)

  if block<1:
    print("ERROR: blocksize has to be positive!")
    sys.exit(1)


  secondary_samples=np.empty(samples, dtype=np.float)

  for sample in range(samples):
    if show_progressbar==1:
      pb.progress_bar(sample, samples)

    primary_samples=[]

    numblocks=int(len(args[0][1])/block)
    end =  block * numblocks

    resampling = np.random.randint(0,numblocks,size=numblocks) 

    for arg in args:
      func_l, vec_l = arg

      # cut vec_in to have a number of columns multiple of "block" and apply "func" 
      data=func_l(vec_l[:end])  

      #block
      block_sum_data=bs.blocksum(data, block)/float(block)

      #sample average
      tmp = np.average([block_sum_data[i] for i in resampling])  

      primary_samples.append(tmp)

    aux=func2(primary_samples)
    secondary_samples[sample]=aux

  ris=np.mean(secondary_samples)
  err=np.std(secondary_samples, ddof=1)

  return ris, err





#***************************
# unit testing

if __name__=="__main__":
  
  print("**********************")
  print("UNIT TESTING")
  print()

  def id(x):
    return x

  def square(x):
    return x*x

  def susc(x):
    return x[0]-x[1]*x[1]

  size=5000
  samples=500

  # gaussian independent data 
  mu=1.0
  sigma=0.2
  test_noauto=np.random.normal(mu, sigma, size)

  # NO AUTOCORRELATION

  # test for primary
  print("Test for primary observables without autocorrelation")
  print("result must be compatible with %f" % mu)

  ris, err = bootstrap_for_primary(id, test_noauto, 1, samples)

  print("average = %f" % ris)
  print("    err = %f" % err)
  print()

  # test for secondary
  print("Test for secondary observables without autocorrelation")
  print("result must be compatible with %f" % (sigma*sigma))

  list0=[square, test_noauto]
  list1=[id, test_noauto]
  ris, err = bootstrap_for_secondary(susc, 1, samples, 1, list0, list1)

  print("average = %f" % ris)
  print("    err = %f" % err)
  print()

  print("**********************")

