#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp
import progressbar as pb

__all__ = ["integrator_yerr"]


def integrator_yerr(x, y, dy, xmin, xmax, samples, spline_order, \
   plot_spline=1, plot_band=1, plot_distribution=1, save_figs=0, show_progressbar=1):
   """
   Compute the integral of the data in the range [xmin, xmax] using a spline
   interpolation of order "spline_order" and using "samples" bootstrap samples to
   evaluate the errors.

   plot_spline: if =1 plot the optimal fit together with the data
   plot_band: if =1 plot the 1std band together with data
   plot_distribtion: if =1 plot the bootstrapped distributions of the parameters
   save_figs: if =1 save the figures in png insted of displaying them
   show_progressbar: if =1 show the progressbar

   return the estimated value of the integral, its error and 
   the bootstrap samples of the integral.
   """

   mask = ((x<=xmax) & (x>=xmin))
   x=x[mask]
   y=y[mask]
   dy=dy[mask]

   band_size=1000

   data_length=len(x)
 
   # array to store the bootstrapped results 
   boot_sample=np.empty(samples, dtype=np.float)

   if plot_band==1:
     x_band=np.linspace(xmin, xmax, band_size)
     boot_band=np.empty((band_size, samples), dtype=np.float)
     
   for i in range(samples): 
     if show_progressbar==1:
       pb.progress_bar(i, samples)

     # bootstrap sample
     booty=y+np.random.normal(0, dy, data_length) 

     ris=1
     err=1

     # spline interpolation
     s = interp.UnivariateSpline(x, booty, k=spline_order, s=0)

     # integrate the spline interpolation
     boot_sample[i]=s.integral(xmin, xmax)

     if plot_band==1:
       boot_band[:,i]=s(x_band)

   # optimal parameters and errors
   ris=np.mean(boot_sample)
   err=np.std(boot_sample, ddof=1)

   if plot_spline==1:
     x_aux=np.linspace(xmin, xmax, 1000)
     s = interp.UnivariateSpline(x, y, k=spline_order, s=0)
     y_aux=s(x_aux)

     plt.figure('Plot of the spline interpolation')
     plt.xlim(0.9*xmin, 1.1*xmax)
     plt.errorbar(x, y, yerr=dy, fmt='ob', ms=5)
     plt.plot(x_aux,y_aux,'r-')
      
     if plot_band==1:
       band_mean=np.mean(boot_band, axis=1)
       band_std=np.std(boot_band, axis=1)
       plt.plot(x_band, band_mean + band_std,'g-')
       plt.plot(x_band, band_mean - band_std,'g-')

     if save_figs==1:
       plt.savefig('spline.png')
     else:
       plt.show()

   if plot_distribution==1:
     plt.figure('Bootstrapped distribution of the integral')
     plt.xlabel('value of the integral')
     plt.ylabel('distribution histogram')
     plt.hist(boot_sample, bins='auto')
     
     if save_figs==1:
       plt.savefig('boot_integral.png')
     else:
       plt.show()

   return ris, err, boot_sample



#***************************
# unit testing

if __name__=="__main__":
  
  print("**********************")
  print("UNIT TESTING")
  print()


  numsamples=1000

  length=10
  x=np.linspace(1, 10, length)
  y=2*x*x+np.random.normal(0, 2, length)
  dy=2*np.ones(length)

  print("Points generated with y=2*x*x + gaussian noise")
  print("The result should be comatible with 666")
  print()

  ris, err, boot_sample = integrator_yerr(x, y, dy, 1, 10, 10000, 2)

  print("  ris = {: f} +- {:f}".format(ris, err) ) 
  print("**********************")

