#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import scipy.stats as stats
import progressbar as pb

__all__ = ["fit_with_yerr", "fit_with_xyerr"]


def _residuals_yerr(params, x, y, dy, func):
   """
   Return the vector of residuals (normalized with the error)
   for func(x, param).
   """

   ris=(y-func(x, params) )/dy
   return ris


def fit_with_yerr(x, y, dy, xmin, xmax, func, params, samples, \
   stop_param=1.0e-15, plot_fit=1, plot_band=1, plot_residuals=1, \
   plot_distribution=1, save_figs=0, show_progressbar=1):
   """
   Perform a fit to data on [xmin, xmax] with the function func(x, param),
   using "samples" bootstrap samples to evaluate the errors.

   stop_param: stopping parameter for the least square regression
   plot_fit: if =1 plot the optimal fit together with the data
   plot_band: if =1 plot the 1std band together with data
   plot_residuals: if =1 plot residuals after convergence
   plot_distribtion: if =1 plot the bootstrapped distributions of the parameters
   save_figs: if =1 save the figures in png insted of displaying them
   show_progressbar: if =1 show the progressbar

   return the optimal vales of the parameters, their errors, 
   the value of chi^2,the number of dof, the p-value
   and the bootstrap samples of the parameters.
   """

   mask = ((x<=xmax) & (x>=xmin))
   x=x[mask]
   y=y[mask]
   dy=dy[mask]

   band_size=1000

   data_length=len(x)
 
   # array to store the bootstrapped results 
   boot_sample=np.empty((len(params), samples), dtype=np.float)

   if plot_band==1:
     x_band=np.linspace(xmin, xmax, band_size)
     boot_band=np.empty((band_size, samples), dtype=np.float)
     
   for i in range(samples): 
     if show_progressbar==1:
       pb.progress_bar(i, samples)

     # bootstrap sample
     booty=y+np.random.normal(0, dy, data_length) 

     # least square regression
     ris = opt.leastsq(_residuals_yerr, params, ftol=stop_param, args=(x, booty, dy, func))
     boot_sample[:,i]=ris[0]

     if plot_band==1:
       boot_band[:,i]=func(x_band, ris[0])

   # optimal parameters and errors
   ris=np.mean(boot_sample, axis=1)
   err=np.std(boot_sample, axis=1, ddof=1)

   # auxilliary stuff
   opt_res=_residuals_yerr(ris, x, y, dy, func)
   chi2=np.sum(opt_res*opt_res)
   dof=data_length - len(params)
   pvalue=1.0 - stats.chi2.cdf(chi2, dof)


   if plot_fit==1:
     x_aux=np.linspace(xmin, xmax, 1000)
     y_aux=func(x_aux, ris)

     plt.figure('Best fit (chi2/dof=%.4f/%d=%f)' % (chi2, dof, chi2/dof))
     plt.xlim(0.9*xmin, 1.1*xmax)
     plt.errorbar(x, y, yerr=dy, fmt='ob', ms=5)
     plt.plot(x_aux,y_aux,'r-')
      
     if plot_band==1:
       band_mean=np.mean(boot_band, axis=1)
       band_std=np.std(boot_band, axis=1)
       plt.plot(x_band, band_mean + band_std,'g-')
       plt.plot(x_band, band_mean - band_std,'g-')

     if save_figs==1:
       plt.savefig('fit.png')
     else:
       plt.show()

   if plot_residuals==1:
     x_aux=np.linspace(xmin, xmax, 1000)
     y_aux=np.ones(len(x_aux))

     plt.figure('Residuals')
     plt.xlim(0.9*xmin, 1.1*xmax)
     plt.errorbar(x, opt_res, yerr=1, fmt='ob', ms=5)
     plt.plot(x_aux, -2*y_aux, 'g:')
     plt.plot(x_aux, -y_aux, 'r--')
     plt.plot(x_aux, 0*y_aux, 'r-')
     plt.plot(x_aux, y_aux, 'r--')
     plt.plot(x_aux, 2*y_aux, 'g:')
     
     if save_figs==1:
       plt.savefig('residuals.png')
     else:
       plt.show()

   if plot_distribution==1:
     for i in range(len(params)):
       plt.figure('Bootstrapped distribution of param[%d]' % i)
       plt.xlabel('param[%d] values' % i)
       plt.ylabel('distribution histogram')
       plt.hist(boot_sample[i], bins='auto')
       
       if save_figs==1:
         plt.savefig('param'+str(i)+'.png')
       else:
         plt.show()

   return ris, err, chi2, dof, pvalue, boot_sample


def _residuals_xyerr(extended_params, x, dx, y, dy, true_param_length, func):
   """
   Return the vector of residuals (normalized with the error)
   for func(x, param).
   
   extended_param[:true_param_length] = real parameters
   extended_param[true_param_length:] = auxilliary parameters like x
   """

   risy=(y-func(extended_params[true_param_length:], extended_params[:true_param_length]) )/dy
   risx=(x-extended_params[true_param_length:])/dx
   return np.append(risy, risx)


def fit_with_xyerr(x, dx, y, dy, xmin, xmax, func, params, samples, \
   stop_param=1.0e-15, plot_fit=1, plot_band=1, plot_residuals=1, \
   plot_distribution=1, save_figs=0, show_progressbar=1):
   """
   Perform a fit to data on [xmin, xmax] with the function func(x, param),
   using "samples" bootstrap samples to evaluate the errors.

   stop_param: stopping parameter for the least square regression
   plot_fit: if =1 plot the optimal fit together with the data
   plot_band: if =1 plot the 1std band together with data
   plot_residuals: if =1 plot residuals after convergence
   plot_distribtion: if =1 plot the bootstrapped distributions of the parameters
   save_figs: if =1 save the figures in png insted of displaying them
   show_progressbar: if =1 show the progress bar

   return the optimal vales of the parameters, their errors, 
   the value of chi^2,the number of dof, the p-value
   and the bootstrap samples of the parameters.
   """

   mask = ((x<=xmax) & (x>=xmin))
   x=x[mask]
   dx=dx[mask]
   y=y[mask]
   dy=dy[mask]

   band_size=1000

   data_length=len(x)

   true_param_length=len(params)

   extended_params=np.append(params, x)
 
   # array to store the bootstrapped results 
   boot_sample=np.empty((len(params), samples), dtype=np.float)

   if plot_band==1:
     x_band=np.linspace(xmin, xmax, band_size)
     boot_band=np.empty((band_size, samples), dtype=np.float)
  
   for i in range(samples):
     if show_progressbar==1: 
       pb.progress_bar(i, samples)

     # bootstrap sample
     bootx=x+np.random.normal(0, dx, data_length) 
     booty=y+np.random.normal(0, dy, data_length) 

     # least square regression
     ris = opt.leastsq(_residuals_xyerr, extended_params, ftol=stop_param, args=(bootx, dx, booty, dy, true_param_length, func))
     boot_sample[:,i]=ris[0][:true_param_length]
     if plot_band==1:
       boot_band[:,i]=func(x_band, ris[0][:true_param_length])

   # optimal parameters and errors
   ris=np.mean(boot_sample, axis=1)
   err=np.std(boot_sample, axis=1, ddof=1)

   # auxilliary stuff
   opt_res=_residuals_yerr(ris, x, y, dy, func)
   chi2=np.sum(opt_res*opt_res)
   dof=data_length - true_param_length
   pvalue=1.0 - stats.chi2.cdf(chi2, dof)

   if plot_fit==1:
     x_aux=np.linspace(xmin, xmax, 1000)
     y_aux=func(x_aux, ris)

     plt.figure('Best fit (chi2/dof=%.4f/%d=%f)' % (chi2, dof, chi2/dof))
     plt.xlim(0.9*xmin, 1.1*xmax)
     plt.errorbar(x, y, xerr=dx, yerr=dy, fmt='ob', ms=5)
     plt.plot(x_aux,y_aux,'r-')

     if plot_band==1:
       band_mean=np.mean(boot_band, axis=1)
       band_std=np.std(boot_band, axis=1)
       plt.plot(x_band, band_mean + band_std,'g-')
       plt.plot(x_band, band_mean - band_std,'g-')

     if save_figs==1:
       plt.savefig('fit.png')
     else:
       plt.show()

   if plot_residuals==1:
     x_aux=np.linspace(xmin, xmax, 1000)
     y_aux=np.ones(len(x_aux))

     plt.figure('Residuals')
     plt.xlim(0.9*xmin, 1.1*xmax)
     plt.errorbar(x, opt_res, xerr=dx, yerr=1, fmt='ob', ms=5)
     plt.plot(x_aux, -2*y_aux, 'g:')
     plt.plot(x_aux, -y_aux, 'r--')
     plt.plot(x_aux, 0*y_aux, 'r-')
     plt.plot(x_aux, y_aux, 'r--')
     plt.plot(x_aux, 2*y_aux, 'g:')

     if save_figs==1:
       plt.savefig('residuals.png')
     else:
       plt.show()

   if plot_distribution==1:
     for i in range(0, len(params)):
       plt.figure('Bootstrapped distribution of param[%d]' % i)
       plt.xlabel('param[%d] values' % i)
       plt.ylabel('distribution histogram')
       plt.hist(boot_sample[i], bins='auto')

       if save_figs==1:
         plt.savefig('param'+str(i)+'.png')
       else:
         plt.show()

   return ris, err, chi2, dof, pvalue, boot_sample


#***************************
# unit testing

if __name__=="__main__":
  
  print("**********************")
  print("UNIT TESTING")
  print()


  numsamples=1000

  length=10
  x=np.linspace(1, 10, length)
  y=2*x*x+np.random.normal(0, 10, length)
  dy=10*np.ones(length)

  print("Points generated with y=2*x*x + gaussian noise")
  print()

  fitparams=np.array([0.1, 1.2], dtype=np.float)
  def parabfit(x, param):
    return param[0]+ x*x*param[1]

  print("Fit of the form param[0]+param[1]*x*x")
  print()

  print("Case with errors only on y")
  print()

  ris, err, chi2, dof, pvalue, boot_sample = fit_with_yerr(x, y, dy, 1, 10, parabfit, fitparams, numsamples)

  print("  chi^2/dof = {:3.3f}/{:d} = {:3.3f}".format(chi2, dof, chi2/dof))
  print("  p-value = %f" % pvalue)
  print()
  print("  a = {: f} +- {:f}".format(ris[0], err[0]) ) 
  print("  b = {: f} +- {:f}".format(ris[1], err[1]) ) 
  print()

  print("Case with errors both on x and y")
  print()

  dx=2.0*np.ones(length)/(length)

  ris, err, chi2, dof, pvalue, boot_sample = fit_with_xyerr(x, dx, y, dy, 1, 10, parabfit, fitparams, numsamples)

  print("  chi^2/dof = {:3.3f}/{:d} = {:3.3f}".format(chi2, dof, chi2/dof))
  print("  p-value = %f" % pvalue)
  print()
  print("  a = {: f} +- {:f}".format(ris[0], err[0]) ) 
  print("  b = {: f} +- {:f}".format(ris[1], err[1]) ) 
  print()
  print("**********************")

