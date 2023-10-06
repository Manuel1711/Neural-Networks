#!/usr/bin/env python3

import time

__all__ = ["progress_bar"]


def progress_bar(iteration, total, barLength=50, disappear=0):
    """
    Show the progress bar of a loop.

    If disappear=0 the progress bar disappear at the ened of the loop.
    """
    percent = int(round((iteration / total) * 100))
    nb_bar_fill = int(round((barLength * percent) / 100))
    bar_fill = '#' * nb_bar_fill
    bar_empty = ' ' * (barLength - nb_bar_fill)
    print("[{0}] {1}%".format(str(bar_fill + bar_empty), percent), end="\r")
    if disappear==0 and iteration == total-1:
      bar_empty = ' ' * (barLength+10)
      print("{0}".format(str(bar_empty)), end="\r")



#***************************
# unit testing

if __name__=="__main__":
  
  print("**********************")
  print("UNIT TESTING")
  print()

  maxval=50

  for i in range(maxval):
    progress_bar(i, maxval)
    time.sleep(0.1)
  print("process completed")

