#!/usr/bin/env python3

import glob
import os.path
import sys

__all__=["concatenate"]


def concatenate(indir, therm, use_line0_as_header=1):
  """
  Concatenate all the files of the "indir" directory whose names differs only
  for what is written after the last occurrence of the character "_",
  removing the first "therm" lines.

  If use_line0_as_header the first line of each group of file 
  is copied irrespective of the therm value
  """

  complete_list=glob.glob("%s/*.dat" % indir)

  # cut the final part of the name
  tmp_list=[]
  for element in complete_list:
    final=element.rfind("_")
    tmp_list.append(element[:final])

  # keep only the different entries of tmp_list
  unique_list=list(set(tmp_list))

  # concatenate the files
  for element in unique_list:

    outname=os.path.basename(element)+".dat"
    with open(outname, 'a') as outf:
      list_to_merge=glob.glob("%s_*" % element )

      if use_line0_as_header==1:
        # for the first occurrence print also the first line
        # that is used as header
        i=0
        with open(list_to_merge[0], 'r') as inf:
          for line in inf:
            if i>therm or i==0:  
              outf.write(line)
            i=i+1
        inf.close()

        for k in range(1, len(list_to_merge) ):
          i=0
          with open(list_to_merge[k], 'r') as inf:
            for line in inf:
              if i>therm:
                outf.write(line)
              i=i+1
          inf.close() 
      else:
        for k in range(0, len(list_to_merge) ):
          i=0
          with open(list_to_merge[k], 'r') as inf:
            for line in inf:
              if i>therm:
                outf.write(line)
              i=i+1
          inf.close() 

    outf.close()

#***************************
if __name__=="__main__":

  try:
    indir =sys.argv[1]
  except:
    print("USE: %s dir_name thermalization\n" % sys.argv[0])
    sys.exit(1)

  try:
    therm =int(sys.argv[2])
  except:
    print("USE: %s dir_name thermalization\n" % sys.argv[0])
    sys.exit(1)

  if not os.path.isdir(indir):
    print("ERROR: file %s does not exists\n" % infile)
    sys.exit(1)

  if therm<0:
    print("ERROR: 'therm' must be a nonnegative integer\n")
    sys.exit(1)

  concatenate(indir, therm)
