import numpy as np
import os
import sys

"""
This script parses a Gaussian09 output file and .rwf file. 
It saves the number of basis function, number of electrons, Fock matrix/Kohn-Sham matrix and overlap matrix in a new folder. 
"""

def lower_triangle2full(elements, nbf):
  matrix = np.zeros((nbf, nbf))
  count1 = 0 # index of row
  count2 = 0 # maximum indices in each column
  for n, element in enumerate(elements):
#    print count1, count2, float(element)
    matrix[count1, count2] = float(element)
    count1+=1
    if count1==count2+1:
      count1=0
      count2+=1
  matrix= matrix+matrix.T - np.diag(np.diag(matrix))
  return matrix

# Determine logfile 
if len(sys.argv) == 1:
    print "Please provide Gaussian09 output file name as input. Example:\n python parse_g09.py mol.log"
    sys.exit(1)

g09file = sys.argv[1]

fileName, fileExtension = os.path.splitext(g09file)


if os.path.exists(fileName):
  print "Directory exists. Potentially overwriting..."
else:
  os.makedirs(fileName)

np.set_printoptions(precision=10, suppress=True) 

# Extract the number of basis functions from G09 .log file
logfile = open(g09file,'r')
data = logfile.read()
first =data.find('There are')
last = first+30
nbf = int(data[first:last].split()[2])
print "Number of basis functions", nbf

first =data.find('alpha electrons')-5
last = first+5
nelec = int(data[first:last].split()[0])*2
print "Number of electrons", nelec

np.savetxt(fileName + '/number_basis_functions.dat',np.array([nbf]),fmt='%d')
np.savetxt(fileName + '/number_electrons.dat',np.array([nelec]),fmt='%d')

#### Run rwfdump ####
import subprocess
rwfdumps = ["514R","536R"] #, "536R"
rwfname = ['overlap','fock_matrix'] #'eig_vectors',
for n, rwf in enumerate(rwfdumps):
  print 'rwfdump %s %s/%s %s'%(fileName+'.rwf', fileName, rwfname[n], rwfdumps[n])
#  quit()
  subprocess.call('rwfdump %s %s/%s.txt %s'%(fileName+'.rwf', fileName, rwfname[n], rwfdumps[n]), shell=True)
  f= open('%s/%s.txt'%(fileName, rwfname[n]),'r')
  data = f.read()
  matrix = np.zeros((nbf, nbf))
  first =len(data)-data[::-1].find('(read left to right):'[::-1])
  matrix_elements =data[first:].replace('D','E').split()
  matrix = lower_triangle2full(matrix_elements, nbf)
  print matrix
  np.save('%s/%s.npy'%(fileName, rwfname[n]), matrix)
  np.savetxt('%s/%s.dat'%(fileName, rwfname[n]), matrix)
  f.close()

# dump eigenvalues
subprocess.call('rwfdump %s %s/%s %s'%(fileName+'.rwf', fileName, 'eig_values.txt', "522R" ), shell=True)
f= open('%s/eig_values.txt'%(fileName),'r')
data = f.read()
first =len(data)-data[::-1].find('(read left to right):'[::-1])
data =data[first:]
data = data.replace('D','E')
eigs = data.split()
eig_values = np.zeros((nbf))
for n in range(nbf):
  eig_values[n] = float(eigs[n])
np.save('%s/%s.npy'%(fileName,'eig_values'), eig_values)
np.savetxt('%s/%s.dat'%(fileName,'eig_values'), eig_values)

