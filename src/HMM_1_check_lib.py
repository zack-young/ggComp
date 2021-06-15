import sys
import os

try:
    from multiprocessing import Pool, Process, Manager
except:
    sys.stderr.write('Error in loading python lib "multiprocessing". Please check the installation of python.\n')
    exit(1)

try:
    import numpy
except:
    sys.stderr.write('Error in loading python lib "numpy". Please try "pip install numpy".\n')
    exit(1)

try:
    from hmmlearn import hmm
except:
    sys.stderr.write('Error in loading python lib "hmmlearn". Please try "pip install hmmlearn".\n')
    exit(1)

try:
    from sklearn.preprocessing import normalize
except:
    sys.stderr.write('Error in loading python lib "sklearn". Please try "pip install sklearn".\n')
    exit(1)

