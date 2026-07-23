# We are computing the matrix W_mm which gives the perturbations 
# in omega as part of the splitting matrix. What we would like to do 
# is compare them against the values for the b*m term of Dahlen and 
# Sailor 1979. 
import numpy as np 
import matplotlib.pyplot as plt 

n = "0"
l = "2"

matrixdir = "./SemiAnalyicalMatrices/"

# Load the matrix: 
# This should be dimensionalised, in units of angular freq.
W = np.diag(np.loadtxt(f"{matrixdir}/Whole_semi_{n}S{l}_{n}S{l}.txt"))

# Since the spacing is linear, we can get beta as 
beta = np.mean(np.abs(W[1:] - W[:-1]))

# The table in DT98 says to multiply by 1e3 for comparison
beta *= 1e3

# Now we want to compare to the b from DS79, but first we need to 
# rescale the W 
print(W)
print(beta)

print("difference: ", beta/14.905, 14.905/beta)
