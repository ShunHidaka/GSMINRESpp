import numpy as np
import scipy as sp

def is_positive_definite(matrix):
	if sp.sparse.isspmatrix(matrix):
		matrix = matrix.toarray()
	try:
		np.linalg.cholesky(matrix)
		return True
	except np.linalg.LinAlgError:
		return False

file_path = input("file name: ")
matrix    = sp.io.mmread(file_path)
if is_positive_definite(matrix):
	print("Positive Definite")
else:
	print("Not Positive Definite")
input("press to end")