##
# \file converter.py
# \brief Python script to convert Matrix Market (.mtx) files into a custom CSR text format.
# \author Shuntaro Hidaka
#
# This script reads a sparse matrix in Matrix Market format (.mtx) and outputs
# a plain-text CSR (Compressed Sparse Row) format suitable for use with the
# GSMINRES++ example programs (e.g., \ref sample2.cpp).
#
# The output format contains:
# - Matrix dimension and number of non-zero entries
# - Row pointer array (size N+1)
# - Column indices
# - Complex values
#
# This simplified format allows fast loading of sparse matrices in educational
# and testing environments without relying on external libraries.
#
# \par Usage:
# \code
#  $ python converter.py A.mtx A.csr
# \endcode
#

import sys
import numpy as np
import scipy as sp

# OPEN INPUT FILE
if len(sys.argv) > 1:
        infile = sys.argv[1]
else:
        infile = input("input MATRIX file name: ")
info   = sp.io.mminfo(infile)
matrix = sp.io.mmread(infile).tocsr()

# OPEN OUTPUT FILE
if len(sys.argv) > 2:
        outfile = sys.argv[2]
else:
        outfile = input("output MATRIX file name: ")
fp = open(outfile, 'w', encoding="utf-8")

# WRITE MATRIX INFORMATION
print(f"# {infile}", file=fp)
print(f"# {info}", file=fp)
row_ptr_size = matrix.indptr.size
col_ind_size = matrix.indices.size
data_size    = matrix.data.size
print(f"{row_ptr_size} {col_ind_size} {data_size}", file=fp)

# WRITE MATRIX DATA IN CSR FORMAT
i, j, k = 0, 0, 0
while True:
        ptr, ind, data = -1, -1, -1
        if i < row_ptr_size:
            ptr = matrix.indptr[i]
        if j < col_ind_size:
            ind = matrix.indices[j]
        if k < data_size:
            data = matrix.data[j]
        if ind == ptr == data == -1:
            break
        if info[4] == 'real':
            print("{:d} {:d} {:.20f} 0.0".format(ptr, ind, data.real), file=fp)
        if info[4] == 'complex':
            print("{:d} {:d} {:.20f} {:.20f}".format(ptr, ind, data.real, data.imag), file=fp)
        i,j,k = i+1,j+1,k+1

# CLOSE FILE
fp.close()                                                                                        
