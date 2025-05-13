# Download matrix archives from the ELSES matrix library
wget http://www.damp.tottori-u.ac.jp/~hoshi/elses_matrix/ELSES_MATRIX_PPE3594_20160426.tgz
wget http://www.damp.tottori-u.ac.jp/~hoshi/elses_matrix/ELSES_MATRIX_VCNT900h_20130501.tgz
wget http://www.damp.tottori-u.ac.jp/~hoshi/elses_matrix/ELSES_MATRIX_VCNT10800h_20130501.tgz
# Extract each archive
tar -xzvf ./ELSES_MATRIX_PPE3594_20160426.tgz
tar -xzvf ./ELSES_MATRIX_VCNT900h_20130501.tgz
tar -xzvf ./ELSES_MATRIX_VCNT10800h_20130501.tgz
# Convert Matrix Market format to CSR format using the provided Python script
python converter.py ./ELSES_MATRIX_PPE3594_20160426/ELSES_MATRIX_PPE3594_20160426_A.mtx PPE3594_A.csr
python converter.py ./ELSES_MATRIX_PPE3594_20160426/ELSES_MATRIX_PPE3594_20160426_B.mtx PPE3594_B.csr
python converter.py ./ELSES_MATRIX_VCNT900h_20130501/ELSES_MATRIX_VCNT900h_A.mtx VCNT900h_A.csr
python converter.py ./ELSES_MATRIX_VCNT900h_20130501/ELSES_MATRIX_VCNT900h_B.mtx VCNT900h_B.csr
python converter.py ./ELSES_MATRIX_VCNT10800h_20130501/ELSES_MATRIX_VCNT10800h_A.mtx VCNT10800h_A.csr
python converter.py ./ELSES_MATRIX_VCNT10800h_20130501/ELSES_MATRIX_VCNT10800h_B.mtx VCNT10800h_B.csr
# Clear .tgz and extract directories
rm -rf ELSES_MATRIX_PPE3594_20160426.tgz ELSES_MATRIX_PPE3594_20160426
rm -rf ELSES_MATRIX_VCNT900h_20130501.tgz ELSES_MATRIX_VCNT900h_20130501
rm -rf ELSES_MATRIX_VCNT10800h_20130501.tgz ELSES_MATRIX_VCNT10800h_20130501
