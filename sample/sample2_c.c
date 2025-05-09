#include "gsminres_c_api.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

extern void zcopy_(const int*, const double _Complex*, const int*, double _Complex*, const int*);
extern void zaxpy_(const int*, const double _Complex*, const double _Complex*, const int*, double _Complex*, const int*);
extern double dznrm2_(const int*, const double _Complex*, const int*);
extern double _Complex zdotc_(const int*, const double _Complex*, const int*, const double _Complex*, const int*);
extern void zscal_(const int*, const double _Complex*, double _Complex*, const int*);

void SpMV(const int *A_row, const int *A_col, const double _Complex *A_ele,
          const double _Complex *x, double _Complex *b, int N);
int CG_method(const int *B_row, const int *B_col, const double _Complex *B_ele,
              double _Complex *x, const double _Complex *b,
              int N, const double tol);
FILE* fopen_mtx(const char *fname, const char *mode,
                int *row_size, int *col_size, int *ele_size);
void read_csr(const char *fname,
              int *N, int *DATASIZE,
              int **row_ptr, int **col_ind, double _Complex **element);

int main(int argc, char *argv[]) {
  int n;
  if (argc < 1) {
    fprintf(stderr, "Usage: %s <MTX_file(A)> <MTX_file(B)>\n", argv[0]);
    exit(1);
  }
  // Read Matrix file (CSR format)
  int *A_row, *A_col, A_nnz; double _Complex *A_ele;
  int *B_row, *B_col, B_nnz; double _Complex *B_ele;
  read_csr(argv[1], &n, &A_nnz, &A_row, &A_col, &A_ele);
  read_csr(argv[2], &n, &B_nnz, &B_row, &B_col, &B_ele);
  //read_csr("../data/DIAB18h_A.csr", &n,&A_nnz,&A_row,&A_col,&A_ele);
  //read_csr("../data/DIAB18h_B.csr", &n,&B_nnz,&B_row,&B_col,&B_ele);
  const size_t N=n, M=10;
  // Prepare sigma
  double _Complex *sigma = (double _Complex *)calloc(M, sizeof(double _Complex));
  for(size_t i=0; i<M; i++) sigma[i] = 0.1 * cexp(2 * M_PI * I * (i-0.5) / M);
  // Prepare RHS vector
  double _Complex *b = (double _Complex *)calloc(N, sizeof(double _Complex)*N);
  for(size_t i=0; i<N; i++) b[i] = 1.0;

  // Allocate vectors
  double _Complex *x = (double _Complex *)calloc(N, sizeof(double _Complex)*M*N);
  double _Complex *w = (double _Complex *)calloc(N, sizeof(double _Complex)*N);
  double _Complex *u = (double _Complex *)calloc(N, sizeof(double _Complex)*N);
  int    itr[M];
  double res[M];

  // Create solver
  gsminres_handle solver = gsminres_create(N, M);

  // Initialize
  if (CG_method(B_row,B_col,B_ele, w, b, n, 1e-13) == 0){
    fprintf(stderr, "# Inner CG failed\n");
    exit(1);
  }
  gsminres_initialize(solver, x, b, w, sigma, 1e-13, N, M);

  // Solve
  for (size_t j=0; j<10000; ++j) {
    SpMV(A_row,A_col,A_ele, w, u, N);
    gsminres_glanczos_pre(solver, u, N);
    if (CG_method(B_row,B_col,B_ele, w, u, n, 1e-13) == 0){
      fprintf(stderr, "# Inner CG failed\n");
      exit(1);
    }
    gsminres_glanczos_pst(solver, w, u, N);
    if (gsminres_update(solver, x, N, M))
      break;
    gsminres_get_residual(solver, res, M);
  }
  gsminres_finalize(solver, itr, res, M);

  // Destroy solver
  gsminres_destroy(solver);

  // Check results
  int ONE=1;
  double _Complex *r_tmp = calloc(N, sizeof(double _Complex));
  double _Complex *tmp = calloc(N, sizeof(double _Complex));
  double _Complex cTMP;
  for(size_t k=0; k<M; k++){
    zcopy_(&n, b, &ONE, r_tmp, &ONE);
    SpMV(A_row,A_col,A_ele, &(x[k*N]), tmp, n);
    cTMP = -1.0;
    zaxpy_(&n, &cTMP, tmp, &ONE, r_tmp, &ONE);
    SpMV(B_row,B_col,B_ele, &(x[k*N]), tmp, n);
    cTMP = -sigma[k];
    zaxpy_(&n, &cTMP, tmp, &ONE, r_tmp, &ONE);
    fprintf(stdout, "%2ld %10.6lf %10.6lf %5d %12.5e %12.5e\n",
            k, creal(sigma[k]), cimag(sigma[k]), itr[k],
            res[k], dznrm2_(&n,r_tmp,&ONE));
  }
  return 0;
}


// Utility functions
void SpMV(const int *A_row, const int *A_col, const double _Complex *A_ele,
	  const double _Complex *x, double _Complex *b, int N)
{
  double _Complex tmp;
#pragma omp parallel for private(tmp)
  for(int i=0; i<N; i++){
    tmp = 0.0 + 0.0*_Complex_I;
    for(int j=A_row[i]; j<A_row[i+1]; j++)
      tmp += A_ele[j]*x[A_col[j]];
    b[i] = tmp;
  }
}

int CG_method(const int *B_row, const int *B_col, const double _Complex *B_ele,
              double _Complex *x, const double _Complex *b, int N, const double tol)
{
  int status = 0, ONE=1;
  double _Complex r[N], p[N], Bp[N];
  double _Complex alpha, beta, rr, rr_old, cTMP;
  for(int i=0; i<N; i++){
    x[i] = 0.0 + 0.0*_Complex_I;
    r[i] = b[i];
    p[i] = r[i];
  }
  rr = zdotc_(&N, r, &ONE, r, &ONE);
  double r0_nrm = dznrm2_(&N, r, &ONE);
  for(int j=0; j<10000; j++){
    SpMV(B_row,B_col,B_ele, p, Bp, N);
    alpha = rr / zdotc_(&N, p, &ONE, Bp, &ONE);
    zaxpy_(&N, &alpha, p, &ONE, x, &ONE);
    alpha = -alpha;
    zaxpy_(&N, &alpha, Bp, &ONE, r, &ONE);
    if(dznrm2_(&N, r, &ONE)/r0_nrm <= tol){
      //fprintf(stderr, "[CG] %d %lf %lf\n", j, dznrm2_(&N, r, &ONE), r0_nrm);
      status = 1;
      break;
    }
    alpha = -alpha;
    rr_old = rr;
    rr = zdotc_(&N, r, &ONE, r, &ONE);
    beta = rr / rr_old;
    zscal_(&N, &beta, p, &ONE);
    cTMP = 1.0 + 0.0*_Complex_I;
    zaxpy_(&N, &cTMP, r, &ONE, p, &ONE);
  }
  return status;
}

FILE *fopen_mtx(const char *fname, const char *mode, int *row_size, int *col_size, int *ele_size)
{
  FILE *fp = NULL;
  fp = fopen(fname, mode);
  if(fp == NULL){
    fprintf(stderr, "Can not open file : %s\n", fname);
    exit(1);
  }
  char chr;
  while( (chr = fgetc(fp)) != EOF && (chr=='%' || chr=='#') )
    while( (chr = fgetc(fp)) != EOF )
      if(chr == '\n') break;
  fseek(fp, -sizeof(char), SEEK_CUR);
  int num, tmp1, tmp2, tmp3;
  num = fscanf(fp, "%d %d %d", &tmp1, &tmp2, &tmp3);
  if(num != 3){ fprintf(stderr, "fopen err\n"); exit(1);}
  if(col_size != NULL) *row_size = tmp1;
  if(row_size != NULL) *col_size = tmp2;
  if(ele_size != NULL) *ele_size = tmp3;
  return fp;
}
void read_csr(const char *fname, int *N, int *DATASIZE,
              int **row_ptr, int **col_ind, double _Complex **element)
{
  FILE *fp;
  int rsize, csize, esize;
  fp = fopen_mtx(fname, "r", &rsize, &csize, &esize);
  *N = rsize-1; *DATASIZE = esize;

  *row_ptr = (int *)calloc(rsize, sizeof(int));
  *col_ind = (int *)calloc(csize, sizeof(int));
  *element = (double _Complex *)calloc(esize, sizeof(double _Complex));

  int row, col, i=0;
  double real, imag=0;
  while( fscanf(fp, "%d %d %lf %lf", &row, &col, &real, &imag) != EOF ){
    if(i < rsize) (*row_ptr)[i] = row;
    if(i < csize) (*col_ind)[i] = col;
    if(i < esize) (*element)[i] = real + imag*I;
    i = i + 1;
  }
  fclose(fp);
}
