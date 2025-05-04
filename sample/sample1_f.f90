program sample1_f
  use gsminres_mod
  use iso_c_binding
  implicit none

  ! Parameters
  integer, parameter :: dp = kind(0.0d0)
  integer(c_size_t) :: n, m
  integer :: i, j, info
  complex(c_double_complex), allocatable :: A(:), B(:), r(:)
  complex(c_double_complex), pointer :: x(:), rhs(:), w(:), u(:), sigma(:)
  real(c_double), allocatable :: res(:)
  integer, allocatable :: itr(:)
  complex(c_double_complex) :: ONE, ZERO
  ! External BLAS/LAPACK routines
  external :: zhpmv, zpptrf, zpptrs
  double precision :: dznrm2
  ! Handle
  type(gsminres_handle) :: solver

  ! BLAS/LAPACK constants
  ONE  = cmplx(1.0d0, 0.0d0, kind=c_double_complex)
  ZERO = cmplx(0.0d0, 0.0d0, kind=c_double_complex)

  ! Read matrices A,B form MTX files
  call load_matrix_from_mm("../data/ELSES_MATRIX_DIAB18h_A.mtx", A, n)
  call load_matrix_from_mm("../data/ELSES_MATRIX_DIAB18h_B.mtx", B, n)

  ! Allocate vectors
  m = 10
  allocate(x(n*m), rhs(n), w(n), u(n), sigma(m), r(n))
  allocate(res(m), itr(m))
  rhs = (1.0d0, 0.0d0)
  do i = 1,m
     sigma(i) = 0.1d0 * exp( cmplx(0.0d0, 2*acos(-1.0d0)*(i-0.5d0) / real(m)) )
  end do

  ! Pre-process: Solve Bu = b
  call zpptrf('U', n, B, info)
  if (info /= 0) stop "zpptrf failed"

  ! Initialize solver
  solver = gsminres_create(n, m);
  call gsminres_initialize(solver, x, rhs, w, sigma, 1.0d-13, n, m)
  w = rhs
  call zpptrs('U', n, 1, B, w, n, info)
  call gsminres_get_residual(solver, res, m)
  ! Solving
  do j = 1, 32
     call zhpmv('U', n, ONE, A, w, 1, ZERO, u, 1)
     call gsminres_glanczos_pre(solver, u, n)
     u = w
     call zpptrs('U', n, 1, B, u, n, info)
     if (info /= 0) stop "zpptrs failed"
     call gsminres_glanczos_pst(solver, w, u, n)
     if (gsminres_update(solver, x, n, m) /= 0) then
        exit
     end if
     call gsminres_get_residual(solver, res, m)
  end do

  ! Finalize solver
  call gsminres_finalize(solver, itr, res, m)
  call gsminres_destroy(solver)
  ! Output Results
  do i = 1, m
     call zhpmv('U', n, ONE,      A, x((i-1)*n+1), 1, ZERO, r, 1)
     call zhpmv('U', n, sigma(i), B, x((i-1)*n+1), 1, ONE,  r, 1)
     r = r - rhs
     write(*, '(I2, 1X, 2F10.6, 1X, I5, 1X, 1P, E12.5, 1X, E12.5)') &
          i, real(sigma(i)), aimag(sigma(i)), itr(i), res(i), dznrm2(n,r,1)
  end do

contains

  subroutine load_matrix_from_mm(fname, A, n)
    character(*), intent(in) ::fname
    complex(c_double_complex), allocatable, intent(out) :: A(:)
    integer(c_size_t), intent(out) :: n
    integer :: i, j, k, row,col,nnz
    complex(c_double_complex), allocatable :: full(:,:)
    real(dp) :: re, im
    character(len=256) :: line
    logical :: sym
    integer :: ios
    open(unit=10, file=fname, status='old', action='read', iostat=ios)
    if(ios /= 0) stop "Cannot open MTX file"
    ! skip comments
    do
       read(10, '(A)', iostat=ios) line
       if (ios /= 0) exit
       if (line(1:1) /= '%') exit
    end do
    read(line,*) n, j, nnz
    allocate(full(n,n)); full=(0.0d0, 0.0d0)
    ! read matrix elements
    do i = 1,nnz
       read(10,*) row, col, re, im
       full(row,col) = cmplx(re, im, kind=dp)
       if (row /= col) full(col,row) = conjg(full(row,col))
    end do
    close(10)
    ! convert packed upper triangular
    allocate(A(n*(n+1)/2))
    k = 0
    do j = 1, n
       do i = 1, j
          k = k + 1
          A(k) = full(i,j)
       end do
    end do
  end subroutine load_matrix_from_mm

end program sample1_f
