module gsminres_mod
  ! Note:
  ! BIND(C) で complex 配列（例: x(:)）を C に inout 引数として渡す場合、
  ! Fortran の仕様により一時配列が使われることがあり、C 側の書き込みが反映されないことがある。
  ! 確実に書き込み結果を反映させるには、c_loc() により実アドレスを取得して c_ptr として渡す必要がある。
  use iso_c_binding, only: c_size_t, c_ptr, c_null_ptr, c_int, c_double, c_double_complex, c_loc
  implicit none

  private

  public :: gsminres_handle
  public :: gsminres_create, gsminres_destroy
  public :: gsminres_initialize
  public :: gsminres_glanczos_pre, gsminres_glanczos_pst
  public :: gsminres_update
  public :: gsminres_finalize
  public :: gsminres_get_residual

  type :: gsminres_handle
     type(c_ptr) :: ref
  end type gsminres_handle

contains

  function gsminres_create(n, m) result(handle)
    integer(c_size_t), intent(in), value :: n, m
    type(gsminres_handle) :: handle
    interface
       function c_gsminres_create(n, m) bind(C, name="gsminres_create")
         import :: c_size_t, c_ptr
         integer(c_size_t), value :: n, m
         type(c_ptr)              :: c_gsminres_create
       end function c_gsminres_create
    end interface
    handle%ref = c_gsminres_create(n, m)
  end function gsminres_create

  subroutine gsminres_destroy(handle)
    type(gsminres_handle), intent(inout) :: handle
    interface
       subroutine c_gsminres_destroy(h) bind(C, name="gsminres_destroy")
         import :: c_ptr
         type(c_ptr), value :: h
       end subroutine c_gsminres_destroy
    end interface
    call c_gsminres_destroy(handle%ref)
    handle%ref = c_null_ptr
  end subroutine gsminres_destroy

  subroutine gsminres_initialize(handle, x, b, w, sigma, threshold, n, m)
    type(gsminres_handle),     intent(in)            :: handle
    complex(c_double_complex), intent(inout), target :: x(*), w(*)
    complex(c_double_complex), intent(in),    target :: b(*), sigma(*)
    real(c_double),            intent(in)            :: threshold
    integer(c_size_t),         intent(in)            :: n, m
    type(c_ptr) :: xp, bp, wp, sigmap
    interface
       subroutine c_gsminres_initialize(h, x, b, w, sigma, threshold, n, m) bind(C, name="gsminres_initialize")
         import :: c_size_t, c_ptr, c_double
         type(c_ptr),        value :: h
         type(c_ptr),        value :: x, b, w, sigma
         real(c_double),     value :: threshold
         integer(c_size_t),  value :: n, m
       end subroutine c_gsminres_initialize
    end interface
    xp = c_loc(x(1)); bp = c_loc(b(1)); wp = c_loc(w(1)); sigmap = c_loc(sigma(1))
    call c_gsminres_initialize(handle%ref, xp, bp, wp, sigmap, threshold, n, m)
  end subroutine gsminres_initialize

  subroutine gsminres_glanczos_pre(handle, u, n)
    type(gsminres_handle),     intent(in)            :: handle
    complex(c_double_complex), intent(inout), target :: u(*)
    integer(c_size_t),         intent(in)            :: n
    type(c_ptr) :: up
    interface
       subroutine c_gsminres_glanczos_pre(h, u, n) bind(C, name="gsminres_glanczos_pre")
         import :: c_size_t, c_ptr
         type(c_ptr),       value :: h
         type(c_ptr),       value :: u
         integer(c_size_t), value :: n
       end subroutine c_gsminres_glanczos_pre
    end interface
    up = c_loc(u(1))
    call c_gsminres_glanczos_pre(handle%ref, up, n)
  end subroutine gsminres_glanczos_pre

  subroutine gsminres_glanczos_pst(handle, w, u, n)
    type(gsminres_handle),     intent(in)            :: handle
    complex(c_double_complex), intent(inout), target :: w(*), u(*)
    integer(c_size_t),         intent(in)            :: n
    type(c_ptr) :: wp, up
    interface
       subroutine c_gsminres_glanczos_pst(h, w, u, n) bind(C, name="gsminres_glanczos_pst")
         import :: c_size_t, c_ptr
         type(c_ptr),       value :: h
         type(c_ptr),       value :: w, u
         integer(c_size_t), value :: n
       end subroutine c_gsminres_glanczos_pst
    end interface
    wp = c_loc(w(1)); up = c_loc(u(1))
    call c_gsminres_glanczos_pst(handle%ref, wp, up, n)
  end subroutine gsminres_glanczos_pst

  function gsminres_update(handle, x, n, m) result(flag)
    type(gsminres_handle),     intent(in)            :: handle
    complex(c_double_complex), intent(inout), target :: x(*)
    integer(c_size_t),         intent(in)            :: n, m
    integer(c_int) :: flag
    type(c_ptr) :: xp
    interface
       function c_gsminres_update(h, x, n, m) bind(C, name="gsminres_update")
         import :: c_size_t, c_ptr, c_int
         type(c_ptr),       value :: h
         type(c_ptr),       value :: x
         integer(c_size_t), value :: n, m
         integer(c_int)           :: c_gsminres_update
       end function c_gsminres_update
    end interface
    xp = c_loc(x(1))
    !write(*,*) 'Before: ', x(1), x(2)
    flag = c_gsminres_update(handle%ref, xp, n ,m)
    !write(*,*) 'After: ', x(1), x(2)
  end function gsminres_update

  subroutine gsminres_finalize(handle, conv_itr, conv_res, m)
    type(gsminres_handle), intent(in)          :: handle
    integer(c_int),        intent(out), target :: conv_itr(*)
    real(c_double),        intent(out), target :: conv_res(*)
    integer(c_size_t),     intent(in)          :: m
    type(c_ptr) :: itrp, resp
    interface
       subroutine c_gsminres_finalize(h, conv_itr, conv_res, m) bind(C, name="gsminres_finalize")
         import :: c_size_t, c_ptr
         type(c_ptr),       value :: h
         type(c_ptr),       value :: conv_itr, conv_res
         integer(c_size_t), value :: m
       end subroutine c_gsminres_finalize
    end interface
    itrp = c_loc(conv_itr(1)); resp = c_loc(conv_res(1))
    call c_gsminres_finalize(handle%ref, itrp, resp, m)
  end subroutine gsminres_finalize

  subroutine gsminres_get_residual(handle, res, m)
    type(gsminres_handle), intent(in)          :: handle
    real(c_double),        intent(out), target :: res(*)
    integer(c_size_t),     intent(in)          :: m
    type(c_ptr) :: resp
    interface
       subroutine c_gsminres_get_residual(h, res, m) bind(C, name="gsminres_get_residual")
         import :: c_size_t, c_ptr
         type(c_ptr),       value :: h
         type(c_ptr),       value :: res
         integer(c_size_t), value :: m
       end subroutine c_gsminres_get_residual
    end interface
    resp = c_loc(res(1))
    call c_gsminres_get_residual(handle%ref, resp, m)
  end subroutine gsminres_get_residual

end module gsminres_mod
