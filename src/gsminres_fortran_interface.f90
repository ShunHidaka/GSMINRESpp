module gsminres_mod
  use iso_c_binding, only: c_size_t, c_ptr, c_null_ptr, c_int, c_double, c_double_complex
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
    integer(c_size_t), value :: n, m
    type(gsminres_handle)    :: handle
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
         import :: c_ptr, c_null_ptr
         type(c_ptr), value :: h
       end subroutine c_gsminres_destroy
    end interface
    call c_gsminres_destroy(handle%ref)
    handle%ref = c_null_ptr
  end subroutine gsminres_destroy

  subroutine gsminres_initialize(handle, x, b, w, sigma, threshold, n, m)
    type(gsminres_handle),     intent(in)        :: handle
    complex(c_double_complex), intent(inout)     :: x(*), w(*)
    complex(c_double_complex), intent(in)        :: b(*), sigma(*)
    real(c_double),            intent(in)        :: threshold
    integer(c_size_t),         intent(in), value :: n, m
    interface
       subroutine c_gsminres_initialize(h, x, b, w, sigma, threshold, n, m) bind(C, name="gsminres_initialize")
         import :: c_size_t, c_ptr, c_double, c_double_complex
         type(c_ptr),              value :: h
         complex(c_double_complex)       :: x(*), b(*), w(*), sigma(*)
         real(c_double),           value :: threshold
         integer(c_size_t),        value :: n, m
       end subroutine c_gsminres_initialize
    end interface
    call c_gsminres_initialize(handle%ref, x, b, w, sigma, threshold, n, m)
  end subroutine gsminres_initialize

  subroutine gsminres_glanczos_pre(handle, u, n)
    type(gsminres_handle),     intent(in)        :: handle
    complex(c_double_complex), intent(inout)     :: u(*)
    integer(c_size_t),         intent(in), value :: n
    interface
       subroutine c_gsminres_glanczos_pre(h, u, n) bind(C, name="gsminres_glanczos_pre")
         import :: c_size_t, c_ptr, c_double_complex
         type(c_ptr),              value :: h
         complex(c_double_complex)       :: u(*)
         integer(c_size_t),        value :: n
       end subroutine c_gsminres_glanczos_pre
    end interface
    call c_gsminres_glanczos_pre(handle%ref, u, n)
  end subroutine gsminres_glanczos_pre

  subroutine gsminres_glanczos_pst(handle, w, u, n)
    type(gsminres_handle),     intent(in)        :: handle
    complex(c_double_complex), intent(inout)     :: w(*), u(*)
    integer(c_size_t),         intent(in), value :: n
    interface
       subroutine c_gsminres_glanczos_pst(h, w, u, n) bind(C, name="gsminres_glanczos_pst")
         import :: c_size_t, c_ptr, c_double_complex
         type(c_ptr),              value :: h
         complex(c_double_complex)       :: w(*), u(*)
         integer(c_size_t),        value :: n
       end subroutine c_gsminres_glanczos_pst
    end interface
    call c_gsminres_glanczos_pst(handle%ref, w, u, n)
  end subroutine gsminres_glanczos_pst

  function gsminres_update(handle, x, n, m) result(flag)
    type(gsminres_handle),     intent(in)        :: handle
    complex(c_double_complex), intent(inout)     :: x(*)
    integer(c_size_t),         intent(in), value :: n, m
    integer(c_int)                               :: flag
    interface
       function c_gsminres_update(h, x, n, m) bind(C, name="gsminres_update")
         import :: c_size_t, c_ptr, c_int, c_double_complex
         type(c_ptr),              value :: h
         complex(c_double_complex)       :: x(*)
         integer(c_size_t),        value :: n, m
         integer(c_int)                  :: c_gsminres_update
       end function c_gsminres_update
    end interface
    print *, "Before update: x(1) = ", x(1)
    flag = c_gsminres_update(handle%ref, x, n ,m)
    print *, "After update:  x(1) = ", x(1)
  end function gsminres_update

  subroutine gsminres_finalize(handle, conv_itr, conv_res, m)
    type(gsminres_handle), intent(in)        :: handle
    integer(c_int),        intent(out)       :: conv_itr(*)
    real(c_double),        intent(out)       :: conv_res(*)
    integer(c_size_t),     intent(in), value :: m
    interface
       subroutine c_gsminres_finalize(h, conv_itr, conv_res, m) bind(C, name="gsminres_finalize")
         import :: c_size_t, c_ptr, c_int, c_double
         type(c_ptr),       value :: h
         integer(c_int)           :: conv_itr(*)
         real(c_double)           :: conv_res(*)
         integer(c_size_t), value :: m
       end subroutine c_gsminres_finalize
    end interface
    call c_gsminres_finalize(handle%ref, conv_itr, conv_res, m)
  end subroutine gsminres_finalize

  subroutine gsminres_get_residual(handle, res, m)
    type(gsminres_handle), intent(in)        :: handle
    real(c_double),        intent(out)       :: res(m)
    integer(c_size_t),     intent(in), value :: m
    interface
       subroutine c_gsminres_get_residual(h, res, m) bind(C, name="gsminres_get_residual")
         import :: c_size_t, c_ptr, c_double
         type(c_ptr),       value :: h
         real(c_double)           :: res(*)
         integer(c_size_t), value :: m
       end subroutine c_gsminres_get_residual
    end interface
    call c_gsminres_get_residual(handle%ref, res, m)
  end subroutine gsminres_get_residual

end module gsminres_mod
