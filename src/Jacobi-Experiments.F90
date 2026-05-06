module Jacobi_Experiments
   use stdlib_linalg_constants, only: dp, ilp, lk
   implicit none(external)
   private

   integer, parameter :: i8 = selected_int_kind(18)  ! 64-bit integer kind
   real(dp), parameter :: eps = epsilon(1.0_dp)
   real(dp), parameter :: tol = sqrt(eps)

   interface
      module function textbook_solver(b, maxiter) result(u)
         implicit none(external)
#if NDIM == 3
         real(dp), intent(in) :: b(:, :, :)
         real(dp), allocatable :: u(:, :, :)
#else
         real(dp), intent(in) :: b(:, :)
         real(dp), allocatable :: u(:, :)
#endif
         integer(ilp), intent(in) :: maxiter
      end function textbook_solver

      module function nocopy_solver(b, maxiter) result(u)
         implicit none(external)
#if NDIM == 3
         real(dp), intent(in) :: b(:, :, :)
         real(dp), allocatable :: u(:, :, :)
#else
         real(dp), intent(in) :: b(:, :)
         real(dp), allocatable :: u(:, :)
#endif
         integer(ilp), intent(in) :: maxiter
      end function nocopy_solver

      module function otf_norm_solver(b, maxiter) result(u)
         implicit none(external)
#if NDIM == 3
         real(dp), intent(in) :: b(:, :, :)
         real(dp), allocatable :: u(:, :, :)
#else
         real(dp), intent(in) :: b(:, :)
         real(dp), allocatable :: u(:, :)
#endif
         integer(ilp), intent(in) :: maxiter
      end function otf_norm_solver

      module function doconcurrent_solver(b, maxiter) result(u)
         implicit none(external)
#if NDIM == 3
         real(dp), intent(in) :: b(:, :, :)
         real(dp), allocatable :: u(:, :, :)
#else
         real(dp), intent(in) :: b(:, :)
         real(dp), allocatable :: u(:, :)
#endif
         integer(ilp), intent(in) :: maxiter
      end function doconcurrent_solver

      module function zorder_solver(n, b, lut, bc, maxiter) result(u)
        implicit none(external)
        integer(ilp), intent(in) :: n
        real(dp), intent(in) :: b(:)
        integer(ilp), intent(in) :: lut(:, :)
        logical(lk), intent(in) :: bc(:)
        integer(ilp), intent(in) :: maxiter
        real(dp), allocatable :: u(:)
      end function zorder_solver
   end interface

   public :: textbook_solver
   public :: nocopy_solver
   public :: otf_norm_solver
   public :: doconcurrent_solver
   public :: zorder_solver

end module Jacobi_Experiments
