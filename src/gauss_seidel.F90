module Gauss_Seidel_Experiments
    use, intrinsic :: iso_fortran_env, only: ilp => int32, dp => real64
    implicit none (type, external)
    private

    real(dp), parameter :: eps = epsilon(1.0_dp)
    real(dp), parameter :: tol = sqrt(eps)

    interface
        module function lexicographic_solver(b, maxiter) result(u)
            implicit none(type, external)
#if NDIM == 3
            real(dp), intent(in) :: b(:, :, :)
            real(dp), allocatable :: u(:, :, :)
#else
            real(dp), intent(in) :: b(:, :)
            real(dp), allocatable :: u(:, :)
#endif
        integer(ilp), intent(in) :: maxiter
        end function lexicographic_solver

        module function peeled_solver(b, maxiter) result(u)
            implicit none(type, external)
#if NDIM == 3
            real(dp), intent(in) :: b(:, :, :)
            real(dp), allocatable :: u(:, :, :)
#else
            real(dp), intent(in) :: b(:, :)
            real(dp), allocatable :: u(:, :)
#endif
        integer(ilp), intent(in) :: maxiter
        end function peeled_solver

    module function double_peeled_solver(b, maxiter) result(u)
            implicit none(type, external)
#if NDIM == 3
            real(dp), intent(in) :: b(:, :, :)
            real(dp), allocatable :: u(:, :, :)
#else
            real(dp), intent(in) :: b(:, :)
            real(dp), allocatable :: u(:, :)
#endif
        integer(ilp), intent(in) :: maxiter
    end function double_peeled_solver

    module function redblack_solver(b, maxiter) result(u)
            implicit none(type, external)
#if NDIM == 3
            real(dp), intent(in) :: b(:, :, :)
            real(dp), allocatable :: u(:, :, :)
#else
            real(dp), intent(in) :: b(:, :)
            real(dp), allocatable :: u(:, :)
#endif
        integer(ilp), intent(in) :: maxiter
    end function redblack_solver
    end interface

    public :: lexicographic_solver
    public :: peeled_solver
    public :: double_peeled_solver
    public :: redblack_solver

end module Gauss_Seidel_Experiments
