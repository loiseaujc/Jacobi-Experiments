submodule(Gauss_Seidel_Experiments) gs_solvers
    implicit none(type, external)
contains
    module procedure lexicographic_solver
    ! Internal variables.
    integer(ilp) :: n, i, j, iteration
#if NDIM == 3
    real(dp), allocatable :: v(:, :, :)
#else
    real(dp), allocatable :: v(:, :)
#endif
    real(dp) :: dx, l2_norm

    ! Sanity check.
    if (any(shape(b) /= size(b, 1))) then
        error stop "Number of points in each direction need to be equal."
    end if

    ! Initialize variables.
    n = size(b, 1)
    dx = 1.0_dp / (n-1)
    l2_norm = 1.0_dp
    iteration = 0

#if NDIM == 3
    allocate(u(n, n, n), v(n, n, n), source=0.0_dp)
#else
    allocate(u(n, n), v(n, n), source=0.0_dp)
#endif

    ! Gauss-Seidel updates.
    do while ((iteration < maxiter) .and. (l2_norm > tol))
        ! Store previous solution for convergence check.
        if (mod(iteration, 1000) == 0) v = u
        ! Gauss-Seidel iteration.
        call lexicographic_kernel(n, u, b, dx)
        ! Compute convergence check.
        if (mod(iteration, 1000) == 0) l2_norm = norm2(u - v)
        ! Update iteration counter.
        iteration = iteration + 1
    end do

    ! Print info.
    print *, "Lexicographic solver :"
    print *, "    - Number of iterations :", iteration
    print *, "    - l2-norm of the error :", norm2(u - v)
    end procedure lexicographic_solver

    !----------------------------------------
    !-----     GAUSS-SEIDEL KERNELS     -----
    !----------------------------------------

    pure subroutine lexicographic_kernel(n, u, b, dx)
        implicit none (type, external)
        integer(ilp), intent(in) :: n
#if NDIM == 3
        real(dp), intent(inout) :: u(:, :, :)
        real(dp), intent(in) :: b(:, :, :), dx
        integer(ilp) :: i, j, k
#else
        real(dp), intent(inout) :: u(:, :)
        real(dp), intent(in) :: b(:, :), dx
        integer(ilp) :: i, j
        do j = 2, n-1
            do i = 2, n-1
                u(i, j) = 0.25_dp*(b(i, j)*dx**2 + (u(i+1, j) + u(i-1, j) &
                                                   + u(i, j+1) + u(i, j-1)))
            end do
        end do
#endif
    end subroutine lexicographic_kernel
end submodule gs_solvers
