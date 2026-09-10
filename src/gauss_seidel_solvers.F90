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

    module procedure peeled_solver
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
        call peeled_kernel(n, u, b, dx)
        ! Compute convergence check.
        if (mod(iteration, 1000) == 0) l2_norm = norm2(u - v)
        ! Update iteration counter.
        iteration = iteration + 1
    end do

    ! Print info.
    print *, "Peeled solver :"
    print *, "    - Number of iterations :", iteration
    print *, "    - l2-norm of the error :", norm2(u - v)
    end procedure peeled_solver

    module procedure double_peeled_solver
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
    call double_peeled_kernel(n, u, b, dx)
        ! Compute convergence check.
        if (mod(iteration, 1000) == 0) l2_norm = norm2(u - v)
        ! Update iteration counter.
        iteration = iteration + 1
    end do

    ! Print info.
    print *, "Double peeled solver :"
    print *, "    - Number of iterations :", iteration
    print *, "    - l2-norm of the error :", norm2(u - v)
    end procedure double_peeled_solver

    module procedure redblack_solver
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
    call redblack_kernel(n, u, b, dx)
        ! Compute convergence check.
        if (mod(iteration, 1000) == 0) l2_norm = norm2(u - v)
        ! Update iteration counter.
        iteration = iteration + 1
    end do

    ! Print info.
    print *, "Red/Black solver :"
    print *, "    - Number of iterations :", iteration
    print *, "    - l2-norm of the error :", norm2(u - v)
    end procedure redblack_solver

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
        real(dp) :: dx2
        dx2 = dx**2
        do j = 2, n-1
            do i = 2, n-1
                u(i, j) = 0.25_dp*(b(i, j)*dx2 + (u(i+1, j) + u(i-1, j) &
                                               + u(i, j+1) + u(i, j-1)))
            end do
        end do
#endif
    end subroutine lexicographic_kernel

    pure subroutine peeled_kernel(n, u, b, dx)
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
        real(dp) :: tmp1, tmp2, dx2
        dx2 = dx**2
        do j = 2, n-1
            tmp1 = dx2*b(2, j) + (u(3, j) + u(2, j+1) + u(2, j-1))
            do i = 2, n-2
                tmp2 = dx2*b(i+1, j) + (u(i+2, j) + u(i+1, j+1) + u(i+1, j-1))
                u(i, j) = 0.25_dp*(u(i-1, j) + tmp1)
                tmp1 = tmp2
            end do
            u(n-1, j) = 0.25_dp*(u(n-2, j) + tmp1)
        end do
#endif
    end subroutine peeled_kernel

    pure subroutine double_peeled_kernel(n, u, b, dx)
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
        real(dp) :: tmp1, tmp2, um, dx2
        real(dp), parameter :: c = 0.25_dp
        real(dp), parameter :: c2 = 0.0625_dp
        dx2 = dx**2
        do j = 2, n-1
            do i = 2, n-2, 2
                um = u(i-1, j)
                tmp1 = dx2*b(i, j) + (u(i+1, j) + u(i, j+1) + u(i, j-1))
                tmp2 = dx2*b(i+1, j) + (u(i+2, j) + u(i+1, j+1) + u(i+1, j-1))
                u(i, j) = c*(um + tmp1)
                u(i+1, j) = c2*um + (c2*tmp1 + c*tmp2)
            end do
            if (mod(n, 2) == 1) then
                tmp1 = dx2*b(n-1, j) + (u(n, j) + u(n-1, j+1) + u(n-1, j-1))
                u(n-1, j) = c*(u(n-2, j) + tmp1)
            end if
        end do
#endif
    end subroutine double_peeled_kernel

    pure subroutine redblack_kernel(n, u, b, dx)
        implicit none (type, external)
        integer(ilp), intent(in) :: n
#if NDIM == 3
        real(dp), intent(inout) :: u(:, :, :)
        real(dp), intent(in) :: b(:, :, :), dx
        integer(ilp) :: i, j, k
#else
        real(dp), intent(inout) :: u(:, :)
        real(dp), intent(in) :: b(:, :), dx
        integer(ilp) :: i, j, k, istart
        real(dp) :: dx2
        dx2 = dx**2
        do k = 0, 1
            ! Red points (k = 0) / Black points (k = 1)
            do concurrent (j=2:n-1) default(none) shared(n, u, b, dx2, k) local(i, istart)
                istart = 2 + mod(k+j, 2)
                do concurrent (i=istart:n-1:2) default(none) shared(u, b, dx2, j)
                    u(i, j) = 0.25_dp*(b(i, j)*dx2 + (u(i+1, j) + u(i-1, j) &
                                                   + u(i, j+1) + u(i, j-1)))
                end do
            end do
        end do
#endif
    end subroutine redblack_kernel
end submodule gs_solvers
