submodule(Jacobi_Experiments) solvers
   implicit none(type, external)
contains
   !----------------------------------
   !-----     JACOBI SOLVERS     -----
   !----------------------------------

   module procedure textbook_solver
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
   dx = 1.0_dp/(n - 1)
   l2_norm = 1.0_dp
   iteration = 0

#if NDIM == 3
   allocate (u(n, n, n), v(n, n, n), source=0.0_dp)
#else
   allocate (u(n, n), v(n, n), source=0.0_dp)
#endif

   ! Jacobi updates.
   do while ((iteration < maxiter) .and. (l2_norm > tol))
      ! Jacobi iteration.
      call textbook_kernel(n, v, u, b, dx)
      ! Compute error norm.
      if (mod(iteration, 1000) == 0) l2_norm = norm2(u - v)
      ! Update variable.
      u = v
      ! Update iteration counter.
      iteration = iteration + 1
   end do

   ! Print info.
   print *, "Textbook solver :"
   print *, "    - Number of iterations :", iteration
   print *, "    - l2-norm of the error :", l2_norm
   end procedure textbook_solver

   module procedure nocopy_solver
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
   l2_norm = 1.0_dp
   iteration = 0
   dx = 1.0_dp/(n - 1)
#if NDIM == 3
   allocate (u(n, n, n), v(n, n, n), source=0.0_dp)
#else
   allocate (u(n, n), v(n, n), source=0.0_dp)
#endif

   ! Jacobi updates.
   do while ((iteration < maxiter) .and. (l2_norm > tol))
      ! Jacobi kernel.
      call textbook_kernel(n, v, u, b, dx)
      ! Update variables.
      call textbook_kernel(n, u, v, b, dx)
      ! Compute error norm.
      if (mod(iteration, 1000) == 0) l2_norm = norm2(u - v)
      ! Update iteration counter.
      iteration = iteration + 2
   end do

   ! Print info.
   print *, "No-copy solver :"
   print *, "    - Number of iterations :", iteration
   print *, "    - l2-norm of the error :", l2_norm
   end procedure nocopy_solver

   module procedure otf_norm_solver
   ! Internal variables.
   integer(ilp) :: n, i, j, iteration
#if NDIM == 3
   real(dp), allocatable :: v(:, :, :)
#else
   real(dp), allocatable :: v(:, :)
#endif
   real(dp) :: dx, l2_norm

   ! Initialize variables
   if (any(shape(b) /= size(b, 1))) then
      error stop "Number of points in each direction need to be equal."
   end if

   ! Initialize variables.
   n = size(b, 1)
   dx = 1.0_dp/(n - 1)
   l2_norm = 1.0_dp
   iteration = 0
#if NDIM == 3
   allocate (u(n, n, n), v(n, n, n), source=0.0_dp)
#else
   allocate (u(n, n), v(n, n), source=0.0_dp)
#endif

   ! Jacobi updates.
   do while ((iteration < maxiter) .and. (l2_norm > tol))
      ! Jacobi iteration.
      call textbook_kernel(n, v, u, b, dx)
      ! Update variable.
      call textbook_kernel(n, u, v, b, dx)
      ! Update iteration counter.
      iteration = iteration + 2
      ! Compute error norm.
      if (mod(iteration, 1000) == 0) l2_norm = norm2(u - v)
   end do

   ! Print info.
   print *, "On-the-fly solver :"
   print *, "    - Number of iterations :", iteration
   print *, "    - l2-norm of the error :", l2_norm
   end procedure otf_norm_solver

   module procedure doconcurrent_solver
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
   dx = 1.0_dp/(n - 1)
   l2_norm = 1.0_dp
   iteration = 0
#if NDIM == 3
   allocate (u(n, n, n), v(n, n, n), source=0.0_dp)
#else
   allocate (u(n, n), v(n, n), source=0.0_dp)
#endif

   ! Jacobi updates.
   do while ((iteration < maxiter) .and. (l2_norm > tol))
      ! Jacobi kernel.
      call doconcurrent_kernel(n, v, u, b, dx)
      ! Update variables.
      call doconcurrent_kernel(n, u, v, b, dx)
      ! Update iteration counter.
      iteration = iteration + 2
      ! Compute error norm.
      if (mod(iteration, 1000) == 0) l2_norm = norm2(u - v)
   end do

   ! Print info.
   print *, "Do-concurrent solver :"
   print *, "    - Number of iterations :", iteration
   print *, "    - l2-norm of the error :", norm2(u - v)
   end procedure doconcurrent_solver

   module procedure zorder_solver
   ! Internal variables.
   integer(ilp) :: iteration, npts
   real(dp), allocatable :: v(:)
   real(dp) :: dx, l2_norm

   ! Sanity check.
#if NDIM == 3
   if (size(b) /= n**3) then
      error stop "Number of points in each direction need to be equal."
   end if
#else
   if (size(b) /= n**2) then
      error stop "Number of points in each direction need to be equal."
   end if
#endif

   ! Initialize variables.
   npts = size(b)
   dx = 1.0_dp/(n - 1)
   l2_norm = 1.0_dp
   iteration = 0
   allocate (u(npts), v(npts), source=0.0_dp)

   ! Jacobi updates.
   do while ((iteration < maxiter) .and. (l2_norm > tol))
      ! Jacobi kernel.
      call zorder_kernel(n, v, u, b, lut, bc, dx)
      ! Update variables.
      call zorder_kernel(n, u, v, b, lut, bc, dx)
      ! Update iteration counter.
      iteration = iteration + 2
      ! Compute error norm.
      if (mod(iteration, 1000) == 0) l2_norm = norm2(u - v)
   end do

   ! Print info.
   print *, "Z-order solver :"
   print *, "    - Number of iterations :", iteration
   print *, "    - l2-norm of the error :", norm2(u - v)
   end procedure zorder_solver

   !----------------------------------
   !-----     JACOBI KERNELS     -----
   !----------------------------------

   pure subroutine textbook_kernel(n, u, v, b, dx)
      implicit none(external)
      integer(ilp), intent(in) :: n
#if NDIM == 3
      real(dp), intent(out) :: u(n, n, n)
      real(dp), intent(in)  :: v(n, n, n), b(n, n, n), dx
      integer(ilp) :: i, j, k
      do k = 2, n - 1
         do j = 2, n - 1
            do i = 2, n - 1
               u(i, j, k) = 1.0_dp/6.0_dp*(b(i, j, k)*dx**2 + (v(i + 1, j, k) + v(i - 1, j, k) &
                                                               + v(i, j + 1, k) + v(i, j - 1, k) &
                                                               + v(i, j, k + 1) + v(i, j, k - 1)))
            end do
         end do
      end do
#else
      real(dp), intent(out) :: u(n, n)
      real(dp), intent(in) :: v(n, n), b(n, n), dx
      integer(ilp) :: i, j
      do j = 2, n - 1
         do i = 2, n - 1
            u(i, j) = 0.25_dp*(b(i, j)*dx**2 + (v(i + 1, j) + v(i - 1, j) &
                                                + v(i, j + 1) + v(i, j - 1)))
         end do
      end do
#endif
   end subroutine textbook_kernel

   pure subroutine doconcurrent_kernel(n, u, v, b, dx)
      implicit none(external)
      integer(ilp), intent(in) :: n
#if NDIM == 3
      real(dp), intent(out) :: u(n, n, n)
      real(dp), intent(in)  :: v(n, n, n), b(n, n, n), dx
      integer(ilp) :: i, j, k
      do concurrent(k=2:n - 1, j=2:n - 1, i=2:n - 1)
         u(i, j, k) = 1.0_dp/6.0_dp*(b(i, j, k)*dx**2 + (v(i + 1, j, k) + v(i - 1, j, k) &
                                                         + v(i, j + 1, k) + v(i, j - 1, k) &
                                                         + v(i, j, k + 1) + v(i, j, k - 1)))
      end do
#else
      real(dp), intent(out) :: u(n, n)
      real(dp), intent(in) :: v(n, n), b(n, n), dx
      integer(ilp) :: i, j
      do concurrent(j=2:n - 1, i=2:n - 1)
         u(i, j) = 0.25_dp*(b(i, j)*dx**2 + (v(i + 1, j) + v(i - 1, j) &
                                             + v(i, j + 1) + v(i, j - 1)))
      end do
#endif
   end subroutine doconcurrent_kernel

   pure subroutine zorder_kernel(n, u, v, b, lut, bc, dx)
      implicit none(external)
      integer(ilp), intent(in) :: n
#if NDIM == 3
      real(dp), intent(out) :: u(0:n**3 - 1)
      real(dp), intent(in) :: v(0:n**3 - 1), b(0:n**3 - 1)
      integer(ilp), intent(in) :: lut(6, 0:n**3 - 1)
      logical, intent(in) :: bc(0:n**3 - 1)
      real(dp), intent(in) :: dx
      integer(ilp) :: z

      do concurrent(z=0:n**3 - 1)
         if (.not. bc(z)) u(z) = 1.0_dp/6.0_dp*(b(z)*dx**2 + sum(v(lut(:, z))))
      end do
#else
      real(dp), intent(out) :: u(0:n**2 - 1)
      real(dp), intent(in) :: v(0:n**2 - 1), b(0:n**2 - 1)
      integer(ilp), intent(in) :: lut(4, 0:n**2 - 1)
      logical, intent(in) :: bc(0:n**2 - 1)
      real(dp), intent(in) :: dx
      integer(ilp) :: z

      do concurrent(z=0:n**2 - 1,.not. bc(z))
         u(z) = 0.25_dp*(b(z)*dx**2 + sum(v(lut(:, z))))
      end do
#endif
   end subroutine zorder_kernel

end submodule solvers
