program main
   use omp_lib, only: omp_get_wtime
   use stdlib_constants, only: pi => pi_dp
   use stdlib_linalg_constants, only: ilp, dp, lk
   use stdlib_math, only: linspace
   use stdlib_io_npy, only: save_npy
   use Jacobi_Experiments, only: textbook_solver, nocopy_solver, &
                                 otf_norm_solver, doconcurrent_solver

   implicit none(type, external)
   integer(ilp), parameter :: n = 128, maxiter = 5000
#if NDIM == 3
   real(dp), allocatable :: b(:, :, :), u(:, :, :)
#else
   real(dp), allocatable :: b(:, :), u(:, :)
#endif
   real(dp), dimension(n) :: x, y, z
   integer(ilp) :: i, j, k
   real(dp) :: start_time, end_time

   !----- Initialize variables -----!
    print *, "---------------"
#if NDIM == 3
   print *, "RUNNING 3D CASE."
   allocate(b(n, n, n), source=0.0_dp)
   x = linspace(0, 1, n)
   y = linspace(0, 1, n)
   z = linspace(0, 1, n)
   do concurrent(i=1:n, j=1:n, k=1:n)
      b(i, j, k) = sin(2*pi*x(i))*sin(2*pi*y(j))*sin(2*pi*z(k))
   end do
#else
    print *, "RUNNING 2D CASE."
   allocate(b(n, n), source=0.0_dp)
   x = linspace(0, 1, n)
   y = linspace(0, 1, n)
   ! Create right-hand side vector.
   do concurrent(i=1:n, j=1:n)
      b(i, j) = sin(2*pi*x(i))*sin(2*pi*y(j))
   end do
#endif
   print *, "  - Number of points per direction :", n
   print *
   b = b/norm2(b)

   !----------------------------------
   !-----     JACOBI SOLVERS     -----
   !----------------------------------

   !> Textbook Jacobi solver.
   start_time = omp_get_wtime()
   u = textbook_solver(b, maxiter)
   end_time = omp_get_wtime()
   print *, "    - Time-to-solution     :", end_time - start_time

   !> No-copy Jacobi solver.
   start_time = omp_get_wtime()
   u = nocopy_solver(b, maxiter)
   end_time = omp_get_wtime()
   print *, "    - Time-to-solution     :", end_time - start_time

   !> doconcurrent Jacobi solver.
   start_time = omp_get_wtime()
   u = doconcurrent_solver(b, maxiter)
   end_time = omp_get_wtime()
   print *, "    - Time-to-solution     :", end_time - start_time

end program main
