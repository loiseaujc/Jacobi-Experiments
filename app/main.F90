program main
   use omp_lib, only: omp_get_wtime
   use stdlib_constants, only: pi => pi_dp
   use stdlib_linalg_constants, only: ilp, dp, lk
   use stdlib_math, only: linspace
   use stdlib_io_npy, only: save_npy
   use morton, only: field2vec, lookup_tables
   use Jacobi_Experiments, only: textbook_solver, nocopy_solver, &
                                 otf_norm_solver, doconcurrent_solver, zorder_solver

   implicit none(type, external)
   integer(ilp), parameter :: n = 512
   integer(ilp), parameter :: maxiter = 5000
   integer(ilp) :: i, j, k
   integer(ilp), allocatable :: lut(:, :)
   logical(lk), allocatable :: bc(:)
#if NDIM == 3
   real(dp), allocatable :: b(:, :, :), u(:, :, :), uref(:, :, :)
#else
   real(dp), allocatable :: b(:, :), u(:, :), uref(:, :)
#endif
   real(dp), allocatable :: bvec(:), uvec(:)
   real(dp) :: start_time, end_time
   real(dp), allocatable :: x(:), y(:), z(:)

   !----- Create right-hand side term -----!
    print *, "---------------"
#if NDIM == 3
   print *, "RUNNING 3D CASE."
   allocate(b(n, n, n), uref(n, n, n), source=0.0_dp)
   x = linspace(0, 1, n)
   y = linspace(0, 1, n)
   z = linspace(0, 1, n)
   do concurrent(k=1:n, j=1:n, i=1:n)
      b(i, j, k) = 12*pi**2 * sin(2*pi*x(i))*sin(2*pi*y(j))*sin(2*pi*z(k))
      uref(i, j, k) = sin(2*pi*x(i))*sin(2*pi*y(j))*sin(2*pi*z(k))
   end do
#else
    print *, "RUNNING 2D CASE."
   allocate(b(n, n), x(n), y(n), uref(n, n), source=0.0_dp)
   x = linspace(0, 1, n)
   y = linspace(0, 1, n)
   ! Create right-hand side vector.
   do concurrent(j=1:n, i=1:n)
      b(i, j) = 8*pi**2 * sin(2*pi*x(i)) * sin(2*pi*y(j))
      uref(i, j) = sin(2*pi*x(i)) * sin(2*pi*y(j))
   end do
#endif
   print *, "    - Number of points per direction :", n
   print *

   !----------------------------------
   !-----     JACOBI SOLVERS     -----
   !----------------------------------

   !> doconcurrent Jacobi solver.
   start_time = omp_get_wtime()
   u = doconcurrent_solver(b, maxiter)
   end_time = omp_get_wtime()
   print *, "    - Time-to-solution     :", end_time - start_time
   print *, "    - Max. pointwise error :", maxval(abs(u - uref))
   print *

   !> Z-fill ordering.
   print *, "Computing look-up table..."
   call lookup_tables(n, lut, bc)
   bvec = field2vec(b)
   print *, "Computation done."
   print *
   start_time = omp_get_wtime()
   uvec = zorder_solver(n, bvec, lut, bc, maxiter) 
   end_time = omp_get_wtime()
   print *, "    - Time-to-solution     :", end_time - start_time
   print *, "    - Max. pointwise error :", maxval(abs(uvec - field2vec(uref)))

end program main
