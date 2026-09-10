program main
   use, intrinsic:: iso_fortran_env, only: ilp => int32, dp => real64
   use omp_lib, only: omp_get_wtime

   use morton, only: field2vec, lookup_tables
   use Jacobi_Experiments, only: jacobi_solver => doconcurrent_solver
   use Gauss_Seidel_Experiments, only: lexicographic_solver, peeled_solver, &
                                       double_peeled_solver, redblack_solver

   implicit none(type, external)
   integer(ilp), parameter :: n = 512
   integer(ilp), parameter :: maxiter = n**2
   real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)

   integer(ilp) :: i, j, k
   integer(ilp), allocatable :: lut(:, :)
   logical, allocatable :: bc(:)
#if NDIM == 3
   real(dp), allocatable :: b(:, :, :), u(:, :, :), uref(:, :, :)
#else
   real(dp), allocatable :: b(:, :), u(:, :), uref(:, :)
#endif
   real(dp), allocatable :: bvec(:), uvec(:)
   real(dp) :: start_time, end_time
   real(dp), allocatable :: x(:), y(:), z(:)
   real(dp) :: scale

   !----- Create right-hand side term -----!
   print *, "---------------"
#if NDIM == 3
   print *, "RUNNING 3D CASE."
   allocate (b(n, n, n), uref(n, n, n), source=0.0_dp)
   x = linspace(0.0_dp, 1.0_dp, n)
   y = linspace(0.0_dp, 1.0_dp, n)
   z = linspace(0.0_dp, 1.0_dp, n)
   do concurrent(k=1:n, j=1:n, i=1:n)
      b(i, j, k) = 12*pi**2*sin(2*pi*x(i))*sin(2*pi*y(j))*sin(2*pi*z(k))
      uref(i, j, k) = sin(2*pi*x(i))*sin(2*pi*y(j))*sin(2*pi*z(k))
   end do
#else
   print *, "RUNNING 2D CASE."
   allocate (b(n, n), x(n), y(n), uref(n, n), source=0.0_dp)
   x = linspace(0.0_dp, 1.0_dp, n)
   y = linspace(0.0_dp, 1.0_dp, n)
   ! Create right-hand side vector.
   do concurrent(j=1:n, i=1:n)
      b(i, j) = 8*pi**2*sin(2*pi*x(i))*sin(2*pi*y(j))
      uref(i, j) = sin(2*pi*x(i))*sin(2*pi*y(j))
   end do
#endif
    scale = norm2(b) * (x(2)-x(1))
    b = b / scale
    uref = uref / scale
   print *, "    - Number of points per direction :", n
   print *

   !----------------------------------
   !-----     JACOBI SOLVERS     -----
   !----------------------------------

   !> doconcurrent Jacobi solver.
   start_time = omp_get_wtime()
   u = jacobi_solver(b, maxiter)
   end_time = omp_get_wtime()
   print *, "    - Time-to-solution     :", end_time - start_time
   print *, "    - Max. pointwise error :", maxval(abs(u - uref))
   print *

   !> Textbook Gauss-Seidel solver.
   start_time = omp_get_wtime()
   u = lexicographic_solver(b, maxiter)
   end_time = omp_get_wtime()
   print *, "    - Time-to-solution     :", end_time - start_time
   print *, "    - Max. pointwise error :", maxval(abs(u - uref))
   print *

   !> Peeled Gauss-Seidel solver.
   start_time = omp_get_wtime()
   u = peeled_solver(b, maxiter)
   end_time = omp_get_wtime()
   print *, "    - Time-to-solution     :", end_time - start_time
   print *, "    - Max. pointwise error :", maxval(abs(u - uref))
   print *

   !> Double peeled Gauss-Seidel solver.
   start_time = omp_get_wtime()
   u = double_peeled_solver(b, maxiter)
   end_time = omp_get_wtime()
   print *, "    - Time-to-solution     :", end_time - start_time
   print *, "    - Max. pointwise error :", maxval(abs(u - uref))
   print *

   !> Red/Black Gauss-Seidel solver.
   start_time = omp_get_wtime()
   u = redblack_solver(b, maxiter)
   end_time = omp_get_wtime()
   print *, "    - Time-to-solution     :", end_time - start_time
   print *, "    - Max. pointwise error :", maxval(abs(u - uref))
   print *

contains

   pure function linspace(start, end, n) result(x)
      real(dp), intent(in) :: start, end
      integer(ilp), intent(in) :: n
      real(dp), allocatable :: x(:)
      real(dp) :: dx
      dx = (end - start)/(n - 1)
      x = [(i*dx, i=0, n - 1)]
   end function linspace

end program main
