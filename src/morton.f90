module morton
   use, intrinsic:: iso_fortran_env, only: ilp => int32, i8 => int64, dp => real64
   implicit none(type, external)
   private
   public :: lookup_tables
   public :: field2vec

contains

   subroutine lookup_tables(n, lut, bc)
      implicit none(type, external)
      integer(ilp), intent(in) :: n
      integer(ilp), allocatable, intent(out) :: lut(:, :)
      logical, allocatable, intent(out) :: bc(:)
      integer(ilp) :: i, j, k, npts
      ! Initialize arrays.
#if NDIM == 3
      integer(i8) :: z
      npts = n**3
      allocate (lut(6, 0:npts - 1), source=0)
#else
      integer(ilp) :: z
      npts = n**2
      allocate (lut(4, 0:npts - 1), source=0)
#endif
      allocate (bc(0:npts - 1), source=.false.)

      ! Loop through every points in the domain.
#if NDIM == 3
      do concurrent(k=0:n - 1, j=0:n - 1, i=0:n - 1)
         ! Get Morton index.
         z = morton3d(int(i, kind=i8), int(j, kind=i8), int(k, kind=i8))

         ! Is point on the boundaries?
         if ((i == 0) .or. (i == n - 1) .or. &
             (j == 0) .or. (j == n - 1) .or. &
             (k == 0) .or. (k == n - 1)) bc(z) = .true.

         ! Compute Morton index of the neighbors.
         lut(1, z) = int(morton3d(int(i + 1, kind=i8), int(j, kind=i8), int(k, kind=i8)), kind=ilp)
         lut(2, z) = int(morton3d(int(i - 1, kind=i8), int(j, kind=i8), int(k, kind=i8)), kind=ilp)
         lut(3, z) = int(morton3d(int(i, kind=i8), int(j + 1, kind=i8), int(k, kind=i8)), kind=ilp)
         lut(4, z) = int(morton3d(int(i, kind=i8), int(j - 1, kind=i8), int(k, kind=i8)), kind=ilp)
         lut(5, z) = int(morton3d(int(i, kind=i8), int(j, kind=i8), int(k + 1, kind=i8)), kind=ilp)
         lut(6, z) = int(morton3d(int(i, kind=i8), int(j, kind=i8), int(k - 1, kind=i8)), kind=ilp)
      end do
#else
      do concurrent(z=0:npts - 1)
         ! Decode Morton index.
         call imc(z, i, j)

         ! Is point on the boundaries?
         if ((i == 0) .or. (i == n - 1) .or. (j == 0) .or. (j == n - 1)) bc(z) = .true.

         ! Compute Morton index of the neighbours.
         lut(1, z) = mc(i + 1, j)
         lut(2, z) = mc(i - 1, j)
         lut(3, z) = mc(i, j + 1)
         lut(4, z) = mc(i, j - 1)
      end do
#endif
   end subroutine lookup_tables

   function field2vec(b) result(out)
#if NDIM == 3
      real(dp), intent(in) :: b(:, :, :)
      integer(ilp) :: z
#else
      real(dp), intent(in) :: b(:, :)
      integer(ilp) :: z
#endif
      real(dp), allocatable :: out(:)
      integer(ilp) :: i, j, k, n

      ! Number of points in the domain.
      n = product(shape(b))
      allocate (out(0:n - 1), source=0.0_dp)

#if NDIM == 3
      do concurrent(k=0:size(b, 3) - 1, j=0:size(b, 2) - 1, i=0:size(b, 1) - 1)
         z = int(morton3d(int(i, kind=i8), int(j, kind=i8), int(k, kind=i8)), kind=ilp)
         out(z) = b(i + 1, j + 1, k + 1)
      end do
#else
      ! Transform to vector representation.
      do concurrent(z=0:n - 1)
         call imc(z, i, j)
         out(z) = b(i + 1, j + 1)
      end do
#endif
   end function field2vec

   pure function part1by2(n) result(r)
      integer(i8), intent(in) :: n
      integer(i8)             :: r

      r = iand(n, int(z'1fffff', i8))          ! keep only 21 bits

      r = iand(ior(r, ishft(r, 32)), int(z'1f00000000ffff', i8))
      r = iand(ior(r, ishft(r, 16)), int(z'1f0000ff0000ff', i8))
      r = iand(ior(r, ishft(r, 8)), int(z'100f00f00f00f00f', i8))
      r = iand(ior(r, ishft(r, 4)), int(z'10c30c30c30c30c3', i8))
      r = iand(ior(r, ishft(r, 2)), int(z'1249249249249249', i8))
   end function part1by2

   ! Compute 3D Morton code from (x,y,z).
   pure function morton3d(x, y, z) result(code)
      integer(i8), intent(in) :: x, y, z
      integer(i8)             :: code

      code = ior(ior(part1by2(x), ishft(part1by2(y), 1)), &
                 ishft(part1by2(z), 2))
   end function morton3d

   pure elemental integer function mc(x, y)
      integer, value :: x, y
      integer :: i
      mc = 0
      do i = 0, bit_size(x)
         block
            integer :: mi, mx, my
            mi = lshift(1, i)
            mx = iand(x, mi)
            my = iand(y, mi)
            mc = ior(mc, lshift(mx, i))
            mc = ior(mc, lshift(my, i + 1))
         end block
      end do
   end function

   pure elemental subroutine imc(m, x, y)
      integer, value :: m
      integer, intent(out) :: x, y
      x = mc_extract(m)
      y = mc_extract(rshift(m, 1))
   end subroutine

   pure elemental integer function mc_extract(m)
      integer, value :: m
      m = iand(m, z'55555555')
      m = iand(ior(m, rshift(m, 1)), z'33333333')
      m = iand(ior(m, rshift(m, 2)), z'0F0F0F0F')
      m = iand(ior(m, rshift(m, 4)), z'00FF00FF')
      m = iand(ior(m, rshift(m, 8)), z'0000FFFF')
      mc_extract = m
   end function
end module morton
