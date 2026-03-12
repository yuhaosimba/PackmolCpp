program test_comparegrad_probe
   use sizes
   use compute_data
   implicit none

   integer, parameter :: n = 30
   integer :: i
   double precision :: x(n)
   external :: comparegrad

   nn = n
   ntotmol = 0
   ntype = 0
   init1 = .true.
   move = .false.

   allocate(gxcar(1, 3))
   gxcar = 0.d0

   do i = 1, n
      x(i) = dble(i)
   end do

   call comparegrad(n, x)
end program test_comparegrad_probe
