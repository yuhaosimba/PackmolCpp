program test_restmol
   use sizes
   use compute_data
   use usegencan
   use test_assertions
   implicit none

   integer :: n
   double precision :: fx, x(6)
   external :: restmol

   nn = 6
   n = 6
   ntype = 1
   ntotmol = 1
   init1 = .false.

   allocate(nmols(1), natoms(1), idfirst(1), comptype(1), compsafe(1), xmol(6), gxcar(1,3))
   nmols(1) = 1
   natoms(1) = 0
   idfirst(1) = 1
   comptype(1) = .true.
   gxcar = 0.d0
   iprint1 = 0
   iprint2 = 0

   x = (/ 1.d0, 2.d0, 3.d0, 4.d0, 5.d0, 6.d0 /)
   fx = -1.d0
   call restmol(1, 0, n, x, fx, .false.)

   call assert_equal_int(n, 6, "restmol should restore the original variable count")
   call assert_equal_int(ntotmol, 1, "restmol should restore ntotmol")
   call assert_equal_logical(init1, .false., "restmol should restore init1")
   call assert_close_real8(fx, 0.d0, 1.d-12, "restmol solve=false should evaluate the empty toy problem")
   call assert_close_real8(x(1), 1.d0, 1.d-12, "restmol should preserve x(1)")
   call assert_close_real8(x(6), 6.d0, 1.d-12, "restmol should preserve x(6)")

   deallocate(nmols, natoms, idfirst, comptype, compsafe, xmol, gxcar)
end program test_restmol
