program test_pgencan
   use sizes
   use compute_data
   use usegencan
   use test_assertions
   implicit none

   integer :: n
   double precision :: fx, x(6)
   external :: pgencan

   nn = 6
   n = 6
   ntype = 1
   ntotmol = 1
   init1 = .true.
   move = .false.
   maxit = 3
   iprint1 = 0
   iprint2 = 0

   allocate(nmols(1), natoms(1), idfirst(1), comptype(1), constrain_rot(1,3), rot_bound(1,3,2))
   allocate(gxcar(1,3), l(n), u(n), g(n), wi(n), wd(8*n))
   nmols(1) = 1
   natoms(1) = 0
   idfirst(1) = 1
   comptype(1) = .true.
   constrain_rot = .false.
   rot_bound = 0.d0
   gxcar = 0.d0

   x = (/ 1.d0, 2.d0, 3.d0, 4.d0, 5.d0, 6.d0 /)
   fx = -1.d0
   call pgencan(n, x, fx)

   call assert_close_real8(fx, 0.d0, 1.d-12, "pgencan should solve the empty toy problem immediately")
   call assert_close_real8(l(1), -1.d20, 1.d10, "pgencan should initialize lower bounds")
   call assert_close_real8(u(1), 1.d20, 1.d10, "pgencan should initialize upper bounds")

   deallocate(nmols, natoms, idfirst, comptype, constrain_rot, rot_bound)
   deallocate(gxcar, l, u, g, wi, wd)
end program test_pgencan
