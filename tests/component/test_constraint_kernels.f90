program test_constraint_kernels
   use compute_data
   use pbc
   use test_assertions
   implicit none

   double precision :: f
   external :: comprest, gparc, gwalls

   scale = 2.d0
   scale2 = 3.d0

   allocate(xcart(2, 3), gxcar(2, 3))
   allocate(nratom(1), iratom(1, 1), ityperest(1), restpars(1, 9))
   allocate(radius(2), fscale(2), use_short_radius(2), fixedatom(2))
   allocate(short_radius(2), short_radius_scale(2), ibmol(2), ibtype(2), latomnext(2))
   allocate(comptype(1))

   xcart = 0.d0
   gxcar = 0.d0
   nratom = 1
   iratom(1, 1) = 1
   ityperest(1) = 3
   restpars = 0.d0
   restpars(1, 1:6) = (/ 0.d0, 0.d0, 0.d0, 1.d0, 1.d0, 1.d0 /)

   xcart(1, :) = (/ -1.d0, 0.5d0, 0.5d0 /)
   call comprest(1, f)
   call assert_close_real8(f, 2.d0, 1.d-12, "comprest should penalize distance outside the box")

   gxcar = 0.d0
   call gwalls(1, 1)
   call assert_close_real8(gxcar(1, 1), -4.d0, 1.d-12, "gwalls should push the atom back toward the box")
   call assert_close_real8(gxcar(1, 2), 0.d0, 1.d-12, "gwalls should keep y unchanged for this fixture")
   call assert_close_real8(gxcar(1, 3), 0.d0, 1.d-12, "gwalls should keep z unchanged for this fixture")

   xcart = 0.d0
   xcart(1, 1) = 0.d0
   xcart(2, 1) = 1.d0
   gxcar = 0.d0
   radius = 1.d0
   fscale = 1.d0
   use_short_radius = .false.
   short_radius = 0.d0
   short_radius_scale = 1.d0
   fixedatom = .false.
   ibmol = (/ 1, 2 /)
   ibtype = (/ 1, 1 /)
   latomnext = 0
   comptype(1) = .true.
   pbc_length = (/ 100.d0, 100.d0, 100.d0 /)

   call gparc(1, 2)
   call assert_true(gxcar(1, 1) > 0.d0, "gparc should push the first atom away from overlap")
   call assert_true(gxcar(2, 1) < 0.d0, "gparc should push the second atom in the opposite direction")
   call assert_close_real8(gxcar(1, 1) + gxcar(2, 1), 0.d0, 1.d-12, "gparc should conserve pairwise force balance")

   deallocate(xcart, gxcar, nratom, iratom, ityperest, restpars)
   deallocate(radius, fscale, use_short_radius, fixedatom)
   deallocate(short_radius, short_radius_scale, ibmol, ibtype, latomnext, comptype)
end program test_constraint_kernels
