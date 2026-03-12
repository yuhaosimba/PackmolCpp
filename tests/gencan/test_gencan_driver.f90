program test_gencan_driver
   use compute_data
   use input, only : precision
   use test_assertions
   implicit none

   integer, parameter :: n = 6, m = 1, maxitngp = 5
   integer :: ind(n), inform, fcnt, gcnt, iter, cgcnt, spgiter, spgfcnt
   integer :: tniter, tnfcnt, tnstpcnt, tnintcnt, tnexgcnt, tnexbcnt
   integer :: tnintfe, tnexgfe, tnexbfe
   double precision :: x(n), g(n), l(n), u(n), d(n), s(n), y(n), w(5 * n)
   double precision :: lambda(m), rho(m), f, gpsupn, gpeucn2, lastgpns(0:maxitngp-1)
   external :: gencan

   ntotmol = 0
   ntype = 0
   init1 = .true.
   move = .false.
   precision = 1.d-6
   allocate(gxcar(1, 3))
   gxcar = 0.d0

   ind = (/ 1, 2, 3, 4, 5, 6 /)
   x = 0.d0
   g = 0.d0
   l = -1.d20
   u = 1.d20
   d = 0.d0
   s = 0.d0
   y = 0.d0
   w = 0.d0
   lambda = 0.d0
   rho = 0.d0
   f = -1.d0
   gpsupn = -1.d0
   gpeucn2 = -1.d0
   lastgpns = 0.d0
   iter = -1
   fcnt = -1
   gcnt = -1
   cgcnt = -1
   spgiter = -1
   spgfcnt = -1
   tniter = -1
   tnfcnt = -1
   tnstpcnt = -1
   tnintcnt = -1
   tnexgcnt = -1
   tnexbcnt = -1
   tnintfe = -1
   tnexgfe = -1
   tnexbfe = -1
   inform = -999

   call gencan(n, x, l, u, m, lambda, rho, 0.d0, 1.d-6, 5, 0.d0, maxitngp, -1.d0, 5, 5, -1.d0, -1, 2, 1.d-6, 1.d-1, 1.d-5, 1.d-4, 5, .false., 2.d0, 2.d0, 4, 10, 0, 1, 1, 0, 5, f, g, gpeucn2, gpsupn, iter, fcnt, gcnt, cgcnt, spgiter, spgfcnt, tniter, tnfcnt, tnstpcnt, tnintcnt, tnexgcnt, tnexbcnt, tnintfe, tnexgfe, tnexbfe, inform, s, y, d, ind, lastgpns, w, 0.9d0, 1.d-2, 1.d10, 1.d-10, 1.d-6, 1.d-4, 0.5d0, 0.1d0, 0.9d0, 1.d-7, 1.d-10, 1.d-10, 1.d-20, 1.d20, 1.d99)

   call assert_equal_int(inform, 0, "gencan should return a non-error status on the empty toy problem")
   call assert_close_real8(f, 0.d0, 1.d-12, "gencan should keep the objective at zero")

   deallocate(gxcar)
end program test_gencan_driver
