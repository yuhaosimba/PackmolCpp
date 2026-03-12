program test_gencan_hessian_diff
   use compute_data
   use test_assertions
   implicit none

   integer, parameter :: n = 4, nind = 2, m = 1
   integer :: ind(nind), inform
   double precision :: x(n), xc(n), g(n), d(n), hd(n), xtmp(n), lambda(m), rho(m)
   external :: evalnaldiff, calchddiff

   ntotmol = 0
   ntype = 0
   init1 = .true.
   move = .false.
   fdist = 0.d0
   frest = 0.d0

   allocate(gxcar(1,3))
   gxcar = 0.d0

   x = (/ 5.d0, 7.d0, 11.d0, 13.d0 /)
   lambda = 0.d0
   rho = 0.d0
   g = -1.d0

   call evalnaldiff(n, x, m, lambda, rho, g, 1.d-6, 1.d-8, inform)
   call assert_equal_int(inform, 0, "evalnaldiff should finish cleanly on the empty toy problem")
   call assert_close_real8(g(1), 0.d0, 1.d-12, "evalnaldiff first component")
   call assert_close_real8(g(2), 0.d0, 1.d-12, "evalnaldiff second component")
   call assert_close_real8(g(3), 0.d0, 1.d-12, "evalnaldiff third component")
   call assert_close_real8(g(4), 0.d0, 1.d-12, "evalnaldiff fourth component")
   call assert_close_real8(x(1), 5.d0, 1.d-12, "evalnaldiff should preserve x(1)")
   call assert_close_real8(x(2), 7.d0, 1.d-12, "evalnaldiff should preserve x(2)")

   ind = (/ 1, 3 /)
   x = (/ 5.d0, 7.d0, -99.d0, -88.d0 /)
   xc = (/ 5.d0, 20.d0, 7.d0, 40.d0 /)
   d = (/ 2.d0, -3.d0, 99.d0, 88.d0 /)
   g = 0.d0
   hd = -1.d0
   xtmp = 0.d0

   call calchddiff(nind, ind, x, d, g, n, xc, m, lambda, rho, 0, hd, xtmp, 1.d-6, 1.d-8, inform)
   call assert_equal_int(inform, 0, "calchddiff with analytical gradients should finish cleanly")
   call assert_close_real8(hd(1), 0.d0, 1.d-12, "calchddiff analytical first reduced component")
   call assert_close_real8(hd(2), 0.d0, 1.d-12, "calchddiff analytical second reduced component")

   hd = -1.d0
   xtmp = 0.d0
   call calchddiff(nind, ind, x, d, g, n, xc, m, lambda, rho, 1, hd, xtmp, 1.d-6, 1.d-8, inform)
   call assert_equal_int(inform, 0, "calchddiff with finite differences should finish cleanly")
   call assert_close_real8(hd(1), 0.d0, 1.d-12, "calchddiff finite-difference first reduced component")
   call assert_close_real8(hd(2), 0.d0, 1.d-12, "calchddiff finite-difference second reduced component")

   deallocate(gxcar)
end program test_gencan_hessian_diff
