program test_gencan_calchd
   use compute_data
   use test_assertions
   implicit none

   integer, parameter :: n = 4, nind = 2, m = 1
   integer :: ind(nind), inform
   double precision :: x(n), xc(n), d(n), g(n), hd(n), xtmp(n), lambda(m), rho(m)
   external :: calchd

   ntotmol = 0
   ntype = 0
   init1 = .true.
   move = .false.

   ind = (/ 1, 3 /)
   x = (/ 5.d0, 7.d0, -99.d0, -88.d0 /)
   xc = (/ 5.d0, 20.d0, 7.d0, 40.d0 /)
   d = (/ 2.d0, -3.d0, 99.d0, 88.d0 /)
   g = (/ 11.d0, 13.d0, 17.d0, 19.d0 /)
   hd = (/ 23.d0, 29.d0, 31.d0, 37.d0 /)
   xtmp = 0.d0
   lambda = 0.d0
   rho = 0.d0
   inform = 0

   call calchd(nind, ind, x, d, g, n, xc, m, lambda, rho, hd, xtmp, 1.d-6, 1.d-8, inform)

   call assert_close_real8(x(1), 5.d0, 1.d-12, "calchd should preserve reduced-space x(1)")
   call assert_close_real8(x(2), 7.d0, 1.d-12, "calchd should preserve reduced-space x(2)")
   call assert_close_real8(d(1), 2.d0, 1.d-12, "calchd should preserve reduced-space d(1)")
   call assert_close_real8(d(2), -3.d0, 1.d-12, "calchd should preserve reduced-space d(2)")
   call assert_close_real8(g(1), 11.d0, 1.d-12, "calchd should preserve reduced-space g(1)")
   call assert_close_real8(g(2), 13.d0, 1.d-12, "calchd should restore g(2) after expand/shrink")
   call assert_close_real8(hd(1), 23.d0, 1.d-12, "calchd should keep the first reduced-space component in the current stub")
   call assert_close_real8(hd(2), 31.d0, 1.d-12, "calchd should keep the second reduced-space component in the current stub")
end program test_gencan_calchd
