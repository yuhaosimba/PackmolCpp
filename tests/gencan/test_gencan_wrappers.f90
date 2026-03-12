program test_gencan_wrappers
   use compute_data
   use input
   use test_assertions
   implicit none

   integer, parameter :: n = 4, nind = 2, m = 1
   integer :: ind(nind), inform
   double precision :: x(n), xc(n), g(n), lambda(m), rho(m), f
   external :: calcf, calcg, calcgdiff

   ntotmol = 0
   ntype = 0
   init1 = .true.
   move = .false.
   fdist = 0.d0
   frest = 0.d0

   allocate(gxcar(1,3))
   gxcar = 0.d0

   ind = (/ 1, 3 /)
   x = (/ 5.d0, 7.d0, -99.d0, -88.d0 /)
   xc = (/ 5.d0, 20.d0, 7.d0, 40.d0 /)
   lambda = 0.d0
   rho = 0.d0
   g = -1.d0

   call calcf(nind, ind, x, n, xc, m, lambda, rho, f, inform)
   call assert_close_real8(f, 0.d0, 1.d-12, "calcf should return zero on the empty toy problem")
   call assert_close_real8(x(1), 5.d0, 1.d-12, "calcf should preserve reduced-space x(1)")
   call assert_close_real8(x(2), 7.d0, 1.d-12, "calcf should preserve reduced-space x(2)")

   call calcg(nind, ind, x, n, xc, m, lambda, rho, g, inform)
   call assert_close_real8(g(1), 0.d0, 1.d-12, "calcg first reduced gradient component")
   call assert_close_real8(g(2), 0.d0, 1.d-12, "calcg second reduced gradient component")

   g = -1.d0
   call calcgdiff(nind, ind, x, n, xc, m, lambda, rho, g, 1.d-6, 1.d-8, inform)
   call assert_close_real8(g(1), 0.d0, 1.d-12, "calcgdiff first reduced gradient component")
   call assert_close_real8(g(2), 0.d0, 1.d-12, "calcgdiff second reduced gradient component")

   deallocate(gxcar)
end program test_gencan_wrappers
