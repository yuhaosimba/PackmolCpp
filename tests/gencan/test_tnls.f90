program test_tnls
   use compute_data
   use test_assertions
   implicit none

   integer, parameter :: n = 6, m = 1
   integer :: ind(n), fcnt, gcnt, intcnt, exgcnt, exbcnt, inform
   double precision :: x(n), g(n), l(n), u(n), d(n), xplus(n), xtmp(n), xbext(n)
   double precision :: lambda(m), rho(m), f, amax
   external :: tnls

   ntotmol = 0
   ntype = 0
   init1 = .true.
   move = .false.
   allocate(gxcar(1, 3))
   gxcar = 0.d0

   ind = (/ 1, 2, 3, 4, 5, 6 /)
   x = (/ 1.d0, 2.d0, 3.d0, 4.d0, 5.d0, 6.d0 /)
   g = 0.d0
   l = -1.d20
   u = 1.d20
   d = 0.d0
   xplus = 0.d0
   xtmp = 0.d0
   xbext = 0.d0
   lambda = 0.d0
   rho = 0.d0
   f = 0.d0
   amax = 1.d0
   fcnt = 0
   gcnt = 0
   intcnt = 0
   exgcnt = 0
   exbcnt = 0

   call tnls(n, ind, n, x, m, lambda, rho, l, u, f, g, d, amax, 1, 1, 2.d0, 2.d0, 4, 5, -1.d0, 5, 0, 0, fcnt, gcnt, intcnt, exgcnt, exbcnt, inform, xplus, xtmp, xbext, 1.d-4, 0.5d0, 0.1d0, 0.9d0, 1.d-7, 1.d-10, 1.d-10, 1.d-20, 1.d20, 1.d99)

   call assert_equal_int(inform, 0, "tnls should accept the zero step on the empty toy problem")
   call assert_close_real8(f, 0.d0, 1.d-12, "tnls should preserve the zero objective")

   deallocate(gxcar)
end program test_tnls
