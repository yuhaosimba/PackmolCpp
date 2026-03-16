program test_tnls_limits
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
   g = (/ 1.d0, -1.d0, 0.5d0, -0.5d0, 0.25d0, -0.25d0 /)
   l = -1.d20
   u = 1.d20
   d = -g
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

   call tnls(n, ind, n, x, m, lambda, rho, l, u, f, g, d, amax, 1, 1, 2.d0, 2.d0, 1, 0, -1.d0, 5, 0, 0, fcnt, gcnt, intcnt, exgcnt, exbcnt, inform, xplus, xtmp, xbext, 1.d-4, 0.5d0, 0.1d0, 0.9d0, 1.d-7, 1.d-10, 1.d-10, 1.d-20, 1.d20, 1.d99)

   call assert_equal_int(inform, 8, "tnls should stop by maxfc budget on a non-zero search direction")
   call assert_true(fcnt >= 0, "tnls should keep a valid function-evaluation counter")

   deallocate(gxcar)
end program test_tnls_limits

