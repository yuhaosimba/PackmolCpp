program test_cg_limits
   use compute_data
   use test_assertions
   implicit none

   integer, parameter :: n = 6, m = 1
   integer :: ind(n), iter, rbdtype, rbdind, inform
   double precision :: x(n), g(n), l(n), u(n), s(n), w(n), y(n), r(n), d(n), sprev(n)
   double precision :: lambda(m), rho(m)
   external :: cg

   ntotmol = 0
   ntype = 0
   init1 = .true.
   move = .false.
   allocate(gxcar(1, 3))
   gxcar = 0.d0

   ind = (/ 1, 2, 3, 4, 5, 6 /)
   x = (/ 1.d0, 2.d0, 3.d0, 4.d0, 5.d0, 6.d0 /)
   g = (/ 1.d0, -1.d0, 2.d0, -2.d0, 0.5d0, -0.5d0 /)
   l = -1.d20
   u = 1.d20
   s = 123.d0
   w = 0.d0
   y = 0.d0
   r = 0.d0
   d = 0.d0
   sprev = 0.d0
   lambda = 0.d0
   rho = 0.d0

   call cg(n, ind, n, x, m, lambda, rho, g, 1.d0, l, u, 1.d-4, 1.d-4, 0, 5, .false., 0, 1, 1, 0, n, s, iter, rbdtype, rbdind, inform, w, y, r, d, sprev, 1.d-6, 1.d-7, 1.d-10, 1.d-10, 1.d-20, 1.d20, 1.d99)

   call assert_true(inform >= 0, "cg should return a non-negative status on bounded non-zero-gradient input")
   call assert_true(any(abs(s + 123.d0) > 1.d-12), "cg should update the step vector from the initialized sentinel values")

   deallocate(gxcar)
end program test_cg_limits
