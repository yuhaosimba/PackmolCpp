program test_spgls
   use compute_data
   use test_assertions
   implicit none

   integer, parameter :: n = 6, m = 1
   integer :: fcnt, inform
   double precision :: x(n), g(n), l(n), u(n), d(n), xtrial(n)
   double precision :: lambda(m), rho(m), f, lamspg
   external :: spgls

   ntotmol = 0
   ntype = 0
   init1 = .true.
   move = .false.
   allocate(gxcar(1, 3))
   gxcar = 0.d0

   x = (/ 1.d0, 2.d0, 3.d0, 4.d0, 5.d0, 6.d0 /)
   g = 0.d0
   l = -1.d20
   u = 1.d20
   d = 0.d0
   xtrial = 0.d0
   lambda = 0.d0
   rho = 0.d0
   f = 0.d0
   fcnt = 0
   lamspg = 1.d0

   call spgls(n, x, m, lambda, rho, f, g, l, u, lamspg, 2.d0, 4, -1.d0, 5, 0, fcnt, inform, xtrial, d, 1.d-4, 0.1d0, 0.9d0, 1.d-7, 1.d-10, 1.d-10, 1.d-20, 1.d20, 1.d99)

   call assert_equal_int(inform, 0, "spgls should accept the zero-gradient trial")
   call assert_close_real8(f, 0.d0, 1.d-12, "spgls should preserve the zero objective")
   call assert_close_real8(x(1), 1.d0, 1.d-12, "spgls should preserve x(1)")
   call assert_close_real8(x(6), 6.d0, 1.d-12, "spgls should preserve x(6)")

   deallocate(gxcar)
end program test_spgls
