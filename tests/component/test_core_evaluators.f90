program test_core_evaluators
   use compute_data
   use test_assertions
   implicit none

   integer, parameter :: n = 6
   double precision :: x(n), g(n), f
   external :: computef, computeg

   ntotmol = 0
   ntype = 0
   init1 = .true.
   move = .false.
   fdist = -1.d0
   frest = -1.d0

   allocate(gxcar(1, 3))
   gxcar = 0.d0

   x = (/ 1.d0, 2.d0, 3.d0, 4.d0, 5.d0, 6.d0 /)
   g = -1.d0

   call computef(n, x, f)
   call assert_close_real8(f, 0.d0, 1.d-12, "computef should return zero on the empty toy problem")
   call assert_close_real8(fdist, 0.d0, 1.d-12, "computef should reset fdist on the empty toy problem")
   call assert_close_real8(frest, 0.d0, 1.d-12, "computef should reset frest on the empty toy problem")
   call assert_close_real8(x(1), 1.d0, 1.d-12, "computef should preserve x(1)")
   call assert_close_real8(x(6), 6.d0, 1.d-12, "computef should preserve x(6)")

   call computeg(n, x, g)
   call assert_close_real8(g(1), 0.d0, 1.d-12, "computeg first component")
   call assert_close_real8(g(2), 0.d0, 1.d-12, "computeg second component")
   call assert_close_real8(g(3), 0.d0, 1.d-12, "computeg third component")
   call assert_close_real8(g(4), 0.d0, 1.d-12, "computeg fourth component")
   call assert_close_real8(g(5), 0.d0, 1.d-12, "computeg fifth component")
   call assert_close_real8(g(6), 0.d0, 1.d-12, "computeg sixth component")

   deallocate(gxcar)
end program test_core_evaluators
