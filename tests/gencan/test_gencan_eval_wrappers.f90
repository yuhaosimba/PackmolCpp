program test_gencan_eval_wrappers
   use compute_data
   use test_assertions
   implicit none

   integer, parameter :: n = 6, m = 1
   integer :: flag
   double precision :: x(n), g(n), lambda(m), rho(m), f
   external :: evalal, evalnal

   ntotmol = 0
   ntype = 0
   init1 = .true.
   move = .false.
   fdist = -1.d0
   frest = -1.d0

   allocate(gxcar(1, 3))
   gxcar = 0.d0

   x = (/ 1.d0, 2.d0, 3.d0, 4.d0, 5.d0, 6.d0 /)
   lambda = 0.d0
   rho = 0.d0
   g = -1.d0

   call evalal(n, x, m, lambda, rho, f, flag)
   call assert_equal_int(flag, 0, "evalal should return flag zero on the empty toy problem")
   call assert_close_real8(f, 0.d0, 1.d-12, "evalal should return zero on the empty toy problem")
   call assert_close_real8(fdist, 0.d0, 1.d-12, "evalal should propagate zero fdist")
   call assert_close_real8(frest, 0.d0, 1.d-12, "evalal should propagate zero frest")

   call evalnal(n, x, m, lambda, rho, g, flag)
   call assert_equal_int(flag, 0, "evalnal should return flag zero on the empty toy problem")
   call assert_close_real8(g(1), 0.d0, 1.d-12, "evalnal first component")
   call assert_close_real8(g(2), 0.d0, 1.d-12, "evalnal second component")
   call assert_close_real8(g(3), 0.d0, 1.d-12, "evalnal third component")
   call assert_close_real8(g(4), 0.d0, 1.d-12, "evalnal fourth component")
   call assert_close_real8(g(5), 0.d0, 1.d-12, "evalnal fifth component")
   call assert_close_real8(g(6), 0.d0, 1.d-12, "evalnal sixth component")

   deallocate(gxcar)
end program test_gencan_eval_wrappers
