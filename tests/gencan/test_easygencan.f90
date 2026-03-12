program test_easygencan
   use compute_data
   use input, only : precision
   use test_assertions
   implicit none

   integer, parameter :: n = 6, m = 1
   integer :: wi(n), iter, fcnt, gcnt, cgcnt, inform, trtype
   double precision :: x(n), l(n), u(n), g(n), wd(8 * n), lambda(m), rho(m), f, gpsupn, delmin
   external :: easygencan

   ntotmol = 0
   ntype = 0
   init1 = .true.
   move = .false.
   precision = 1.d-6
   allocate(gxcar(1, 3))
   gxcar = 0.d0

   x = 0.d0
   l = -1.d20
   u = 1.d20
   g = 0.d0
   wd = 0.d0
   wi = 0
   lambda = 0.d0
   rho = 0.d0
   f = -1.d0
   gpsupn = -1.d0
   iter = -1
   fcnt = -1
   gcnt = -1
   cgcnt = -1
   inform = -999
   delmin = 1.d-2
   trtype = 1

   call easygencan(n, x, l, u, m, lambda, rho, 1.d-6, 5, 5, trtype, 0, 5, f, g, gpsupn, iter, fcnt, gcnt, cgcnt, inform, wi, wd, delmin)

   call assert_equal_int(inform, 0, "easygencan should return a non-error status on the empty toy problem")
   call assert_close_real8(f, 0.d0, 1.d-12, "easygencan should keep the objective at zero")

   deallocate(gxcar)
end program test_easygencan
