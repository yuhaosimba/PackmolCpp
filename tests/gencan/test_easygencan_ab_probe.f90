program test_easygencan_ab_probe
   use iso_c_binding, only : c_int
   use compute_data
   use input, only : precision
   implicit none

   interface
      integer(c_int) function packmol_gencan_impl_mode_c() bind(C, name="packmol_gencan_impl_mode_c")
         import :: c_int
      end function packmol_gencan_impl_mode_c
   end interface

   integer, parameter :: n = 6, m = 1
   integer :: wi(n), iter, fcnt, gcnt, cgcnt, inform, trtype, i, mode_id
   double precision :: x(n), l(n), u(n), g(n), wd(8 * n), lambda(m), rho(m), f, gpsupn, delmin
   double precision :: xsum, gnorm2
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
   trtype = 0

   call easygencan(n, x, l, u, m, lambda, rho, 1.d-6, 5, 5, trtype, 0, 5, &
      f, g, gpsupn, iter, fcnt, gcnt, cgcnt, inform, wi, wd, delmin)

   mode_id = packmol_gencan_impl_mode_c()
   xsum = 0.d0
   gnorm2 = 0.d0
   do i = 1, n
      xsum = xsum + x(i)
      gnorm2 = gnorm2 + g(i) * g(i)
   end do

   write(*,'(A,1X,I0,1X,A,1X,I0,1X,A,1X,I0,1X,A,1X,I0,1X,A,1X,I0,1X,A,1X,ES24.16,1X,A,1X,ES24.16,1X,A,1X,ES24.16,1X,A,1X,ES24.16)') &
      'mode', mode_id, 'inform', inform, 'iter', iter, 'fcnt', fcnt, 'gcnt', gcnt, &
      'f', f, 'gpsupn', gpsupn, 'xsum', xsum, 'gnorm2', gnorm2
end program test_easygencan_ab_probe
