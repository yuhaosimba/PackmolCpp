program test_tnls_interior_beta_extrap_ab_probe
   use iso_c_binding, only : c_int
   use compute_data
   implicit none

   interface
      integer(c_int) function packmol_gencan_impl_mode_c() bind(C, name="packmol_gencan_impl_mode_c")
         import :: c_int
      end function packmol_gencan_impl_mode_c
   end interface

   integer, parameter :: n = 6, m = 1
   integer :: ind(n), i, mode_id
   integer :: fcnt, gcnt, intcnt, exgcnt, exbcnt, inform
   double precision :: x(n), g(n), l(n), u(n), d(n), xplus(n), xtmp(n), xbext(n)
   double precision :: lambda(m), rho(m), f, amax, xsum, gnorm2
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
   d = -g
   l = -1.d20
   u = 1.d20
   xplus = 0.d0
   xtmp = 0.d0
   xbext = 0.d0
   lambda = 0.d0
   rho = 0.d0
   f = 1.d0
   amax = 2.d0
   fcnt = 0
   gcnt = 0
   intcnt = 0
   exgcnt = 0
   exbcnt = 0
   inform = -999

   call tnls(n, ind, n, x, m, lambda, rho, l, u, f, g, d, amax, 1, 1, 2.d0, 2.d0, 1, 2, -1.d0, 20, 0, 0, &
      fcnt, gcnt, intcnt, exgcnt, exbcnt, inform, xplus, xtmp, xbext, 1.d-4, -1.d0, 0.1d0, 0.9d0, &
      1.d-7, 1.d-10, 1.d-10, 1.d-20, 1.d20, 1.d99)

   xsum = 0.d0
   gnorm2 = 0.d0
   do i = 1, n
      xsum = xsum + x(i)
      gnorm2 = gnorm2 + g(i) * g(i)
   end do
   mode_id = packmol_gencan_impl_mode_c()

   write(*,'(A,1X,I0,1X,A,1X,I0,1X,A,1X,I0,1X,A,1X,I0,1X,A,1X,I0,1X,A,1X,I0,1X,A,1X,I0,1X,A,1X,ES24.16,1X,A,1X,ES24.16,1X,A,1X,ES24.16)') &
      'mode', mode_id, 'inform', inform, 'fcnt', fcnt, 'gcnt', gcnt, 'intcnt', intcnt, 'exgcnt', exgcnt, &
      'exbcnt', exbcnt, 'f', f, 'xsum', xsum, 'gnorm2', gnorm2

   deallocate(gxcar)
end program test_tnls_interior_beta_extrap_ab_probe
