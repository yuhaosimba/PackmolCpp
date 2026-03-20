program test_spgls_ab_probe
   use iso_c_binding, only : c_int
   use compute_data
   implicit none

   interface
      integer(c_int) function packmol_gencan_impl_mode_c() bind(C, name="packmol_gencan_impl_mode_c")
         import :: c_int
      end function packmol_gencan_impl_mode_c
   end interface

   integer, parameter :: n = 6, m = 1
   integer :: i, mode_id, fcnt, inform
   double precision :: x(n), g(n), l(n), u(n), d(n), xtrial(n)
   double precision :: lambda(m), rho(m), f, lamspg, xsum, dnorm2
   external :: spgls

   ntotmol = 0
   ntype = 0
   init1 = .true.
   move = .false.
   allocate(gxcar(1, 3))
   gxcar = 0.d0

   x = (/ 1.d0, 2.d0, 3.d0, 4.d0, 5.d0, 6.d0 /)
   g = (/ 3.d0, -2.d0, 1.d0, 0.5d0, -4.d0, 2.5d0 /)
   l = -1.d20
   u = 1.d20
   d = 0.d0
   xtrial = 0.d0
   lambda = 0.d0
   rho = 0.d0
   f = 0.d0
   fcnt = 0
   inform = -999
   lamspg = 1.d0

   call spgls(n, x, m, lambda, rho, f, g, l, u, lamspg, 2.d0, 1, -1.d0, 1, 0, fcnt, inform, xtrial, d, &
      1.d-4, 0.1d0, 0.9d0, 1.d-7, 1.d-10, 1.d-10, 1.d-20, 1.d20, 1.d99)

   xsum = 0.d0
   dnorm2 = 0.d0
   do i = 1, n
      xsum = xsum + x(i)
      dnorm2 = dnorm2 + d(i) * d(i)
   end do
   mode_id = packmol_gencan_impl_mode_c()

   write(*,'(A,1X,I0,1X,A,1X,I0,1X,A,1X,I0,1X,A,1X,ES24.16,1X,A,1X,ES24.16,1X,A,1X,ES24.16)') &
      'mode', mode_id, 'inform', inform, 'fcnt', fcnt, 'f', f, 'xsum', xsum, 'dnorm2', dnorm2

   deallocate(gxcar)
end program test_spgls_ab_probe
