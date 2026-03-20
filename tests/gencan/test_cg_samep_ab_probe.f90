program test_cg_samep_ab_probe
   use iso_c_binding, only : c_int
   use compute_data
   implicit none

   interface
      integer(c_int) function packmol_gencan_impl_mode_c() bind(C, name="packmol_gencan_impl_mode_c")
         import :: c_int
      end function packmol_gencan_impl_mode_c
   end interface

   integer, parameter :: n = 6, m = 1
   integer :: ind(n), i, mode_id, iter, rbdtype, rbdind, inform
   double precision :: x(n), g(n), l(n), u(n), s(n), w(n), y(n), r(n), d(n), sprev(n)
   double precision :: lambda(m), rho(m), snorm2
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
   s = -1.d0
   w = 0.d0
   y = 0.d0
   r = 0.d0
   d = 0.d0
   sprev = 0.d0
   lambda = 0.d0
   rho = 0.d0

   call cg(n, ind, n, x, m, lambda, rho, g, 1.d0, l, u, 1.d-4, 1.d-4, 5, 20, .false., 0, 1, 1, 0, n, s, iter, &
      rbdtype, rbdind, inform, w, y, r, d, sprev, 1.d-6, 1.d-7, 1.d0, 1.d6, 1.d-20, 1.d20, 1.d99)

   snorm2 = 0.d0
   do i = 1, n
      snorm2 = snorm2 + s(i) * s(i)
   end do
   mode_id = packmol_gencan_impl_mode_c()

   write(*,'(A,1X,I0,1X,A,1X,I0,1X,A,1X,I0,1X,A,1X,I0,1X,A,1X,I0,1X,A,1X,ES24.16)') &
      'mode', mode_id, 'inform', inform, 'iter', iter, 'rbdtype', rbdtype, 'rbdind', rbdind, 'snorm2', snorm2

   deallocate(gxcar)
end program test_cg_samep_ab_probe
