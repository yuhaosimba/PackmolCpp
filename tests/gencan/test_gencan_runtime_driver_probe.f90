program test_gencan_runtime_driver_probe
   use iso_c_binding, only : c_int
   use compute_data
   use test_assertions
   use test_runtime_setup
   implicit none

   interface
      integer(c_int) function packmol_gencan_impl_mode_c() bind(C, name="packmol_gencan_impl_mode_c")
         import :: c_int
      end function packmol_gencan_impl_mode_c
   end interface

   integer, parameter :: maxitngp = 5
   integer :: n, m, i, mode_id
   integer :: inform, fcnt, gcnt, iter, cgcnt, spgiter, spgfcnt
   integer :: tniter, tnfcnt, tnstpcnt, tnintcnt, tnexgcnt, tnexbcnt
   integer :: tnintfe, tnexgfe, tnexbfe
   integer, allocatable :: ind(:)
   double precision :: f, gpsupn, gpeucn2, xsum, gnorm2
   double precision, allocatable :: x(:), g(:), l(:), u(:), d(:), s(:), y(:), w(:)
   double precision, allocatable :: lambda(:), rho(:), lastgpns(:)
   character(len=1000) :: input_path
   external :: initial, gencan

   call get_command_argument(1, input_path)
   call load_runtime_problem(trim(input_path), n, x)
   call initial(n, x)

   m = max(1, ntype)
   allocate(ind(n), g(n), l(n), u(n), d(n), s(n), y(n), w(5 * n), &
      lambda(m), rho(m), lastgpns(0:maxitngp-1))

   do i = 1, n
      ind(i) = i
   end do
   g = 0.d0
   l = -1.d20
   u = 1.d20
   d = 0.d0
   s = 0.d0
   y = 0.d0
   w = 0.d0
   lambda = 0.d0
   rho = 0.d0
   f = -1.d0
   gpsupn = -1.d0
   gpeucn2 = -1.d0
   lastgpns = 0.d0
   iter = -1
   fcnt = -1
   gcnt = -1
   cgcnt = -1
   spgiter = -1
   spgfcnt = -1
   tniter = -1
   tnfcnt = -1
   tnstpcnt = -1
   tnintcnt = -1
   tnexgcnt = -1
   tnexbcnt = -1
   tnintfe = -1
   tnexgfe = -1
   tnexbfe = -1
   inform = -999

   call gencan(n, x, l, u, m, lambda, rho, 0.d0, 1.d-6, 5, 0.d0, maxitngp, &
      -1.d0, 20, 500, -1.d0, -1, 2, 1.d-6, 1.d-1, 1.d-5, 1.d-4, 5, .false., &
      2.d0, 2.d0, 4, 10, 0, 1, 1, 0, 20, f, g, gpeucn2, gpsupn, iter, fcnt, &
      gcnt, cgcnt, spgiter, spgfcnt, tniter, tnfcnt, tnstpcnt, tnintcnt, &
      tnexgcnt, tnexbcnt, tnintfe, tnexgfe, tnexbfe, inform, s, y, d, ind, &
      lastgpns, w, 0.9d0, 1.d-2, 1.d10, 1.d-10, 1.d-6, 1.d-4, 0.5d0, 0.1d0, &
      0.9d0, 1.d-7, 1.d-10, 1.d-10, 1.d-20, 1.d20, 1.d99)

   mode_id = packmol_gencan_impl_mode_c()
   xsum = sum(x)
   gnorm2 = sum(g * g)

   write(*,'(A,1X,I0,1X,A,1X,I0,1X,A,1X,ES24.16,1X,A,1X,ES24.16,1X,A,1X,ES24.16,1X,A,1X,I0,1X,A,1X,I0,1X,A,1X,I0,1X,A,1X,I0)') &
      'mode', mode_id, 'inform', inform, 'f', f, 'gpsupn', gpsupn, 'xsum', xsum, &
      'iter', iter, 'fcnt', fcnt, 'gcnt', gcnt, 'cgcnt', cgcnt
   write(*,'(A,1X,ES24.16)') 'gnorm2', gnorm2

   call assert_true(inform >= 0, "runtime gencan driver probe should finish without internal error")
   call assert_true(f == f, "runtime gencan driver probe should return finite objective")
   call assert_true(gpsupn == gpsupn, "runtime gencan driver probe should return finite gpsupn")
end program test_gencan_runtime_driver_probe
