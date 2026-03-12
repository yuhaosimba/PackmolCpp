program test_discret
   use compute_data
   use test_assertions
   implicit none

   integer, parameter :: n = 4
   double precision :: x(n), gcomp
   external :: discret

   ntotmol = 0
   ntype = 0
   init1 = .true.
   move = .false.

   x = (/ 1.d0, 2.d0, 3.d0, 4.d0 /)
   call discret(2, n, x, gcomp, 1.d-4)

   call assert_close_real8(gcomp, 0.d0, 1.d-12, "discret should be zero for the empty toy problem")
   call assert_close_real8(x(1), 1.d0, 1.d-12, "discret should preserve x(1)")
   call assert_close_real8(x(2), 2.d0, 1.d-12, "discret should preserve x(2)")
   call assert_close_real8(x(3), 3.d0, 1.d-12, "discret should preserve x(3)")
   call assert_close_real8(x(4), 4.d0, 1.d-12, "discret should preserve x(4)")
end program test_discret
