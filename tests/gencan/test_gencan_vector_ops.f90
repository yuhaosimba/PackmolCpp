program test_gencan_vector_ops
   use test_assertions
   implicit none

   integer :: ind(3)
   double precision :: values(5)
   external :: shrink, expand, evalhd

   ind = (/ 1, 3, 5 /)
   values = (/ 10.d0, 20.d0, 30.d0, 40.d0, 50.d0 /)
   call shrink(3, ind, 5, values)
   call assert_close_real8(values(1), 10.d0, 1.d-12, "shrink first value")
   call assert_close_real8(values(2), 30.d0, 1.d-12, "shrink second value")
   call assert_close_real8(values(3), 50.d0, 1.d-12, "shrink third value")

   call expand(3, ind, 5, values)
   call assert_close_real8(values(1), 10.d0, 1.d-12, "expand restores first value")
   call assert_close_real8(values(3), 30.d0, 1.d-12, "expand restores third value")
   call assert_close_real8(values(5), 50.d0, 1.d-12, "expand restores fifth value")

   call evalhd(5)
end program test_gencan_vector_ops
