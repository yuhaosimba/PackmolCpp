program test_flashsort
   use test_assertions
   implicit none

   double precision :: values(5)
   integer :: classes(2), ind(5)
   external :: flash1

   values = (/ 5.d0, 2.d0, 3.d0, 1.d0, 4.d0 /)
   call flash1(values, 5, classes, 2, ind)

   call assert_close_real8(values(1), 1.d0, 1.d-12, "flashsort sorted first element")
   call assert_close_real8(values(5), 5.d0, 1.d-12, "flashsort sorted last element")
   call assert_equal_int(ind(1), 4, "flashsort preserves index map")
   call assert_equal_int(ind(5), 1, "flashsort preserves index map for last element")
end program test_flashsort
