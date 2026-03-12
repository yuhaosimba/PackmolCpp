program test_pbc
   use pbc
   use test_assertions
   implicit none

   call assert_close_real8(v_in_box(11.d0, 0.d0, 10.d0), 1.d0, 1.d-12, "v_in_box wraps coordinates")
   call assert_close_real8(delta_vector(9.d0, 1.d0, 10.d0), -2.d0, 1.d-12, "delta_vector minimum image")
   call assert_equal_int(cell_ind(0, 4), 4, "cell_ind wraps left")
   call assert_equal_int(cell_ind(5, 4), 1, "cell_ind wraps right")
end program test_pbc
