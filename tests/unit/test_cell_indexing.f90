program test_cell_indexing
   use compute_data, only : cell_length, ncells
   use pbc, only : pbc_min, pbc_length
   use cell_indexing
   use test_assertions
   implicit none

   integer :: cell(3), back(3)
   double precision :: x(3)

   ncells = (/ 4, 5, 6 /)
   cell_length = (/ 2.d0, 2.d0, 2.d0 /)
   pbc_min = (/ 0.d0, 0.d0, 0.d0 /)
   pbc_length = (/ 8.d0, 10.d0, 12.d0 /)

   x = (/ 3.5d0, 4.1d0, 5.9d0 /)
   call setcell(x, cell)
   call assert_equal_int(cell(1), 2, "setcell x")
   call assert_equal_int(cell(2), 3, "setcell y")
   call assert_equal_int(cell(3), 3, "setcell z")

   call icell_to_cell(index_cell(cell, ncells), ncells, back)
   call assert_equal_int(back(1), cell(1), "icell_to_cell x")
   call assert_equal_int(back(2), cell(2), "icell_to_cell y")
   call assert_equal_int(back(3), cell(3), "icell_to_cell z")
end program test_cell_indexing
