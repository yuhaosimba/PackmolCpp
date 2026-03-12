program test_resetcells
   use compute_data
   use input, only : fix, ntype_with_fixed
   use pbc, only : pbc_min, pbc_length
   use test_assertions
   implicit none

   external :: resetcells

   ncells = (/ 2, 2, 2 /)
   ntype = 0
   ntype_with_fixed = 1
   ntotat = 1
   nfixedat = 1
   fix = .true.
   pbc_min = 0.d0
   pbc_length = (/ 2.d0, 2.d0, 2.d0 /)
   cell_length = (/ 1.d0, 1.d0, 1.d0 /)

   allocate(natoms(1), latomfirst(2,2,2), lcellnext(8), empty_cell(2,2,2), latomnext(1), xcart(1,3))
   natoms(1) = 1
   latomfirst = -7
   lcellnext = -3
   latomnext = -5
   empty_cell = .false.
   lcellfirst = -9
   xcart(1,:) = (/ 0.25d0, 0.25d0, 0.25d0 /)

   call resetcells()

   call assert_equal_int(lcellfirst, 1, "resetcells should register the occupied cell")
   call assert_equal_int(latomfirst(1,1,1), 1, "resetcells should insert the fixed atom")
   call assert_equal_int(latomnext(1), 0, "resetcells should terminate the linked list")
   call assert_equal_logical(empty_cell(1,1,1), .false., "resetcells should mark the occupied cell")
   call assert_equal_logical(empty_cell(2,2,2), .true., "resetcells should clear unrelated cells")

   deallocate(natoms, latomfirst, lcellnext, empty_cell, latomnext, xcart)
end program test_resetcells
