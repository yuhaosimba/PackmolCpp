program test_setsizes
   use sizes
   use compute_data
   use input
   use usegencan
   use test_assertions
   implicit none

   character(len=strl) :: path
   integer :: argc
   external :: setsizes

   argc = command_argument_count()
   call assert_equal_int(argc, 1, "test_setsizes expects the input path")
   call get_command_argument(1, path)
   input_file_name = trim(path)

   call setsizes()

   call assert_equal_int(ntype, 1, "setsizes should find one structure")
   call assert_equal_int(natoms(1), 2, "setsizes should read two atoms from the fixture structure")
   call assert_equal_int(nmols(1), 1, "setsizes should read one molecule")
   call assert_equal_int(ntotat, 2, "setsizes should compute total atoms")
   call assert_equal_int(ntotmol, 1, "setsizes should compute total molecules")
   call assert_equal_int(nn, 6, "setsizes should set six optimization variables for one molecule")
   call assert_true(allocated(keyword), "setsizes should allocate keyword table")
   call assert_true(allocated(xcart), "setsizes should allocate core coordinate arrays")
   call assert_true(allocated(l), "setsizes should allocate gencan bounds")
end program test_setsizes
