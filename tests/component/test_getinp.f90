program test_getinp
   use sizes
   use compute_data
   use input
   use usegencan
   use test_assertions
   implicit none

   character(len=strl) :: path
   integer :: argc
   external :: setsizes, getinp

   argc = command_argument_count()
   call assert_equal_int(argc, 1, "test_getinp expects the input path")
   call get_command_argument(1, path)
   input_file_name = trim(path)

   call setsizes()
   call getinp()

   call assert_close_real8(dism, 2.d0, 1.d-12, "getinp should parse tolerance")
   call assert_equal_int(seed, 1234567, "getinp should keep the default seed")
   call assert_equal_int(maxit, 20, "getinp should keep the default maxit for the minimal fixture")
   call assert_equal_int(ntotmol, 1, "getinp should keep the one-molecule fixture")
   call assert_equal_int(nn, 6, "getinp should preserve the six optimization variables")
   call assert_equal_int(resnumbers(1), 0, "getinp should infer residue numbering for the minimal structure")
end program test_getinp
