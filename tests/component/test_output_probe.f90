program test_output_probe
   use sizes
   use compute_data
   use input
   use test_assertions
   implicit none

   character(len=strl) :: input_path, output_path
   double precision, allocatable :: x(:)
   integer :: argc, i
   external :: setsizes, getinp, output

   argc = command_argument_count()
   call assert_equal_int(argc, 2, "test_output_probe expects input and output paths")
   call get_command_argument(1, input_path)
   call get_command_argument(2, output_path)
   input_file_name = trim(input_path)

   call setsizes()
   call getinp()
   fix = .false.
   ntype_with_fixed = ntype
   do i = 1, ntype
      input_itype(i) = i
   end do
   fixedoninput(1:ntype) = .false.

   allocate(x(nn))
   x = 0.d0
   call output(nn, x, trim(output_path))
end program test_output_probe
