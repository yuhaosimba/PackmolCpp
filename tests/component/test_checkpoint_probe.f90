program test_checkpoint_probe
   use sizes
   use compute_data
   use input
   use test_assertions
   implicit none

   character(len=strl) :: input_path
   double precision, allocatable :: x(:)
   integer :: argc, i
   external :: setsizes, getinp, checkpoint

   argc = command_argument_count()
   call assert_equal_int(argc, 1, "test_checkpoint_probe expects the input path")
   call get_command_argument(1, input_path)
   input_file_name = trim(input_path)

   call setsizes()
   call getinp()
   fix = .false.
   ntype_with_fixed = ntype
   do i = 1, ntype
      input_itype(i) = i
   end do
   fixedoninput(1:ntype) = .false.
   init1 = .true.
   move = .false.
   nratom = 0
   iratom = 0
   xyzout = "output.pdb"

   allocate(x(nn))
   x = 0.d0
   call checkpoint(nn, x)
end program test_checkpoint_probe
