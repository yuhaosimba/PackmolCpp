program test_initial
   use compute_data
   use input
   use swaptypemod
   use test_assertions
   use test_runtime_setup
   implicit none

   integer :: n
   double precision, allocatable :: x(:)
   character(len=1000) :: input_path
   external :: initial

   call get_command_argument(1, input_path)
   call load_runtime_problem(trim(input_path), n, x)

   call initial(n, x)

   call assert_true(allocated(latomfirst), "initial should allocate cell linked-list storage")
   call assert_true(allocated(empty_cell), "initial should allocate empty-cell flags")
   call assert_true(all(ncells >= 1), "initial should compute at least one cell per axis")
   call assert_true(any(dabs(x) > 0.d0), "initial should build a non-trivial trial point")
   call assert_equal_logical(init1, .false., "initial should return with init1 disabled")
end program test_initial
