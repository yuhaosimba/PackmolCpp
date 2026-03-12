program test_movebad
   use compute_data
   use input
   use test_assertions
   use test_runtime_setup
   implicit none

   integer :: n, i
   double precision :: fx
   double precision, allocatable :: x(:), before(:)
   character(len=1000) :: input_path
   logical :: movebadprint
   external :: initial, movebad, computef, init_random_number

   call get_command_argument(1, input_path)
   call load_runtime_problem(trim(input_path), n, x)
   call initial(n, x)

   do i = 1, ntotat
      radius_ini(i) = radius(i)
   end do

   x(1:3) = (/ 1.d0, 1.d0, 1.d0 /)
   x(4:6) = (/ 1.d0, 1.d0, 1.d0 /)
   x(7:12) = 0.d0

   movefrac = 1.d0
   movebadrandom = .true.
   maxmove(1) = nmols(1)
   call init_random_number(12345)

   allocate(before(n))
   before = x
   call computef(n, x, fx)
   movebadprint = .false.
   call movebad(n, x, fx, movebadprint)

   call assert_true(any(dabs(x(1:6) - before(1:6)) > 1.d-12), "movebad should modify at least one molecule center")
   call assert_true(fx == fx, "movebad should keep fx finite")
end program test_movebad
