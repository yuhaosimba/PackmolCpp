program test_gencan_cpp_bridge
   use gencan_cpp_bridge
   use test_assertions
   implicit none

   integer, parameter :: n = 3
   integer :: inform
   double precision :: x(n), fx

   x = (/ 1.d0, 2.d0, 3.d0 /)
   fx = -1.d0
   inform = -999

   call gencan_cpp_probe(n, x, fx, inform)

   call assert_equal_int(inform, 0, "C++ probe bridge should return success")
   call assert_close_real8(fx, 0.d0, 1.d-12, "C++ probe bridge should set fx_out to zero in Phase B skeleton")
end program test_gencan_cpp_bridge

