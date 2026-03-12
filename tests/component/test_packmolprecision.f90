program test_packmolprecision
   use compute_data
   use input
   use test_assertions
   implicit none

   integer, parameter :: n = 4
   double precision :: x(n)
   logical, external :: packmolprecision

   ntotmol = 0
   ntype = 0
   init1 = .true.
   move = .false.

   x = (/ 1.d0, 2.d0, 3.d0, 4.d0 /)

   precision = 1.d-6
   call assert_true(packmolprecision(n, x), "packmolprecision should accept the empty toy problem")

   precision = -1.d0
   call assert_true(.not. packmolprecision(n, x), "packmolprecision should reject when precision is negative")
end program test_packmolprecision
