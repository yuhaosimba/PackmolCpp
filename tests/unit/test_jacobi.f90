program test_jacobi
   use test_assertions
   implicit none

   integer, parameter :: n = 2, np = 2
   integer :: nrot
   double precision :: a(np, n), d(n), v(np, n)
   external :: jacobi

   a = reshape((/ 2.d0, 1.d0, 1.d0, 2.d0 /), shape(a))
   nrot = 10
   call jacobi(a, n, np, d, v, nrot)

   call assert_close_real8(d(1), 1.d0, 1.d-10, "jacobi first eigenvalue")
   call assert_close_real8(d(2), 3.d0, 1.d-10, "jacobi second eigenvalue")
   call assert_close_real8(abs(v(1, 1)), 1.d0 / sqrt(2.d0), 1.d-10, "jacobi eigenvector normalization")
   call assert_true(nrot >= 1, "jacobi should perform at least one rotation")
end program test_jacobi
