program test_geometry
   use test_assertions
   implicit none

   double precision :: v1(3), v2(3), v3(3), x(3), xcm(3), xref(3)
   external :: eulerrmat, compcart, eulerfixed

   call eulerrmat(0.d0, 0.d0, 0.d0, v1, v2, v3)
   call assert_close_real8(v1(1), 1.d0, 1.d-12, "eulerrmat identity v1")
   call assert_close_real8(v2(2), 1.d0, 1.d-12, "eulerrmat identity v2")
   call assert_close_real8(v3(3), 1.d0, 1.d-12, "eulerrmat identity v3")

   xcm = (/ 1.d0, 2.d0, 3.d0 /)
   xref = (/ 4.d0, 5.d0, 6.d0 /)
   call compcart(x, xcm, xref, v1, v2, v3)
   call assert_close_real8(x(1), 5.d0, 1.d-12, "compcart x")
   call assert_close_real8(x(2), 7.d0, 1.d-12, "compcart y")
   call assert_close_real8(x(3), 9.d0, 1.d-12, "compcart z")

   call eulerfixed(0.d0, 0.d0, 0.d0, v1, v2, v3)
   call assert_close_real8(v1(1), 1.d0, 1.d-12, "eulerfixed identity v1")
   call assert_close_real8(v2(2), 1.d0, 1.d-12, "eulerfixed identity v2")
   call assert_close_real8(v3(3), 1.d0, 1.d-12, "eulerfixed identity v3")
end program test_geometry
