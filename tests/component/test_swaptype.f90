program test_swaptype
   use sizes
   use compute_data
   use input
   use swaptypemod
   use test_assertions
   implicit none

   integer :: n
   double precision :: x(12)
   external :: swaptype

   nn = 12
   n = 12
   ntype = 2
   ntotmol = 2
   nloop_all = 50
   nloop = nloop_all

   allocate(nmols(2), comptype(2), nloop_type(2), xfull(nn))
   nmols = (/ 1, 1 /)
   comptype = .true.
   nloop_type = (/ 7, 9 /)

   x = (/ 1.d0, 2.d0, 3.d0, 4.d0, 5.d0, 6.d0, 7.d0, 8.d0, 9.d0, 10.d0, 11.d0, 12.d0 /)

   call swaptype(n, x, 1, 0)
   call swaptype(n, x, 2, 1)
   call assert_equal_int(n, 6, "swaptype action=1 should reduce the variable count to one molecule")
   call assert_equal_int(ntotmol, 1, "swaptype action=1 should reduce ntotmol to one molecule")
   call assert_equal_int(nloop, 9, "swaptype action=1 should load the per-type loop count")
   call assert_close_real8(x(1), 4.d0, 1.d-12, "swaptype action=1 should expose molecule-2 center x")
   call assert_close_real8(x(2), 5.d0, 1.d-12, "swaptype action=1 should expose molecule-2 center y")
   call assert_close_real8(x(3), 6.d0, 1.d-12, "swaptype action=1 should expose molecule-2 center z")
   call assert_close_real8(x(4), 10.d0, 1.d-12, "swaptype action=1 should expose molecule-2 angle beta")
   call assert_close_real8(x(5), 11.d0, 1.d-12, "swaptype action=1 should expose molecule-2 angle gamma")
   call assert_close_real8(x(6), 12.d0, 1.d-12, "swaptype action=1 should expose molecule-2 angle theta")

   x(1:6) = x(1:6) + 100.d0
   call swaptype(n, x, 2, 2)
   call swaptype(n, x, 2, 3)
   call assert_equal_int(n, 12, "swaptype action=3 should restore the full variable count")
   call assert_equal_int(ntotmol, 2, "swaptype action=3 should restore ntotmol")
   call assert_equal_int(nloop, 50, "swaptype action=3 should restore the global nloop")
   call assert_close_real8(x(4), 104.d0, 1.d-12, "swaptype action=2 should save the updated molecule-2 center x")
   call assert_close_real8(x(10), 110.d0, 1.d-12, "swaptype action=2 should save the updated molecule-2 beta angle")

   deallocate(nmols, comptype, nloop_type, xfull)
end program test_swaptype
