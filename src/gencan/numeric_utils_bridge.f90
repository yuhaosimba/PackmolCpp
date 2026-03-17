subroutine packmol_norm2s_fortran_c(n, x, val) bind(C, name="packmol_norm2s_fortran_c")
   use iso_c_binding, only : c_int, c_double
   implicit none

   integer(c_int), intent(in) :: n
   real(c_double), intent(in) :: x(*)
   real(c_double), intent(out) :: val

   external :: norm2s
   double precision :: norm2s

   val = norm2s(n, x)
end subroutine packmol_norm2s_fortran_c

subroutine packmol_dot_fortran_c(n, x, y, val) bind(C, name="packmol_dot_fortran_c")
   use iso_c_binding, only : c_int, c_double
   implicit none

   integer(c_int), intent(in) :: n
   real(c_double), intent(in) :: x(*), y(*)
   real(c_double), intent(out) :: val
   integer :: i

   val = 0.0d0
   do i = 1, n
      val = val + x(i) * y(i)
   end do
end subroutine packmol_dot_fortran_c
