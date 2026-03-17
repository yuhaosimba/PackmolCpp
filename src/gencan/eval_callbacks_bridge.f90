subroutine packmol_evalal_fortran_c(n, x, m, lambda, rho, f, flag) bind(C, name="packmol_evalal_fortran_c")
   use iso_c_binding, only : c_int, c_double
   implicit none

   integer(c_int), intent(in) :: n, m
   real(c_double), intent(in) :: x(*), lambda(*), rho(*)
   real(c_double), intent(out) :: f
   integer(c_int), intent(out) :: flag

   external :: evalal

   call evalal(n, x, m, lambda, rho, f, flag)
end subroutine packmol_evalal_fortran_c
