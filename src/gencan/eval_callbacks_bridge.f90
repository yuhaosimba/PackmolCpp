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

subroutine packmol_evalnal_fortran_c(n, x, m, lambda, rho, g, flag) bind(C, name="packmol_evalnal_fortran_c")
   use iso_c_binding, only : c_int, c_double
   implicit none

   integer(c_int), intent(in) :: n, m
   real(c_double), intent(inout) :: x(*)
   real(c_double), intent(in) :: lambda(*), rho(*)
   real(c_double), intent(out) :: g(*)
   integer(c_int), intent(out) :: flag

   external :: evalnal

   call evalnal(n, x, m, lambda, rho, g, flag)
end subroutine packmol_evalnal_fortran_c

subroutine packmol_evalnaldiff_fortran_c(n, x, m, lambda, rho, g, sterel, steabs, flag) bind(C, name="packmol_evalnaldiff_fortran_c")
   use iso_c_binding, only : c_int, c_double
   implicit none

   integer(c_int), intent(in) :: n, m
   real(c_double), intent(inout) :: x(*)
   real(c_double), intent(in) :: lambda(*), rho(*)
   real(c_double), intent(out) :: g(*)
   real(c_double), intent(in) :: sterel, steabs
   integer(c_int), intent(out) :: flag

   external :: evalnaldiff

   call evalnaldiff(n, x, m, lambda, rho, g, sterel, steabs, flag)
end subroutine packmol_evalnaldiff_fortran_c
