subroutine packmol_calcf_fortran_c(nind, ind, x, n, xc, m, lambda, rho, f, inform) bind(C, name="packmol_calcf_fortran_c")
   use iso_c_binding, only : c_int, c_double
   implicit none

   integer(c_int), intent(in) :: nind, n, m
   integer(c_int), intent(in) :: ind(*)
   real(c_double), intent(inout) :: x(*)
   real(c_double), intent(in) :: xc(*), lambda(*), rho(*)
   real(c_double), intent(out) :: f
   integer(c_int), intent(out) :: inform

   external :: calcf

   call calcf(nind, ind, x, n, xc, m, lambda, rho, f, inform)
end subroutine packmol_calcf_fortran_c

subroutine packmol_calcg_fortran_c(nind, ind, x, n, xc, m, lambda, rho, g, inform) bind(C, name="packmol_calcg_fortran_c")
   use iso_c_binding, only : c_int, c_double
   implicit none

   integer(c_int), intent(in) :: nind, n, m
   integer(c_int), intent(in) :: ind(*)
   real(c_double), intent(inout) :: x(*), g(*)
   real(c_double), intent(in) :: xc(*), lambda(*), rho(*)
   integer(c_int), intent(out) :: inform

   external :: calcg

   call calcg(nind, ind, x, n, xc, m, lambda, rho, g, inform)
end subroutine packmol_calcg_fortran_c

subroutine packmol_calcgdiff_fortran_c(nind, ind, x, n, xc, m, lambda, rho, g, sterel, steabs, inform) &
   bind(C, name="packmol_calcgdiff_fortran_c")
   use iso_c_binding, only : c_int, c_double
   implicit none

   integer(c_int), intent(in) :: nind, n, m
   integer(c_int), intent(in) :: ind(*)
   real(c_double), intent(inout) :: x(*), g(*)
   real(c_double), intent(in) :: xc(*), lambda(*), rho(*), sterel, steabs
   integer(c_int), intent(inout) :: inform

   external :: calcgdiff

   call calcgdiff(nind, ind, x, n, xc, m, lambda, rho, g, sterel, steabs, inform)
end subroutine packmol_calcgdiff_fortran_c

subroutine packmol_calchd_fortran_c(nind, ind, x, d, g, n, xc, m, lambda, rho, hd, xtmp, sterel, steabs, inform) &
   bind(C, name="packmol_calchd_fortran_c")
   use iso_c_binding, only : c_int, c_double
   implicit none

   integer(c_int), intent(in) :: nind, n, m
   integer(c_int), intent(in) :: ind(*)
   real(c_double), intent(inout) :: x(*), d(*), g(*), hd(*), xtmp(*)
   real(c_double), intent(in) :: xc(*), lambda(*), rho(*), sterel, steabs
   integer(c_int), intent(inout) :: inform

   external :: calchd

   call calchd(nind, ind, x, d, g, n, xc, m, lambda, rho, hd, xtmp, sterel, steabs, inform)
end subroutine packmol_calchd_fortran_c

subroutine packmol_calchddiff_fortran_c(nind, ind, x, d, g, n, xc, m, lambda, rho, gtype, hd, xtmp, sterel, steabs, inform) &
   bind(C, name="packmol_calchddiff_fortran_c")
   use iso_c_binding, only : c_int, c_double
   implicit none

   integer(c_int), intent(in) :: nind, n, m, gtype
   integer(c_int), intent(in) :: ind(*)
   real(c_double), intent(inout) :: x(*), d(*), g(*), hd(*), xtmp(*)
   real(c_double), intent(in) :: xc(*), lambda(*), rho(*), sterel, steabs
   integer(c_int), intent(out) :: inform

   external :: calchddiff

   call calchddiff(nind, ind, x, d, g, n, xc, m, lambda, rho, gtype, hd, xtmp, sterel, steabs, inform)
end subroutine packmol_calchddiff_fortran_c
