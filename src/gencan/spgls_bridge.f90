subroutine spgls(n, x, m, lambda, rho, f, g, l, u, lamspg, nint, mininterp, &
   fmin, maxfc, iprint, fcnt, inform, xtrial, d, gamma, sigma1, sigma2, &
   sterel, steabs, epsrel, epsabs, infrel, infabs)
   use iso_c_binding, only : c_int, c_double
   implicit none

   integer, intent(in) :: n, m, mininterp, maxfc, iprint
   integer, intent(inout) :: fcnt
   integer, intent(out) :: inform
   double precision, intent(inout) :: x(n), f
   double precision, intent(in) :: lambda(m), rho(m), g(n), l(n), u(n)
   double precision, intent(in) :: lamspg, nint, fmin, gamma, sigma1, sigma2
   double precision, intent(in) :: sterel, steabs, epsrel, epsabs, infrel, infabs
   double precision, intent(inout) :: xtrial(n), d(n)

   interface
      subroutine packmol_gencan_spgls_bridge(n, x, m, lambda, rho, f, g, l, u, lamspg, nint, mininterp, &
         fmin, maxfc, iprint, fcnt, inform, xtrial, d, gamma, sigma1, sigma2, sterel, steabs, epsrel, epsabs, &
         infrel, infabs) bind(C, name="packmol_gencan_spgls_bridge")
         import :: c_int, c_double
         integer(c_int), intent(in) :: n, m, mininterp, maxfc, iprint
         integer(c_int), intent(inout) :: fcnt
         integer(c_int), intent(out) :: inform
         real(c_double), intent(inout) :: x(*), f
         real(c_double), intent(in) :: lambda(*), rho(*), g(*), l(*), u(*)
         real(c_double), intent(in) :: lamspg, nint, fmin, gamma, sigma1, sigma2
         real(c_double), intent(in) :: sterel, steabs, epsrel, epsabs, infrel, infabs
         real(c_double), intent(inout) :: xtrial(*), d(*)
      end subroutine packmol_gencan_spgls_bridge
   end interface

   call packmol_gencan_spgls_bridge(n, x, m, lambda, rho, f, g, l, u, lamspg, nint, mininterp, &
      fmin, maxfc, iprint, fcnt, inform, xtrial, d, gamma, sigma1, sigma2, sterel, steabs, epsrel, epsabs, &
      infrel, infabs)
end subroutine spgls

subroutine packmol_spgls_fortran_c(n, x, m, lambda, rho, f, g, l, u, lamspg, nint, mininterp, &
   fmin, maxfc, iprint, fcnt, inform, xtrial, d, gamma, sigma1, sigma2, sterel, steabs, epsrel, epsabs, &
   infrel, infabs) bind(C, name="packmol_spgls_fortran_c")
   use iso_c_binding, only : c_int, c_double
   implicit none

   integer(c_int), intent(in) :: n, m, mininterp, maxfc, iprint
   integer(c_int), intent(inout) :: fcnt
   integer(c_int), intent(out) :: inform
   real(c_double), intent(inout) :: x(*), f
   real(c_double), intent(in) :: lambda(*), rho(*), g(*), l(*), u(*)
   real(c_double), intent(in) :: lamspg, nint, fmin, gamma, sigma1, sigma2
   real(c_double), intent(in) :: sterel, steabs, epsrel, epsabs, infrel, infabs
   real(c_double), intent(inout) :: xtrial(*), d(*)

   external :: spglsf

   call spglsf(n, x, m, lambda, rho, f, g, l, u, lamspg, nint, mininterp, &
      fmin, maxfc, iprint, fcnt, inform, xtrial, d, gamma, sigma1, sigma2, sterel, steabs, epsrel, epsabs, &
      infrel, infabs)
end subroutine packmol_spgls_fortran_c
