subroutine tnls(nind, ind, n, x, m, lambda, rho, l, u, f, g, d, amax, rbdtype, &
   rbdind, nint, next, mininterp, maxextrap, fmin, maxfc, gtype, iprint, fcnt, &
   gcnt, intcnt, exgcnt, exbcnt, inform, xplus, xtmp, xbext, gamma, beta, &
   sigma1, sigma2, sterel, steabs, epsrel, epsabs, infrel, infabs)
   use iso_c_binding, only : c_int, c_double
   implicit none

   integer, intent(in) :: nind, n, m, rbdtype, rbdind, mininterp, maxextrap, maxfc, gtype, iprint
   integer, intent(inout) :: fcnt, gcnt, intcnt, exgcnt, exbcnt
   integer, intent(out) :: inform
   integer, intent(in) :: ind(nind)
   double precision, intent(inout) :: x(n), f, g(n)
   double precision, intent(in) :: lambda(m), rho(m), l(n), u(n), d(n), amax, nint, next
   double precision, intent(in) :: fmin, gamma, beta, sigma1, sigma2, sterel, steabs, epsrel, epsabs, infrel, infabs
   double precision, intent(inout) :: xplus(n), xtmp(n), xbext(n)

   interface
      subroutine packmol_gencan_tnls_bridge(nind, ind, n, x, m, lambda, rho, l, u, f, g, d, amax, rbdtype, &
         rbdind, nint, next, mininterp, maxextrap, fmin, maxfc, gtype, iprint, fcnt, gcnt, intcnt, exgcnt, &
         exbcnt, inform, xplus, xtmp, xbext, gamma, beta, sigma1, sigma2, sterel, steabs, epsrel, epsabs, &
         infrel, infabs) bind(C, name="packmol_gencan_tnls_bridge")
         import :: c_int, c_double
         integer(c_int), intent(in) :: nind, n, m, rbdtype, rbdind, mininterp, maxextrap, maxfc, gtype, iprint
         integer(c_int), intent(inout) :: fcnt, gcnt, intcnt, exgcnt, exbcnt
         integer(c_int), intent(out) :: inform
         integer(c_int), intent(in) :: ind(*)
         real(c_double), intent(inout) :: x(*), f, g(*)
         real(c_double), intent(in) :: lambda(*), rho(*), l(*), u(*), d(*), amax, nint, next
         real(c_double), intent(in) :: fmin, gamma, beta, sigma1, sigma2, sterel, steabs, epsrel, epsabs, infrel, infabs
         real(c_double), intent(inout) :: xplus(*), xtmp(*), xbext(*)
      end subroutine packmol_gencan_tnls_bridge
   end interface

   call packmol_gencan_tnls_bridge(nind, ind, n, x, m, lambda, rho, l, u, f, g, d, amax, rbdtype, &
      rbdind, nint, next, mininterp, maxextrap, fmin, maxfc, gtype, iprint, fcnt, gcnt, intcnt, exgcnt, &
      exbcnt, inform, xplus, xtmp, xbext, gamma, beta, sigma1, sigma2, sterel, steabs, epsrel, epsabs, &
      infrel, infabs)
end subroutine tnls

subroutine packmol_tnls_fortran_c(nind, ind, n, x, m, lambda, rho, l, u, f, g, d, amax, rbdtype, &
   rbdind, nint, next, mininterp, maxextrap, fmin, maxfc, gtype, iprint, fcnt, gcnt, intcnt, exgcnt, &
   exbcnt, inform, xplus, xtmp, xbext, gamma, beta, sigma1, sigma2, sterel, steabs, epsrel, epsabs, &
   infrel, infabs) bind(C, name="packmol_tnls_fortran_c")
   use iso_c_binding, only : c_int, c_double
   implicit none

   integer(c_int), intent(in) :: nind, n, m, rbdtype, rbdind, mininterp, maxextrap, maxfc, gtype, iprint
   integer(c_int), intent(inout) :: fcnt, gcnt, intcnt, exgcnt, exbcnt
   integer(c_int), intent(out) :: inform
   integer(c_int), intent(in) :: ind(*)
   real(c_double), intent(inout) :: x(*), f, g(*)
   real(c_double), intent(in) :: lambda(*), rho(*), l(*), u(*), d(*), amax, nint, next
   real(c_double), intent(in) :: fmin, gamma, beta, sigma1, sigma2, sterel, steabs, epsrel, epsabs, infrel, infabs
   real(c_double), intent(inout) :: xplus(*), xtmp(*), xbext(*)

   external :: tnlsf

   call tnlsf(nind, ind, n, x, m, lambda, rho, l, u, f, g, d, amax, rbdtype, rbdind, nint, next, &
      mininterp, maxextrap, fmin, maxfc, gtype, iprint, fcnt, gcnt, intcnt, exgcnt, exbcnt, inform, xplus, &
      xtmp, xbext, gamma, beta, sigma1, sigma2, sterel, steabs, epsrel, epsabs, infrel, infabs)
end subroutine packmol_tnls_fortran_c
