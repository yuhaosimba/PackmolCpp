module gencan_cpp_bridge
   use iso_c_binding
   implicit none
   private

   public :: gencan_cpp_probe

   interface
      integer(c_int) function packmol_gencan_cpp_probe_c(n, x, fx_out) bind(C, name="packmol_gencan_cpp_probe")
         import :: c_int, c_double
         integer(c_int), value :: n
         real(c_double), intent(in) :: x(*)
         real(c_double), intent(out) :: fx_out
      end function packmol_gencan_cpp_probe_c
   end interface

contains

subroutine gencan_cpp_probe(n, x, fx_out, inform)
   integer, intent(in) :: n
   double precision, intent(in) :: x(n)
   double precision, intent(out) :: fx_out
   integer, intent(out) :: inform

   inform = packmol_gencan_cpp_probe_c(n, x, fx_out)
end subroutine gencan_cpp_probe

end module gencan_cpp_bridge

