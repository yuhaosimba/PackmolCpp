module test_runtime_setup
   use sizes
   use compute_data
   use input
   use swaptypemod
   implicit none
   private
   public :: load_runtime_problem

contains

subroutine load_runtime_problem(input_path, n, x)
   character(len=*), intent(in) :: input_path
   integer, intent(out) :: n
   double precision, allocatable, intent(out) :: x(:)

   integer :: i, iline, icart, irest, itype, imol, iatom
   external :: setsizes, getinp, cenmass

   input_file_name = trim(input_path)
   call setsizes()
   call getinp()
   call cenmass()

   fix = .false.
   ntype_with_fixed = ntype
   nfixedat = 0
   do i = 1, ntype
      input_itype(i) = i
      fixedoninput(i) = .false.
   end do

   ntotmol = 0
   do itype = 1, ntype
      ntotmol = ntotmol + nmols(itype)
   end do
   n = ntotmol * 6

   if (.not. allocated(xfull)) allocate(xfull(nn))
   allocate(x(n))
   x = 0.d0

   do i = 1, ntotat
      radius(i) = dism / 2.d0
      fscale(i) = 1.d0
      if (use_short_tol) then
         use_short_radius(i) = .true.
      else
         use_short_radius(i) = .false.
      end if
      short_radius(i) = short_tol_dist / 2.d0
      short_radius_scale(i) = short_tol_scale
   end do

   radmax = 0.d0
   do i = 1, ntotat
      radmax = dmax1(radmax, 2.d0 * radius(i))
   end do

   icart = 0
   do itype = 1, ntype
      do imol = 1, nmols(itype)
         do iatom = 1, natoms(itype)
            icart = icart + 1
            nratom(icart) = 0
            iratom(icart, :) = 0
            do iline = linestrut(itype, 1), linestrut(itype, 2)
               if (keyword(iline, 1) == 'inside' .or. &
                   keyword(iline, 1) == 'outside' .or. &
                   keyword(iline, 1) == 'over' .or. &
                   keyword(iline, 1) == 'above' .or. &
                   keyword(iline, 1) == 'below') then
                  nratom(icart) = nratom(icart) + 1
                  do irest = 1, nrest
                     if (irestline(irest) == iline) then
                        iratom(icart, nratom(icart)) = irest
                     end if
                  end do
               end if
            end do
         end do
      end do
   end do

   constrain_rot = .false.
   rot_bound = 0.d0
end subroutine load_runtime_problem

end module test_runtime_setup
