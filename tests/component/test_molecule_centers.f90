program test_molecule_centers
   use sizes
   use compute_data
   use input
   use test_assertions
   implicit none

   external :: tobar, cenmass

   maxkeywords = 4
   ntype = 1
   nlines = 0

   allocate(natoms(1), idfirst(1), coor(2,3), amass(2), keyword(1,maxkeywords), linestrut(1,2))
   natoms(1) = 2
   idfirst(1) = 1
   linestrut = 0
   keyword = ""

   coor = 0.d0
   coor(1,1) = 1.d0
   coor(2,1) = 3.d0
   call tobar()
   call assert_close_real8(coor(1,1), -1.d0, 1.d-12, "tobar should center first atom")
   call assert_close_real8(coor(2,1), 1.d0, 1.d-12, "tobar should center second atom")

   coor = 0.d0
   coor(1,1) = 0.d0
   coor(2,1) = 4.d0
   amass(1) = 1.d0
   amass(2) = 3.d0
   call cenmass()
   call assert_close_real8(coor(1,1), -3.d0, 1.d-12, "cenmass should subtract mass center from first atom")
   call assert_close_real8(coor(2,1), 1.d0, 1.d-12, "cenmass should subtract mass center from second atom")

   deallocate(natoms, idfirst, coor, amass, keyword, linestrut)
end program test_molecule_centers
