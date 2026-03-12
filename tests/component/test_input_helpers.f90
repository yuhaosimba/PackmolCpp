program test_input_helpers
   use sizes
   use input, only : maxcon
   use test_assertions
   implicit none

   integer :: nres
   character(len=1) :: chain
   character(len=strl) :: record, pdbfile
   character(len=64) :: xyzfile
   external :: chainc, clear, setrnum, setcon

   call chainc(0, chain)
   call assert_char_equal(chain, " ", "chainc zero should map to blank")
   call chainc(1, chain)
   call assert_char_equal(chain, "A", "chainc one should map to A")
   call chainc(36, chain)
   call assert_char_equal(chain, "0", "chainc thirty-six should map to 0")
   call chainc(99, chain)
   call assert_char_equal(chain, "#", "chainc overflow should map to #")

   record = "abc"
   call clear(record)
   call assert_equal_int(len_trim(record), 0, "clear should blank the full record")

   pdbfile = "test_setrnum_input.pdb"
   open(10, file=trim(pdbfile), status="replace")
   write(10, "(a)") "ATOM      1  H   MOL A   1      11.104  13.207   2.100"
   write(10, "(a)") "ATOM      2  O   MOL A   1      12.104  13.207   2.100"
   close(10)
   call setrnum(pdbfile, nres)
   call assert_equal_int(nres, 1, "setrnum should keep a single residue grouped")

   open(10, file=trim(pdbfile), status="replace")
   write(10, "(a)") "ATOM      1  H   MOL A   1      11.104  13.207   2.100"
   write(10, "(a)") "ATOM      2  O   MOL A   2      12.104  13.207   2.100"
   close(10)
   call setrnum(pdbfile, nres)
   call assert_equal_int(nres, 2, "setrnum should detect multiple residues")

   allocate(maxcon(2))
   xyzfile = "test_setcon_input.xyz"
   open(20, file=trim(xyzfile), status="replace")
   write(20, "(a)") "2"
   write(20, "(a)") "1 C 0.0 0.0 0.0 2"
   write(20, "(a)") "2 H 0.0 0.0 1.0 1"
   close(20)
   call setcon(xyzfile, 1)
   call assert_equal_int(maxcon(1), 1, "setcon should count one connection for atom 1")
   call assert_equal_int(maxcon(2), 1, "setcon should count one connection for atom 2")
   deallocate(maxcon)
end program test_input_helpers
