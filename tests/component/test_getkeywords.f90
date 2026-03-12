program test_getkeywords
   use sizes
   use input, only : nlines, inputfile, keyword, forbidden_char
   use test_assertions
   implicit none

   external :: getkeywords
   character(len=strl) :: quoted_line

   maxkeywords = 4
   nlines = 2
   allocate(inputfile(nlines), keyword(nlines, maxkeywords))

   inputfile(1) = "output output.pdb"
   quoted_line = 'structure "foo'
   quoted_line = trim(quoted_line) // forbidden_char // 'bar"'
   inputfile(2) = trim(quoted_line)

   call getkeywords()

   call assert_char_equal(keyword(1,1), "output", "getkeywords first keyword")
   call assert_char_equal(keyword(1,2), "output.pdb", "getkeywords first value")
   call assert_char_equal(keyword(2,1), "structure", "getkeywords quoted keyword")
   call assert_char_equal(keyword(2,2), "foo bar", "getkeywords should restore forbidden_char as space")

   deallocate(inputfile, keyword)
end program test_getkeywords
