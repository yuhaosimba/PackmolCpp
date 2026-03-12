program test_strings
   use sizes
   use input, only : forbidden_char
   use test_assertions
   implicit none

   integer :: strlength
   logical :: empty_char
   character(len=strl) :: alltospace
   external :: parse_spaces

   character(len=strl) :: s
   character(len=strl) :: normalized

   s = "abc" // achar(9) // "   "
   call assert_equal_int(strlength(s), 3, "strlength should ignore tabs and trailing spaces")
   call assert_equal_logical(empty_char(" "), .true., "space should be empty")
   call assert_equal_logical(empty_char(achar(9)), .true., "tab should be empty")
   call assert_equal_logical(empty_char("x"), .false., "letter should not be empty")

   s = "a" // achar(9) // "b"
   normalized = alltospace(s)
   call assert_char_equal(adjustl(normalized(1:3)), "a b", "alltospace should normalize whitespace")

   s = 'keyword "has spaces" tail'
   call parse_spaces(s)
   call assert_true(index(s, forbidden_char) > 0, "parse_spaces should protect quoted spaces")

   s = "escaped\ path"
   call parse_spaces(s)
   call assert_true(index(s, forbidden_char) == 0, "parse_spaces currently leaves escaped spaces untouched")
end program test_strings
