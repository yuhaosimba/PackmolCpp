program test_output_helpers
   use input, only : hexadecimal_indices
   use test_assertions
   implicit none

   character(len=4) :: i4hex
   character(len=5) :: i5hex

   hexadecimal_indices = .false.
   call assert_char_equal(i4hex(42), "  42", "i4hex decimal output")
   call assert_char_equal(i5hex(42), "   42", "i5hex decimal output")

   hexadecimal_indices = .true.
   call assert_char_equal(adjustl(i4hex(255)), "FF", "i4hex hex output")
   call assert_char_equal(adjustl(i5hex(255)), "FF", "i5hex hex output")
end program test_output_helpers
