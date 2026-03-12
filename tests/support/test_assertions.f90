module test_assertions
   use, intrinsic :: iso_fortran_env, only : error_unit
   implicit none
   private
   public :: assert_true, assert_equal_int, assert_equal_logical, assert_close_real8
   public :: assert_char_equal

contains

subroutine fail_assert(message)
   character(len=*), intent(in) :: message
   write(error_unit, "(a)") trim(message)
   error stop 1
end subroutine fail_assert

subroutine assert_true(condition, message)
   logical, intent(in) :: condition
   character(len=*), intent(in) :: message
   if (.not. condition) call fail_assert(message)
end subroutine assert_true

subroutine assert_equal_int(actual, expected, message)
   integer, intent(in) :: actual, expected
   character(len=*), intent(in) :: message
   if (actual /= expected) call fail_assert(message)
end subroutine assert_equal_int

subroutine assert_equal_logical(actual, expected, message)
   logical, intent(in) :: actual, expected
   character(len=*), intent(in) :: message
   if (actual .neqv. expected) call fail_assert(message)
end subroutine assert_equal_logical

subroutine assert_close_real8(actual, expected, tolerance, message)
   double precision, intent(in) :: actual, expected, tolerance
   character(len=*), intent(in) :: message
   if (abs(actual - expected) > tolerance) call fail_assert(message)
end subroutine assert_close_real8

subroutine assert_char_equal(actual, expected, message)
   character(len=*), intent(in) :: actual, expected
   character(len=*), intent(in) :: message
   if (trim(actual) /= trim(expected)) call fail_assert(message)
end subroutine assert_char_equal

end module test_assertions
