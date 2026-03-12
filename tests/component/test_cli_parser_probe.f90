program test_cli_parser_probe
   use cli_parser
   use input, only : input_file_name, xyzout
   use test_assertions
   implicit none

   call parse_command_line_args()
   call assert_true(len_trim(input_file_name) > 0, "input file should be captured")
   call assert_true(trim(input_file_name) == "fixtures/input.inp", "input file path should match")

   if (command_argument_count() == 4) then
      call assert_true(trim(xyzout) == "fixtures/output.pdb", "output file path should match")
   else
      call assert_true(len_trim(xyzout) == 0, "output file should stay empty when omitted")
   end if
end program test_cli_parser_probe
