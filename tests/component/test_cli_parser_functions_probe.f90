program test_cli_parser_functions_probe
   use cli_parser, only : get_filename, parse_command
   use test_assertions
   implicit none

   character(len=1000) :: input_name, output_name

   call assert_equal_int(parse_command("-i"), 0, "parse_command should detect the input flag")
   call assert_equal_int(parse_command("-o"), 1, "parse_command should detect the output flag")

   call get_filename(1, input_name)
   call get_filename(3, output_name)

   call assert_char_equal(input_name, "fixtures/input.inp", "get_filename should read the input filename")
   call assert_char_equal(output_name, "fixtures/output.pdb", "get_filename should read the output filename")
end program test_cli_parser_functions_probe
