program test_failopen_probe
   use sizes
   implicit none
   external :: failopen
   character(len=strl) :: record

   record = "missing-file.pdb"
   call failopen(record)
end program test_failopen_probe
