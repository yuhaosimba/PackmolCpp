program test_write_connect
   use sizes
   use input
   use test_assertions
   implicit none

   character(len=strl) :: line
   integer :: ioerr
   external :: write_connect

   allocate(maxcon(3), nconnect(3, 3))
   maxcon = 0
   nconnect = 0
   hexadecimal_indices = .false.

   maxcon(1) = 2
   nconnect(1, 1:2) = (/ 2, 3 /)

   open(10, file="test_write_connect.tmp", status="replace")
   call write_connect(10, 0, 1, 1)
   close(10)

   open(11, file="test_write_connect.tmp", status="old")
   read(11, "(a)", iostat=ioerr) line
   close(11, status="delete")

   call assert_equal_int(ioerr, 0, "write_connect should produce one output line")
   call assert_true(index(line, "CONECT") > 0, "write_connect should emit a CONECT record")
   call assert_true(index(line, "2") > 0, "write_connect should include the first neighbour")
   call assert_true(index(line, "3") > 0, "write_connect should include the second neighbour")

   deallocate(maxcon, nconnect)
end program test_write_connect
