program test_output_messages_probe
   use compute_data
   use input
   implicit none

   external :: title, printbar, writesuccess

   ntype = 1
   allocate(input_itype(1))
   input_itype(1) = 7

   call title()
   call printbar(2, 0, 3)
   call printbar(2, 1, 3)
   call printbar(2, 3, 3)
   call writesuccess(1, 0.125d0, 0.25d0, 0.5d0)
   call writesuccess(2, 0.125d0, 0.25d0, 0.5d0)
end program test_output_messages_probe
