program test_random
   use test_assertions
   implicit none

   double precision :: rnd, a, b, c, d
   external :: init_random_number, seed_from_time
   integer :: seed

   call init_random_number(17)
   a = rnd()
   b = rnd()
   call init_random_number(17)
   c = rnd()
   d = rnd()

   call assert_close_real8(a, c, 1.d-12, "seeded random sequence first value")
   call assert_close_real8(b, d, 1.d-12, "seeded random sequence second value")

   call seed_from_time(seed)
   call assert_true(seed > 0, "seed_from_time should return a positive seed")
end program test_random
