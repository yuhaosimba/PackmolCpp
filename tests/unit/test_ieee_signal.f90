program test_ieee_signal
   use test_assertions
   implicit none

   double precision :: gpsupn, acgeps, bcgeps, cgepsf, cgepsi, cggpnf
   double precision :: kappa, gpeucn2, gpeucn20, epsgpen2, epsgpsn, cgeps, gpsupn0
   integer :: cgmaxit
   external :: gp_ieee_signal1, gp_ieee_signal2

   gpsupn = 1.d-2
   cgepsi = 1.d-1
   cgepsf = 1.d-6
   cggpnf = 1.d-8
   call gp_ieee_signal1(gpsupn, acgeps, bcgeps, cgepsf, cgepsi, cggpnf)
   call assert_true(acgeps > 0.d0, "gp_ieee_signal1 should compute a positive slope")
   call assert_close_real8(acgeps * log10(gpsupn) + bcgeps, log10(cgepsi), 1.d-12, &
      "gp_ieee_signal1 should satisfy the interpolation formula")

   gpeucn2 = 1.d-4
   gpeucn20 = 1.d0
   epsgpen2 = 1.d-8
   epsgpsn = 1.d-6
   gpsupn0 = 1.d0
   call gp_ieee_signal2(cgmaxit, 12, .false., 0, 1, kappa, gpeucn2, gpeucn20, &
      epsgpen2, epsgpsn, cgeps, acgeps, bcgeps, cgepsf, cgepsi, gpsupn, gpsupn0)
   call assert_true(cgmaxit >= 1 .and. cgmaxit <= 12, "gp_ieee_signal2 cgmaxit should be bounded")
   call assert_true(cgeps >= cgepsf .and. cgeps <= cgepsi, "gp_ieee_signal2 cgeps should be clipped")
end program test_ieee_signal
