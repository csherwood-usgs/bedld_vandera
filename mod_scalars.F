      MODULE mod_scalars
     
      USE mod_kinds

      real(r8) :: Eradius = 6371315.0_r8      ! m
      real(r8) :: StefBo = 5.78E-8_r8         ! Watts/m2/K4 (Match MU)
      real(r8) :: emmiss = 0.97_r8            ! non_dimensional
      real(r8) :: rho0 = 1025.0_r8
      real(r8) :: g = 9.81_r8                 ! m/s2
      real(r8) :: vonKar = 0.41_r8            ! non-dimensional
      real(r8) :: nu = 1.36E-6_r8            ! non-dimensional
      real(r8), parameter :: pi = 3.14159265358979323846_r8
      real(r8), parameter :: deg2rad = pi / 180.0_r8
      real(r8), parameter :: rad2deg = 180.0_r8 / pi
      real(r8), parameter :: day2sec = 86400.0_r8
      real(r8), parameter :: sec2day = 1.0_r8 / 86400.0_r8

      END MODULE mod_scalars
