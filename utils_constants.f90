      module utils_constants
       
        public
  
      ! Kind helpers
      integer , parameter :: dp = kind(1.0D0)
      integer , parameter :: sp = kind(1.0)
      integer , parameter :: r8 = kind(1.0D0)
      integer , parameter :: r4 = kind(1.0)
      integer , parameter :: shortint = 2
      integer , parameter :: normint = 4
      integer , parameter :: longint = 8
  
      ! Numbers 1-10
      real(dp) , parameter :: d_0 = 0.0_r8
      real(dp) , parameter :: d_1 = 1.0_r8
      real(dp) , parameter :: d_2 = 2.0_r8
      real(dp) , parameter :: d_3 = 3.0_r8
      real(dp) , parameter :: d_4 = 4.0_r8
      real(dp) , parameter :: d_5 = 5.0_r8
      real(dp) , parameter :: d_6 = 6.0_r8
      real(dp) , parameter :: d_7 = 7.0_r8
      real(dp) , parameter :: d_8 = 8.0_r8
      real(dp) , parameter :: d_9 = 9.0_r8
      real(dp) , parameter :: d_10 = 10.0_r8
      real(dp) , parameter :: d_16 = 16.0_r8
      real(dp) , parameter :: d_32 = 32.0_r8
  
      ! Simple Fractions
      real(dp) , parameter :: d_half = d_1/d_2
  
      real(dp) , parameter :: d_1q2 = d_1/d_2
      real(dp) , parameter :: d_3q2 = d_3/d_2
      real(dp) , parameter :: d_5q2 = d_5/d_2
      real(dp) , parameter :: d_7q2 = d_7/d_2
      real(dp) , parameter :: d_9q2 = d_9/d_2
      real(dp) , parameter :: d_1q3 = d_1/d_3
      real(dp) , parameter :: d_2q3 = d_2/d_3
      real(dp) , parameter :: d_4q3 = d_4/d_3
      real(dp) , parameter :: d_5q3 = d_5/d_3
      real(dp) , parameter :: d_1q4 = d_1/d_4
      real(dp) , parameter :: d_3q4 = d_3/d_4
      real(dp) , parameter :: d_5q4 = d_5/d_4
      real(dp) , parameter :: d_7q4 = d_7/d_4
      real(dp) , parameter :: d_1q5 = d_1/d_5

      end module utils_constants

