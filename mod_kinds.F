      module mod_kinds
       
        public
  
      ! Kind helpers
  
       integer, parameter :: i1b= selected_int_kind(2)        !  1-byte
       integer, parameter :: i2b= selected_int_kind(4)        !  2-byte
       integer, parameter :: i4b= selected_int_kind(9)        !  4-byte
       integer, parameter :: c8 = selected_real_kind(6,30)    ! 32-bit
       integer, parameter :: r4 = selected_real_kind(6,30)    ! 32-bit
       integer, parameter :: r8 = selected_real_kind(12,300)  ! 64-bit


      end module mod_kinds

