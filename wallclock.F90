module wallclock_mod
  implicit none
  INTERFACE
     FUNCTION wallclock()
       integer, parameter :: dp=kind(1.0d0)
       REAL(dp) :: wallclock
     END FUNCTION wallclock
  END INTERFACE

end module wallclock_mod
