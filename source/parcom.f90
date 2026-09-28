module parcom
      use iso_fortran_env, only : wp => real64
      implicit none
      private :: wp
      public

      real(wp) glob_par  (20)     
      real(wp) expscal (4,10,86)

      real(wp) ener_par1 (10,86)
      real(wp) ener_par2 (10,86)
      real(wp) ener_par3 (10,86)
      real(wp) ener_par4 (10,86)
      real(wp) ener_par5 (10,86)
      real(wp) ener_par6 (10,86)

      real(wp) shell_xi  (13,86)
      real(wp) shell_cnf1(13,86)
      real(wp) shell_cnf2(13,86)
      real(wp) shell_cnf3(13,86)
      real(wp) shell_cnf4(10,86)
      real(wp) shell_resp(13,86,2)

!     Optional experimental radial response.  It is deliberately not stored in
!     atompara: zero is the production/default model.  The command-line flag
!     is used only by guarded Yb model-development screens.
      real(wp) :: yb_d7_charge_response = 0.0_wp

!     real(wp),parameter :: mull_loew14 = 0.1666666667_wp ! Mulliken-Loewdin mixing factor (0=M, 0.5=L)
!     real(wp),parameter :: mull_loew14 = 0.2500000000_wp ! Mulliken-Loewdin mixing factor (0=M, 0.5=L)
      real(wp),parameter :: mull_loew14 = 0.3333333333_wp ! Mulliken-Loewdin mixing factor (0=M, 0.5=L)
                                                          ! was 1/4 for quite some time

end module parcom       
