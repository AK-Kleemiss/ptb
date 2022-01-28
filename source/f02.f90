REAL*8 FUNCTION F02(ARG)
      IMPLICIT REAL*8 (A-H,O-Z)    
      DATA PI/3.14159265358979D0/

      IF(ARG.LT.1.0D-10) GOTO 1

      F02=0.50D0*DSQRT(PI/ARG)*DERF(DSQRT(ARG))
      RETURN

    1 F02=1.0D0-0.33333333333333D0*ARG
      RETURN

END 
