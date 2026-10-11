! Driver of the standalone NL6 check (check_nl6.py): W3SNL6 of model/src/w3snl6md.F90 on spectra E(f,theta) from a file,
! called as W3SRCE calls it, on the grid given on the command line. SIG as ww3_grid builds it from FR1 and XFR in single
! precision (W3GRID); deep-water CG; A = E * CONX and S(E) = S / CONX, CONX = CG / (2 pi SIG) as in W3SNL1. The depth,
! 10 km, makes the WAM depth factor that W3SNL6 applies equal to 1, so that the operator alone is checked.
! usage: nl6_driver in.bin out.bin nb fr1 xfr nk nth
!   in.bin  nb spectra E, REAL*4, WAVEWATCH III order (direction fastest), stream access; out.bin: S(E), same layout,
!           followed by the diagonal term D (dS/dA = dS(E)/dE per bin)
PROGRAM NL6_DRIVER
  USE CONSTANTS, ONLY: TPI, TPIINV, GRAV
  USE W3GDATMD,  ONLY: NK, NTH, NSPEC, SIG
  USE W3SNL6MD,  ONLY: W3SNL6
  IMPLICIT NONE
  REAL, PARAMETER     :: DEPTH = 1.E4
  CHARACTER(LEN=1024) :: FIN, FOUT, ARG
  INTEGER             :: NB, IB, IK, ITH, ISP, IU
  REAL                :: SIGMA, FR1, XFR
  REAL, ALLOCATABLE   :: E(:,:), SOUT(:,:), DOUT(:,:), A(:), S(:), D(:), CG(:), CONX(:)
  !
  CALL GET_COMMAND_ARGUMENT ( 1, FIN )
  CALL GET_COMMAND_ARGUMENT ( 2, FOUT )
  CALL GET_COMMAND_ARGUMENT ( 3, ARG ); READ (ARG,*) NB
  CALL GET_COMMAND_ARGUMENT ( 4, ARG ); READ (ARG,*) FR1
  CALL GET_COMMAND_ARGUMENT ( 5, ARG ); READ (ARG,*) XFR
  CALL GET_COMMAND_ARGUMENT ( 6, ARG ); READ (ARG,*) NK
  CALL GET_COMMAND_ARGUMENT ( 7, ARG ); READ (ARG,*) NTH
  NSPEC = NK * NTH
  ALLOCATE ( SIG(0:NK+1), CG(NK), CONX(NK), A(NSPEC), S(NSPEC), D(NSPEC) )
  SIGMA = FR1 * TPI / XFR**2                             ! as W3GRID
  DO IK=0, NK+1
    SIGMA   = SIGMA * XFR
    SIG(IK) = SIGMA
  END DO
  DO IK=1, NK
    CG(IK)   = 0.5 * GRAV / SIG(IK)                      ! deep water
    CONX(IK) = TPIINV / SIG(IK) * CG(IK)
  END DO
  !
  ALLOCATE ( E(NSPEC,NB), SOUT(NSPEC,NB), DOUT(NSPEC,NB) )
  OPEN ( NEWUNIT=IU, FILE=TRIM(FIN), ACCESS='STREAM', FORM='UNFORMATTED', STATUS='OLD', ACTION='READ' )
  READ (IU) E
  CLOSE ( IU )
  DO IB=1, NB
    DO IK=1, NK
      DO ITH=1, NTH
        ISP    = ITH + (IK-1)*NTH
        A(ISP) = E(ISP,IB) * CONX(IK)
      END DO
    END DO
    CALL W3SNL6 ( A, CG, DEPTH, S, D )
    DO IK=1, NK
      DO ITH=1, NTH
        ISP           = ITH + (IK-1)*NTH
        SOUT(ISP,IB)  = S(ISP) / CONX(IK)
        DOUT(ISP,IB)  = D(ISP)
      END DO
    END DO
  END DO
  OPEN ( NEWUNIT=IU, FILE=TRIM(FOUT), ACCESS='STREAM', FORM='UNFORMATTED', STATUS='REPLACE', ACTION='WRITE' )
  WRITE (IU) SOUT, DOUT
  CLOSE ( IU )
END PROGRAM NL6_DRIVER
