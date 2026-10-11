! Stand-ins for the WAVEWATCH III modules that model/src/w3snl6md.F90 uses (standalone check, check_nl6.py): only the
! entities it reads, with the values W3NOUT sets (units 6, one process). CONSTANTS is the model's own constants.F90.
MODULE W3GDATMD
  INTEGER          :: NK, NTH, NSPEC
  REAL, POINTER    :: SIG(:)                               ! SIG(0:NK+1), as in WAVEWATCH III
END MODULE W3GDATMD

MODULE W3ODATMD
  INTEGER          :: NDSO = 6, NDSE = 6, NDST = 6, IAPROC = 1, NAPERR = 1, NAPOUT = 1
END MODULE W3ODATMD

MODULE W3SERVMD
CONTAINS
  SUBROUTINE EXTCDE ( IEXIT, UNIT, MSG, FILE, LINE, COMM ) ! as W3SERVMD without MPI: optional message, then exit
    INTEGER, INTENT(IN)                :: IEXIT
    INTEGER, INTENT(IN), OPTIONAL      :: UNIT, LINE, COMM
    CHARACTER(*), INTENT(IN), OPTIONAL :: MSG, FILE
    IF ( PRESENT(MSG) ) WRITE (0,'(A)') 'WW3 ERROR: '//TRIM(MSG)
    CALL EXIT ( IEXIT )
  END SUBROUTINE EXTCDE
END MODULE W3SERVMD
