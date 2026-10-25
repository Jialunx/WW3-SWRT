!> @file
!> @brief Nonlinear 4-wave interactions by sparse WRT (sWRT), NL6 slot.
!>
!> @date   25-Oct-2026
!>

#include "w3macros.h"
!/ ------------------------------------------------------------------- /
!>
!> @brief sWRT: the WRT pair sum restricted to K mirror-pair units of
!>  its quadrature table, one gain per unit (WAM depth factor).
!>
MODULE W3SNL6MD
  !/
  !/                  +-----------------------------------+
  !/                  | WAVEWATCH III                     |
  !/                  |      sWRT (sparse WRT) as NL6     |
  !/                  |                        FORTRAN 90 |
  !/                  | Last update :         25-Oct-2026 |
  !/                  +-----------------------------------+
  !/
  !/    24-Sep-2026 : Origination, replaces the ONNX surrogate of the
  !/                  WW3-ML-Snl source in this copy.
  !/    24-Sep-2026 : Diagonal term as in Xnl (was D = 0).
  !/    07-Oct-2026 : Vectorized loop (FASTLOOP) by default; the
  !/                  original loop if WW3_SWRT_KERNEL = reference.
  !/                  Same transfer and diagonal term to round-off.
  !/    09-Oct-2026 : WAM depth scaling of S as in Xnl (q_dscale,
  !/                  IQTYPE = 2); diagonal term unscaled as in Xnl.
  !/                  DEPTH enters only this factor.
  !/    20-Oct-2026 : Rows grouped by resonance curve (LOCUSLOOP) by
  !/                  default: the rows of one curve share N1 and N3,
  !/                  so two sums per curve give the transfer and the
  !/                  diagonal term; anchors up to NKLIM (0: all).
  !/                  WW3_SWRT_KERNEL = vector selects FASTLOOP,
  !/                  = reference the original loop.
  !/    25-Oct-2026 : W3SRCE sets NKLIM = 0, so that S_nl is computed
  !/                  at all frequencies, as with DIA and WRT.
  !/
  !  1. Purpose :
  !
  !     Nonlinear interactions by sparse WRT (sWRT): the WRT pair sum
  !     (Xnl, as NL2) restricted to K mirror-pair units (a row of WRT's
  !     own quadrature table plus its mirror row), one gain per unit,
  !     times the WAM depth factor of Xnl. The transfer is that of the
  !     Python operator of tools/swrt.py, for one spectrum, on
  !     E(f,theta) = A / CONX with CONX = CG / (2 pi SIG) (as W3SNL1);
  !     S = CONX * dE/dt. D = dS/dA per bin as Xnl (NL2) computes it
  !     for WRT (q_t13v4, diagk1 / diagk3): every term differentiated
  !     with respect to its own N1 and N3, N2 and N4 held fixed (as
  !     tools/swrt.py). D = 0 if the environment variable WW3_SWRT_DIAG
  !     is 0 (for sensitivity tests).
  !
  !  2. Operator file (tools/export_operator.py) :
  !
  !     Plain text: magic line WW3-SWRT-OPERATOR, format version (1),
  !     'nk na nr K', source name, then blocks of a name line followed
  !     by the values: sigma(nk), di, dj, i2, f2, j2, g2, i4, f4, j4,
  !     g4 (nr each), w, fa, fp (nk,nr), conv(nk), qfrac, kpow(nk).
  !     Read at the first call of W3SNL6 from the file named by the
  !     environment variable WW3_SWRT_OPERATOR, default
  !     'swrt_operator.ww3op' in the run directory. The model must have
  !     NK = nk, NTH = na and SIG = sigma to a relative 1.E-5.
  !
  !  3. Subroutines :
  !
  !      Name      Type  Scope    Description
  !     ----------------------------------------------------------------
  !      W3SNL6    Subr. Public   sWRT source term and diagonal.
  !      W3SNLGQM  Subr. Public   Same, for the IQTYPE <= 0 call sites.
  !      INSNL1    Subr. Public   Called by W3IOGR, nothing to set up.
  !      INSNLGQM  Subr. Public   Called by W3IOGR, nothing to set up.
  !      SWRTRD    Subr. Private  Read and check the operator file.
  !     ----------------------------------------------------------------
  !
  !  4. Remarks :
  !
  !     DEPTH enters only the WAM depth factor of S. The SNL1 namelist
  !     of ww3_grid (DIA parameters) is still read, but not used.
  !
  !/ ------------------------------------------------------------------- /
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: W3SNL6, W3SNLGQM, INSNL1, INSNLGQM, NKLIM
  !/
  !/ Operator: the row arrays of the kernel, read once by SWRTRD
  !/
  LOGICAL, SAVE              :: OPREAD = .FALSE.
  INTEGER, SAVE              :: NR = 0
  INTEGER, ALLOCATABLE, SAVE :: DI(:), DJ(:), I2(:), J2(:), I4(:), J4(:)
  REAL, ALLOCATABLE, SAVE    :: F2(:), G2(:), F4(:), G4(:), CONV(:),   &
                                KPOW(:), W(:,:), FA(:,:), FP(:,:)
  REAL, SAVE                 :: QFRAC = 0.
  LOGICAL, SAVE              :: DODIAG = .TRUE.
  LOGICAL, SAVE              :: FASTK  = .FALSE.
  !/ Rows grouped by resonance curve (same di, dj): curve G holds rows
  !/ GS(G) to GS(G+1)-1; NKLIM = last anchor bin (0: all), set by W3SRCE
  LOGICAL, SAVE              :: LOCUSK = .FALSE.
  INTEGER, SAVE              :: NG = 0, NKLIM = 0
  INTEGER, ALLOCATABLE, SAVE :: GS(:)
  !/
CONTAINS
  !/ ------------------------------------------------------------------- /
!>
!> @brief sWRT nonlinear interactions for one action spectrum.
!>
!> @param[in]  A      Action spectrum A(ISP), ISP = ITH + (IK-1)*NTH.
!> @param[in]  CG     Group velocities (dimension NK).
!> @param[in]  DEPTH  Water depth (WAM depth factor of S).
!> @param[out] S      Source term, action units per second.
!> @param[out] D      Diagonal term of the derivative, dS/dA per bin.
!>
  SUBROUTINE W3SNL6 ( A, CG, DEPTH, S, D )
    USE CONSTANTS, ONLY: TPIINV
    USE W3GDATMD,  ONLY: NK, NTH, NSPEC, SIG
    IMPLICIT NONE
    REAL, INTENT(IN)        :: A(NSPEC), CG(NK), DEPTH
    REAL, INTENT(OUT)       :: S(NSPEC), D(NSPEC)
    INTEGER                 :: IK, ITH, ISP, R, I1, J1, I3, J3,       &
                               IA2, IB2, IA4, IB4, JA2, JB2, JA4, JB4
    REAL                    :: CONX(NK), N(NK,NTH), SE(NK,NTH), NMAX, &
                               QMIN, N1, N2, N3, N4, T, FU2, FU4,     &
                               DE(NK,NTH)
    !
    ! 0.  Operator, read at the first call
    !
    IF ( .NOT. OPREAD ) THEN
!$OMP CRITICAL (W3SNL6_READ)
      IF ( .NOT. OPREAD ) CALL SWRTRD
!$OMP END CRITICAL (W3SNL6_READ)
    END IF
    !
    ! 1.  Action A -> E(f,theta) = A / CONX -> N(kx,ky) = E * conv
    !
    DO IK=1, NK
      CONX(IK) = TPIINV / SIG(IK) * CG(IK)
      DO ITH=1, NTH
        ISP       = ITH + (IK-1)*NTH
        N(IK,ITH) = ( A(ISP) / CONX(IK) ) * CONV(IK)
      END DO
    END DO
    NMAX = MAXVAL(N)
    SE   = 0.
    DE   = 0.
    !
    ! 2.  Rows x anchor bins (i1,j1)
    !
    IF ( LOCUSK ) THEN
      CALL LOCUSLOOP
    ELSE IF ( FASTK ) THEN
      CALL FASTLOOP
    ELSE
    DO R=1, NR
      DO I1=1, NK
        I3 = I1 + DI(R)
        IF ( I3 .GT. NK ) CYCLE
        IF ( W(I1,R) .EQ. 0. ) CYCLE
        QMIN = QFRAC * NMAX * KPOW(I3)
        CALL CORNER ( I1+I2(R), F2(R), IA2, IB2, FU2 )
        CALL CORNER ( I1+I4(R), F4(R), IA4, IB4, FU4 )
        DO J1=1, NTH
          J3  = MODULO(J1-1+DJ(R),NTH) + 1
          N1  = N(I1,J1)
          N3  = N(I3,J3)
          IF ( N1.LT.QMIN .AND. N3.LT.QMIN ) CYCLE
          JA2 = MODULO(J1-1+J2(R),NTH) + 1
          JB2 = MODULO(JA2,NTH) + 1
          JA4 = MODULO(J1-1+J4(R),NTH) + 1
          JB4 = MODULO(JA4,NTH) + 1
          N2  = (1.-FU2) * ( (1.-G2(R))*N(IA2,JA2) + G2(R)*N(IA2,JB2) ) &
              +     FU2  * ( (1.-G2(R))*N(IB2,JA2) + G2(R)*N(IB2,JB2) )
          N4  = (1.-FU4) * ( (1.-G4(R))*N(IA4,JA4) + G4(R)*N(IA4,JB4) ) &
              +     FU4  * ( (1.-G4(R))*N(IB4,JA4) + G4(R)*N(IB4,JB4) )
          T   = W(I1,R) * ( N1*N3*(N4-N2) + N2*N4*(N3-N1) )
          SE(I1,J1) = SE(I1,J1) + T*FA(I1,R)
          SE(I3,J3) = SE(I3,J3) - T*FP(I1,R)
          DE(I1,J1) = DE(I1,J1) + W(I1,R) *                           &
                      ( N3*(N4-N2) - N2*N4 ) * FA(I1,R)
          DE(I3,J3) = DE(I3,J3) - W(I1,R) *                           &
                      ( N1*(N4-N2) + N2*N4 ) * FP(I1,R)
        END DO
      END DO
    END DO
    END IF
    !
    ! 2b. WAM depth scaling of the transfer, as Xnl does for IQTYPE = 2
    !
    SE = SE * DFAC ( A, CG, DEPTH )
    !
    ! 3.  Back to action units; D = dS/dA = dSE/dE = dSE/dN * CONV
    !
    IF ( .NOT. DODIAG ) DE = 0.
    DO IK=1, NK
      DO ITH=1, NTH
        ISP    = ITH + (IK-1)*NTH
        S(ISP) = SE(IK,ITH) * CONX(IK)
        D(ISP) = DE(IK,ITH) * CONV(IK)
      END DO
    END DO
    !
  CONTAINS
    !/ ------------------------------------------------------------------- /
    SUBROUTINE CORNER ( I0, FIN, IA, IB, FU )
      ! Frequency corners as swrt.Rows._corners: lower corner clamped to
      ! the grid with fraction 0 when it is off the grid, upper = IA+1.
      INTEGER, INTENT(IN)  :: I0
      REAL,    INTENT(IN)  :: FIN
      INTEGER, INTENT(OUT) :: IA, IB
      REAL,    INTENT(OUT) :: FU
      IA = I0
      FU = FIN
      IF ( IA.LT.1 .OR. IA.GT.NK ) THEN
        FU = 0.
        IA = MIN ( MAX ( IA, 1 ), NK )
      END IF
      IB = MIN ( IA+1, NK )
    END SUBROUTINE CORNER
    !/ ------------------------------------------------------------------- /
    REAL FUNCTION DFAC ( A, CG, DEPTH )
      ! WAM depth factor of Xnl (q_dscale): mean wavenumber of the action
      ! sigma spectrum A/CG with weights k**(-1/2), steps of z_steps,
      ! wavenumbers of z_wnumb (Hunt), kd = max(0.5, 0.75 km d).
      USE CONSTANTS, ONLY: GRAV
      REAL, INTENT(IN) :: A(NSPEC), CG(NK), DEPTH
      REAL    :: DSG(NK), SUM0, SUMK, Y, XX, KK, KMS, KD, DNN
      INTEGER :: IK, ITH
      DO IK=2, NK-1
        DSG(IK) = 0.5 * ( SIG(IK+1) - SIG(IK-1) )
      END DO
      DSG(1)  = DSG(2) * DSG(2) / DSG(3)
      DSG(NK) = DSG(NK-1) * DSG(NK-1) / DSG(NK-2)
      SUM0 = 0.
      SUMK = 0.
      DO IK=1, NK
        Y  = SIG(IK)**2 / GRAV * DEPTH
        XX = Y*(Y+1./(1.+Y*(0.66667+Y*(0.35550+Y*(0.16084+Y*(0.06320+Y* &
             (0.02174+Y*(0.00654+Y*(0.00171+Y*(0.00039+Y*0.00011))))))))))
        KK = SQRT(XX) / DEPTH
        DO ITH=1, NTH
          DNN  = A(ITH+(IK-1)*NTH) / CG(IK) * DSG(IK)
          SUM0 = SUM0 + DNN
          SUMK = SUMK + DNN / SQRT(KK)
        END DO
      END DO
      IF ( SUM0 .GT. 0. .AND. DEPTH .GT. 0. ) THEN
        KMS  = ( SUM0 / SUMK )**2
        KD   = MAX ( 0.5, 0.75*KMS*DEPTH )
        DFAC = 1. + 5.5/KD * (1.-5./6.*KD) * EXP(-5./4.*KD)
      ELSE
        DFAC = 1.
      END IF
    END FUNCTION DFAC
    !/ ------------------------------------------------------------------- /
    SUBROUTINE FASTLOOP
      ! The loop of section 2 reorganized for SIMD (the default):
      ! direction-first arrays, a padded copy of N (one period on each side)
      ! instead of MODULO, the density filter as a 0/1 mask; the padded sums
      ! are folded back at the end. Same transfer and diagonal term.
      REAL    :: NX(1-NTH:2*NTH,NK), SX(1-NTH:2*NTH,NK), DX(1-NTH:2*NTH,NK)
      REAL    :: QM, M, X1, X2, X3, X4, TT, A2, B2, A4, B4, WR, WA, WP,   &
                 U2, U4
      INTEGER :: RR, K1, K3, JJ, O3, O2, O4, KA2, KB2, KA4, KB4
      DO K1=1, NK
        NX(1:NTH,K1)       = N(K1,:)
        NX(1-NTH:0,K1)     = N(K1,:)
        NX(NTH+1:2*NTH,K1) = N(K1,:)
      END DO
      SX = 0.
      DX = 0.
      DO RR=1, NR
        O3 = DJ(RR)
        O2 = J2(RR)
        O4 = J4(RR)
        B2 = G2(RR)
        A2 = 1. - B2
        B4 = G4(RR)
        A4 = 1. - B4
        DO K1=1, NK
          K3 = K1 + DI(RR)
          IF ( K3 .GT. NK ) CYCLE
          WR = W(K1,RR)
          IF ( WR .EQ. 0. ) CYCLE
          QM = QFRAC * NMAX * KPOW(K3)
          WA = WR * FA(K1,RR)
          WP = WR * FP(K1,RR)
          CALL CORNER ( K1+I2(RR), F2(RR), KA2, KB2, U2 )
          CALL CORNER ( K1+I4(RR), F4(RR), KA4, KB4, U4 )
          DO JJ=1, NTH
            X1 = NX(JJ,K1)
            X3 = NX(JJ+O3,K3)
            M  = MERGE ( 1., 0., X1.GE.QM .OR. X3.GE.QM )
            X2 = (1.-U2) * ( A2*NX(JJ+O2,KA2) + B2*NX(JJ+O2+1,KA2) )     &
               +     U2  * ( A2*NX(JJ+O2,KB2) + B2*NX(JJ+O2+1,KB2) )
            X4 = (1.-U4) * ( A4*NX(JJ+O4,KA4) + B4*NX(JJ+O4+1,KA4) )     &
               +     U4  * ( A4*NX(JJ+O4,KB4) + B4*NX(JJ+O4+1,KB4) )
            TT = M * ( X1*X3*(X4-X2) + X2*X4*(X3-X1) )
            SX(JJ,K1)    = SX(JJ,K1)    + TT * WA
            SX(JJ+O3,K3) = SX(JJ+O3,K3) - TT * WP
            DX(JJ,K1)    = DX(JJ,K1)    + M * ( X3*(X4-X2) - X2*X4 ) * WA
            DX(JJ+O3,K3) = DX(JJ+O3,K3) - M * ( X1*(X4-X2) + X2*X4 ) * WP
          END DO
        END DO
      END DO
      DO K1=1, NK
        SE(K1,:) = SX(1:NTH,K1) + SX(1-NTH:0,K1) + SX(NTH+1:2*NTH,K1)
        DE(K1,:) = DX(1:NTH,K1) + DX(1-NTH:0,K1) + DX(NTH+1:2*NTH,K1)
      END DO
    END SUBROUTINE FASTLOOP
    !/ ------------------------------------------------------------------- /
    SUBROUTINE LOCUSLOOP
      ! FASTLOOP with the rows taken curve by curve (the default). All rows
      ! of a resonance curve share k1 and k3, hence N1 and N3, so the sums
      ! A = sum W (N4 - N2) and B = sum W N2 N4 over its rows give the
      ! transfer T = N1 N3 A + (N3 - N1) B and the diagonal terms
      ! dT/dN1 = N3 A - B and dT/dN3 = N1 A + B, deposited once per curve.
      ! A curve is skipped at an anchor when the frequency rows of k1 and k3
      ! are both under the density threshold (every term is masked there).
      ! Anchors run up to NKLIM when W3SRCE has set it (WAVEWATCH III
      ! replaces the spectrum above its last prognostic bin).
      REAL    :: NX(1-NTH:2*NTH,NK), SX(1-NTH:2*NTH,NK), DX(1-NTH:2*NTH,NK)
      REAL    :: X1(NTH), X3(NTH), MM(NTH), AA(NTH), BB(NTH), AS1(NTH), &
                 AD1(NTH), RMX(NK)
      REAL    :: QM, X2, X4, A2, B2, A4, B4, WR, WA, WP, U2, U4, TS
      INTEGER :: G, RR, R0, K1, K3, JJ, O3, O2, O4, KA2, KB2, KA4, KB4,   &
                 NU, KMX
      DO K1=1, NK
        NX(1:NTH,K1)       = N(K1,:)
        NX(1-NTH:0,K1)     = N(K1,:)
        NX(NTH+1:2*NTH,K1) = N(K1,:)
        RMX(K1)            = MAXVAL(N(K1,:))
      END DO
      SX = 0.
      DX = 0.
      KMX = NK
      IF ( NKLIM .GT. 0 ) KMX = MIN ( NK, NKLIM )
      DO K1=1, KMX
        X1(1:NTH) = NX(1:NTH,K1)
        AS1 = 0.
        AD1 = 0.
        DO G=1, NG
          R0 = GS(G)
          O3 = DJ(R0)
          K3 = K1 + DI(R0)
          IF ( K3 .GT. NK ) CYCLE
          QM = QFRAC * NMAX * KPOW(K3)
          IF ( RMX(K1).LT.QM .AND. RMX(K3).LT.QM ) CYCLE
          DO JJ=1, NTH
            X3(JJ) = NX(JJ+O3,K3)
            MM(JJ) = MERGE ( 1., 0., X1(JJ).GE.QM .OR. X3(JJ).GE.QM )
          END DO
          AA = 0.
          BB = 0.
          NU = 0
          DO RR=GS(G), GS(G+1)-1
            WR = W(K1,RR)
            IF ( WR .EQ. 0. ) CYCLE
            NU = NU + 1
            O2 = J2(RR)
            O4 = J4(RR)
            B2 = G2(RR)
            A2 = 1. - B2
            B4 = G4(RR)
            A4 = 1. - B4
            CALL CORNER ( K1+I2(RR), F2(RR), KA2, KB2, U2 )
            CALL CORNER ( K1+I4(RR), F4(RR), KA4, KB4, U4 )
            DO JJ=1, NTH
              X2 = (1.-U2) * ( A2*NX(JJ+O2,KA2) + B2*NX(JJ+O2+1,KA2) )   &
                 +     U2  * ( A2*NX(JJ+O2,KB2) + B2*NX(JJ+O2+1,KB2) )
              X4 = (1.-U4) * ( A4*NX(JJ+O4,KA4) + B4*NX(JJ+O4+1,KA4) )   &
                 +     U4  * ( A4*NX(JJ+O4,KB4) + B4*NX(JJ+O4+1,KB4) )
              AA(JJ) = AA(JJ) + WR * ( X4 - X2 )
              BB(JJ) = BB(JJ) + WR * X2 * X4
            END DO
          END DO
          IF ( NU .EQ. 0 ) CYCLE
          WA = FA(K1,R0)
          WP = FP(K1,R0)
          DO JJ=1, NTH
            TS = MM(JJ) * ( X1(JJ)*X3(JJ)*AA(JJ) + (X3(JJ)-X1(JJ))*BB(JJ) )
            AS1(JJ) = AS1(JJ) + TS * WA
            SX(JJ+O3,K3) = SX(JJ+O3,K3) - TS * WP
          END DO
          IF ( DODIAG ) THEN
            DO JJ=1, NTH
              AD1(JJ) = AD1(JJ) + MM(JJ) * ( X3(JJ)*AA(JJ) - BB(JJ) ) * WA
              DX(JJ+O3,K3) = DX(JJ+O3,K3) - MM(JJ) * ( X1(JJ)*AA(JJ) + BB(JJ) ) * WP
            END DO
          END IF
        END DO
        SX(1:NTH,K1) = SX(1:NTH,K1) + AS1
        DX(1:NTH,K1) = DX(1:NTH,K1) + AD1
      END DO
      DO K1=1, NK
        SE(K1,:) = SX(1:NTH,K1) + SX(1-NTH:0,K1) + SX(NTH+1:2*NTH,K1)
        DE(K1,:) = DX(1:NTH,K1) + DX(1-NTH:0,K1) + DX(NTH+1:2*NTH,K1)
      END DO
    END SUBROUTINE LOCUSLOOP
    !/
  END SUBROUTINE W3SNL6
  !/ ------------------------------------------------------------------- /
!>
!> @brief sWRT for the call sites of the GQM option (IQTYPE <= 0).
!>
!> @param[in]  A      Action spectrum.
!> @param[in]  CG     Group velocities.
!> @param[in]  WN     Wavenumbers (not used).
!> @param[in]  DEPTH  Water depth (WAM depth factor of S).
!> @param[out] S      Source term.
!> @param[out] D      Diagonal term.
!>
  SUBROUTINE W3SNLGQM ( A, CG, WN, DEPTH, S, D )
    USE W3GDATMD, ONLY: NK, NSPEC
    IMPLICIT NONE
    REAL, INTENT(IN)        :: A(NSPEC), CG(NK), WN(NK), DEPTH
    REAL, INTENT(OUT)       :: S(NSPEC), D(NSPEC)
    CALL W3SNL6 ( A, CG, DEPTH, S, D )
  END SUBROUTINE W3SNLGQM
  !/ ------------------------------------------------------------------- /
!>
!> @brief Called by W3IOGR when a model definition file is read. The
!>  operator is read at the first call of W3SNL6 instead, so that
!>  programs that do not compute S_nl do not need the file.
!>
!> @param[in] IMOD  Model number.
!>
  SUBROUTINE INSNL1 ( IMOD )
    IMPLICIT NONE
    INTEGER, INTENT(IN)     :: IMOD
  END SUBROUTINE INSNL1
  !/ ------------------------------------------------------------------- /
!>
!> @brief Called by W3IOGR for IQTYPE <= 0; nothing to set up.
!>
  SUBROUTINE INSNLGQM
    IMPLICIT NONE
  END SUBROUTINE INSNLGQM
  !/ ------------------------------------------------------------------- /
!>
!> @brief Read the operator file and check it against the model grid.
!>
  SUBROUTINE SWRTRD
    USE W3GDATMD, ONLY: NK, NTH, SIG
    USE W3ODATMD, ONLY: NDSE, NDSO, IAPROC, NAPERR, NAPOUT
    USE W3SERVMD, ONLY: EXTCDE
    IMPLICIT NONE
    CHARACTER(LEN=1024)     :: FNAME
    CHARACTER(LEN=256)      :: LINE, SRCNAM, MSG
    INTEGER                 :: NDSOP, IERR, ILEN, IVER, FNK, FNA, KU, &
                               IK, IKW
    CHARACTER(LEN=16)       :: DGVAL
    REAL                    :: RDIF, RMAX
    REAL, ALLOCATABLE       :: FSIG(:)
    !
    ! 1.  File name and header
    !
    CALL GET_ENVIRONMENT_VARIABLE ( 'WW3_SWRT_OPERATOR', FNAME, ILEN, IERR )
    IF ( IERR .EQ. -1 ) CALL SWRTER ( 'WW3_SWRT_OPERATOR is too long' )
    IF ( IERR .NE. 0 .OR. ILEN .EQ. 0 ) FNAME = 'swrt_operator.ww3op'
    !
    OPEN ( NEWUNIT=NDSOP, FILE=TRIM(FNAME), STATUS='OLD',             &
           ACTION='READ', FORM='FORMATTED', IOSTAT=IERR )
    IF ( IERR .NE. 0 ) CALL SWRTER ( 'cannot open the file ' //       &
         '(environment variable WW3_SWRT_OPERATOR, default ' //       &
         'swrt_operator.ww3op in the run directory)' )
    READ (NDSOP,'(A)',IOSTAT=IERR) LINE
    IF ( IERR.NE.0 .OR. TRIM(ADJUSTL(LINE)).NE.'WW3-SWRT-OPERATOR' )  &
         CALL SWRTER ( 'first line is not WW3-SWRT-OPERATOR' )
    READ (NDSOP,*,IOSTAT=IERR) IVER
    IF ( IERR.NE.0 .OR. IVER.NE.1 )                                   &
         CALL SWRTER ( 'unknown format version (this code reads 1)' )
    READ (NDSOP,*,IOSTAT=IERR) FNK, FNA, NR, KU
    IF ( IERR.NE.0 .OR. MIN(FNK,FNA,NR).LT.1 )                        &
         CALL SWRTER ( 'cannot read the line nk na nr K' )
    READ (NDSOP,'(A)',IOSTAT=IERR) SRCNAM
    IF ( IERR .NE. 0 ) CALL SWRTER ( 'cannot read the source name' )
    !
    ! 2.  Spectral grid: NK, NTH and SIG
    !
    IF ( FNK.NE.NK .OR. FNA.NE.NTH ) THEN
      WRITE (MSG,'(A,I0,A,I0,A,I0,A,I0)') 'spectral grid mismatch: '//&
           'model NK = ', NK, ', NTH = ', NTH, '; operator nk = ', FNK, &
           ', na = ', FNA
      CALL SWRTER ( MSG )
    END IF
    ALLOCATE ( FSIG(NK) )
    CALL RDREAL ( 'sigma', FSIG, NK )
    RMAX = 0.
    IKW  = 1
    DO IK=1, NK
      RDIF = ABS(SIG(IK)-FSIG(IK)) / MAX(ABS(FSIG(IK)),TINY(1.))
      IF ( RDIF .GT. RMAX ) THEN
        RMAX = RDIF
        IKW  = IK
      END IF
    END DO
    IF ( RMAX .GT. 1.E-5 ) THEN
      WRITE (MSG,'(A,I0,A,ES14.7,A,ES14.7,A,ES9.2,A)')                &
           'frequency mismatch: SIG(', IKW, ') = ', SIG(IKW),           &
           ' rad/s, operator sigma = ', FSIG(IKW), ' (relative ',       &
           RMAX, ' > 1.E-5)'
      CALL SWRTER ( MSG )
    END IF
    !
    ! 3.  Row arrays
    !
    ALLOCATE ( DI(NR), DJ(NR), I2(NR), F2(NR), J2(NR), G2(NR), I4(NR), &
               F4(NR), J4(NR), G4(NR), W(NK,NR), FA(NK,NR), FP(NK,NR), &
               CONV(NK), KPOW(NK) )
    CALL RDINT  ( 'di', DI, NR )
    CALL RDINT  ( 'dj', DJ, NR )
    CALL RDINT  ( 'i2', I2, NR )
    CALL RDREAL ( 'f2', F2, NR )
    CALL RDINT  ( 'j2', J2, NR )
    CALL RDREAL ( 'g2', G2, NR )
    CALL RDINT  ( 'i4', I4, NR )
    CALL RDREAL ( 'f4', F4, NR )
    CALL RDINT  ( 'j4', J4, NR )
    CALL RDREAL ( 'g4', G4, NR )
    CALL RDREAL ( 'w' , W , NK*NR )
    CALL RDREAL ( 'fa', FA, NK*NR )
    CALL RDREAL ( 'fp', FP, NK*NR )
    CALL RDREAL ( 'conv', CONV, NK )
    CALL RDLAB  ( 'qfrac' )
    READ (NDSOP,*,IOSTAT=IERR) QFRAC
    IF ( IERR .NE. 0 ) CALL SWRTER ( 'cannot read the value of qfrac' )
    CALL RDREAL ( 'kpow', KPOW, NK )
    CLOSE ( NDSOP )
    IF ( MINVAL(DI) .LT. 0 ) CALL SWRTER ( 'negative offset di' )
    !
    ! 4.  Diagonal term on (default) or off (WW3_SWRT_DIAG = 0)
    !
    CALL GET_ENVIRONMENT_VARIABLE ( 'WW3_SWRT_DIAG', DGVAL, ILEN, IERR )
    DODIAG = .NOT. ( IERR.EQ.0 .AND. TRIM(ADJUSTL(DGVAL)).EQ.'0' )
    !
    ! 5.  Rows grouped by resonance curve (default), the vectorized row
    !     loop (WW3_SWRT_KERNEL = vector) or the original loop
    !     (WW3_SWRT_KERNEL = reference, for checks)
    !
    CALL GET_ENVIRONMENT_VARIABLE ( 'WW3_SWRT_KERNEL', DGVAL, ILEN, IERR )
    FASTK  = .NOT. ( IERR.EQ.0 .AND. TRIM(ADJUSTL(DGVAL)).EQ.'reference' )
    LOCUSK = FASTK .AND. .NOT. ( IERR.EQ.0 .AND. TRIM(ADJUSTL(DGVAL)).EQ.'vector' )
    IF ( FASTK .AND. ( MINVAL(J2).LT.-NTH .OR. MAXVAL(J2).GT.NTH-1 .OR. &
         MINVAL(J4).LT.-NTH .OR. MAXVAL(J4).GT.NTH-1 .OR.               &
         MAXVAL(ABS(DJ)).GT.NTH ) )                                      &
         CALL SWRTER ( 'direction offsets too large for the fast loop' )
    IF ( IAPROC.EQ.NAPOUT .AND. .NOT.FASTK ) WRITE (NDSO,'(A)')         &
         '  W3SNL6 : original loop (WW3_SWRT_KERNEL = reference)'
    IF ( LOCUSK ) THEN
      ! curves = runs of rows with the same (di, dj); the operator files
      ! list the rows of a curve together
      ALLOCATE ( GS(NR+1) )
      NG    = 1
      GS(1) = 1
      DO IK=2, NR
        IF ( DI(IK).NE.DI(IK-1) .OR. DJ(IK).NE.DJ(IK-1) ) THEN
          NG     = NG + 1
          GS(NG) = IK
        END IF
      END DO
      GS(NG+1) = NR + 1
      IF ( IAPROC.EQ.NAPOUT ) WRITE (NDSO,'(A,I0,A)')                   &
           '  W3SNL6 : rows grouped by resonance curve, ', NG, ' curves'
    ELSE IF ( IAPROC.EQ.NAPOUT .AND. FASTK ) THEN
      WRITE (NDSO,'(A)') '  W3SNL6 : vectorized row loop (WW3_SWRT_KERNEL = vector)'
    END IF
    !
    IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,900) TRIM(FNAME),            &
         TRIM(SRCNAM), KU, NR, MERGE('on ','off',DODIAG)
    OPREAD = .TRUE.
    !
900 FORMAT (/'  W3SNL6 : sWRT operator ',A/                            &
             '           (',A,', K = ',I0,' units, ',I0,' rows,',      &
             ' diagonal term ',A,')'/)
    !
  CONTAINS
    !/ ------------------------------------------------------------------- /
    SUBROUTINE SWRTER ( TXT )
      ! Error message (file name and TXT) and stop.
      CHARACTER(LEN=*), INTENT(IN) :: TXT
      IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,'(/A/2A/2A/)')             &
           ' *** WAVEWATCH III ERROR IN W3SNL6 (sWRT) :',               &
           '     OPERATOR FILE : ', TRIM(FNAME), '     ', TRIM(TXT)
      CALL EXTCDE ( 1 )
    END SUBROUTINE SWRTER
    !/
    SUBROUTINE RDLAB ( LABEL )
      ! Name line of the next block.
      CHARACTER(LEN=*), INTENT(IN) :: LABEL
      READ (NDSOP,'(A)',IOSTAT=IERR) LINE
      IF ( IERR.NE.0 .OR. TRIM(ADJUSTL(LINE)).NE.LABEL )              &
           CALL SWRTER ( 'expected block '//LABEL//', found: '//      &
           TRIM(LINE) )
    END SUBROUTINE RDLAB
    !/
    SUBROUTINE RDINT ( LABEL, X, NV )
      CHARACTER(LEN=*), INTENT(IN) :: LABEL
      INTEGER, INTENT(IN)          :: NV
      INTEGER, INTENT(OUT)         :: X(NV)
      CALL RDLAB ( LABEL )
      READ (NDSOP,*,IOSTAT=IERR) X
      IF ( IERR .NE. 0 ) CALL SWRTER ( 'cannot read block '//LABEL )
    END SUBROUTINE RDINT
    !/
    SUBROUTINE RDREAL ( LABEL, X, NV )
      CHARACTER(LEN=*), INTENT(IN) :: LABEL
      INTEGER, INTENT(IN)          :: NV
      REAL, INTENT(OUT)            :: X(NV)
      CALL RDLAB ( LABEL )
      READ (NDSOP,*,IOSTAT=IERR) X
      IF ( IERR .NE. 0 ) CALL SWRTER ( 'cannot read block '//LABEL )
    END SUBROUTINE RDREAL
    !/
  END SUBROUTINE SWRTRD
  !/
  !/ End of module W3SNL6MD -------------------------------------------- /
  !/
END MODULE W3SNL6MD
