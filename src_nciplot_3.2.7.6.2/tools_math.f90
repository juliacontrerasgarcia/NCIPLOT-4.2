! Copyright (c) 2013 Alberto Otero de la Roza <aoterodelaroza@ucmerced.edu>,
! Julia Conteras-Garcia <julia.contreras.garcia@gmail.com>, 
! Erin R. Johnson <ejohnson29@ucmerced.edu>, and Weitao Yang
! <weitao.yang@duke.edu>
!
! nciplot is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

! tools_math: Mathetical routines


module tools_math
  
  private

  public :: rs,nrspline,nrsplint

contains
  
  SUBROUTINE RS(NM,N,A,W,MATZ,Z,FV1,FV2,IERR)
!
      INTEGER N,NM,IERR,MATZ
      REAL*8 A(NM,N),W(N),Z(NM,N),FV1(N),FV2(N)
!
!     THIS SUBROUTINE CALLS THE RECOMMENDED SEQUENCE OF
!     SUBROUTINES FROM THE EIGENSYSTEM SUBROUTINE PACKAGE (EISPACK)
!     TO FIND THE EIGENVALUES AND EIGENVECTORS (IF DESIRED)
!     OF A REAL SYMMETRIC MATRIX.
!
!     ON INPUT:
!
!        NM  MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL
!        ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!        DIMENSION STATEMENT;
!
!        N  IS THE ORDER OF THE MATRIX  A;
!
!        A  CONTAINS THE REAL SYMMETRIC MATRIX;
!
!        MATZ  IS AN INTEGER VARIABLE SET EQUAL TO ZERO IF
!        ONLY EIGENVALUES ARE DESIRED;  OTHERWISE IT IS SET TO
!        ANY NON-ZERO INTEGER FOR BOTH EIGENVALUES AND EIGENVECTORS.
!
!     ON OUTPUT:
!
!        W  CONTAINS THE EIGENVALUES IN ASCENDING ORDER;
!
!        Z  CONTAINS THE EIGENVECTORS IF MATZ IS NOT ZERO;
!
!        IERR  IS AN INTEGER OUTPUT VARIABLE SET EQUAL TO AN
!        ERROR COMPLETION CODE DESCRIBED IN SECTION 2B OF THE
!        DOCUMENTATION.  THE NORMAL COMPLETION CODE IS ZERO;
!
!        FV1  AND  FV2  ARE TEMPORARY STORAGE ARRAYS.
!
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!
!     ------------------------------------------------------------------
!
      IF (N .LE. NM) GO TO 10
      IERR = 10 * N
      GO TO 50
!
   10 IF (MATZ .NE. 0) GO TO 20
!     :::::::::: FIND EIGENVALUES ONLY ::::::::::
      CALL  TRED1(NM,N,A,W,FV1,FV2)
      CALL  TQLRAT(N,W,FV2,IERR)
      GO TO 50
!     :::::::::: FIND BOTH EIGENVALUES AND EIGENVECTORS ::::::::::
   20 CALL  TRED2(NM,N,A,W,FV1,Z)
      CALL  TQL2(NM,N,W,FV1,Z,IERR)
   50 RETURN
!     :::::::::: LAST CARD OF RS ::::::::::
 END subroutine RS
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!
      SUBROUTINE TRED1(NM,N,A,D,E,E2)
!
      INTEGER I,J,K,L,N,II,NM,JP1
      REAL*8 A(NM,N),D(N),E(N),E2(N)
      REAL*8 F,G,H,SCALE
      REAL*8 DSQRT,DABS,DSIGN
!
!     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED1,
!     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
!
!     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX
!     TO A SYMMETRIC TRIDIAGONAL MATRIX USING
!     ORTHOGONAL SIMILARITY TRANSFORMATIONS.
!
!     ON INPUT:
!
!        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!          DIMENSION STATEMENT;
!
!        N IS THE ORDER OF THE MATRIX;
!
!        A CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE
!          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.
!
!     ON OUTPUT:
!
!        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANS-
!          FORMATIONS USED IN THE REDUCTION IN ITS STRICT LOWER
!          TRIANGLE.  THE FULL UPPER TRIANGLE OF A IS UNALTERED;
!
!        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX;
!
!        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
!          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO;
!
!        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
!          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.
!
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!
!     ------------------------------------------------------------------
!
      DO 100 I = 1, N
  100 D(I) = A(I,I)
!     :::::::::: FOR I=N STEP -1 UNTIL 1 DO -- ::::::::::
      DO 300 II = 1, N
         I = N + 1 - II
         L = I - 1
         H = 0.0D0
         SCALE = 0.0D0
         IF (L .LT. 1) GO TO 130
!     :::::::::: SCALE ROW (ALGOL TOL THEN NOT NEEDED) ::::::::::
         DO 120 K = 1, L
  120    SCALE = SCALE + DABS(A(I,K))
!
         IF (SCALE .NE. 0.0D0) GO TO 140
  130    E(I) = 0.0D0
         E2(I) = 0.0D0
         GO TO 290
!
  140    DO 150 K = 1, L
            A(I,K) = A(I,K) / SCALE
            H = H + A(I,K) * A(I,K)
  150    CONTINUE
!
         E2(I) = SCALE * SCALE * H
         F = A(I,L)
         G = -DSIGN(DSQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         A(I,L) = F - G
         IF (L .EQ. 1) GO TO 270
         F = 0.0D0
!
         DO 240 J = 1, L
            G = 0.0D0
!     :::::::::: FORM ELEMENT OF A*U ::::::::::
            DO 180 K = 1, J
  180       G = G + A(J,K) * A(I,K)
!
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
!
            DO 200 K = JP1, L
  200       G = G + A(K,J) * A(I,K)
!     :::::::::: FORM ELEMENT OF P ::::::::::
  220       E(J) = G / H
            F = F + E(J) * A(I,J)
  240    CONTINUE
!
         H = F / (H + H)
!     :::::::::: FORM REDUCED A ::::::::::
         DO 260 J = 1, L
            F = A(I,J)
            G = E(J) - H * F
            E(J) = G
!
            DO 260 K = 1, J
               A(J,K) = A(J,K) - F * E(K) - G * A(I,K)
  260    CONTINUE
!
  270    DO 280 K = 1, L
  280    A(I,K) = SCALE * A(I,K)
!
  290    H = D(I)
         D(I) = A(I,I)
         A(I,I) = H
  300 CONTINUE
!
      RETURN
!     :::::::::: LAST CARD OF TRED1 ::::::::::
  END subroutine tred1
!
!     ------------------------------------------------------------------
!
  SUBROUTINE TRED2(NM,N,A,D,E,Z)
!
      INTEGER I,J,K,L,N,II,NM,JP1
      REAL*8 A(NM,N),D(N),E(N),Z(NM,N)
      REAL*8 F,G,H,HH,SCALE
      REAL*8 DSQRT,DABS,DSIGN
!
!     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED2,
!     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
!
!     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX TO A
!     SYMMETRIC TRIDIAGONAL MATRIX USING AND ACCUMULATING
!     ORTHOGONAL SIMILARITY TRANSFORMATIONS.
!
!     ON INPUT:
!
!        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!          DIMENSION STATEMENT;
!
!        N IS THE ORDER OF THE MATRIX;
!
!        A CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE
!          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.
!
!     ON OUTPUT:
!
!        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX;
!
!        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
!          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO;
!
!        Z CONTAINS THE ORTHOGONAL TRANSFORMATION MATRIX
!          PRODUCED IN THE REDUCTION;
!
!        A AND Z MAY COINCIDE.  IF DISTINCT, A IS UNALTERED.
!
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!
!     ------------------------------------------------------------------
!
      DO 100 I = 1, N
!
         DO 100 J = 1, I
            Z(I,J) = A(I,J)
  100 CONTINUE
!
      IF (N .EQ. 1) GO TO 320
!     :::::::::: FOR I=N STEP -1 UNTIL 2 DO -- ::::::::::
      DO 300 II = 2, N
         I = N + 2 - II
         L = I - 1
         H = 0.0D0
         SCALE = 0.0D0
         IF (L .LT. 2) GO TO 130
!     :::::::::: SCALE ROW (ALGOL TOL THEN NOT NEEDED) ::::::::::
         DO 120 K = 1, L
  120    SCALE = SCALE + DABS(Z(I,K))
!
         IF (SCALE .NE. 0.0D0) GO TO 140
  130    E(I) = Z(I,L)
         GO TO 290
!
  140    DO 150 K = 1, L
            Z(I,K) = Z(I,K) / SCALE
            H = H + Z(I,K) * Z(I,K)
  150    CONTINUE
!
         F = Z(I,L)
         G = -DSIGN(DSQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         Z(I,L) = F - G
         F = 0.0D0
!
         DO 240 J = 1, L
            Z(J,I) = Z(I,J) / H
            G = 0.0D0
!     :::::::::: FORM ELEMENT OF A*U ::::::::::
            DO 180 K = 1, J
  180       G = G + Z(J,K) * Z(I,K)
!
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
!
            DO 200 K = JP1, L
  200       G = G + Z(K,J) * Z(I,K)
!     :::::::::: FORM ELEMENT OF P ::::::::::
  220       E(J) = G / H
            F = F + E(J) * Z(I,J)
  240    CONTINUE
!
         HH = F / (H + H)
!     :::::::::: FORM REDUCED A ::::::::::
         DO 260 J = 1, L
            F = Z(I,J)
            G = E(J) - HH * F
            E(J) = G
!
            DO 260 K = 1, J
               Z(J,K) = Z(J,K) - F * E(K) - G * Z(I,K)
  260    CONTINUE
!
  290    D(I) = H
  300 CONTINUE
!
  320 D(1) = 0.0D0
      E(1) = 0.0D0
!     :::::::::: ACCUMULATION OF TRANSFORMATION MATRICES ::::::::::
      DO 500 I = 1, N
         L = I - 1
         IF (D(I) .EQ. 0.0D0) GO TO 380
!
         DO 360 J = 1, L
            G = 0.0D0
!
            DO 340 K = 1, L
  340       G = G + Z(I,K) * Z(K,J)
!
            DO 360 K = 1, L
               Z(K,J) = Z(K,J) - G * Z(K,I)
  360    CONTINUE
!
  380    D(I) = Z(I,I)
         Z(I,I) = 1.0D0
         IF (L .LT. 1) GO TO 500
!
         DO 400 J = 1, L
            Z(I,J) = 0.0D0
            Z(J,I) = 0.0D0
  400    CONTINUE
!
  500 CONTINUE
!
      RETURN
!     :::::::::: LAST CARD OF TRED2 ::::::::::
  END subroutine tred2
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!
  SUBROUTINE TQLRAT(N,D,E2,IERR)
!
      INTEGER I,J,L,M,N,II,L1,MML,IERR
      REAL*8 D(N),E2(N)
      REAL*8 B,C,F,G,H,P,R,S,MACHEP
      REAL*8 DSQRT,DABS,DSIGN
!
!     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQLRAT,
!     ALGORITHM 464, COMM. ACM 16, 689(1973) BY REINSCH.
!
!     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC
!     TRIDIAGONAL MATRIX BY THE RATIONAL QL METHOD.
!
!     ON INPUT:
!
!        N IS THE ORDER OF THE MATRIX;
!
!        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX;
!
!        E2 CONTAINS THE SQUARES OF THE SUBDIAGONAL ELEMENTS OF THE
!          INPUT MATRIX IN ITS LAST N-1 POSITIONS.  E2(1) IS ARBITRARY.
!
!      ON OUTPUT:
!
!        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
!          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND
!          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE
!          THE SMALLEST EIGENVALUES;
!
!        E2 HAS BEEN DESTROYED;
!
!        IERR IS SET TO
!          ZERO       FOR NORMAL RETURN,
!          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
!                     DETERMINED AFTER 30 ITERATIONS.
!
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!
!     ------------------------------------------------------------------
!
!     :::::::::: MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
!                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
      DATA MACHEP/3.33D-16/
!
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
!
      DO 100 I = 2, N
  100 E2(I-1) = E2(I)
!
      F = 0.0D0
      B = 0.0D0
      E2(N) = 0.0D0
!
      DO 290 L = 1, N
         J = 0
         H = MACHEP * (DABS(D(L)) + DSQRT(E2(L)))
         IF (B .GT. H) GO TO 105
         B = H
         C = B * B
!     :::::::::: LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT ::::::::::
  105    DO 110 M = L, N
            IF (E2(M) .LE. C) GO TO 120
!     :::::::::: E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
!                THROUGH THE BOTTOM OF THE LOOP ::::::::::
  110    CONTINUE
!
  120    IF (M .EQ. L) GO TO 210
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
!     :::::::::: FORM SHIFT ::::::::::
         L1 = L + 1
         S = DSQRT(E2(L))
         G = D(L)
         P = (D(L1) - G) / (2.0D0 * S)
         R = DSQRT(P*P+1.0D0)
         D(L) = S / (P + DSIGN(R,P))
         H = G - D(L)
!
         DO 140 I = L1, N
  140    D(I) = D(I) - H
!
         F = F + H
!     :::::::::: RATIONAL QL TRANSFORMATION ::::::::::
         G = D(M)
         IF (G .EQ. 0.0D0) G = B
         H = G
         S = 0.0D0
         MML = M - L
!     :::::::::: FOR I=M-1 STEP -1 UNTIL L DO -- ::::::::::
         DO 200 II = 1, MML
            I = M - II
            P = G * H
            R = P + E2(I)
            E2(I+1) = S * R
            S = E2(I) / R
            D(I+1) = H + S * (H + D(I))
            G = D(I) - E2(I) / G
            IF (G .EQ. 0.0D0) G = B
            H = G * P / R
  200    CONTINUE
!
         E2(L) = S * G
         D(L) = H
!     :::::::::: GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST ::::::::::
         IF (H .EQ. 0.0D0) GO TO 210
         IF (DABS(E2(L)) .LE. DABS(C/H)) GO TO 210
         E2(L) = H * E2(L)
         IF (E2(L) .NE. 0.0D0) GO TO 130
  210    P = D(L) + F
!     :::::::::: ORDER EIGENVALUES ::::::::::
         IF (L .EQ. 1) GO TO 250
!     :::::::::: FOR I=L STEP -1 UNTIL 2 DO -- ::::::::::
         DO 230 II = 2, L
            I = L + 2 - II
            IF (P .GE. D(I-1)) GO TO 270
            D(I) = D(I-1)
  230    CONTINUE
!
  250    I = 1
  270    D(I) = P
  290 CONTINUE
!
      GO TO 1001
!     :::::::::: SET ERROR -- NO CONVERGENCE TO AN
!                EIGENVALUE AFTER 30 ITERATIONS ::::::::::
 1000 IERR = L
 1001 RETURN
!     :::::::::: LAST CARD OF TQLRAT ::::::::::
  END subroutine tqlrat
!
!     ------------------------------------------------------------------
!
  SUBROUTINE TQL2(NM,N,D,E,Z,IERR)
!
      INTEGER I,J,K,L,M,N,II,L1,NM,MML,IERR
      REAL*8 D(N),E(N),Z(NM,N)
      REAL*8 B,C,F,G,H,P,R,S,MACHEP
      REAL*8 DSQRT,DABS,DSIGN
!
!     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQL2,
!     NUM. MATH. 11, 293-306(1968) BY BOWDLER, MARTIN, REINSCH, AND
!     WILKINSON.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).
!
!     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
!     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE QL METHOD.
!     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO
!     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS
!     FULL MATRIX TO TRIDIAGONAL FORM.
!
!     ON INPUT:
!
!        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!          DIMENSION STATEMENT;
!
!        N IS THE ORDER OF THE MATRIX;
!
!        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX;
!
!        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
!          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY;
!
!        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
!          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS
!          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN
!          THE IDENTITY MATRIX.
!
!      ON OUTPUT:
!
!        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
!          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT
!          UNORDERED FOR INDICES 1,2,...,IERR-1;
!
!        E HAS BEEN DESTROYED;
!
!        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC
!          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,
!          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED
!          EIGENVALUES;
!
!        IERR IS SET TO
!          ZERO       FOR NORMAL RETURN,
!          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
!                     DETERMINED AFTER 30 ITERATIONS.
!
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!
!     ------------------------------------------------------------------
!
!     :::::::::: MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
!                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
      DATA MACHEP/3.33D-16/
!
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
!
      DO 100 I = 2, N
  100 E(I-1) = E(I)
!
      F = 0.0D0
      B = 0.0D0
      E(N) = 0.0D0
!
      DO 240 L = 1, N
         J = 0
         H = MACHEP * (DABS(D(L)) + DABS(E(L)))
         IF (B .LT. H) B = H
!     :::::::::: LOOK FOR SMALL SUB-DIAGONAL ELEMENT ::::::::::
         DO 110 M = L, N
            IF (DABS(E(M)) .LE. B) GO TO 120
!     :::::::::: E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
!                THROUGH THE BOTTOM OF THE LOOP ::::::::::
  110    CONTINUE
!
  120    IF (M .EQ. L) GO TO 220
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
!     :::::::::: FORM SHIFT ::::::::::
         L1 = L + 1
         G = D(L)
         P = (D(L1) - G) / (2.0D0 * E(L))
         R = DSQRT(P*P+1.0D0)
         D(L) = E(L) / (P + DSIGN(R,P))
         H = G - D(L)
!
         DO 140 I = L1, N
  140    D(I) = D(I) - H
!
         F = F + H
!     :::::::::: QL TRANSFORMATION ::::::::::
         P = D(M)
         C = 1.0D0
         S = 0.0D0
         MML = M - L
!     :::::::::: FOR I=M-1 STEP -1 UNTIL L DO -- ::::::::::
         DO 200 II = 1, MML
            I = M - II
            G = C * E(I)
            H = C * P
            IF (DABS(P) .LT. DABS(E(I))) GO TO 150
            C = E(I) / P
            R = DSQRT(C*C+1.0D0)
            E(I+1) = S * P * R
            S = C / R
            C = 1.0D0 / R
            GO TO 160
  150       C = P / E(I)
            R = DSQRT(C*C+1.0D0)
            E(I+1) = S * E(I) * R
            S = 1.0D0 / R
            C = C * S
  160       P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
!     :::::::::: FORM VECTOR ::::::::::
            DO 180 K = 1, N
               H = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * H
               Z(K,I) = C * Z(K,I) - S * H
  180       CONTINUE
!
  200    CONTINUE
!
         E(L) = S * P
         D(L) = C * P
         IF (DABS(E(L)) .GT. B) GO TO 130
  220    D(L) = D(L) + F
  240 CONTINUE
!     :::::::::: ORDER EIGENVALUES AND EIGENVECTORS ::::::::::
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
!
         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE
!
         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P
!
         DO 280 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  280    CONTINUE
!
  300 CONTINUE
!
      GO TO 1001
!     :::::::::: SET ERROR -- NO CONVERGENCE TO AN
!                EIGENVALUE AFTER 30 ITERATIONS ::::::::::
 1000 IERR = L
 1001 RETURN
!     :::::::::: LAST CARD OF TQL2 ::::::::::
  END subroutine tql2

!     ------------------------------------------------------------------

  SUBROUTINE NRSPLINE(X,Y,N,YP1,YPN,Y2)

! CUBIC SPLINE INTERPOLATION TAKEN FROM NUMERICAL RECIPLE

        INTEGER N,I,K
        PARAMETER(NMAX=500)
        REAL*8 YP1,YPN,X(N),Y(N), Y2(N),P,QN,SIG,UN,U(NMAX)
        IF(YP1.GT..99E30) THEN
          Y2(1)=0D0
          U(1)=0D0
        ELSE
          Y2(1)=-0.5D0
          U(1)=(3D0/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
        ENDIF
        DO I=2,N-1
           SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
           P=SIG*Y2(I-1)+2
           Y2(I)=(SIG-1D0)/P
           U(I)=(6D0*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1)) &
           /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
        ENDDO
        IF(YPN.GT..99E30) THEN
          QN=0D0
          UN=0D0
        ELSE
          QN=0.5D0
          UN=(3D0/(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
        ENDIF
        Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1D0)
        DO K=N-1,1,-1
           Y2(K)=Y2(K+1)+U(K)
        ENDDO
        RETURN
        END subroutine nrspline

!     ------------------------------------------------------------------

        SUBROUTINE NRSPLINT(XA,YA,Y2A,N,X,Y)

! CUBIC SPLINE INTERPOLATION TAKEN FROM NUMERICAL RECIPLE

        INTEGER N,K,KHI,KLO
        REAL*8 X,Y,XA(N),YA(N),Y2A(N),A,B,H
        KLO=1
        KHI=N
1       IF(KHI-KLO.GT.1) THEN
          K=(KHI+KLO)/2
          IF(XA(K).GT.X) THEN
            KHI=K
          ELSE
            KLO=K
          ENDIF
        GOTO 1
        ENDIF
        H=XA(KHI)-XA(KLO)
        A=(XA(KHI)-X)/H
        B=(X-XA(KLO))/H
        Y=A*YA(KLO)+B*YA(KHI)+((A**3-A)*Y2A(KLO) &
     +    (B**3-B)*Y2A(KHI))*(H**2)/6D0
        RETURN
        END subroutine nrsplint

end module tools_math
