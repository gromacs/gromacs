      SUBROUTINE FQL77(N,X,D,E,NMAX)
C
C     MATRIX DIAGONALIZATION ROUTINE
C
C     CODE DERIVED FROM EISPACK
C     UPDATE : FMP 16.2.86 (CONVERSION TO FORTRAN 77)
C     UPDATE : FMP 09.3.89 (get rid of double precision)
C
C
C     INPUT
C
C     N           ORDER OF THE MATRIX, NOT LARGER THAN NMAX
C     X(N,N)      THE REAL SYMMETRIC MATRIX TO BE DIAGONALIZED AS A
C                 TWO-INDEXED SQUARE ARRAY OF WHICH ONLY THE LOWER LEFT
C                 TRIANGLE IS USED AND NEED BE SUPPLIED
C     NMAX        MAXIMUM VALUE OF N AS DECLARED IN THE CALLING ROUTINE
C                 DIMENSION STATEMENT
C     E           WORKING SPACE OF LENGTH NMAX
C
C     OUTPUT
C
C     D(N)        COMPUTED EIGENVALUES IN DESCENDING ORDER 
C     X(N,N)      EIGENVECTORS (IN COLUMNS) IN THE SAME ORDER AS THE
C                 EIGENVALUES
C

      IMPLICIT NONE
      INTEGER*4 N,NMAX,NI,I,J,K,L,J1
      REAL      X(NMAX,NMAX), D(NMAX), E(NMAX)
      REAL      EPS,TOL,F,G,H,S,B,P,R,C,ABSP
      EPS=7.E-14
      TOL=1.E-30
C
C     HANDLE SPECIAL CASE                                               
C
      IF (N.NE.1) GOTO 10
      D(1)=X(1,1)
      X(1,1)=1.0
      RETURN
C
C     HOUSEHOLDER REDUCTION TO TRIDIAGONAL FORM                         
C
   10 DO 150 NI=2,N
      I=N+2-NI
      L=I-2
      H=0.0                                                           
      G=X(I,I-1)
      IF (L.le.0)goto 140
      DO 30 K=1,L
   30 H=H+X(I,K)**2                                                     
      S=H+G*G                                                           
      IF(S.GE.TOL) GOTO 50
      H=0.0                                                           
      GOTO 140
   50 IF (H.le.0.0)goto 140
      L=L+1
      F=G
      G=SQRT(S)
      IF (F.le.0.0)goto 70
   70 G=-G
      H=S-F*G                                                           
      X(I,I-1)=F-G
      F=0.0
      DO 110 J=1,L
      X(J,I)=X(I,J)/H
      S=0.0
      DO 80 K=1,J
   80 S=S+X(J,K)*X(I,K)
      J1=J+1
      IF (J1.GT.L) GOTO 100
      DO 90 K=J1,L
   90 S=S+X(K,J)*X(I,K)
  100 E(J)=S/H
  110 F=F+S*X(J,I)
      F=F/(H+H)
      DO 120 J=1,L
  120 E(J)=E(J)-F*X(I,J)
      DO 130 J=1,L
      F=X(I,J)
      S=E(J)
      DO 130 K=1,J
  130 X(J,K)=X(J,K)-F*E(K)-X(I,K)*S
  140 D(I)=H
  150 E(I-1)=G
C
C     ACCUMULATION OF TRANSFORMATION MATRIX AND INTERMEDIATE D VECTOR
C
      D(1)=X(1,1)
      X(1,1)=1.0
      DO 220 I=2,N
      L=I-1
      IF (D(I).le.0.0)goto 200
      DO 190 J=1,L
      S=0.0
      DO 180 K=1,L
  180 S=S+X(I,K)*X(K,J)
      DO 190 K=1,L
  190 X(K,J)=X(K,J)-S*X(K,I)
  200 D(I)=X(I,I)
      X(I,I)=1.0
      DO 220 J=1,L
      X(I,J)=0.0
  220 X(J,I)=0.0
C
C     QL ITERATES
C
      B=0.0
      F=0.0
      E(N)=0.0
      DO 340 L=1,N
      H=EPS*(ABS(D(L))+ABS(E(L)))                                     
      IF (H.GT.B) B=H                                                   
      DO 240 J=L,N
      IF(ABS(E(J)).LE.B) GOTO 250
  240 CONTINUE
  250 IF(J.EQ.L) GOTO 340
  260 G=D(L)
      P=(D(L+1)-G)*0.5/E(L)
      R=SQRT(P*P+1.0)
      IF(P.ge.0.0)goto 280
      P=P-R
      GOTO 290
  280 P=P+R
  290 D(L)=E(L)/P
      H=G-D(L)                                                       
      K=L+1
      DO 300 I=K,N
  300 D(I)=D(I)-H                                                       
      F=F+H                                                             
      P=D(J)
      C=1.0
      S=0.0
      J1=J-1
      DO 330 NI=L,J1
      I=L+J1-NI
      G=C*E(I)
      H=C*P
      IF(ABS(P).LT.ABS(E(I))) GOTO 310
      C=E(I)/P
      R=SQRT(C*C+1.0)
      E(I+1)=S*P*R
      S=C/R
      C=1.0/R
      GOTO 320
  310 C=P/E(I)
      R=SQRT(C*C+1.0)
      E(I+1)=S*E(I)*R
      S=1.0/R
      C=C/R
  320 P=C*D(I)-S*G
      D(I+1)=H+S*(C*G+S*D(I))
      DO 330 K=1,N
      H=X(K,I+1)
      X(K,I+1)=X(K,I)*S+H*C
  330 X(K,I)=X(K,I)*C-H*S
      E(L)=S*P
      D(L)=C*P
      IF(ABS(E(L)).GT.B) GOTO 260
  340 D(L)=D(L)+F

C
C**** PUT EIGENVALUES AND EIGENVECTORS IN 
C**** DESIRED ASCENDING ORDER
C

      NI = N-1

      DO I = 1,NI

         K    = I
         P    = D(I)
         ABSP = ABS(D(I))
         J1   = I+1

         DO J = J1,N
            IF(ABS(D(J)) .LT. ABSP) THEN
               K    = J
               P    = D(J)
               ABSP = ABS(D(J))
            ENDIF
         ENDDO

         IF (K. NE. I) THEN
            D(K) = D(I)
            D(I) = P
            DO J = 1,N
               P      = X(J,I)
               X(J,I) = X(J,K)
               X(J,K) = P
            ENDDO
         ENDIF

      ENDDO

      RETURN

C
C     FMP
C     LAST BUT NOT LEAST I HAVE TO CONFESS THAT THERE IS AN ORIGINAL    
C     G. BINSCH REMARK ON THIS ROUTINE: 'QL IS PURPORTED TO BE ONE OF
C     THE FASTEST AND MOST COMPACT ROUTINES OF ITS KIND PRESENTLY KNOWN
C     TO MANKIND.'
C
      END


