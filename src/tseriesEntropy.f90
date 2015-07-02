! ## ***************************************************************************************************
!  #  The author of this software is
!  #
!  #  Simone Giannerini, Copyright (c) 2007 -
!  #
!  #  Permission to use, copy, modify, and distribute this software for any
!  #  purpose without fee is hereby granted, provided that this entire notice
!  #  is included in all copies of any software which is or includes a copy
!  #  or modification of this software and in all copies of the supporting
!  #  documentation for such software.
!  #
!  #  This program is free software; you can redistribute it and/or modify
!  #  it under the terms of the GNU General Public License as published by
!  #  the Free Software Foundation; either version 2 of the License, or
!  #  (at your option) any later version.
!  #
!  #  This program is distributed in the hope that it will be useful,
!  #  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  #  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  #  GNU General Public License for more details.
!  #
!  #  A copy of the GNU General Public License is available at
!  #  http://www.r-project.org/Licenses/
!
! ## ***************************************************************************************************

!! *************************************************************************************************
!! SUBROUTINE LIST FOR USE WITH THE R PACKAGE tseriesEntropy
!! Simone Giannerini 2007 -
!! *************************************************************************************************
!! *************************************************************************************************
!! SSUNI   : COMPUTES SRHO ON INTEGER/CATEGORICAL SERIES
!! SSUNI2  : THE SAME AS SSUNI BUT ASSUMING STATIONARITY (estimates univariate frequencies on the whole data set)
!! SSUNIB  : SSUNI AND SSUNI2 -- PERMUTATION VERSION
!! SSBIV   : COMPUTES SRHO ON BIVARIATE INTEGER/CATEGORICAL SERIES (CROSS ENTROPY)
!! SSBIV2  : THE SAME AS SSBIV BUT ASSUMING STATIONARITY
!! SSBIVB  : SSBIV AND SSBIV2-- PERMUTATION VERSION

!! ************* INTERNALS THAT DO NOT WORK IN R ***************************************************

!! TABF    : ONE WAY FREQUENCY TABLE FOR CATEGORICAL DATA
!! TABFD   : TWO WAY FREQUENCY TABLE FOR CATEGORICAL DATA
!! TABFD2  : AS TABFD BUT WITH MARGINS GIVEN AS INPUT
!! [ ... ]
!! *************************************************************************************************
!! *************************************************************************************************
!! *************************************************************************************************

MODULE SHARED_DATA
    IMPLICIT NONE
    REAL(8),PARAMETER :: M_PI         = 3.141592653589793238462643383280
    REAL(8),PARAMETER :: M_1_SQRT_2PI = 0.398942280401432677939946059934 ! 1/(sqrt(2*pi))

    CONTAINS
!! ****************************************************************************
    FUNCTION PERM(N)
        ! GIVES A RANDOM PERMUTATION OF THE FIRST N INTEGERS

        INTEGER,INTENT(IN) :: N
        INTEGER,DIMENSION(N) :: PERM
        INTEGER,DIMENSION(N) :: B
        REAL(8) :: u(N)
        INTEGER :: i,j,k

        B   =  (/(i,i=1,N)/)
        PERM = 0  ! vettore degli n numeri estratti
        CALL RANDOM_NUMBER(u)
        do i = N,1,-1
            j = int(u(i)*i)+1   ! numero casuale compreso fra 1 e i
            PERM(i) = B(j)
            k = B(j)
            B(j)=B(i)
            B(i)=k
        enddo
    END FUNCTION PERM
!! *************************************************************************************************

SUBROUTINE DNORMF(X,N,OUT)
    IMPLICIT NONE
    INTEGER:: N
    REAL(8) :: X(N),OUT(N)
    OUT = EXP(-.5*X**2)*M_1_SQRT_2PI
END SUBROUTINE DNORMF
!! *************************************************************************************************

SUBROUTINE SRHOBIVA(TX,TY,NUNI,NX,NY,T,NBIV,S,nor)
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NX,NY,NUNI,NBIV,TX(NX,2),TY(NY,2),T(NX,NY),nor
    REAL(8),INTENT(OUT):: S

    REAL(8)  :: PX(NX),PY(NY),P(NX,NY),smax
    INTEGER :: ix,iy
    S = 0.
    PX = DBLE(TX(:,2))/NUNI
    PY = DBLE(TY(:,2))/NUNI
    P  = DBLE(T)/NBIV
    DO ix = 1,nx
        DO iy = 1,ny
            S = S + (DSQRT(P(ix,iy))-DSQRT(PX(ix)*PY(iy)))**2;
        ENDDO
    ENDDO
    S = S/2
    IF(nor>0) THEN
        smax = MAX(1 - SUM(PX**1.5), 1 - SUM(PY**1.5))
        S=S/smax
    ENDIF
!    call dblepr("S", -1, S, 1)
!    call intpr("nor_SB", -1, nor, 1)
END SUBROUTINE SRHOBIVA
!! *************************************************************************************************

SUBROUTINE TABFD(X,Y,N,TX,TY,T)
!!DEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:'tabfd' :: tabfd
!! ******************************************************************************
!! TWO WAY frequency table for integer data
!! ******************************************************************************
!! ** INPUT *********************************************************************
!! X    : integer/categorical series
!! Y    : integer/categorical series
!! N    : length of X and Y
!! ** OUTPUT ********************************************************************
!! TX   : matrix with 2 columns containing values and counts for X (rows)
!! TY   : matrix with 2 columns containing values and counts for Y (columns)
!! T    : matrix containing two-way counts for X and Y
!! ******************************************************************************
!! Simone Giannerini April 2007
!! ******************************************************************************
    IMPLICIT NONE

    INTEGER,INTENT(IN):: N,X(N),Y(N)
    INTEGER:: nx,ny,ix,iy
    INTEGER,allocatable,INTENT(OUT):: TX(:,:),TY(:,:),T(:,:)
    CALL TABF(X,N,TX)
    CALL TABF(Y,N,TY)
    nx = SIZE(TX,1)
    ny = SIZE(TY,1)
    ALLOCATE(T(nx,ny))
    DO ix=1,nx
        DO iy=1,ny
            T(ix,iy)=COUNT(X==TX(ix,1).AND.Y==TY(iy,1))
        ENDDO
    ENDDO
END SUBROUTINE TABFD
!! *************************************************************************************************
SUBROUTINE TABFD2(X,Y,N,TX,TY,NX,NY,T)

!! ******************************************************************************
!! TWO WAY frequency table for integer data **** MARGINS ARE TAKEN AS INPUT *****
!! ******************************************************************************
!! ** INPUT *********************************************************************
!! X    : integer/categorical series
!! Y    : integer/categorical series
!! N    : length of X and Y
!! TX   : matrix with 2 columns containing values and counts for X (rows)
!! TY   : matrix with 2 columns containing values and counts for Y (columns)
!! NX   : number of distinct values of X: NX = size(TX,1)
!! NY   : number of distinct values of Y: NY = size(TY,1)
!! ** OUTPUT ********************************************************************
!! T    : matrix containing two-way counts for X and Y
!! ******************************************************************************
!! Simone Giannerini April 2007
!! ******************************************************************************

    IMPLICIT NONE
    INTEGER,INTENT(IN):: N,NX,NY,X(N),Y(N),TX(NX,2),TY(NY,2)
    INTEGER:: ix,iy
    INTEGER,INTENT(OUT):: T(NX,NY)
    DO ix=1,NX
        DO iy=1,NY
            T(ix,iy)=COUNT(X==TX(ix,1).AND.Y==TY(iy,1))
        ENDDO
    ENDDO
END SUBROUTINE TABFD2
!! *************************************************************************************************

SUBROUTINE TABF(X,N,T)
!! ******************************************************************************
!! One dimensional frequency table for integer data
!! ******************************************************************************
!! ** INPUT *********************************************************************
!! X    : integer/categorical series
!! N    : length of X
!! ** OUTPUT ********************************************************************
!! T    : matrix with 2 columns containing values and counts
!! ******************************************************************************
!! Simone Giannerini April 2007
!! ******************************************************************************

    IMPLICIT NONE

    INTEGER,INTENT(IN):: N,X(N)
    INTEGER:: DUM(N,2),YES(N),IND(N),i,cnt
    INTEGER,allocatable,INTENT(OUT):: T(:,:)
    DUM=-9999
    YES=0
    cnt=0
    DO i=1,N
        IF(YES(i)==0) THEN
            cnt=cnt+1
            IND=0
            DUM(cnt,1) = X(i)
            DUM(cnt,2) = COUNT(X==X(i))
            WHERE(X==X(i)) IND=(ABS(X)+1)/(ABS(X)+1)
            YES=YES+IND
        ENDIF
    ENDDO
    IF(cnt>0) THEN
        ALLOCATE(T(cnt,2))
        T(:,:)=DUM(1:cnt,:)
    ENDIF
END SUBROUTINE TABF
!! *************************************************************************************************
END MODULE SHARED_DATA
!! *************************************************************************************************

SUBROUTINE MBBOOT(X,N,B,l,M)

    IMPLICIT NONE
    INTEGER,INTENT(IN):: N,B,l
    REAL(8),INTENT(IN):: X(N)
    REAL(8),INTENT(OUT):: M(N,B)

    INTEGER:: nblocks,ind(N),indmat(N-l+1,l),i
    INTEGER,allocatable:: indblock(:),indx(:,:)
    REAL(8),allocatable:: M2(:,:)
    nblocks = int(N/l) + 1
    ALLOCATE(indblock(nblocks*B),indx(nblocks*l*B,1),M2(nblocks*l, B))
    ind      = (/(i,i=1,N)/)
    indmat   = EMBED(ind,1,l) ! embedding of the indices
    indblock = BOOT(N-l+1,nblocks*B) ! indices of the blocks
    indx     = RESHAPE(TRANSPOSE(indmat(indblock,:)), (/nblocks*l*B, 1/)) ! vectorized indices of the selected blocks
    M2       = RESHAPE(X(indx(:,1)),(/nblocks*l, B/))
    M        = M2(1:N,:)
CONTAINS
    FUNCTION BOOT(N,size)
        ! GIVES A RESAMPLE WITH REPLACEMENT OF LENGTH size OF THE FIRST N INTEGERS
        INTEGER,INTENT(IN) :: N,size
        INTEGER,DIMENSION(size) :: BOOT
        REAL(8) :: u(size)
        CALL init_random_seed()
        CALL RANDOM_NUMBER(u)
        BOOT = int(u*N)+1   ! numero casuale compreso fra 1 e N
    END FUNCTION BOOT

    SUBROUTINE init_random_seed()
        INTEGER :: i, n, clock
        INTEGER, DIMENSION(:), ALLOCATABLE :: seed

        CALL RANDOM_SEED(size = n)
        ALLOCATE(seed(n))

        CALL SYSTEM_CLOCK(COUNT=clock)

        seed = clock + 37 * (/ (i - 1, i = 1, n) /)
        CALL RANDOM_SEED(PUT = seed)

        DEALLOCATE(seed)
    END SUBROUTINE

    FUNCTION EMBED(X,j,d)

    !	TIME DELAY EMBEDDING
    !        X : time series
    !        j : time delay
    !        d : embedding dimension
        implicit none
        integer,intent(IN):: j,d
        integer,dimension(:),intent(IN) :: X
        integer,dimension(size(X,1)-(d-1)*j,d):: EMBED
     	
        integer :: npun,h
        npun = size(X,1)-(d-1)*j
        do h=0,d-1
            EMBED(:,h+1)=X(1+h*j:npun+h*j)
        enddo

    END FUNCTION EMBED

END SUBROUTINE MBBOOT

!! *************************************************************************************************
SUBROUTINE SSBIV(X,Y,N,nlag,S,nor)
!!DEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:'ssbiv' :: ssbiv

!! **************************************************************************
!! Implementation of the entropy-based dependence measure
!! proposed by Granger et al. (2004) Journal of Time Series Analysis.
!!
!! Deals with integer/categorical time series
!!
!! BIVARIATE VERSION
!! It is the extension of Srho to the case of two series, includes the cross-entropy
!!
!! INPUT:   X,Y     : integer/categorical eries of length N
!!          nlag    : number of lags computed
!!          nor     : 1 - normalized version
!! OUTPUT:  Srho(k) where k = -nlag,...,0,...,nlag
!!          i.e. Srho[1]        contains Srho between X_{t} and Y_{t-nlag}
!!               Srho[nlag+1]   contains Srho between X_{t} and Y_{t}
!!               Srho[2*nlag+1] contains Srho between X_{t} and Y_{t+nlag}
!! **************************************************************************
!! Simone Giannerini DECEMBER 2007
!!
    USE SHARED_DATA
    IMPLICIT NONE
    INTEGER,intent(IN):: N,nlag,X(N),Y(N),nor
    REAL(8),intent(OUT):: S(2*nlag+1)

    INTEGER,allocatable:: TX(:,:),TY(:,:),T(:,:)
    REAL(8)  :: dum
    INTEGER :: k,nx,ny

    S = 999
    CALL TABFD(X,Y,N,TX,TY,T)
        nx = SIZE(TX,1)
        ny = SIZE(TY,1)
        CALL SRHOBIVA(TX,TY,N,NX,NY,T,N,dum,nor)
        S(nlag+1) = dum;
    DO k = 1,nlag
        CALL TABFD(X(1:(N-k)),Y((k+1):N),N-k,TX,TY,T)
        nx = SIZE(TX,1)
        ny = SIZE(TY,1)

        CALL SRHOBIVA(TX,TY,N-k,NX,NY,T,N-k,dum,nor)
        S(nlag+1+k) =   dum;

        CALL TABFD(X((k+1):N),Y(1:(N-k)),N-k,TX,TY,T)
        nx = SIZE(TX,1)
        ny = SIZE(TY,1)
        CALL SRHOBIVA(TX,TY,N-k,NX,NY,T,N-k,dum,nor)
        S(nlag+1-k) =   dum;
    ENDDO
END SUBROUTINE SSBIV

!! *************************************************************************************************
SUBROUTINE SSBIVB(X,Y,N,nlag,B,S,M,STAT,nor)
!!DEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:'ssbivb' :: ssbivb

!! ******************************************************************************
!! PERMUTATION VERSION OF SSBIV
!! ** INPUT *********************************************************************
!! X    : integer/categorical series
!! Y    : integer/categorical series
!! N    : length of X
!! nlag : number of lags computed
!! B    : number of bootstrap replications
!! STAT : 1 - Stationary version 0 - Non-stationary version
!! nor  : 1 - normalized version
!! ** OUTPUT ********************************************************************
!! S    : Srho computed on the original series
!! M    : nlag by B matrix containing Srho computed on the replications
!! ******************************************************************************
!! Simone Giannerini April 2007
!! ******************************************************************************
    USE SHARED_DATA
    IMPLICIT NONE

INTERFACE
    SUBROUTINE SSBIV(X,Y,N,nlag,S,nor)
        INTEGER,INTENT(IN):: N,nlag,X(N),Y(N),nor
        REAL(8),intent(OUT):: S(2*nlag+1)
    END SUBROUTINE SSBIV
    SUBROUTINE SSBIV2(X,Y,N,nlag,S,nor)
        INTEGER,INTENT(IN):: N,nlag,X(N),Y(N),nor
        REAL(8),intent(OUT):: S(2*nlag+1)
    END SUBROUTINE SSBIV2
END INTERFACE

    INTEGER,INTENT(IN):: N,nlag,X(N),Y(N),B,STAT,nor
    REAL(8),intent(OUT):: S(2*nlag+1),M(2*nlag+1,B)

    REAL(8) :: dum(2*nlag+1)
    INTEGER:: ind(N),XB(N),YB(N)
    INTEGER:: i

    S    = 0.;
    M    = 0.;
    dum  = 0.;
    ind  = 0 ;
    XB   = 0 ;
    YB   = 0 ;
!    CALL intpr ("nor", -1, nor, 1)
    IF (STAT==0) THEN
        CALL SSBIV(X,Y,N,nlag,S,nor);
        do i = 1,B
            ind = PERM(N)
            XB  = X(ind)
            ind = PERM(N)
            YB  = Y(ind)
            CALL SSBIV(XB,YB,N,nlag,dum,nor);
            M(:,i) = dum;
        enddo
    ELSE
        CALL SSBIV2(X,Y,N,nlag,S,nor);
        do i = 1,B
            ind = PERM(N)
            XB  = X(ind)
            ind = PERM(N)
            YB  = Y(ind)
            CALL SSBIV2(XB,YB,N,nlag,dum,nor);
            M(:,i) = dum;
        enddo
    ENDIF

END SUBROUTINE SSBIVB
!! *************************************************************************************************

SUBROUTINE SSBIV2(X,Y,N,nlag,S,nor)
!!DEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:'ssbiv2' :: ssbiv2
!! **************************************************************************
!! Implementation of the entropy-based dependence measure
!! proposed by Granger et al. (2004) Journal of Time Series Analysis.
!!
!! Deals with integer/categorical time series
!!
!! BIVARIATE VERSION ** SAME AS SSBIV BUT ASSUMES STATIONARITY **************
!! ** ESTIMATES PROBABILITIES WITH RELATIVE FREQUENCIES ON THE WHOLE SERIES *
!! It is the extension of Srho to the case of two series, includes the cross-entropy
!!
!! INPUT:   X,Y     : integer/categorical eries of length N
!!          nlag    : number of lags computed
!!          nor     : 1 - normalized version
!! OUTPUT:  Srho(k) where k = -nlag,...,0,...,nlag
!!          i.e. Srho[1]        contains Srho between X_{t} and Y_{t-nlag}
!!               Srho[nlag+1]   contains Srho between X_{t} and Y_{t}
!!               Srho[2*nlag+1] contains Srho between X_{t} and Y_{t+nlag}
!! **************************************************************************
!! Simone Giannerini DECEMBER 2007
!!
    USE SHARED_DATA
    IMPLICIT NONE
    INTEGER,intent(IN):: N,nlag,X(N),Y(N),nor
    REAL(8),intent(OUT)::S(2*nlag+1)

    INTEGER,allocatable:: TX(:,:),TY(:,:),T(:,:)
    REAL(8)  :: dum
    INTEGER :: k,NX,NY

    S = 999
    CALL TABF(X,N,TX)
    CALL TABF(Y,N,TY)
    NX = SIZE(TX,1)
    NY = SIZE(TY,1)
    ALLOCATE(T(NX,NY))
    T = -999
    CALL TABFD2(X,Y,N,TX,TY,NX,NY,T)
    CALL SRHOBIVA(TX,TY,N,NX,NY,T,N,dum,nor)
        S(nlag+1) = dum;
    DO k = 1,nlag
        CALL TABFD2(X(1:(N-k)),Y((k+1):N),N-k,TX,TY,NX,NY,T)
        CALL SRHOBIVA(TX,TY,N,NX,NY,T,N-k,dum,nor)
        S(nlag+1+k) =   dum;

        CALL TABFD2(X((k+1):N),Y(1:(N-k)),N-k,TX,TY,NX,NY,T)
        CALL SRHOBIVA(TX,TY,N,NX,NY,T,N-k,dum,nor)
        S(nlag+1-k) =   dum;
    ENDDO
END SUBROUTINE SSBIV2

!! *************************************************************************************************

SUBROUTINE SSUNIB(X,N,nlag,B,S,M,STAT,nor)
!!DEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:'ssunib' :: ssunib
!! ******************************************************************************
!! PERMUTATION VERSION OF SSUNI
!! ** INPUT *********************************************************************
!! X    : integer/categorical series
!! N    : length of X
!! nlag : number of lags computed
!! B    : number of bootstrap replications
!! STAT : 1 - Stationary version 0 - Non-stationary version
!! nor  : 1 - normalized version
!! ** OUTPUT ********************************************************************
!! S    : Srho computed on the original series
!! M    : nlag by B matrix containing Srho computed on the replications
!! ******************************************************************************
!! Simone Giannerini April 2007
!! ******************************************************************************

    USE SHARED_DATA
    IMPLICIT NONE
INTERFACE
    SUBROUTINE SSUNI(X,N,nlag,S,nor)
    INTEGER,intent(IN):: N,nlag,X(N),nor
    REAL(8),intent(OUT)::S(nlag)
    END SUBROUTINE SSUNI

    SUBROUTINE SSUNI2(X,N,nlag,S,nor)
    INTEGER,intent(IN):: N,nlag,X(N),nor
    REAL(8),intent(OUT)::S(nlag)
    END SUBROUTINE SSUNI2
END INTERFACE

    INTEGER,intent(IN):: N,nlag,X(N),B,STAT,nor
    REAL(8),intent(OUT):: S(nlag),M(nlag,B)
    REAL(8):: dum(nlag)
    INTEGER:: ind(N),XB(N)
    INTEGER:: i

    S    = 0.;
    M    = 0.;
    dum  = 0.;
    ind  = 0 ;
    XB   = 0 ;
    IF (STAT==0) THEN
        CALL SSUNI(X,N,nlag,S,nor);
        do i = 1,B
            ind = PERM(N)
            XB  = X(ind)
            CALL SSUNI(XB,N,nlag,dum,nor);
            M(:,i) = dum;
        enddo
    ELSE
        CALL SSUNI2(X,N,nlag,S,nor);
        do i = 1,B
            ind = PERM(N)
            XB  = X(ind)
            CALL SSUNI2(XB,N,nlag,dum,nor);
            M(:,i) = dum;
        enddo
    ENDIF
END SUBROUTINE SSUNIB

!! *************************************************************************************************

SUBROUTINE SSUNI(X,N,nlag,S,nor)
!!DEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:'ssuni' :: ssuni
!! **************************************************************************
!! Implementation of the entropy-based dependence measure
!! proposed by Granger et al. (2004) Journal of Time Series Analysis.
!!
!! Deals with integer/categorical time series
!! **************************************************************************
!! Simone Giannerini March 2007
!!
    USE SHARED_DATA
    IMPLICIT NONE
    INTEGER,intent(IN):: N,nlag,X(N),nor
    REAL(8),intent(OUT):: S(nlag)

    INTEGER,allocatable:: TX(:,:),TY(:,:),T(:,:)
    REAL(8)  :: dum
    INTEGER :: k,nx,ny

    S = 0
    DO k = 1,nlag
        CALL TABFD(X(1:(N-k)),X((k+1):N),N-k,TX,TY,T)
        nx   = SIZE(TX,1)
        ny   = SIZE(TY,1)
        CALL SRHOBIVA(TX,TY,N-k,nx,ny,T,N-k,dum,nor)
        S(k) = dum
    ENDDO
END SUBROUTINE SSUNI
!! *************************************************************************************************

SUBROUTINE SSUNI2(X,N,nlag,S,nor)
!DEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:'ssuni2' :: ssuni2
!! **************************************************************************
!! Implementation of the entropy-based dependence measure
!! proposed by Granger et al. (2004) Journal of Time Series Analysis.
!!
!! Deals with integer/categorical time series
!! **************************************************************************
!! DIFFERS FROM SSUNI BECAUSE STATIONARITY IS ASSUMED
!! MARGINAL DISTRIBUTIONS ARE COMPUTED ONLY ONCE ON X(N)
!! **************************************************************************
!! Simone Giannerini 2007
!!
    USE SHARED_DATA
    IMPLICIT NONE

    INTEGER,intent(IN):: N,nlag,X(N),nor
    REAL(8),intent(OUT):: S(nlag)

    INTEGER,allocatable:: TX(:,:),T(:,:)
    REAL(8) :: dum,smax
    INTEGER :: k,NX

    CALL TABF(X,N,TX)
    NX = SIZE(TX,1)
    ALLOCATE(T(NX,NX))
    DO k = 1,nlag
        CALL TABFD2(X(1:(N-k)),X((k+1):N),N-k,TX,TX,NX,NX,T)
        CALL SRHOBIVA(TX,TX,N,NX,NX,T,N-k,dum,nor)
        S(k) = dum
    ENDDO
END SUBROUTINE SSUNI2
!! *************************************************************************************************

SUBROUTINE SRhointegrand(X,x1,x2,N,h1,h2,h1biv,h2biv,SINT)
!DEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:'srhointegrand' :: srhointegrand
    USE SHARED_DATA
    INTEGER,INTENT(IN) :: N
    REAL(8),INTENT(IN) :: X(2),x1(N),x2(N),h1,h2,h1biv,h2biv
    REAL(8),INTENT(OUT):: SINT
    REAL(8) :: x1eval(N),x2eval(N),fx1,fx2,fx12,DN1(N),DN2(N)
!   !# input X is evaluation point for x1 and x2, a 2x1 vector
    x1eval(:) = X(1)
    x2eval(:) = X(2)
!# x1 and x2 are the data vectors
!# Compute the marginal densities
    CALL DNORMF((x1eval-x1)/h1,N,DN1)
    fx1 = sum(DN1)/(N*h1)
    CALL DNORMF((x2eval-x2)/h2,N,DN2)
    fx2 = sum(DN2)/(N*h2)
!# Compute the bivariate density
    CALL DNORMF((x1eval-x1)/h1biv,N,DN1)
    CALL DNORMF((x2eval-x2)/h2biv,N,DN2)
    fx12 = (sum(DN1*DN2)/ (N*h1biv*h2biv))
!# Return the integrand
    SINT = (sqrt(fx12)-sqrt(fx1)*sqrt(fx2))**2
END SUBROUTINE SRhointegrand
!! *************************************************************************************************

SUBROUTINE kdenestmlcv(X,N,h,F,DMACH)
!DEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:'kdenestmlcv' :: kdenestmlcv

    USE SHARED_DATA
    IMPLICIT NONE
    INTEGER:: N,i
    REAL(8) :: X(N),h,F
    REAL(8) :: fhat(N),LF(N),Xeval(N),d1(1),d2(1),DN(N)
    REAL(8) :: DMACH(4)
    d1 = 0
    CALL DNORMF(d1,1,d2)
    DO i=1,N
        Xeval(:) = X(i)
        CALL DNORMF((Xeval-X)/h,N,DN)
        fhat(i)=SUM(DN)-d2(1)
    ENDDO
    fhat = fhat/((N-1)*h)
    IF (h>0) THEN
        WHERE(fhat > 0.0)
            LF = LOG(fhat)
        ELSEWHERE
            LF = LOG(DMACH(3))
        END WHERE
        F = -SUM(LF)/N
    ELSE
        F = DMACH(4)
    ENDIF
END SUBROUTINE kdenestmlcv
!! *************************************************************************************************

SUBROUTINE kdenestmlcvb(X,Y,N,h,F,DMACH)
!DEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:'kdenestmlcvb' :: kdenestmlcvb
    USE SHARED_DATA
    IMPLICIT NONE

    INTEGER:: N,i
    REAL(8) :: X(N),Y(N),h(2),F
    REAL(8) :: fhat(N),LF(N),Xeval(N),Yeval(N),d1(1),d2(1),DN1(N),DN2(N)
    REAL(8) :: DMACH(4)
    d1 = 0
    CALL DNORMF(d1,1,d2)
    DO i=1,N
        Xeval(:) = X(i)
        Yeval(:) = Y(i)
        CALL DNORMF((Xeval-X)/h(1),N,DN1)
        CALL DNORMF((Yeval-Y)/h(2),N,DN2)
        fhat(i)= SUM(DN1*DN2) - d2(1)**2
    ENDDO
    fhat = fhat/((N-1)*h(1)*h(2))
    IF ((h(1)>0).AND.(h(2)>0))  THEN
        WHERE(fhat > 0.0)
            LF = LOG(fhat)
        ELSEWHERE
            LF = LOG(DMACH(3))
        END WHERE
        F = -SUM(LF)/N
    ELSE
        F = DMACH(4)
    ENDIF
END SUBROUTINE kdenestmlcvb
!! *************************************************************************************************

SUBROUTINE SRhosum(X,x1,x2,N,h1,h2,h1biv,h2biv,S)
!DEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:'srhosum' :: srhosum
    USE SHARED_DATA
    IMPLICIT NONE

    INTEGER:: N,i
    REAL(8) :: X(2),x1(N),x2(N),h1,h2,h1biv,h2biv,S
!   REAL(8) :: x1eval(N),x2eval(N),fx1,fx2,fx12
    S = 0
    DO i=1,N
        X(1) = x1(i)
        X(2) = x2(i)
        S = S + (1 - Srho_eval(X,x1,x2,N,h1,h2,h1biv,h2biv))**2
    ENDDO
    S = 0.5*S/N
CONTAINS
    FUNCTION Srho_eval(X,x1,x2,N,h1,h2,h1biv,h2biv)
        INTEGER:: N
        REAL(8):: X(2),x1(N),x2(N),h1,h2,h1biv,h2biv,Srho_eval
        REAL(8):: x1eval(N),x2eval(N),fx1,fx2,fx12,DN1(N),DN2(N)
        x1eval(:) = X(1)
        x2eval(:) = X(2)
        CALL DNORMF((x1eval-x1)/h1,N,DN1)
        CALL DNORMF((x2eval-x2)/h2,N,DN2)
        fx1 = SUM(DN1)/(N*h1)
        fx2 = SUM(DN2)/(N*h2)
        CALL DNORMF((x1eval-x1)/h1biv,N,DN1)
        CALL DNORMF((x2eval-x2)/h2biv,N,DN2)
        fx12 = (SUM(DN1*DN2))/ (N*h1biv*h2biv)
        Srho_eval = SQRT(fx1*fx2/fx12)
    END FUNCTION Srho_eval
END SUBROUTINE SRhosum
!! *************************************************************************************************

SUBROUTINE SURROGATEACF(A,N,nlag,Te,RT,eps,nsuccmax,nmax,nsurr,che,S)
!DEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:'surrogateacf' :: surrogateacf

!   Algoritmo di generazione dei surrogati mediante
!   annealing simulato (Giugno 1999, rivisto Gennaio 2000, Maggio 2005, AGOSTO 2005)
!   AGOSTO 2005 VERSIONE PER R
!   DICEMBRE 2008 VERSIONE PER tseriesEntropy
!
!   implicit none
    USE SHARED_DATA

    REAL(8)  :: RT,Te,T1,eps,DeltaC,p,cost1,cost2,x(3)
    INTEGER :: nsucc,nsuccmax,nmax,N,nlag,h,k,che,nsurr,n1,n2,nrep
    INTEGER :: ind(N)

    REAL(8) :: Me,Var
    REAL(8) :: A(N),As(3*N+1),Corig(nlag),Csurr(nlag),Csurr2(nlag),S(N,nsurr)

!     Input  Parameters --------------------------------------------------------
!
!      N                   ! length of the series
!      nlag                ! minimization w.r.t to the firts nlag lags
!      Te                  ! initial temperature
!      RT                  ! reduction facttor for Te
!      eps                 ! target tolerance
!      nsuccmax            ! max number of successes after which Te is lowered
!      nmax                ! max number of iterations after which Te is lowered
!      nsurr               ! number of surrogates
!      che                 ! check
!
    nsuccmax = nsuccmax*N;
    nmax = nmax*N;
    che = che*2*N;

    S = 0
    T1=Te
    nsucc=0
    cost1=10000
    As = 0
    call Meva(A,N,Me,Var)
    A = A-Me
    As(N+2:2*N+1)=A(1:N)

    call ACFsurr(A,N,nlag,Corig,Me,Var) ! funzione di costo originale
!  Correction for robustness 14/09/2014
        Corig = Corig*1.05
!  Correction for robustness 14/09/2014
    do  nrep=1,nsurr
        Te = T1                     ! reinizializza la temperatura per ogni surrogato
        ind = PERM(N)               ! permutazione casuale iniziale della serie
        As(N+2:2*N+1) = A(ind);     ! permutazione casuale iniziale della serie
        call ACFsurr(As(N+2:2*N+1),N,nlag,Csurr,Me,Var)   ! funzione di costo del surrogato

        cost1=Costo(Corig,Csurr,nlag)
        CALL RANDOM_seed()
        k=0                                     ! k     : numero delle iterazioni globali
        do while(cost1>=eps)                    !
            h=0                                 ! h     : numero delle iterazioni per ogni serie
            nsucc=0                             ! nsucc : numero di successi per ogni serie
            do while(nsucc<nsuccmax)
                CALL RANDOM_NUMBER(x)
                n1=int(N*x(1))+1
                n2=int(N*x(2))+1
                do while(n1==n2)            !
                    CALL RANDOM_NUMBER(x(2))!
                    n2=int(N*x(2))+1        !       Check per numeri uguali
                enddo                       !

                call ACFswap(Csurr,As,N,n1,n2,nlag,Csurr2,Var)
                cost2 = Costo(Corig,Csurr2,nlag)

                DeltaC = cost2-cost1

                if (DeltaC>=0) then
                    p=exp((-DeltaC)/Te)
                    if(x(3)<=p) then
                        nsucc=nsucc+1
                        cost1=cost2
                        Csurr=Csurr2
                    else

                    call swap(As(n2+N+1),As(n1+N+1))

                    endif
                else
                    nsucc=nsucc+1
                    cost1=cost2
                    Csurr=Csurr2
                endif
                h=h+1; k=k+1
                if (h>nmax) then
                    nsucc=nsuccmax+1
                endif
                if (k>=che) then
                    nsucc=nsuccmax+1
                    ! STARTS AGAIN
                    Te = T1
                    k = 0
                    ind = PERM(N)
                    As(N+2:2*N+1) = A(ind);
                    call ACFsurr(As(N+2:2*N+1),N,nlag,Csurr,Me,Var)   ! funzione di costo del surrogato
                    cost1=Costo(Corig,Csurr,nlag)
                    CALL RANDOM_seed()
                endif
            enddo       !while(nsucc<nsuccmax)
            Te=Te*RT
            if (h<=nmax) then
            endif
        enddo       ! while(cost1>=eps)

        S(:,nrep) = As(N+2:2*N+1)
        call ACFsurr(As(N+2:2*N+1),N,nlag,Csurr2,Me,Var)
        cost2=Costo(Corig,Csurr2,nlag)
    end do          ! nrep=1,nsurr

CONTAINS

!! *************************************************************************************************

SUBROUTINE SWAP(x,y)
    REAL(8),intent(INOUT):: x,y
    REAL(8) :: tmp
    tmp = x; x = y; y = tmp
end SUBROUTINE SWAP
!! *************************************************************************************************

REAL(8) FUNCTION COSTO(C1,C2,N)
    implicit none
    integer,intent(in):: N
    REAL(8),dimension(N),intent(in):: C1,C2
    COSTO = maxval(abs((C1-C2)))
    return
end FUNCTION Costo
!! *************************************************************************************************

SUBROUTINE MEVA(X,N,Me,Var)
!   Calcola media e varianza (corretta) di una serie X
    INTEGER,intent(in) :: N
    REAL(8),intent(in) :: X(N)
    REAL(8),intent(out):: Me,Var
    REAL(8) :: S2(N),S

    S   = SUM(X)
    Me  = S/N              ! Media
    S2  = (X - Me)**2
    Var = SUM(S2)/(N-1)    ! Varianza campionaria corretta
end SUBROUTINE MEVA
!! *************************************************************************************************

SUBROUTINE ACFsurr(IN,N,nlag,Rho,Me,Var)

!   Gennaio 2000    Calcola la funzione di Autocorrelazione (Rho)
!   per i primi nlag tempi su una serie (IN) lunga N, di media Me e varianza Var.

    REAL(8),dimension(N),intent(in):: IN
    REAL(8),dimension(nlag),intent(out):: Rho
    REAL(8) :: INd(N) !,S
    REAL(8) :: Me,Var !,S2,IN2(N)

    ! write(*,*) Me,Var**.5,np,' Media',' sigma',' NN'
    INd=IN  !-Me        ! togliere il commento per sottrarre la media
    DO  K=1,nlag
    Rho(K)=DOT_PRODUCT(INd(K+1:N),INd(1:N-K))
    enddo
    Rho=Rho/(N*Var)

end SUBROUTINE ACFsurr

!! *************************************************************************************************

SUBROUTINE ACFswap(Rho1,As,N,n1,n2,nlag,Rho2,Var)

!   25/06/01    Parte dalla ACF (Rho1 calcolata su nlag tempi) di una serie in input (As,lunga N)
!   e calcola la ACF (Rho2) della serie dopo avere scambiato i due elementi di posti n1 e n2.
!   N.B. As ha N zeri prima e dopo
!
    integer,intent(in) :: n1,n2,nlag,N
    REAL(8),intent(in) :: Var
    REAL(8),dimension(nlag),intent(in) :: Rho1
    REAL(8),dimension(3*N+1),intent(inout) :: As
    REAL(8),dimension(nlag),intent(out):: Rho2
    REAL(8),dimension(nlag,4):: b

    b=0
    b(:,1)=As(n1+N+1)*As(n1+N+2:n1+nlag+N+1:1)
    b(:,2)=As(n2+N+1)*As(n2+N+2:n2+nlag+N+1:1)
    b(:,3)=As(n1+N+1)*As(n1+N:n1-nlag+N+1:-1)
    b(:,4)=As(n2+N+1)*As(n2+N:n2-nlag+N+1:-1)

    Rho2=Rho1*(N*Var)-b(:,1)-b(:,2)-b(:,3)-b(:,4)

    call swap(As(n1+N+1),As(n2+N+1))

    b=0
    b(:,1)=As(n1+N+1)*As(n1+N+2:n1+nlag+N+1:1)
    b(:,2)=As(n2+N+1)*As(n2+N+2:n2+nlag+N+1:1)
    b(:,3)=As(n1+N+1)*As(n1+N:n1-nlag+N+1:-1)
    b(:,4)=As(n2+N+1)*As(n2+N:n2-nlag+N+1:-1)

    Rho2=Rho2+b(:,1)+b(:,2)+b(:,3)+b(:,4)
    Rho2=Rho2/(N*Var)
end SUBROUTINE ACFswap
!! *************************************************************************************************
END SUBROUTINE  SURROGATEACF
!! *************************************************************************************************
