MODULE Problem
!
! This MODULE CONTAINS routines for initializing the problem and for driving
! the Multigrid routines in the MODULE AFASRoutines.  In particular, the
! routines in AFASRoutines USE RelaxGrid#D, Operator#D, Source#D, and
! SourceCorrection#D where # = 2 and 3.
!
! Templates for a biological growth model.
!
CONTAINS
!
SUBROUTINE SetProb
USE NodeInfoDef
USE ProblemDef
IMPLICIT NONE
!
INTEGER:: ierror
NAMELIST/problemdata/czero, gamma, lambda, eps
!
OPEN(UNIT=75,FILE='problemdata.dat',STATUS='OLD',ACTION='READ',IOSTAT=ierror)
IF(ierror/=0) THEN
  PRINT *,'Error opening input file problemdata.dat. Program stop.'
  STOP
END IF
READ(75,NML=problemdata)
CLOSE(75)
OPEN(UNIT=76,FILE='output.dat',STATUS='OLD',ACTION='WRITE',FORM='FORMATTED', &
     POSITION='APPEND')
WRITE(76,NML=problemdata)
CLOSE(76)
!
eps2 = eps*eps
!
END SUBROUTINE SetProb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE AfterRun
USE NodeInfoDef
USE ProblemDef
IMPLICIT NONE
!
END SUBROUTINE AfterRun
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE QInit2D(q,mx,nrvars,h,xlower)
USE NodeInfoDef
USE GridUtilities, ONLY: ULap2D
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(1:,1:,1:), INTENT(OUT):: q
INTEGER, DIMENSION(1:2), INTENT(IN):: mx
INTEGER, INTENT(IN):: nrvars
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), DIMENSION(1:2), INTENT(IN):: xlower
!
INTEGER:: i, j
REAL(KIND=r8):: r1, r2, r3, r4, tmp, x, y, z
!
q(:,:,1) = -1.0_r8
q(:,:,2) =  1.0_r8
!
DO i = 1, mx(1)
  x = (REAL(i,KIND=r8)-0.5_r8)*h+xlower(1)
  DO j = 1, mx(2)
    y = (REAL(j,KIND=r8)-0.5_r8)*h+xlower(2)
!
    IF(2.6 <= x .AND. x <= 3.8 .AND. 3.0 <=y .AND. y <= 3.4) q(i,j,1) = 1.0_r8
!
  END DO
END DO
!
END SUBROUTINE QInit2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE QInit3D(q,mx,nrvars,h,xlower)
USE NodeInfoDef
USE GridUtilities, ONLY: ULap3D
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(1:,1:,1:,1:), INTENT(OUT):: q
INTEGER, DIMENSION(1:3), INTENT(IN):: mx
INTEGER, INTENT(IN):: nrvars
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), DIMENSION(1:3), INTENT(IN):: xlower
!
INTEGER:: i, j, k
REAL(KIND=r8):: x, y, z
!
DO i = 1, mx(1)
  x = (REAL(i,KIND=r8)-0.5_r8)*h+xlower(1)
  DO j = 1, mx(2)
    y = (REAL(j,KIND=r8)-0.5_r8)*h+xlower(2)
    DO k = 1, mx(3)
      z = (REAL(k,KIND=r8)-0.5_r8)*h+xlower(3)
!
      q(i,j,k,:) = 0.0_r8
!
    END DO
  END DO
END DO
!
END SUBROUTINE QInit3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE AfterStep2D(q,qc,mx,nrvars,h,xlower,level)
USE NodeInfoDef
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: q
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: qc
INTEGER, DIMENSION(1:2), INTENT(IN):: mx
INTEGER, INTENT(IN):: nrvars
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), DIMENSION(1:2), INTENT(IN):: xlower
INTEGER, INTENT(IN):: level
!
INTEGER:: tmp
INTEGER, DIMENSION(1:2):: cmx
REAL(KIND=r8):: ch, ch2, h2
!
cmx(1:2) = mx(1:2)/2
h2 = h*h; ch = h*2.0_r8; ch2 = ch*ch
!
IF(level==0) THEN
  integralresult(1) &
    = integralresult(1)+h2*SUM(q(1:mx(1),1:mx(2),1))
ELSE
  integralresult(1) &
    = integralresult(1)+h2*SUM(q( 1: mx(1),1: mx(2),1)) &
    -                  ch2*SUM(qc(1:cmx(1),1:cmx(2),1))
END IF
!
END SUBROUTINE AfterStep2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE AfterStep3D(q,qc,mx,nrvars,h,xlower,level)
USE NodeInfoDef
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:,0:,1:), INTENT(IN):: q
REAL(KIND=r8), DIMENSION(0:,0:,0:,1:), INTENT(IN):: qc
INTEGER, DIMENSION(1:3), INTENT(IN):: mx
INTEGER, INTENT(IN):: nrvars
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), DIMENSION(1:3), INTENT(IN):: xlower
INTEGER, INTENT(IN):: level
!
END SUBROUTINE AfterStep3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    
SUBROUTINE SetAux(info)
USE NodeInfoDef
USE ProblemDef
IMPLICIT NONE
!
TYPE(nodeinfo):: info
!
INTEGER:: i, j, k
REAL(KIND=r8):: time, x, y, z
!
!time = info%time
!
END SUBROUTINE SetAux
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    
SUBROUTINE SetSrc(info)
USE NodeInfoDef
USE ProblemDef
IMPLICIT NONE
!    
TYPE(nodeinfo):: info
!
INTEGER:: i, j, k
REAL(KIND=r8):: time, x, y, z
!
!time = info%time
!
END SUBROUTINE SetSrc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE SetErrFlagsUser2D(qrte,qnew,errorflags,mx,cmx,nrvars,h,xlower,level)
USE NodeInfoDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(1:,1:,1:), INTENT(IN):: qrte
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: qnew
INTEGER, DIMENSION(1:,1:), INTENT(OUT):: errorflags
INTEGER, DIMENSION(1:2), INTENT(IN):: mx
INTEGER, DIMENSION(1:2), INTENT(IN):: cmx
INTEGER, INTENT(IN):: nrvars
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), DIMENSION(1:2), INTENT(IN):: xlower
INTEGER, INTENT(IN):: level
!
errorflags(1:mx(1),1:mx(2)) = 0
!
! Error tagging based on the second differences of the phase variable:
WHERE(  ABS(       qnew(2:mx(1)+1,1:mx(2)  ,1)  &
      -     2.0_r8*qnew(1:mx(1)  ,1:mx(2)  ,1)  &
      +            qnew(0:mx(1)-1,1:mx(2)  ,1)) &
      + ABS(       qnew(1:mx(1)  ,2:mx(2)+1,1)  &
      -     2.0_r8*qnew(1:mx(1)  ,1:mx(2)  ,1)  &
      +            qnew(1:mx(1)  ,0:mx(2)-1,1)) > 2.0_r8*qtolerance(level))
  errorflags(1:mx(1),1:mx(2)) = 1
END WHERE
!
END SUBROUTINE SetErrFlagsUser2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE SetErrFlagsUser3D(qrte,qnew,errorflags,mx,cmx,nrvars,h,xlower,level)
USE NodeInfoDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(1:,1:,1:,1:), INTENT(IN):: qrte
REAL(KIND=r8), DIMENSION(0:,0:,0:,1:), INTENT(IN):: qnew
INTEGER, DIMENSION(1:,1:,1:), INTENT(OUT):: errorflags
INTEGER, DIMENSION(1:3), INTENT(IN):: mx
INTEGER, DIMENSION(1:3), INTENT(IN):: cmx
INTEGER, INTENT(IN):: nrvars
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), DIMENSION(1:3), INTENT(IN):: xlower
INTEGER, INTENT(IN):: level
!
errorflags(1:mx(1),1:mx(2),1:mx(3)) = 0
!
! Error tagging based on the second differences of the phase variable:
WHERE(  ABS(       qnew(2:mx(1)+1,1:mx(2)  ,1:mx(3)  ,1)  &
      -     2.0_r8*qnew(1:mx(1)  ,1:mx(2)  ,1:mx(3)  ,1)  &
      +            qnew(0:mx(1)-1,1:mx(2)  ,1:mx(3)  ,1)) &
      + ABS(       qnew(1:mx(1)  ,2:mx(2)+1,1:mx(3)  ,1)  &
      -     2.0_r8*qnew(1:mx(1)  ,1:mx(2)  ,1:mx(3)  ,1)  &
      +            qnew(1:mx(1)  ,0:mx(2)-1,1:mx(3)  ,1)) &
      + ABS(       qnew(1:mx(1)  ,1:mx(2)  ,2:mx(3)+1,1)  &
      -     2.0_r8*qnew(1:mx(1)  ,1:mx(2)  ,1:mx(3)  ,1)  &
      +            qnew(1:mx(1)  ,1:mx(2)  ,0:mx(3)-1,1)) &
      > 3.0_r8*qtolerance(level))
  errorflags(1:mx(1),1:mx(2),1:mx(3)) = 1
END WHERE
!
END SUBROUTINE SetErrFlagsUser3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! User-Supplied AFAS Multigrid Routines:
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! 2D Multigrid Routines:
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE RelaxGrid2D(qnew,qold,aux,f,mx,nrvars,maux,h,xlower,ipass)
USE NodeInfoDef
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN OUT):: qnew
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: qold
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: aux
REAL(KIND=r8), DIMENSION(1:,1:,1:), INTENT(IN):: f
INTEGER, DIMENSION(1:2), INTENT(IN):: mx
INTEGER, INTENT(IN):: nrvars
INTEGER, INTENT(IN):: maux
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), DIMENSION(1:2), INTENT(IN):: xlower
INTEGER, INTENT(IN):: ipass
!
INTEGER:: i, j
REAL(KIND=r8):: h2, tmp1, tmp2, tmp3, tmp4
!
h2 = h*h
!
tmp1 = dt*eps2/h2
tmp3 = 1.0_r8/h2
!
DO j = 1, mx(2)
  DO i = 1+MODULO(j+ipass,2), mx(1), 2
!
    tmp2 = 4.0_r8*tmp1+1.0_r8+dt*Sqr(qnew(i,j,1))
!
    qnew(i,j,1) = (f(i,j,1)+dt*Growth(qold(i,j,1))*(qnew(i,j,2)-czero) &
                + tmp1*(qnew(i-1,j,1)+qnew(i+1,j,1) &
         +              qnew(i,j-1,1)+qnew(i,j+1,1)))/tmp2
!
    tmp4 = 4.0_r8*tmp3+Uptake(qold(i,j,1))
!
    qnew(i,j,2) = (f(i,j,2)+tmp3*(qnew(i+1,j,2)+qnew(i-1,j,2) &
         +                        qnew(i,j+1,2)+qnew(i,j-1,2)))/tmp4
!
  END DO
END DO
!
END SUBROUTINE RelaxGrid2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE UpdateAuxVcycle2D(qnew,qold,aux,ll,mx,nrvars,maux,mbc,h,xlower)
USE NodeInfoDef
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(ll:,ll:,1:), INTENT(IN):: qnew
REAL(KIND=r8), DIMENSION(ll:,ll:,1:), INTENT(IN):: qold
REAL(KIND=r8), DIMENSION(ll:,ll:,1:), INTENT(IN OUT):: aux
INTEGER, INTENT(IN):: ll
INTEGER, DIMENSION(1:2), INTENT(IN):: mx
INTEGER, INTENT(IN):: nrvars
INTEGER, INTENT(IN):: maux
INTEGER, INTENT(IN):: mbc
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), DIMENSION(1:2), INTENT(IN):: xlower
!
!aux(0:mx(1)+1,0:mx(2)+1,1) = 
!
END SUBROUTINE UpdateAuxVcycle2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION Operator2D(qnew,qold,aux,mx,nrvars,maux,h,xlower) RESULT(opresult)
USE NodeInfoDef
USE GridUtilities, ONLY: ULap2D
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: qnew
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: qold
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: aux
INTEGER, DIMENSION(1:2), INTENT(IN):: mx
INTEGER, INTENT(IN):: nrvars
INTEGER, INTENT(IN):: maux
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), DIMENSION(1:2), INTENT(IN):: xlower
REAL(KIND=r8), DIMENSION(1:mx(1),1:mx(2),1:nrvars):: opresult
!
REAL(KIND=r8):: h2, tmp1, tmp2, tmp3
!
h2 = h*h
tmp1 = dt*eps2/h2
tmp2 = dt*lambda*0.5_r8
tmp3 = 1.0_r8/h2
!
opresult(1:mx(1),1:mx(2),1) = qnew(1:mx(1),1:mx(2),1) &
                            - tmp1*ULap2D(qnew(0:mx(1)+1,0:mx(2)+1,1)) &
                            + dt*Cube(  qnew(1:mx(1),1:mx(2),1)) &
                            - dt*Growth(qold(1:mx(1),1:mx(2),1)) &
                            * (qnew(1:mx(1),1:mx(2),2)-czero)
!
opresult(1:mx(1),1:mx(2),2) = -tmp3*ULap2D(qnew(0:mx(1)+1,0:mx(2)+1,2)) &
                            + Uptake(qold(1:mx(1),1:mx(2),1)) &
                            * qnew(1:mx(1),1:mx(2),2)
!
END FUNCTION Operator2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION Source2D(qnew,qold,aux,ll,mx,nrvars,maux,mbc,h,xlower) &
         RESULT(sourceresult)
USE NodeInfoDef
USE GridUtilities, ONLY: ULap2D
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(ll:,ll:,1:), INTENT(IN):: qnew
REAL(KIND=r8), DIMENSION(ll:,ll:,1:), INTENT(IN):: qold
REAL(KIND=r8), DIMENSION(ll:,ll:,1:), INTENT(IN):: aux
INTEGER, INTENT(IN):: ll
INTEGER, DIMENSION(1:2), INTENT(IN):: mx
INTEGER, INTENT(IN):: nrvars
INTEGER, INTENT(IN):: maux
INTEGER, INTENT(IN):: mbc
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), DIMENSION(1:2), INTENT(IN):: xlower
REAL(KIND=r8), DIMENSION(1:mx(1),1:mx(2),1:nrvars):: sourceresult
!
sourceresult(1:mx(1),1:mx(2),1) = (1.0_r8+dt)*qold(1:mx(1),1:mx(2),1)
!
sourceresult(1:mx(1),1:mx(2),2) = 0.0_r8
!
END FUNCTION Source2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION SourceUpdate2D(qnew,qold,aux,ll,mx,nrvars,maux,mbc,h,xlower) &
         RESULT(updateresult)
USE NodeInfoDef
USE GridUtilities, ONLY: UDiv2D, ULap2D
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(ll:,ll:,1:), INTENT(IN):: qnew
REAL(KIND=r8), DIMENSION(ll:,ll:,1:), INTENT(IN):: qold
REAL(KIND=r8), DIMENSION(ll:,ll:,1:), INTENT(IN):: aux
INTEGER, INTENT(IN):: ll
INTEGER, DIMENSION(1:2), INTENT(IN):: mx
INTEGER, INTENT(IN):: nrvars
INTEGER, INTENT(IN):: maux
INTEGER, INTENT(IN):: mbc
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), DIMENSION(1:2), INTENT(IN):: xlower
REAL(KIND=r8), DIMENSION(1:mx(1),1:mx(2),1:nrvars):: updateresult
!
updateresult(1:mx(1),1:mx(2),1) = 0.0_r8
!
updateresult(1:mx(1),1:mx(2),2) = 0.0_r8
!
END FUNCTION SourceUpdate2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION ULapMob2D(c,a) RESULT(ulapmobresult)
USE NodeInfoDef
USE GridUtilities, ONLY: UDiv2D
IMPLICIT NONE
!
! Level independent, 2D UNDIVIDED laplacian operator for non-constant mobity.
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: c
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: a
REAL(KIND=r8), DIMENSION(1:SIZE(a,1)-2,1:SIZE(a,2)-2):: ulapmobresult
!
INTEGER, DIMENSION(1:2):: mx
REAL(KIND=r8), DIMENSION(0:SIZE(a,1)-2,1:SIZE(a,2)-2):: f1
REAL(KIND=r8), DIMENSION(1:SIZE(a,1)-2,0:SIZE(a,2)-2):: f2
!
mx(1) = SIZE(a,1)-2; mx(2) = SIZE(a,2)-2
!
! Calculate the UNDIVIDED 2D flux function:
f1(0:mx(1),1:mx(2)) &
  = Mob(0.5_r8*(c(1:mx(1)+1,1:mx(2))+c(0:mx(1),1:mx(2)))) &
  *            (a(1:mx(1)+1,1:mx(2))-a(0:mx(1),1:mx(2)))
!
f2(1:mx(1),0:mx(2)) &
  = Mob(0.5_r8*(c(1:mx(1),1:mx(2)+1)+c(1:mx(1),0:mx(2)))) &
  *            (a(1:mx(1),1:mx(2)+1)-a(1:mx(1),0:mx(2)))
!
! Calculate the UNDIVIDED divergence of the flux:
ulapmobresult(1:mx(1),1:mx(2)) = UDiv2D(f1(0:mx(1),1:mx(2)), &
                                        f2(1:mx(1),0:mx(2)))
!
END FUNCTION ULapMob2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! 3D Multigrid Routines:
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE RelaxGrid3D(qnew,qold,aux,f,mx,nrvars,maux,h,xlower,ipass)
USE NodeInfoDef
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:,0:,1:), INTENT(IN OUT):: qnew
REAL(KIND=r8), DIMENSION(0:,0:,0:,1:), INTENT(IN):: qold
REAL(KIND=r8), DIMENSION(0:,0:,0:,1:), INTENT(IN):: aux
REAL(KIND=r8), DIMENSION(1:,1:,1:,1:), INTENT(IN):: f
INTEGER, DIMENSION(1:3), INTENT(IN):: mx
INTEGER, INTENT(IN):: nrvars
INTEGER, INTENT(IN):: maux
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), DIMENSION(1:3), INTENT(IN):: xlower
INTEGER, INTENT(IN):: ipass
!
INTEGER:: i, j, k
!
! These loops must always be structured this way, i.e., so that i,j,k = 1 is a
! black cell. (The first cell relaxed is always i = 2, j,k = 1, since ipass = 1
! is called first.)
!
DO k = 1, mx(3)
  DO j = 1, mx(2)
    DO i = 1+MODULO(k+j+ipass,2), mx(1), 2
!
      qnew(i,j,k,1) = 0.0_r8
      qnew(i,j,k,2) = 0.0_r8
!
    END DO
  END DO
END DO
!
END SUBROUTINE RelaxGrid3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE UpdateAuxVcycle3D(qnew,qold,aux,ll,mx,nrvars,maux,mbc,h,xlower)
USE NodeInfoDef
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN):: qnew
REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN):: qold
REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN OUT):: aux
INTEGER, INTENT(IN):: ll
INTEGER, DIMENSION(1:3), INTENT(IN):: mx
INTEGER, INTENT(IN):: nrvars
INTEGER, INTENT(IN):: maux
INTEGER, INTENT(IN):: mbc
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), DIMENSION(1:3), INTENT(IN):: xlower
!
!aux(0:mx(1)+1,0:mx(2)+1,0:mx(3)+1,1) = 
!
END SUBROUTINE UpdateAuxVcycle3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION Operator3D(qnew,qold,aux,mx,nrvars,maux,h,xlower) RESULT(opresult)
USE NodeInfoDef
USE GridUtilities, ONLY: ULap3D
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:,0:,1:), INTENT(IN):: qnew
REAL(KIND=r8), DIMENSION(0:,0:,0:,1:), INTENT(IN):: qold
REAL(KIND=r8), DIMENSION(0:,0:,0:,1:), INTENT(IN):: aux
INTEGER, DIMENSION(1:3), INTENT(IN):: mx
INTEGER, INTENT(IN):: nrvars
INTEGER, INTENT(IN):: maux
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), DIMENSION(1:3), INTENT(IN):: xlower
REAL(KIND=r8), DIMENSION(1:mx(1),1:mx(2),1:mx(3),1:nrvars):: opresult
!
opresult(1:mx(1),1:mx(2),1:mx(3),1) = 0.0_r8
!
opresult(1:mx(1),1:mx(2),1:mx(3),2) = 0.0_r8
!
END FUNCTION Operator3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION Source3D(qnew,qold,aux,ll,mx,nrvars,maux,mbc,h,xlower) &
         RESULT(sourceresult)
USE NodeInfoDef
USE GridUtilities, ONLY: ULap3D
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN):: qnew
REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN):: qold
REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN):: aux
INTEGER, INTENT(IN):: ll
INTEGER, DIMENSION(1:3), INTENT(IN):: mx
INTEGER, INTENT(IN):: nrvars
INTEGER, INTENT(IN):: maux
INTEGER, INTENT(IN):: mbc
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), DIMENSION(1:3), INTENT(IN):: xlower
REAL(KIND=r8), DIMENSION(1:mx(1),1:mx(2),1:mx(3),1:nrvars):: sourceresult
!
sourceresult(1:mx(1),1:mx(2),1:mx(3),1) = 0.0_r8
!
sourceresult(1:mx(1),1:mx(2),1:mx(3),2) = 0.0_r8
!
END FUNCTION Source3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION SourceUpdate3D(qnew,qold,aux,ll,mx,nrvars,maux,mbc,h,xlower) &
         RESULT(updateresult)
USE NodeInfoDef
USE GridUtilities, ONLY: UDiv3D, ULap3D
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN):: qnew
REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN):: qold
REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN):: aux
INTEGER, INTENT(IN):: ll
INTEGER, DIMENSION(1:3), INTENT(IN):: mx
INTEGER, INTENT(IN):: nrvars
INTEGER, INTENT(IN):: maux
INTEGER, INTENT(IN):: mbc
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), DIMENSION(1:3), INTENT(IN):: xlower
REAL(KIND=r8), DIMENSION(1:mx(1),1:mx(2),1:mx(3),1:nrvars):: updateresult
!
updateresult(1:mx(1),1:mx(2),1:mx(3),1) = 0.0_r8
!
updateresult(1:mx(1),1:mx(2),1:mx(3),2) = 0.0_r8
!
END FUNCTION SourceUpdate3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION ULapMob3D(c,a) RESULT(ulapmobresult)
USE NodeInfoDef
USE GridUtilities, ONLY: UDiv3D
IMPLICIT NONE
!
! Level independent 3D UNDIVIDED laplacian operator for non-constant mobity.
!
REAL(KIND=r8), DIMENSION(0:,0:,0:), INTENT(IN):: c
REAL(KIND=r8), DIMENSION(0:,0:,0:), INTENT(IN):: a
REAL(KIND=r8), DIMENSION(1:SIZE(a,1)-2, &
                         1:SIZE(a,2)-2, &
                         1:SIZE(a,3)-2):: ulapmobresult
!
INTEGER, DIMENSION(1:3):: mx
REAL(KIND=r8), DIMENSION(0:SIZE(a,1)-2,1:SIZE(a,2)-2,1:SIZE(a,3)-2):: f1
REAL(KIND=r8), DIMENSION(1:SIZE(a,1)-2,0:SIZE(a,2)-2,1:SIZE(a,3)-2):: f2
REAL(KIND=r8), DIMENSION(1:SIZE(a,1)-2,1:SIZE(a,2)-2,0:SIZE(a,3)-2):: f3
!
mx(1) = SIZE(a,1)-2; mx(2) = SIZE(a,2)-2; mx(3) = SIZE(a,3)-2
!
! Calculate the 3D flux function:
f1(0:mx(1),1:mx(2),1:mx(3)) &
  = Mob(0.5_r8*(c(1:mx(1)+1,1:mx(2),1:mx(3))+c(0:mx(1),1:mx(2),1:mx(3)))) &
  *              (a(1:mx(1)+1,1:mx(2),1:mx(3))-a(0:mx(1),1:mx(2),1:mx(3)))
!
f2(1:mx(1),0:mx(2),1:mx(3)) &
  = Mob(0.5_r8*(c(1:mx(1),1:mx(2)+1,1:mx(3))+c(1:mx(1),0:mx(2),1:mx(3)))) &
  *              (a(1:mx(1),1:mx(2)+1,1:mx(3))-a(1:mx(1),0:mx(2),1:mx(3)))
!
f3(1:mx(1),1:mx(2),0:mx(3)) &
  = Mob(0.5_r8*(c(1:mx(1),1:mx(2),1:mx(3)+1)+c(1:mx(1),1:mx(2),0:mx(3)))) &
  *              (a(1:mx(1),1:mx(2),1:mx(3)+1)-a(1:mx(1),1:mx(2),0:mx(3)))
!
! Calculate the divergence of the flux:
ulapmobresult(1:mx(1),1:mx(2),1:mx(3)) &
    = UDiv3D(f1(0:mx(1),1:mx(2),1:mx(3)), &
             f2(1:mx(1),0:mx(2),1:mx(3)), &
             f3(1:mx(1),1:mx(2),0:mx(3)))
!
END FUNCTION ULapMob3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! User supplied physical boundary conditions:
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE UserBC2D(q,ll,mx,nrvars,mbc,bcnow)
USE NodeinfoDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(ll:,ll:,1:), INTENT(IN OUT):: q
INTEGER, INTENT(IN):: ll
INTEGER, DIMENSION(1:2), INTENT(IN):: mx
INTEGER, INTENT(IN):: nrvars
INTEGER, INTENT(IN):: mbc
INTEGER, INTENT(IN):: bcnow
!
INTEGER:: ibc
INTEGER, DIMENSION(1:2):: ul
!
ul(1:2) = mx(1:2)+mbc
!
SELECT CASE(bcnow)
  CASE(1)
!
! bndry 1:
    DO ibc = 1, mbc
      q(1-ibc,ll:ul(2),1) =        q(ibc,ll:ul(2),1)
      q(1-ibc,ll:ul(2),2) = 2.0_r8-q(ibc,ll:ul(2),2)
    END DO
!
  CASE(2)
!
! bndry 2:
    DO ibc = 1, mbc
      q(mx(1)+ibc,ll:ul(2),1) =        q(mx(1)-ibc+1,ll:ul(2),1)
      q(mx(1)+ibc,ll:ul(2),2) = 2.0_r8-q(mx(1)-ibc+1,ll:ul(2),2)
    END DO
!
  CASE(3)
!
! bndry 3:
    DO ibc = 1, mbc
      q(ll:ul(1),1-ibc,1) =        q(ll:ul(1),ibc,1)
      q(ll:ul(1),1-ibc,2) = 2.0_r8-q(ll:ul(1),ibc,2)
    END DO
!
  CASE(4)
!
! bndry 4:
    DO ibc = 1, mbc
      q(ll:ul(1),mx(2)+ibc,1) =        q(ll:ul(1),mx(2)-ibc+1,1)
      q(ll:ul(1),mx(2)+ibc,2) = 2.0_r8-q(ll:ul(1),mx(2)-ibc+1,2)
    END DO          
END SELECT
!
END SUBROUTINE UserBC2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE UserBC3D(q,ll,mx,nrvars,mbc,bcnow)
USE NodeinfoDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN OUT):: q
INTEGER, INTENT(IN):: ll
INTEGER, DIMENSION(1:3), INTENT(IN):: mx
INTEGER, INTENT(IN):: nrvars
INTEGER, INTENT(IN):: mbc
INTEGER, INTENT(IN):: bcnow
!
INTEGER:: ibc
INTEGER, DIMENSION(1:3):: ul
!
ul(1:3) = mx(1:3)+mbc
!
SELECT CASE(bcnow)
  CASE(1)
!
! Dirichlet on bndry 1:
    DO ibc = 1, mbc
      q(1-ibc,ll:ul(2),ll:ul(3),:) = 2.0_r8-q(ibc,ll:ul(2),ll:ul(3),:)
    END DO
!
  CASE(2)
!
! Dirichlet on bndry 2:
    DO ibc = 1, mbc
      q(mx(1)+ibc,ll:ul(2),ll:ul(3),:) = 2.0_r8-q(mx(1)-ibc+1,ll:ul(2),ll:ul(3),:)
    END DO
!
  CASE(3)
!
! Dirichlet on bndry 3:
    DO ibc = 1, mbc
      q(ll:ul(1),1-ibc,ll:ul(3),:) = 2.0_r8-q(ll:ul(1),ibc,ll:ul(3),:)
    END DO
!
  CASE(4)
!
! Dirichlet on bndry 4:
    DO ibc = 1, mbc
      q(ll:ul(1),mx(2)+ibc,ll:ul(3),:) = 2.0_r8-q(ll:ul(1),mx(2)-ibc+1,ll:ul(3),:)
    END DO
!
  CASE(5)
!
! Dirichlet on bndry 5:
    DO ibc = 1, mbc
      q(ll:ul(1),ll:ul(2),1-ibc,:) = 2.0_r8-q(ll:ul(1),ll:ul(2),ibc,:)
    END DO
!
  CASE(6)
!
! Dirichlet on bndry 6:
    DO ibc = 1, mbc
      q(ll:ul(1),ll:ul(2),mx(3)+ibc,:) = 2.0_r8-q(ll:ul(1),ll:ul(2),mx(3)-ibc+1,:)
    END DO
END SELECT
!
END SUBROUTINE UserBC3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Dimensionally Invariant Multigrid Routines:
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
ELEMENTAL FUNCTION Growth(p) RESULT(grthresult)
USE NodeInfoDef
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), INTENT(IN):: p
REAL(KIND=r8):: grthresult
!
grthresult = gamma*0.5_r8*(p+1.0_r8)
!
END FUNCTION Growth
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
ELEMENTAL FUNCTION Uptake(p) RESULT(uptkresult)
USE NodeInfoDef
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), INTENT(IN):: p
REAL(KIND=r8):: uptkresult
!
uptkresult = lambda*0.5_r8*(p+1.0_r8)
!
END FUNCTION Uptake
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
ELEMENTAL FUNCTION Mob(c) RESULT(mobresult)
USE NodeInfoDef
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), INTENT(IN):: c
REAL(KIND=r8):: mobresult
!
mobresult = 1.0_r8
!
END FUNCTION Mob
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
ELEMENTAL FUNCTION Cube(c) RESULT(cuberesult)
USE NodeInfoDef
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), INTENT(IN):: c
REAL(KIND=r8):: cuberesult
!
cuberesult = c*c*c
!
END FUNCTION Cube
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
ELEMENTAL FUNCTION Sqr(c) RESULT(sqrresult)
USE NodeInfoDef
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), INTENT(IN):: c
REAL(KIND=r8):: sqrresult
!
sqrresult = c*c
!
END FUNCTION Sqr
END MODULE Problem
