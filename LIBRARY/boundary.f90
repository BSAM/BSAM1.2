!==========================================================================
! BSAM 1.2: Block-Structured Adaptive Multigrid Solver
!==========================================================================
!
! WiseSoft: Innovators, Brothers.
!
! (c) Copyright Steven M. Wise, 2006
! Department of Mathematics
! University of California at Irvine
! swise@math.uci.edu
!
! (c) Copyright Steven M. Wise, 2007
! Department of Mathematics
! University of Tennessee
! swise@math.utk.edu
!
! (c) Copyright Steven M. Wise, 2015
! Department of Mathematics
! University of Tennessee
! swise@math.utk.edu
!
! Portions of the code
!
! (c) Copyright Sorin Mitran, 2002
! Department of Mathematics
! University of North Carolina at Chapel Hill
! mitran@amath.unc.edu
!
! -----------------------------------------------------------------------
! This software is made available for research and instructional use only.
! You may copy and use this software without charge for these
! non-commercial purposes, provided that the copyright notice and
! associated text is reproduced on all copies. For all other uses,
! including distribution of modified versions, please contact the authors.
!
! Commercial use is strictly forbidden without permission.
!
! This software is made available "as is" without any assurance that it
! will work for your purposes. The software may in fact have defects,
! so use the software at your own risk.
!
! -----------------------------------------------------------------------
! File:             boundary.f90
! Purpose:          Ghost-cell interpolation module.
! Contains:
! Revision History: Ver. 1.0 Oct. 2006 Steven Wise
! Revision History: Ver. 1.1 May. 2007 Steven Wise
! Revision History: Ver. 1.2 Jul. 2015 Steven Wise
! -----------------------------------------------------------------------
MODULE Boundary
USE NodeInfoDef, ONLY: r8
IMPLICIT NONE
!
REAL(KIND=r8), PARAMETER:: a1 = -01.0_r8/05.0_r8
REAL(KIND=r8), PARAMETER:: a2 =  02.0_r8/03.0_r8
REAL(KIND=r8), PARAMETER:: a3 =  08.0_r8/15.0_r8
REAL(KIND=r8), PARAMETER:: b1 =  05.0_r8/32.0_r8
REAL(KIND=r8), PARAMETER:: b2 =  15.0_r8/16.0_r8
REAL(KIND=r8), PARAMETER:: b3 = -03.0_r8/32.0_r8
REAL(KIND=r8), PARAMETER:: c1 =  09.0_r8/16.0_r8
REAL(KIND=r8), PARAMETER:: c2 =  03.0_r8/16.0_r8
REAL(KIND=r8), PARAMETER:: c3 =  03.0_r8/16.0_r8
REAL(KIND=r8), PARAMETER:: c4 =  01.0_r8/16.0_r8
!
CONTAINS
!
SUBROUTINE SetGhost(level,ipass)
USE NodeInfoDef
USE TreeOps, ONLY: ApplyOnLevel, ApplyOnLevelPairs
IMPLICIT NONE
!
INTEGER, INTENT(IN):: level
INTEGER, INTENT(IN):: ipass
!
TYPE(funcparam):: dummy
!
dummy%iswitch = ipass
!
! 0. Interpolate between coarse and fine grids:
IF(level>rootlevel) CALL ApplyOnLevel(level,InterpolateCoarseFine,dummy)
!
! 1. Transfer boundary conditions among grids on same level:
IF(level>rootlevel) CALL ApplyOnLevelPairs(level,TransferBC,dummy)
!
IF(periodicboundaryconditions) THEN
!
! 2. Apply periodic boundary conditions on same grid:
  CALL ApplyOnLevel(level,SelfPeriodicBC,dummy)
!
! 3. Apply periodic boundary conditions on different grids:
  IF(level>rootlevel) CALL ApplyOnLevelPairs(level,TransferPeriodicBC,dummy)
END IF
!
! 4. Apply physical boundary conditions:
CALL ApplyOnLevel(level,SetPhysicalBC,dummy)
!
END SUBROUTINE SetGhost
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION InterpolateCoarseFine(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok, GetParentInfo
IMPLICIT NONE
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
TYPE(nodeinfo), POINTER:: parent
INTEGER:: ierror, ipass, nrvars
INTEGER, DIMENSION(1:maxdims):: mx, pmx
INTEGER, DIMENSION(1:2*maxdims):: mthbc
INTEGER, DIMENSION(1:maxdims,1:2):: mb
!
InterpolateCoarseFine = err_ok
!
IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
ipass = dummy%iswitch
ierror = GetParentInfo(parent)
nrvars = info%nrvars
mx = 1; mx(1:ndims) = info%mx(1:ndims)
pmx = 1; pmx(1:ndims) = parent%mx(1:ndims)
mb = 1; mb(1:ndims,1:2) = info%mbounds(1:ndims,1:2)
mthbc = 1; mthbc(1:2*ndims) = info%mthbc(1:2*ndims)
!
SELECT CASE(ndims)
  CASE(2)
    CALL Interpolate2D(  info%q(0: mx(1)+1,0: mx(2)+1,1,1:nrvars), &
                       parent%q(0:pmx(1)+1,0:pmx(2)+1,1,1:nrvars), &
                       mx(1:2),nrvars,mb(1:2,1:2),mthbc(1:4),ipass)
  CASE(3)
    CALL Interpolate3D(  info%q(0: mx(1)+1,0: mx(2)+1,0: mx(3)+1,1:nrvars), &
                       parent%q(0:pmx(1)+1,0:pmx(2)+1,0:pmx(3)+1,1:nrvars), &
                       mx(1:3),nrvars,mb(1:3,1:2),mthbc(1:6),ipass)
  CASE DEFAULT
    PRINT *, 'BSAM 1.2, InterpolateCoarseFine: only n=2,3 supported.'
    STOP
END SELECT
!
END FUNCTION InterpolateCoarseFine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE Interpolate2D(qf,qc,mx,nrvars,mb,mthbc,ipass)
USE NodeInfoDef, ONLY: r8
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN OUT):: qf
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: qc
INTEGER, DIMENSION(1:2), INTENT(IN):: mx
INTEGER, INTENT(IN):: nrvars
INTEGER, DIMENSION(1:2,1:2), INTENT(IN):: mb
INTEGER, DIMENSION(1:4), INTENT(IN):: mthbc
INTEGER, INTENT(IN):: ipass
!
! Corners:
!
! Red corners:
IF(ipass==1 .OR. ipass==0) THEN
  qf(mx(1)+1,mx(2)+1,1:nrvars) &
    = a1*qf(mx(1)-1,mx(2)-1,1:nrvars) &
    + a2*qf(mx(1)  ,mx(2)  ,1:nrvars) &
    + a3*qc(mb(1,2)+1,mb(2,2)+1,1:nrvars)
!
  qf(0,0,1:nrvars) &
    = a1*qf(2,2,1:nrvars) &
    + a2*qf(1,1,1:nrvars) &
    + a3*qc(mb(1,1)-1,mb(2,1)-1,1:nrvars)
END IF
!
! Black corners:
IF(ipass==2 .OR. ipass==0) THEN
  qf(mx(1)+1,0,1:nrvars) &
    = a1*qf(mx(1)-1,2,1:nrvars) &
    + a2*qf(mx(1)  ,1,1:nrvars) &
    + a3*qc(mb(1,2)+1,mb(2,1)-1,1:nrvars)
!
  qf(0,mx(2)+1,1:nrvars) &
    = a1*qf(2,mx(2)-1,1:nrvars) &
    + a2*qf(1,mx(2)  ,1:nrvars) &
    + a3*qc(mb(1,1)-1,mb(2,2)+1,1:nrvars)
END IF
!
! Faces:
!
SELECT CASE(mthbc(1))
CASE(2,999) ! Periodic or Internal boundary.
!
! Face #1, Red:
  IF(ipass==1 .OR. ipass==0) THEN
    qf(0,2:mx(2)-2:2,1:nrvars) &
      = a1*qf(2,4:mx(2)  :2,1:nrvars) &
      + a2*qf(1,3:mx(2)-1:2,1:nrvars) &
      + a3*qc(mb(1,1)-1,mb(2,1):mb(2,2)-1,1:nrvars)
!
    qf(0,mx(2),1:nrvars) &
      = a1*qf(2,mx(2),1:nrvars) &
      + a2*qf(1,mx(2),1:nrvars) &
      + a3*(b1*qc(mb(1,1)-1,mb(2,2)+1,1:nrvars) &
      +     b2*qc(mb(1,1)-1,mb(2,2)  ,1:nrvars) &
      +     b3*qc(mb(1,1)-1,mb(2,2)-1,1:nrvars))
  END IF
!
! Face #1, Black:
  IF(ipass==2 .OR. ipass==0) THEN
    qf(0,3:mx(2)-1:2,1:nrvars) &
      = a1*qf(2,1:mx(2)-3:2,1:nrvars) &
      + a2*qf(1,2:mx(2)-2:2,1:nrvars) &
      + a3*qc(mb(1,1)-1,mb(2,1)+1:mb(2,2),1:nrvars)
!
    qf(0,1,1:nrvars) &
      =          a1*qf(2,1,1:nrvars) &
      + a2*qf(1,1,1:nrvars) &
      + a3*(b1*qc(mb(1,1)-1,mb(2,1)-1,1:nrvars) &
      +     b2*qc(mb(1,1)-1,mb(2,1)  ,1:nrvars) &
      +     b3*qc(mb(1,1)-1,mb(2,1)+1,1:nrvars))
  END IF
!
CASE DEFAULT
  CONTINUE ! Physical boundary.
END SELECT
!
SELECT CASE(mthbc(2))
CASE(2,999) ! Periodic or Internal boundary.
!
! Face #2, Red:
  IF(ipass==1 .OR. ipass==0) THEN
    qf(mx(1)+1,3:mx(2)-1:2,1:nrvars) &
      = a1*qf(mx(1)-1,1:mx(2)-3:2,1:nrvars) &
      + a2*qf(mx(1)  ,2:mx(2)-2:2,1:nrvars) &
      + a3*qc(mb(1,2)+1,mb(2,1)+1:mb(2,2):1,1:nrvars)
!
    qf(mx(1)+1,1,1:nrvars) &
      = a1*qf(mx(1)-1,1,1:nrvars) &
      + a2*qf(mx(1)  ,1,1:nrvars) &
      + a3*(b1*qc(mb(1,2)+1,mb(2,1)-1,1:nrvars) &
      +     b2*qc(mb(1,2)+1,mb(2,1)  ,1:nrvars) &
      +     b3*qc(mb(1,2)+1,mb(2,1)+1,1:nrvars))
  END IF
!
! Face #2, Black:
  IF(ipass==2 .OR. ipass==0) THEN
    qf(mx(1)+1,2:mx(2)-2:2,1:nrvars) &
      = a1*qf(mx(1)-1,4:mx(2)  :2,1:nrvars) &
      + a2*qf(mx(1)  ,3:mx(2)-1:2,1:nrvars) &
      + a3*qc(mb(1,2)+1,mb(2,1):mb(2,2)-1:1,1:nrvars)
!
    qf(mx(1)+1,mx(2),1:nrvars) &
      = a1*qf(mx(1)-1,mx(2),1:nrvars) &
      + a2*qf(mx(1)  ,mx(2),1:nrvars) &
      + a3*(b1*qc(mb(1,2)+1,mb(2,2)+1,1:nrvars) &
      +     b2*qc(mb(1,2)+1,mb(2,2)  ,1:nrvars) &
      +     b3*qc(mb(1,2)+1,mb(2,2)-1,1:nrvars))
  END IF
!
CASE DEFAULT
  CONTINUE ! Physical boundary.
END SELECT
!
SELECT CASE(mthbc(3))
CASE(2,999) ! Periodic or Internal boundary.
!
! Face #3, Red:
  IF(ipass==1 .OR. ipass==0) THEN
    qf(2:mx(1)-2:2,0,1:nrvars) &
      = a1*qf(4:mx(1)  :2,2,1:nrvars) &
      + a2*qf(3:mx(1)-1:2,1,1:nrvars) &
      + a3*qc(mb(1,1):mb(1,2)-1,mb(2,1)-1,1:nrvars)
!
    qf(mx(1),0,1:nrvars) &
      = a1*qf(mx(1),2,1:nrvars) &
      + a2*qf(mx(1),1,1:nrvars) &
      + a3*(b1*qc(mb(1,2)+1,mb(2,1)-1,1:nrvars) &
      +     b2*qc(mb(1,2)  ,mb(2,1)-1,1:nrvars) &
      +     b3*qc(mb(1,2)-1,mb(2,1)-1,1:nrvars))
  END IF
!
! Face #3, Black:
  IF(ipass==2 .OR. ipass==0) THEN
    qf(3:mx(1)-1:2,0,1:nrvars) &
      = a1*qf(1:mx(1)-3:2,2,1:nrvars) &
      + a2*qf(2:mx(1)-2:2,1,1:nrvars) &
      + a3*qc(mb(1,1)+1:mb(1,2),mb(2,1)-1,1:nrvars)
!
    qf(1,0,1:nrvars) &
      = a1*qf(1,2,1:nrvars) &
      + a2*qf(1,1,1:nrvars) &
      + a3*(b1*qc(mb(1,1)-1,mb(2,1)-1,1:nrvars) &
      +     b2*qc(mb(1,1)  ,mb(2,1)-1,1:nrvars) &
      +     b3*qc(mb(1,1)+1,mb(2,1)-1,1:nrvars))
  END IF
!
CASE DEFAULT
  CONTINUE
END SELECT
!
SELECT CASE(mthbc(4))
CASE(2,999) ! Periodic or Internal boundary.
!
! Face #4, Red:
  IF(ipass==1 .OR. ipass==0) THEN
    qf(3:mx(1)-1:2,mx(2)+1,1:nrvars) &
      = a1*qf(1:mx(1)-3:2,mx(2)-1,1:nrvars) &
      + a2*qf(2:mx(1)-2:2,mx(2)  ,1:nrvars) &
      + a3*qc(mb(1,1)+1:mb(1,2),mb(2,2)+1,1:nrvars)
!
    qf(1,mx(2)+1,1:nrvars) &
      = a1*qf(1,mx(2)-1,1:nrvars) &
      + a2*qf(1,mx(2)  ,1:nrvars) &
      + a3*(b1*qc(mb(1,1)-1,mb(2,2)+1,1:nrvars) &
      +     b2*qc(mb(1,1)  ,mb(2,2)+1,1:nrvars) &
      +     b3*qc(mb(1,1)+1,mb(2,2)+1,1:nrvars))
  END IF
!
! Face #4, Black:
  IF(ipass==2 .OR. ipass==0) THEN
    qf(2:mx(1)-2:2,mx(2)+1,1:nrvars) &
      = a1*qf(4:mx(1)  :2,mx(2)-1,1:nrvars) &
      + a2*qf(3:mx(1)-1:2,mx(2)  ,1:nrvars) &
      + a3*qc(mb(1,1):mb(1,2)-1,mb(2,2)+1,1:nrvars)
!
    qf(mx(1),mx(2)+1,1:nrvars) &
      =          a1*qf(mx(1),mx(2)-1,1:nrvars) &
      + a2*qf(mx(1),mx(2)  ,1:nrvars) &
      + a3*(b1*qc(mb(1,2)+1,mb(2,2)+1,1:nrvars) &
      +     b2*qc(mb(1,2)  ,mb(2,2)+1,1:nrvars) &
      +     b3*qc(mb(1,2)-1,mb(2,2)+1,1:nrvars))
  END IF
!
CASE DEFAULT
  CONTINUE ! Physical boundary.
END SELECT
!
END SUBROUTINE Interpolate2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION TransferBC(grid1,grid2,dummy)
USE NodeinfoDef
USE TreeOps, ONLY: err_ok
IMPLICIT NONE
!
! Transfer boundary values among grids on same level
!
TYPE(nodeinfo):: grid1, grid2
TYPE(funcparam):: dummy
!
TransferBC=err_ok
!
! Check for inactive grids awaiting garbage collection
IF(            grid1%tobedeleted .OR.        grid2%tobedeleted  &
   .OR. (.NOT. grid1%activegrid) .OR. (.NOT. grid2%activegrid)) RETURN
!
! Look for overlap of ghost cell regions of grid1 and grid2.  Order doesn't
! matter, as transfer goes both ways:
CALL GhostOverlap(grid1,grid2)
!
END FUNCTION TransferBC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE GhostOverlap(grid1,grid2,grid2offset)
USE NodeinfoDef
IMPLICIT NONE
!
TYPE(nodeinfo):: grid1
TYPE(nodeinfo):: grid2
INTEGER, DIMENSION(1:maxdims), INTENT(IN), OPTIONAL:: grid2offset
!
! Assumes ghost layer of size mbc:
!
INTEGER:: grid1bdry, grid2bdry, m, mbc, n, nrvars
INTEGER, DIMENSION(1:maxdims,1:2):: mg1, mg2, ml1, ml2, ovg, ovg1, ovg2
!
mbc = grid2%mbc; nrvars = grid1%nrvars
!
mg1 = 1
mg1(1:ndims,1:2) = grid1%mglobal(1:ndims,1:2)
!
mg2 = 1
IF(PRESENT(grid2offset)) THEN
  DO n = 1, ndims
    mg2(n,1:2) = grid2%mglobal(n,1:2)+grid2offset(n)
  END DO
ELSE
  mg2(1:ndims,1:2) = grid2%mglobal(1:ndims,1:2)
END IF
!
! Calculate overlapping region, taking into account the extra ghost layer(s)
! with the addition of +/-mbc:
ovg = 1
DO n = 1, ndims
  ovg(n,1) = MAX(mg1(n,1)-mbc,mg2(n,1)-mbc)
  ovg(n,2) = MIN(mg1(n,2)+mbc,mg2(n,2)+mbc)
END DO
!
! Check for the empty set:
IF(ANY((ovg(1:ndims,2)-ovg(1:ndims,1))<0)) RETURN
!
! Nonempty; begin transfer:
!
! ovg1 is the global overlap of ovg and mg2;
! ovg2 is the global overlap of ovg and mg1;
! These must have nonempty overlap:
ovg1 = 1; ovg2 = 1
DO n = 1, ndims
  ovg1(n,1) = MAX(ovg(n,1),mg2(n,1))
  ovg1(n,2) = MIN(ovg(n,2),mg2(n,2))
!
  ovg2(n,1) = MAX(ovg(n,1),mg1(n,1))
  ovg2(n,2) = MIN(ovg(n,2),mg1(n,2))
END DO
!
! Convert to local coordinates by subtracting off the global offsets:
ml1 = 1; ml2 = 1
ml1(1:ndims,1) = ovg1(1:ndims,1)-mg1(1:ndims,1)+1 ! Grid 1 local
ml1(1:ndims,2) = ovg1(1:ndims,2)-mg1(1:ndims,1)+1 ! Grid 1 local
ml2(1:ndims,1) = ovg1(1:ndims,1)-mg2(1:ndims,1)+1 ! Grid 2 local
ml2(1:ndims,2) = ovg1(1:ndims,2)-mg2(1:ndims,1)+1 ! Grid 2 local
!
! Transfer:
  grid1%q(ml1(1,1):ml1(1,2),ml1(2,1):ml1(2,2),ml1(3,1):ml1(3,2),1:nrvars) &
= grid2%q(ml2(1,1):ml2(1,2),ml2(2,1):ml2(2,2),ml2(3,1):ml2(3,2),1:nrvars)
!
! Convert to local coordinates by subtracting off the global offsets:
ml1(1:ndims,1) = ovg2(1:ndims,1)-mg1(1:ndims,1)+1 ! Grid 1 local
ml1(1:ndims,2) = ovg2(1:ndims,2)-mg1(1:ndims,1)+1 ! Grid 1 local
ml2(1:ndims,1) = ovg2(1:ndims,1)-mg2(1:ndims,1)+1 ! Grid 2 local
ml2(1:ndims,2) = ovg2(1:ndims,2)-mg2(1:ndims,1)+1 ! Grid 2 local
!
! Transfer:
  grid2%q(ml2(1,1):ml2(1,2),ml2(2,1):ml2(2,2),ml2(3,1):ml2(3,2),1:nrvars) &
= grid1%q(ml1(1,1):ml1(1,2),ml1(2,1):ml1(2,2),ml1(3,1):ml1(3,2),1:nrvars)
!
! If the mesh is constructed, no need to procede:
! This code will be used in BSAM 2.0:
IF(meshbuildcomplete) RETURN
!
! In the mesh-build phase, record the level numbers of the neighbors:
!
! Recalculate overlapping region, taking into account the only 1 ghost layer:
!
mg1(1:ndims,1) = (mg1(1:ndims,1)+1)/2
mg1(1:ndims,2) =  mg1(1:ndims,2)   /2
mg2(1:ndims,1) = (mg2(1:ndims,1)+1)/2
mg2(1:ndims,2) =  mg2(1:ndims,2)   /2
!
ovg = 1
DO n = 1, ndims
  ovg(n,1) = MAX(mg1(n,1)-1,mg2(n,1)-1)
  ovg(n,2) = MIN(mg1(n,2)+1,mg2(n,2)+1)
END DO
!
ovg1 = 1; ovg2 = 1
DO n = 1, ndims
  ovg1(n,1) = MAX(ovg(n,1),mg2(n,1))
  ovg1(n,2) = MIN(ovg(n,2),mg2(n,2))
!
  ovg2(n,1) = MAX(ovg(n,1),mg1(n,1))
  ovg2(n,2) = MIN(ovg(n,2),mg1(n,2))
END DO
!
ml1 = 1; ml2 = 1
ml1(1:ndims,1) = ovg1(1:ndims,1)-mg1(1:ndims,1)+1 ! Grid 1 local
ml1(1:ndims,2) = ovg1(1:ndims,2)-mg1(1:ndims,1)+1 ! Grid 1 local
ml2(1:ndims,1) = ovg1(1:ndims,1)-mg2(1:ndims,1)+1 ! Grid 2 local
ml2(1:ndims,2) = ovg1(1:ndims,2)-mg2(1:ndims,1)+1 ! Grid 2 local
!
  grid1%levellandscape(ml1(1,1):ml1(1,2),ml1(2,1):ml1(2,2),ml1(3,1):ml1(3,2)) &
= grid2%levellandscape(ml2(1,1):ml2(1,2),ml2(2,1):ml2(2,2),ml2(3,1):ml2(3,2))
!
ml1(1:ndims,1) = ovg2(1:ndims,1)-mg1(1:ndims,1)+1 ! Grid 1 local
ml1(1:ndims,2) = ovg2(1:ndims,2)-mg1(1:ndims,1)+1 ! Grid 1 local
ml2(1:ndims,1) = ovg2(1:ndims,1)-mg2(1:ndims,1)+1 ! Grid 2 local
ml2(1:ndims,2) = ovg2(1:ndims,2)-mg2(1:ndims,1)+1 ! Grid 2 local
!
  grid2%levellandscape(ml2(1,1):ml2(1,2),ml2(2,1):ml2(2,2),ml2(3,1):ml2(3,2)) &
= grid1%levellandscape(ml1(1,1):ml1(1,2),ml1(2,1):ml1(2,2),ml1(3,1):ml1(3,2))
!
END SUBROUTINE GhostOverlap
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION SelfPeriodicBC(info,timestepparam)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok
IMPLICIT NONE
!
TYPE(nodeinfo):: info
TYPE(funcparam):: timestepparam
!
INTEGER:: ibc, level, mbc, n, nrvars
INTEGER, DIMENSION(1:maxdims):: cmx, mx
!
SelfPeriodicBC = err_OK
!
level = info%level
mbc = info%mbc
nrvars = info%nrvars
!
mx = 1; cmx = 1
mx(1:ndims) = info%mx(1:ndims)
cmx(1:ndims) = mx(1:ndims)/2
!
! Check for inactive grid awaiting garbage collection
IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
! Note: qpo is the quasi-periodic offset in q along dimension 1.
!
SELECT CASE(ndims)
!
! Two dimensions:
  CASE (2)
    IF(info%mthbc(1)==2 .AND. info%mthbc(2)==2) THEN
      DO ibc = 1, mbc
        DO n = 1, nrvars
            info%q(mx(1)+  ibc,1-mbc:mx(2)+mbc,1,n) &
          = info%q(        ibc,1-mbc:mx(2)+mbc,1,n)+qpo(n)
            info%q(      1-ibc,1-mbc:mx(2)+mbc,1,n) &
          = info%q(mx(1)+1-ibc,1-mbc:mx(2)+mbc,1,n)-qpo(n)
        END DO
      END DO
!
      IF(.NOT. meshbuildcomplete) THEN
        info%levellandscape(       0,1:cmx(2),1) = level
        info%levellandscape(cmx(1)+1,1:cmx(2),1) = level
      END IF
    END IF
!
    IF(info%mthbc(3)==2 .AND. info%mthbc(4)==2) THEN
      DO ibc = 1, mbc
          info%q(1-mbc:mx(1)+mbc,mx(2)+  ibc,1,1:nrvars) &
        = info%q(1-mbc:mx(1)+mbc,        ibc,1,1:nrvars)
          info%q(1-mbc:mx(1)+mbc,      1-ibc,1,1:nrvars) &
        = info%q(1-mbc:mx(1)+mbc,mx(2)+1-ibc,1,1:nrvars)
      END DO
!
      IF(.NOT. meshbuildcomplete) THEN
        info%levellandscape(1:cmx(1),       0,1) = level
        info%levellandscape(1:cmx(1),cmx(2)+1,1) = level
      END IF
    END IF
!
! If the whole domain is refined then it is possible that the grid touches
! corner-to-corner. We record this unlikely event:
    IF(info%mthbc(1)==2 .AND. info%mthbc(2)==2 .AND. &
       info%mthbc(3)==2 .AND. info%mthbc(4)==2 .AND. &
       (.NOT. meshbuildcomplete)) THEN
      info%levellandscape(       0,       0,1) = level
      info%levellandscape(cmx(1)+1,       0,1) = level
      info%levellandscape(       0,cmx(2)+1,1) = level
      info%levellandscape(cmx(1)+1,cmx(2)+1,1) = level
    END IF
!
! Three dimensions:
  CASE (3)
    IF(info%mthbc(1)==2 .AND. info%mthbc(2)==2) THEN
      DO ibc = 1, mbc
        DO n = 1, nrvars
          info%q(mx(1)+ibc,:,:,n) = info%q(        ibc,:,:,n)+qpo(n)
          info%q(    1-ibc,:,:,n) = info%q(mx(1)+1-ibc,:,:,n)-qpo(n)
        END DO
      END DO
    END IF
    IF(info%mthbc(3)==2 .AND. info%mthbc(4)==2) THEN
      DO ibc = 1, mbc
        info%q(:,mx(2)+ibc,:,1:nrvars) = info%q(:,        ibc,:,1:nrvars)
        info%q(:,    1-ibc,:,1:nrvars) = info%q(:,mx(2)+1-ibc,:,1:nrvars)
      END DO
    END IF
    IF(info%mthbc(5)==2 .AND. info%mthbc(6)==2) THEN
      DO ibc = 1, mbc
        info%q(:,:,mx(3)+ibc,1:nrvars) = info%q(:,:,        ibc,1:nrvars)
        info%q(:,:,    1-ibc,1:nrvars) = info%q(:,:,mx(3)+1-ibc,1:nrvars)
      END DO
    END IF
  CASE DEFAULT
    PRINT *,'SelfPeriodicBC: Sorry. Code only works up to ndims=3.'
    STOP
END SELECT
END FUNCTION SelfPeriodicBC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION TransferPeriodicBC(grid1,grid2,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok
IMPLICIT NONE
!
! Transfer periodic boundary values among grids on same level:
!
TYPE(nodeinfo):: grid1, grid2
TYPE(funcparam):: dummy
!
INTEGER:: offset, polarity
INTEGER, DIMENSION(1:maxdims):: grid2offset
!
TransferPeriodicBC = err_ok
!
! Check for inactive grids awaiting garbage collection:
IF(            grid1%tobedeleted .OR.        grid2%tobedeleted  &
   .OR. (.NOT. grid1%activegrid) .OR. (.NOT. grid2%activegrid)) RETURN
!
! Look for periodic overlap of ghost cell regions of grid1 and grid2:
!
CALL GetPeriodicOffsets(grid2%level,grid2%nrvars)
!
! Apply +/- poffsets to grid2:
DO polarity = -1, 1, 2
  DO offset = 1, nperiodicoffsets
    grid2offset(1:maxdims) = polarity*poffset(1:maxdims,offset)
    CALL GhostOverlap(grid1,grid2,grid2offset)
  END DO
END DO
!
END FUNCTION TransferPeriodicBC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION SetPhysicalBC(info,dummy)
USE NodeinfoDef
USE TreeOps, ONLY: err_ok
IMPLICIT NONE
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
INTEGER:: ipass
!
ipass = dummy%iswitch ! Red-black switch when used:
!
SetPhysicalBC = err_ok
!
! Check for inactive grids awaiting garbage collection
IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
SELECT CASE(ndims)
  CASE (2)
    CALL SetPhysicalBC2D(info,ipass)      
  CASE (3)
    CALL SetPhysicalBC3D(info,ipass)
  CASE DEFAULT
    PRINT *, 'SetPhysicalBC: Only ndims = 2,3 are supported.'
END SELECT
!
END FUNCTION SetPhysicalBC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE SetPhysicalBC2D(info,ipass)
USE NodeinfoDef
USE Problem, ONLY: UserBC2D
IMPLICIT NONE
!
TYPE(nodeinfo):: info
INTEGER, INTENT(IN):: ipass

INTEGER:: ll, level, mbc, nrvars
INTEGER, DIMENSION(1:2):: cmx, mx, ul
INTEGER, DIMENSION(1:4):: mthbc
!
mbc = info%mbc; nrvars = info%nrvars; level = info%level
ll = 1-mbc
mx(1:2) = info%mx(1:2)
cmx(1:2) = mx(1:2)/2
ul(1:2) = mx(1:2)+mbc
mthbc(1:4) = info%mthbc(1:4)
!
1001 FORMAT('SetPhysicalBC2D: User boundary condition ', i4, ' on boundary ', i1, //, &
            '         not coded in SetPhysicalBC2D.')
!
! Left Boundary:
SELECT CASE(mthbc(1))
  CASE(2,999) ! Periodic (2), internal (999) boundary condition. Coded below:
    CONTINUE
  CASE(10)    ! User boundary condition. Call Problem module:
    CALL UserBC2D(info%q(ll:ul(1),ll:ul(2),1,1:nrvars),ll,mx,nrvars,mbc,1)
    IF(.NOT. meshbuildcomplete) info%levellandscape(0,0:cmx(2)+1,1) = level
  CASE DEFAULT
    PRINT 1001, mthbc(1), 1
    STOP
END SELECT
!
! Right Boundary:
SELECT CASE(mthbc(2))
  CASE(2,999) ! Periodic (2), internal (999) boundary condition. Coded below:
    CONTINUE
  CASE(10)    ! User boundary condition. Call Problem module:
    CALL UserBC2D(info%q(ll:ul(1),ll:ul(2),1,1:nrvars),ll,mx,nrvars,mbc,2)
    IF(.NOT. meshbuildcomplete) info%levellandscape(cmx(1)+1,0:cmx(2)+1,1) = level
  CASE DEFAULT
    PRINT 1001, mthbc(2), 2
    STOP
END SELECT
!
! Bottom Boundary:
SELECT CASE(mthbc(3))
  CASE(2,999) ! Periodic (2), internal (999) boundary condition. Coded below:
    CONTINUE
  CASE(10)    ! User boundary condition. Call Problem module:
    CALL UserBC2D(info%q(ll:ul(1),ll:ul(2),1,1:nrvars),ll,mx,nrvars,mbc,3)
    IF(.NOT. meshbuildcomplete) info%levellandscape(0:cmx(1)+1,0,1) = level
  CASE DEFAULT
    PRINT 1001, mthbc(3), 3
    STOP
END SELECT
!
! Top Boundary:
SELECT CASE(mthbc(4))
  CASE(2,999) ! Periodic (2), internal (999) boundary condition. Coded below:
    CONTINUE
  CASE(10)    ! User boundary condition. Call Problem module:
    CALL UserBC2D(info%q(ll:ul(1),ll:ul(2),1,1:nrvars),ll,mx,nrvars,mbc,4)
    IF(.NOT. meshbuildcomplete) info%levellandscape(0:cmx(1)+1,cmx(2)+1,1) = level
  CASE DEFAULT
    PRINT 1001, mthbc(4), 4
    STOP
END SELECT
!
END SUBROUTINE SetPhysicalBC2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE SetPhysicalBC3D(info,ipass)
USE NodeinfoDef
USE Problem, ONLY: UserBC3D
IMPLICIT NONE
!
TYPE(nodeinfo):: info
INTEGER, INTENT(IN):: ipass

INTEGER:: ll, mbc, nrvars
INTEGER, DIMENSION(1:3):: mx, ul
INTEGER, DIMENSION(1:6):: mthbc
!
mbc = info%mbc; nrvars = info%nrvars
ll = 1-mbc
mx(1:3) = info%mx(1:3)
ul(1:3) = mx(1:3)+mbc
mthbc(1:6) = info%mthbc(1:6)
!
1001 FORMAT('SetPhysicalBC3D: User boundary condition ', i4, ' on boundary ', i1, //, &
            '         not coded in SetPhysicalBC3D.')
!
! Left Boundary along dimension 1:
SELECT CASE(mthbc(1))
  CASE(2,999) ! Periodic, internal boundary condition. Coded below:
    CONTINUE
  CASE(10)    ! User boundary condition. Call Problem module:
    CALL UserBC3D(info%q(ll:ul(1),ll:ul(2),ll:ul(3),1:nrvars),ll,mx,nrvars,mbc,1)
  CASE DEFAULT
    PRINT 1001, mthbc(1), 1
    STOP
END SELECT
!
! Right Boundary along dimension 1:
SELECT CASE(mthbc(2))
  CASE(2,999) ! Periodic, internal boundary condition. Coded below:
    CONTINUE
  CASE(10)    ! User boundary condition. Call Problem module:
    CALL UserBC3D(info%q(ll:ul(1),ll:ul(2),ll:ul(3),1:nrvars),ll,mx,nrvars,mbc,2)
  CASE DEFAULT
    PRINT 1001, mthbc(2), 2
    STOP
END SELECT
!
! Left Boundary along dimension 2:
SELECT CASE(mthbc(3))
  CASE(2,999) ! Periodic, internal boundary condition. Coded below:
    CONTINUE
  CASE(10)    ! User boundary condition. Call Problem module:
    CALL UserBC3D(info%q(ll:ul(1),ll:ul(2),ll:ul(3),1:nrvars),ll,mx,nrvars,mbc,3)
  CASE DEFAULT
    PRINT 1001, mthbc(3), 3
    STOP
END SELECT
!
! Right Boundary along dimension 2:
SELECT CASE(mthbc(4))
  CASE(2,999) ! Periodic, internal boundary condition. Coded below:
    CONTINUE
  CASE(10)    ! User boundary condition. Call Problem module:
    CALL UserBC3D(info%q(ll:ul(1),ll:ul(2),ll:ul(3),1:nrvars),ll,mx,nrvars,mbc,4)
  CASE DEFAULT
    PRINT 1001, mthbc(4), 4
    STOP
END SELECT
!
! Left Boundary along dimension 3:
SELECT CASE(mthbc(5))
  CASE(2,999) ! Periodic, internal boundary condition. Coded below:
    CONTINUE
  CASE(10)    ! User boundary condition. Call Problem module:
    CALL UserBC3D(info%q(ll:ul(1),ll:ul(2),ll:ul(3),1:nrvars),ll,mx,nrvars,mbc,5)
  CASE DEFAULT
    PRINT 1001, mthbc(5), 5
    STOP
END SELECT
!
! Right Boundary along dimension 3:
SELECT CASE(mthbc(6))
  CASE(2,999) ! Periodic, internal boundary condition. Coded below:
    CONTINUE
  CASE(10)    ! User boundary condition. Call Problem module:
    CALL UserBC3D(info%q(ll:ul(1),ll:ul(2),ll:ul(3),1:nrvars),ll,mx,nrvars,mbc,6)
  CASE DEFAULT
    PRINT 1001, mthbc(6), 6
    STOP
END SELECT
!
END SUBROUTINE SetPhysicalBC3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE PeriodicSetup(rootinfo)
USE NodeInfoDef
USE BSAMStorage, ONLY: AllocPeriodicBCStorage
IMPLICIT NONE
!
TYPE(nodeinfo):: rootinfo
!
INTEGER:: np, nrvars
INTEGER, DIMENSION(1:2*maxdims):: mthbc
!
mthbc = 10; mthbc(1:2*ndims) = rootinfo%mthbc(1:2*ndims)
rootinfo%mthbc = mthbc
nrvars = rootinfo%nrvars
!
periodicboundaryconditions = .TRUE.
!
IF(mthbc(1)==2) THEN
  IF(mthbc(3)==2) THEN
    IF(mthbc(5)==2) THEN
      nperiodicoffsets = 13
      np = nperiodicoffsets
      CALL AllocPeriodicBCStorage(np,nrvars)
      periodicoffsetindex(1:np) = (/1,2,3,4,5,6,7,8,9,10,11,12,13/)
    ELSE
      nperiodicoffsets = 4
      np = nperiodicoffsets
      CALL AllocPeriodicBCStorage(np,nrvars)
      periodicoffsetindex(1:np) = (/1,2,4,5/)
    END IF
  ELSE
    IF(mthbc(5)==2) THEN
      nperiodicoffsets = 4
      np = nperiodicoffsets
      CALL AllocPeriodicBCStorage(np,nrvars)
      periodicoffsetindex(1:np) = (/1,3,8,9/)
    ELSE
      nperiodicoffsets = 1
      np = nperiodicoffsets
      CALL AllocPeriodicBCStorage(np,nrvars)
      periodicoffsetindex(1:np) = (/1/)
    END IF
  END IF
ELSE
  IF(mthbc(3)==2) THEN
    IF(mthbc(5)==2) THEN
      nperiodicoffsets = 4
      np = nperiodicoffsets
      CALL AllocPeriodicBCStorage(np,nrvars)
      periodicoffsetindex(1:np) = (/2,3,6,7/)
    ELSE
      nperiodicoffsets = 1
      np = nperiodicoffsets
      CALL AllocPeriodicBCStorage(np,nrvars)
      periodicoffsetindex(1:np) = (/2/)
    END IF
  ELSE
    IF(mthbc(5)==2) THEN
      nperiodicoffsets = 1
      np = nperiodicoffsets
      CALL AllocPeriodicBCStorage(np,nrvars)
      periodicoffsetindex(1:np) = (/3/)
    ELSE
      periodicboundaryconditions = .FALSE.
    END IF
  END IF
END IF
!
END SUBROUTINE PeriodicSetup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE GetPeriodicOffsets(level,nrvars)
USE NodeInfoDef
IMPLICIT NONE
!
INTEGER, INTENT(IN):: level
INTEGER, OPTIONAL:: nrvars
!
INTEGER:: offset
!
DO offset = 1, nperiodicoffsets
  SELECT CASE(periodicoffsetindex(offset))
  CASE(1)
    poffset(1,offset) =  mxmax(level,1)
    poffset(2,offset) =  0
    poffset(3,offset) =  0
    IF(PRESENT(nrvars)) qoffset(offset,1:nrvars) =  qpo(1:nrvars)
  CASE(2)
    poffset(1,offset) =  0
    poffset(2,offset) =  mxmax(level,2)
    poffset(3,offset) =  0
    IF(PRESENT(nrvars)) qoffset(offset,1:nrvars) =  0.0_r8
  CASE(3)
    poffset(1,offset) =  0
    poffset(2,offset) =  0
    poffset(3,offset) =  mxmax(level,3)
    IF(PRESENT(nrvars)) qoffset(offset,1:nrvars) =  0.0_r8
  CASE(4)
    poffset(1,offset) =  mxmax(level,1) 
    poffset(2,offset) =  mxmax(level,2)
    poffset(3,offset) =  0
    IF(PRESENT(nrvars)) qoffset(offset,1:nrvars) =  qpo(1:nrvars)
  CASE(5)
    poffset(1,offset) = -mxmax(level,1) 
    poffset(2,offset) =  mxmax(level,2)
    poffset(3,offset) =  0
    IF(PRESENT(nrvars)) qoffset(offset,1:nrvars) = -qpo(1:nrvars)
  CASE(6)
    poffset(1,offset) =  0
    poffset(2,offset) =  mxmax(level,2) 
    poffset(3,offset) =  mxmax(level,3)
    IF(PRESENT(nrvars)) qoffset(offset,1:nrvars) =  0.0_r8
  CASE(7)
    poffset(1,offset) =  0
    poffset(2,offset) =  mxmax(level,2) 
    poffset(3,offset) = -mxmax(level,3)
    IF(PRESENT(nrvars)) qoffset(offset,1:nrvars) =  0.0_r8
  CASE(8)
    poffset(1,offset) =  mxmax(level,1)
    poffset(2,offset) =  0
    poffset(3,offset) =  mxmax(level,3)
    IF(PRESENT(nrvars)) qoffset(offset,1:nrvars) =  qpo(1:nrvars)
  CASE(9)
    poffset(1,offset) = -mxmax(level,1) 
    poffset(2,offset) =  0
    poffset(3,offset) =  mxmax(level,3)
    IF(PRESENT(nrvars)) qoffset(offset,1:nrvars) = -qpo(1:nrvars)
  CASE(10)
    poffset(1,offset) =  mxmax(level,1) 
    poffset(2,offset) =  mxmax(level,2)
    poffset(3,offset) =  mxmax(level,3)
    IF(PRESENT(nrvars)) qoffset(offset,1:nrvars) =  qpo(1:nrvars)
  CASE(11)
    poffset(1,offset) =  mxmax(level,1) 
    poffset(2,offset) =  mxmax(level,2)
    poffset(3,offset) = -mxmax(level,3)
    IF(PRESENT(nrvars)) qoffset(offset,1:nrvars) =  qpo(1:nrvars)
  CASE(12)
    poffset(1,offset) = -mxmax(level,1) 
    poffset(2,offset) =  mxmax(level,2)
    poffset(3,offset) =  mxmax(level,3)
    IF(PRESENT(nrvars)) qoffset(offset,1:nrvars) = -qpo(1:nrvars)
  CASE(13)
    poffset(1,offset) =  mxmax(level,1) 
    poffset(2,offset) = -mxmax(level,2)
    poffset(3,offset) =  mxmax(level,3)
    IF(PRESENT(nrvars)) qoffset(offset,1:nrvars) =  qpo(1:nrvars)
  END SELECT
END DO
!
END SUBROUTINE GetPeriodicOffsets
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION GetCoarseGhostPoints(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: GetParentInfo, err_ok
IMPLICIT NONE
!
! Dimensionally invariant routine to fill coarse grid ghost points.
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
TYPE(nodeinfo), POINTER:: parent
INTEGER:: ierror, nrvars
INTEGER, DIMENSION(1:maxdims):: cmx, ih, il, jh, jl, mx
INTEGER, DIMENSION(1:maxdims,1:2):: mb
!
GetCoarseGhostPoints = err_ok
!
IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
nrvars = info%nrvars
mx = 1; mx(1:ndims) = info%mx(1:ndims)
cmx = 1; cmx(1:ndims) = mx(1:ndims)/2
il = 1; ih = 1; il(1:ndims) = 0; ih(1:ndims) = cmx(1:ndims)+1
mb = 1; mb(1:ndims,1:2) = info%mbounds(1:ndims,1:2)
jl = 1; jh = 1; jl(1:ndims) = mb(1:ndims,1)-1; jh(1:ndims) = mb(1:ndims,2)+1
!
ierror = GetParentInfo(parent)
!
! Face #1:
      info%qc(il(1),il(2):ih(2),il(3):ih(3),1:nrvars) &
  = parent%q( jl(1),jl(2):jh(2),jl(3):jh(3),1:nrvars)
!
! Face #2:
      info%qc(ih(1),il(2):ih(2),il(3):ih(3),1:nrvars) &
  = parent%q( jh(1),jl(2):jh(2),jl(3):jh(3),1:nrvars)
!
IF(ndims>=2) THEN
!
! Face #3:
        info%qc(     1:cmx(1)  ,il(2),il(3):ih(3),1:nrvars) &
    = parent%q(mb(1,1): mb(1,2),jl(2),jl(3):jh(3),1:nrvars)
!
! Face #4:
        info%qc(     1:cmx(1)  ,ih(2),il(3):ih(3),1:nrvars) &
    = parent%q(mb(1,1): mb(1,2),jh(2),jl(3):jh(3),1:nrvars)
!  
  IF(ndims>=3) THEN
!
! Face #5:
          info%qc(     1:cmx(1) ,       1:cmx(2)  ,il(3),1:nrvars) &
      = parent%q(mb(1,1): mb(1,2),mb(2,1): mb(2,2),jl(3),1:nrvars)
!
! Face #6:
          info%qc(     1:cmx(1) ,       1:cmx(2)  ,ih(3),1:nrvars) &
      = parent%q(mb(1,1): mb(1,2),mb(2,1): mb(2,2),jh(3),1:nrvars)
  END IF
END IF
!
END FUNCTION GetCoarseGhostPoints
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE GetFaceIndex(mx,mb,dim,parity,fi)
IMPLICIT NONE
!
INTEGER, DIMENSION(1:3), INTENT(IN)::mx
INTEGER, DIMENSION(1:3,1:2), INTENT(IN)::mb
INTEGER, INTENT(IN):: dim
INTEGER, INTENT(IN):: parity
INTEGER, DIMENSION(1:5), INTENT(OUT):: fi
!
SELECT CASE(parity)
  CASE(1)
    fi(1) = 0
    fi(2) = 2
    fi(3) = 1
    fi(4) = mb(dim,parity)-1
    fi(5) = mb(dim,parity)
  CASE(2)
    fi(1) = mx(dim)+1
    fi(2) = mx(dim)-1
    fi(3) = mx(dim)
    fi(4) = mb(dim,parity)+1
    fi(5) = mb(dim,parity)
  CASE DEFAULT
    PRINT *, 'GetFaceIndex: Parity error'
    STOP
END SELECT
!
END SUBROUTINE GetFaceIndex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE Interpolate3D(qf,qc,mx,nrvars,mb,mthbc,ipass)
USE NodeInfoDef, ONLY: r8
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:,0:,1:), INTENT(IN OUT):: qf
REAL(KIND=r8), DIMENSION(0:,0:,0:,1:), INTENT(IN):: qc
INTEGER, DIMENSION(1:3), INTENT(IN):: mx
INTEGER, INTENT(IN):: nrvars
INTEGER, DIMENSION(1:3,1:2), INTENT(IN):: mb
INTEGER, DIMENSION(1:6), INTENT(IN):: mthbc
INTEGER, INTENT(IN):: ipass
!
LOGICAL, DIMENSION(1:6):: cycle_face
INTEGER:: fc, parity_1, parity_2, parity_3
INTEGER, DIMENSION(1:3,1:5):: fi
!
! ipass = 1: red (or odd-odd and even-even) squares should be updated,
! ipass = 2: black (or even-odd and odd-even) squares should be updated.
! ipass = 0: both red and black squares should be updated.
!
! Faces that are on physical boundaries have mthbc(face) equal to a number 
! other than 2 and 999. If the face is on a physical boundary its ghost points
! aren't found by interpolation. 
! interpol
!
! (1) Interpolate the 8 corners of each box.  We only interpolate on the
!     physical boundaries if ipass = 0, 2:
parity_1c_loop: DO parity_1 = 1, 2
  IF(mthbc(parity_1)/=2 .AND. mthbc(parity_1)/=999 .AND. ipass==1) &
    CYCLE parity_1c_loop
  CALL GetFaceIndex(mx,mb,1,parity_1,fi(1,1:5))
  parity_2c_loop: DO parity_2 = 1, 2
    IF(mthbc(parity_2+2)/=2 .AND. mthbc(parity_2+2)/=999 .AND. ipass==1) &
      CYCLE parity_2c_loop
    CALL GetFaceIndex(mx,mb,2,parity_2,fi(2,1:5))
    parity_3c_loop: DO parity_3 = 1, 2
      IF(mthbc(parity_3+4)/=2 .AND. mthbc(parity_3+4)/=999 .AND. ipass==1) &
        CYCLE parity_3c_loop
      CALL GetFaceIndex(mx,mb,3,parity_3,fi(3,1:5))
           qf(fi(1,1),fi(2,1),fi(3,1),1:nrvars) = &
        a1*qf(fi(1,2),fi(2,2),fi(3,2),1:nrvars) &
      + a2*qf(fi(1,3),fi(2,3),fi(3,3),1:nrvars) &
      + a3*qc(fi(1,4),fi(2,4),fi(3,4),1:nrvars)
    END DO parity_3c_loop
  END DO parity_2c_loop
END DO parity_1c_loop
!
! (2) Interpolate the 12 edges of each box.  We only interpolate on the
!     physical boundaries if ipass = 0, 2:
!
! Edges 1--4 (coordinates 1 and 3 fixed):
parity_1e1_loop: DO parity_1 = 1, 2
  IF(mthbc(parity_1)/=2 .AND. mthbc(parity_1)/=999 .AND. ipass==1) &
    CYCLE parity_1e1_loop
  CALL GetFaceIndex(mx,mb,1,parity_1,fi(1,1:5))
  parity_3e1_loop: DO parity_3 = 1, 2
    IF(mthbc(parity_3+4)/=2 .AND. mthbc(parity_3+4)/=99 .AND. ipass==1) &
      CYCLE parity_3e1_loop
    CALL GetFaceIndex(mx,mb,3,parity_3,fi(3,1:5))
!
         qf(fi(1,1),      3:mx(2)-1:2,fi(3,1),1:nrvars) = &
      a1*qf(fi(1,2),      1:mx(2)-3:2,fi(3,2),1:nrvars) &
    + a2*qf(fi(1,3),      2:mx(2)-2:2,fi(3,3),1:nrvars) &
    + a3*qc(fi(1,4),mb(2,1)+1:mb(2,2),fi(3,4),1:nrvars)
!
             qf(fi(1,1),        1,fi(3,1),1:nrvars) = &
      a1*    qf(fi(1,2),        1,fi(3,2),1:nrvars) &
    + a2*    qf(fi(1,3),        1,fi(3,3),1:nrvars) &
    + a3*(b1*qc(fi(1,4),mb(2,1)-1,fi(3,4),1:nrvars) &
    +     b2*qc(fi(1,4),mb(2,1)  ,fi(3,4),1:nrvars) &
    +     b3*qc(fi(1,4),mb(2,1)+1,fi(3,4),1:nrvars))
!
         qf(fi(1,1),      2:mx(2)-2:2,fi(3,1),1:nrvars) = &
      a1*qf(fi(1,2),      4:mx(2)  :2,fi(3,2),1:nrvars) &
    + a2*qf(fi(1,3),      3:mx(2)-1:2,fi(3,3),1:nrvars) &
    + a3*qc(fi(1,4),mb(2,1):mb(2,2)-1,fi(3,4),1:nrvars)
!
             qf(fi(1,1),    mx(2),fi(3,1),1:nrvars) = &
      a1*    qf(fi(1,2),    mx(2),fi(3,2),1:nrvars) &
    + a2*    qf(fi(1,3),    mx(2),fi(3,3),1:nrvars) &
    + a3*(b1*qc(fi(1,4),mb(2,2)+1,fi(3,4),1:nrvars) &
    +     b2*qc(fi(1,4),mb(2,2)  ,fi(3,4),1:nrvars) &
    +     b3*qc(fi(1,4),mb(2,2)-1,fi(3,4),1:nrvars))
!
  END DO parity_3e1_loop
END DO parity_1e1_loop
!
! Edges 5--8 (coordinates 1 and 2 fixed):
parity_1e2_loop: DO parity_1 = 1, 2
  IF(mthbc(parity_1)/=2 .AND. mthbc(parity_1)/= 999 .AND. ipass==1) &
    CYCLE parity_1e2_loop
  CALL GetFaceIndex(mx,mb,1,parity_1,fi(1,1:5))
  parity_2e2_loop: DO parity_2 = 1, 2
    IF(mthbc(parity_2+2)/=2 .AND. mthbc(parity_2+2)/=999 .AND. ipass==1) &
      CYCLE parity_2e2_loop
    CALL GetFaceIndex(mx,mb,2,parity_2,fi(2,1:5))
!
         qf(fi(1,1),fi(2,1),      3:mx(3)-1:2,1:nrvars) = &
      a1*qf(fi(1,2),fi(2,2),      1:mx(3)-3:2,1:nrvars) &
    + a2*qf(fi(1,3),fi(2,3),      2:mx(3)-2:2,1:nrvars) &
    + a3*qc(fi(1,4),fi(2,4),mb(3,1)+1:mb(3,2),1:nrvars)
!
             qf(fi(1,1),fi(2,1),        1,1:nrvars) = &
      a1*    qf(fi(1,2),fi(2,2),        1,1:nrvars) &
    + a2*    qf(fi(1,3),fi(2,3),        1,1:nrvars) &
    + a3*(b1*qc(fi(1,4),fi(2,4),mb(3,1)-1,1:nrvars) &
    +     b2*qc(fi(1,4),fi(2,4),mb(3,1)  ,1:nrvars) &
    +     b3*qc(fi(1,4),fi(2,4),mb(3,1)+1,1:nrvars))
!
         qf(fi(1,1),fi(2,1),      2:mx(3)-2:2,1:nrvars) = &
      a1*qf(fi(1,2),fi(2,2),      4:mx(3)  :2,1:nrvars) &
    + a2*qf(fi(1,3),fi(2,3),      3:mx(3)-1:2,1:nrvars) &
    + a3*qc(fi(1,4),fi(2,4),mb(3,1):mb(3,2)-1,1:nrvars)
!
             qf(fi(1,1),fi(2,1),    mx(3),1:nrvars) = &
      a1*    qf(fi(1,2),fi(2,2),    mx(3),1:nrvars) &
    + a2*    qf(fi(1,3),fi(2,3),    mx(3),1:nrvars) &
    + a3*(b1*qc(fi(1,4),fi(2,4),mb(3,2)+1,1:nrvars) &
    +     b2*qc(fi(1,4),fi(2,4),mb(3,2)  ,1:nrvars) &
    +     b3*qc(fi(1,4),fi(2,4),mb(3,2)-1,1:nrvars))
!
  END DO parity_2e2_loop
END DO parity_1e2_loop
!
! Edges 9--12 (coordinates 2 and 3 fixed):
parity_2e3_loop: DO parity_2 = 1, 2
  IF(mthbc(parity_2+2)/=2 .AND. mthbc(parity_2+2)/=999 .AND. ipass==1) &
    CYCLE parity_2e3_loop
  CALL GetFaceIndex(mx,mb,2,parity_2,fi(2,1:5))
  parity_3e3_loop: DO parity_3 = 1, 2
    IF(mthbc(parity_3+4)/=2 .AND. mthbc(parity_3+4)/=999 .AND. ipass==1) &
      CYCLE parity_3e3_loop
    CALL GetFaceIndex(mx,mb,3,parity_3,fi(3,1:5))
!
         qf(      3:mx(1)-1:2,fi(2,1),fi(3,1),1:nrvars) = &
      a1*qf(      1:mx(1)-3:2,fi(2,2),fi(3,2),1:nrvars) &
    + a2*qf(      2:mx(1)-2:2,fi(2,3),fi(3,3),1:nrvars) &
    + a3*qc(mb(1,1)+1:mb(1,2),fi(2,4),fi(3,4),1:nrvars)
!
             qf(        1,fi(2,1),fi(3,1),1:nrvars) = &
      a1*    qf(        1,fi(2,2),fi(3,2),1:nrvars) &
    + a2*    qf(        1,fi(2,3),fi(3,3),1:nrvars) &
    + a3*(b1*qc(mb(1,1)-1,fi(2,4),fi(3,4),1:nrvars) &
    +     b2*qc(mb(1,1)  ,fi(2,4),fi(3,4),1:nrvars) &
    +     b3*qc(mb(1,1)+1,fi(2,4),fi(3,4),1:nrvars))
!
         qf(      2:mx(1)-2:2,fi(2,1),fi(3,1),1:nrvars) = &
      a1*qf(      4:mx(1)  :2,fi(2,2),fi(3,2),1:nrvars) &
    + a2*qf(      3:mx(1)-1:2,fi(2,3),fi(3,3),1:nrvars) &
    + a3*qc(mb(1,1):mb(1,2)-1,fi(2,4),fi(3,4),1:nrvars)
!
             qf(    mx(1),fi(2,1),fi(3,1),1:nrvars) = &
      a1*    qf(    mx(1),fi(2,2),fi(3,2),1:nrvars) &
    + a2*    qf(    mx(1),fi(2,3),fi(3,3),1:nrvars) &
    + a3*(b1*qc(mb(1,2)+1,fi(2,4),fi(3,4),1:nrvars) &
    +     b2*qc(mb(1,2)  ,fi(2,4),fi(3,4),1:nrvars) &
    +     b3*qc(mb(1,2)-1,fi(2,4),fi(3,4),1:nrvars))
!
  END DO parity_3e3_loop
END DO parity_2e3_loop
!
! (3) Interpolate the 6 faces:
!
! Only interpolate periodic or internal interfaces:
DO fc = 1, 6
  SELECT CASE(mthbc(fc))
  CASE(2,999)
    cycle_face(fc) = .FALSE.
  CASE DEFAULT
    cycle_face(fc) = .TRUE.
  END SELECT
END DO
!
! ipass = 0: Red and Black
! ipass = 1: Red
! ipass = 2: Black
! parity = 1: Yellow, Blue = Red
!             Green, Red   = Black
! parity = 2: Green, Red   = Red
!             Yellow, Blue = Black
!
! Faces 1 and 2 (coordinate 1 fixed):
parity_1_loop: DO parity_1 = 1, 2
  IF(cycle_face(parity_1)) CYCLE parity_1_loop
  CALL GetFaceIndex(mx,mb,1,parity_1,fi(1,1:5))
!
  IF(ipass==0 .OR. MODULO(parity_1+ipass,2)==0) THEN
!
! Yellow:
         qf(fi(1,1),      2:mx(2)-2:2,      2:mx(3)-2:2,1:nrvars) = &
      a1*qf(fi(1,2),      4:mx(2)  :2,      4:mx(3)  :2,1:nrvars) &
    + a2*qf(fi(1,3),      3:mx(2)-1:2,      3:mx(3)-1:2,1:nrvars) &
    + a3*qc(fi(1,4),mb(2,1):mb(2,2)-1,mb(3,1):mb(3,2)-1,1:nrvars)
!
             qf(fi(1,1),    mx(2),      2:mx(3)-2:2,1:nrvars) = &
      a1*    qf(fi(1,2),    mx(2),      4:mx(3)  :2,1:nrvars) &
    + a2*    qf(fi(1,3),    mx(2),      3:mx(3)-1:2,1:nrvars) &
    + a3*(b1*qc(fi(1,4),mb(2,2)+1,mb(3,1):mb(3,2)-1,1:nrvars) &
    +     b2*qc(fi(1,4),mb(2,2)  ,mb(3,1):mb(3,2)-1,1:nrvars) &
    +     b3*qc(fi(1,4),mb(2,2)-1,mb(3,1):mb(3,2)-1,1:nrvars))
!
             qf(fi(1,1),      2:mx(2)-2:2,    mx(3),1:nrvars) = &
      a1*    qf(fi(1,2),      4:mx(2)  :2,    mx(3),1:nrvars) &
    + a2*    qf(fi(1,3),      3:mx(2)-1:2,    mx(3),1:nrvars) &
    + a3*(b1*qc(fi(1,4),mb(2,1):mb(2,2)-1,mb(3,2)+1,1:nrvars) &
    +     b2*qc(fi(1,4),mb(2,1):mb(2,2)-1,mb(3,2)  ,1:nrvars) &
    +     b3*qc(fi(1,4),mb(2,1):mb(2,2)-1,mb(3,2)-1,1:nrvars))
!
! Blue:
         qf(fi(1,1),      3:mx(2)-1:2,      3:mx(3)-1:2,1:nrvars) = &
      a1*qf(fi(1,2),      1:mx(2)-3:2,      1:mx(3)-3:2,1:nrvars) &
    + a2*qf(fi(1,3),      2:mx(2)-2:2,      2:mx(3)-2:2,1:nrvars) &
    + a3*qc(fi(1,4),mb(2,1)+1:mb(2,2),mb(3,1)+1:mb(3,2),1:nrvars)
!
             qf(fi(1,1),      3:mx(2)-1:2,        1,1:nrvars) = &
      a1*    qf(fi(1,2),      1:mx(2)-3:2,        1,1:nrvars) &
    + a2*    qf(fi(1,3),      2:mx(2)-2:2,        1,1:nrvars) &
    + a3*(b1*qc(fi(1,4),mb(2,1)+1:mb(2,2),mb(3,1)-1,1:nrvars) &
    +     b2*qc(fi(1,4),mb(2,1)+1:mb(2,2),mb(3,1)  ,1:nrvars) &
    +     b3*qc(fi(1,4),mb(2,1)+1:mb(2,2),mb(3,1)+1,1:nrvars))
!
             qf(fi(1,1),        1,      3:mx(3)-1:2,1:nrvars) = &
      a1*    qf(fi(1,2),        1,      1:mx(3)-3:2,1:nrvars) &
    + a2*    qf(fi(1,3),        1,      2:mx(3)-2:2,1:nrvars) &
    + a3*(b1*qc(fi(1,4),mb(2,1)-1,mb(3,1)+1:mb(3,2),1:nrvars) &
    +     b2*qc(fi(1,4),mb(2,1)  ,mb(3,1)+1:mb(3,2),1:nrvars) &
    +     b3*qc(fi(1,4),mb(2,1)+1,mb(3,1)+1:mb(3,2),1:nrvars))
  END IF
!
  IF(ipass==0 .OR. MODULO(parity_1+ipass,2)==1) THEN
!
! Green:
         qf(fi(1,1),      3:mx(2)-1:2,      2:mx(3)-2:2,1:nrvars) = &
      a1*qf(fi(1,2),      1:mx(2)-3:2,      4:mx(3)  :2,1:nrvars) &
    + a2*qf(fi(1,3),      2:mx(2)-2:2,      3:mx(3)-1:2,1:nrvars) &
    + a3*qc(fi(1,4),mb(2,1)+1:mb(2,2),mb(3,1):mb(3,2)-1,1:nrvars)
!
             qf(fi(1,1),        1,      2:mx(3)-2:2,1:nrvars) = &
      a1*    qf(fi(1,2),        1,      4:mx(3)  :2,1:nrvars) &
    + a2*    qf(fi(1,3),        1,      3:mx(3)-1:2,1:nrvars) &
    + a3*(b1*qc(fi(1,4),mb(2,1)-1,mb(3,1):mb(3,2)-1,1:nrvars) &
    +     b2*qc(fi(1,4),mb(2,1)  ,mb(3,1):mb(3,2)-1,1:nrvars) &
    +     b3*qc(fi(1,4),mb(2,1)+1,mb(3,1):mb(3,2)-1,1:nrvars))
!
             qf(fi(1,1),      3:mx(2)-1:2,    mx(3),1:nrvars) = &
      a1*    qf(fi(1,2),      1:mx(2)-3:2,    mx(3),1:nrvars) &
    + a2*    qf(fi(1,3),      2:mx(2)-2:2,    mx(3),1:nrvars) &
    + a3*(b1*qc(fi(1,4),mb(2,1)+1:mb(2,2),mb(3,2)+1,1:nrvars) &
    +     b2*qc(fi(1,4),mb(2,1)+1:mb(2,2),mb(3,2)  ,1:nrvars) &
    +     b3*qc(fi(1,4),mb(2,1)+1:mb(2,2),mb(3,2)-1,1:nrvars))
!
! Red:
         qf(fi(1,1),      2:mx(2)-2:2,      3:mx(3)-1:2,1:nrvars) = &
      a1*qf(fi(1,2),      4:mx(2)  :2,      1:mx(3)-3:2,1:nrvars) &
    + a2*qf(fi(1,3),      3:mx(2)-1:2,      2:mx(3)-2:2,1:nrvars) &
    + a3*qc(fi(1,4),mb(2,1):mb(2,2)-1,mb(3,1)+1:mb(3,2),1:nrvars)
!
             qf(fi(1,1),      2:mx(2)-2:2,        1,1:nrvars) = &
      a1*    qf(fi(1,2),      4:mx(2)  :2,        1,1:nrvars) &
    + a2*    qf(fi(1,3),      3:mx(2)-1:2,        1,1:nrvars) &
    + a3*(b1*qc(fi(1,4),mb(2,1):mb(2,2)-1,mb(3,1)-1,1:nrvars) &
    +     b2*qc(fi(1,4),mb(2,1):mb(2,2)-1,mb(3,1)  ,1:nrvars) &
    +     b3*qc(fi(1,4),mb(2,1):mb(2,2)-1,mb(3,1)+1,1:nrvars))
!
             qf(fi(1,1),    mx(2),      3:mx(3)-1:2,1:nrvars) = &
      a1*    qf(fi(1,2),    mx(2),      1:mx(3)-3:2,1:nrvars) &
    + a2*    qf(fi(1,3),    mx(2),      2:mx(3)-2:2,1:nrvars) &
    + a3*(b1*qc(fi(1,4),mb(2,2)+1,mb(3,1)+1:mb(3,2),1:nrvars) &
    +     b2*qc(fi(1,4),mb(2,2)  ,mb(3,1)+1:mb(3,2),1:nrvars) &
    +     b3*qc(fi(1,4),mb(2,2)-1,mb(3,1)+1:mb(3,2),1:nrvars))
  END IF
!
! Corner grey cells:
  DO parity_2 = 1, 2
    CALL GetFaceIndex(mx,mb,2,parity_2,fi(2,1:5))
    DO parity_3 = 1, 2
      CALL GetFaceIndex(mx,mb,3,parity_3,fi(3,1:5))
               qf(fi(1,1),fi(2,3),fi(3,3),1:nrvars) = &
        a1*    qf(fi(1,2),fi(2,3),fi(3,3),1:nrvars) &
      + a2*    qf(fi(1,3),fi(2,3),fi(3,3),1:nrvars) &
      + a3*(c1*qc(fi(1,4),fi(2,5),fi(3,5),1:nrvars) &
      +     c2*qc(fi(1,4),fi(2,4),fi(3,5),1:nrvars) &
      +     c3*qc(fi(1,4),fi(2,5),fi(3,4),1:nrvars) &
      +     c4*qc(fi(1,4),fi(2,4),fi(3,4),1:nrvars))
    END DO
  END DO
END DO parity_1_loop
!
! Faces 3 and 4 (coordinate 2 fixed):
parity_2_loop: DO parity_2 = 1, 2
  IF(cycle_face(2+parity_2)) CYCLE parity_2_loop
  CALL GetFaceIndex(mx,mb,2,parity_2,fi(2,1:5))
!
  IF(ipass==0 .OR. MODULO(parity_2+ipass,2)==0) THEN
!
! Yellow:
         qf(      2:mx(1)-2:2,fi(2,1),      2:mx(3)-2:2,1:nrvars) = &
      a1*qf(      4:mx(1)  :2,fi(2,2),      4:mx(3)  :2,1:nrvars) &
    + a2*qf(      3:mx(1)-1:2,fi(2,3),      3:mx(3)-1:2,1:nrvars) &
    + a3*qc(mb(1,1):mb(1,2)-1,fi(2,4),mb(3,1):mb(3,2)-1,1:nrvars)
!
             qf(      2:mx(1)-2:2,fi(2,1),    mx(3),1:nrvars) = &
      a1*    qf(      4:mx(1)  :2,fi(2,2),    mx(3),1:nrvars) &
    + a2*    qf(      3:mx(1)-1:2,fi(2,3),    mx(3),1:nrvars) &
    + a3*(b1*qc(mb(1,1):mb(1,2)-1,fi(2,4),mb(3,2)+1,1:nrvars) &
    +     b2*qc(mb(1,1):mb(1,2)-1,fi(2,4),mb(3,2)  ,1:nrvars) &
    +     b3*qc(mb(1,1):mb(1,2)-1,fi(2,4),mb(3,2)-1,1:nrvars))
!
             qf(    mx(1),fi(2,1),      2:mx(3)-2:2,1:nrvars) = &
      a1*    qf(    mx(1),fi(2,2),      4:mx(3)  :2,1:nrvars) &
    + a2*    qf(    mx(1),fi(2,3),      3:mx(3)-1:2,1:nrvars) &
    + a3*(b1*qc(mb(1,2)+1,fi(2,4),mb(3,1):mb(3,2)-1,1:nrvars) &
    +     b2*qc(mb(1,2)  ,fi(2,4),mb(3,1):mb(3,2)-1,1:nrvars) &
    +     b3*qc(mb(1,2)-1,fi(2,4),mb(3,1):mb(3,2)-1,1:nrvars))
!
! Blue:
         qf(      3:mx(1)-1:2,fi(2,1),      3:mx(3)-1:2,1:nrvars) = &
      a1*qf(      1:mx(1)-3:2,fi(2,2),      1:mx(3)-3:2,1:nrvars) &
    + a2*qf(      2:mx(1)-2:2,fi(2,3),      2:mx(3)-2:2,1:nrvars) &
    + a3*qc(mb(1,1)+1:mb(1,2),fi(2,4),mb(3,1)+1:mb(3,2),1:nrvars)
!
             qf(        1,fi(2,1),      3:mx(3)-1:2,1:nrvars) = &
      a1*    qf(        1,fi(2,2),      1:mx(3)-3:2,1:nrvars) &
    + a2*    qf(        1,fi(2,3),      2:mx(3)-2:2,1:nrvars) &
    + a3*(b1*qc(mb(1,1)-1,fi(2,4),mb(3,1)+1:mb(3,2),1:nrvars) &
    +     b2*qc(mb(1,1)  ,fi(2,4),mb(3,1)+1:mb(3,2),1:nrvars) &
    +     b3*qc(mb(1,1)+1,fi(2,4),mb(3,1)+1:mb(3,2),1:nrvars))
!
             qf(      3:mx(1)-1:2,fi(2,1),        1,1:nrvars) = &
      a1*    qf(      1:mx(1)-3:2,fi(2,2),        1,1:nrvars) &
    + a2*    qf(      2:mx(1)-2:2,fi(2,3),        1,1:nrvars) &
    + a3*(b1*qc(mb(1,1)+1:mb(1,2),fi(2,4),mb(3,1)-1,1:nrvars) &
    +     b2*qc(mb(1,1)+1:mb(1,2),fi(2,4),mb(3,1)  ,1:nrvars) &
    +     b3*qc(mb(1,1)+1:mb(1,2),fi(2,4),mb(3,1)+1,1:nrvars))
  END IF
!
  IF(ipass==0 .OR. MODULO(parity_2+ipass,2)==1) THEN
!
! Green:
         qf(      3:mx(1)-1:2,fi(2,1),      2:mx(3)-2:2,1:nrvars) = &
      a1*qf(      1:mx(1)-3:2,fi(2,2),      4:mx(3)  :2,1:nrvars) &
    + a2*qf(      2:mx(1)-2:2,fi(2,3),      3:mx(3)-1:2,1:nrvars) &
    + a3*qc(mb(1,1)+1:mb(1,2),fi(2,4),mb(3,1):mb(3,2)-1,1:nrvars)
!
             qf(      2:mx(1)-2:2,fi(2,1),        1,1:nrvars) = &
      a1*    qf(      4:mx(1)  :2,fi(2,2),        1,1:nrvars) &
    + a2*    qf(      3:mx(1)-1:2,fi(2,3),        1,1:nrvars) &
    + a3*(b1*qc(mb(1,1):mb(1,2)-1,fi(2,4),mb(3,1)-1,1:nrvars) &
    +     b2*qc(mb(1,1):mb(1,2)-1,fi(2,4),mb(3,1)  ,1:nrvars) &
    +     b3*qc(mb(1,1):mb(1,2)-1,fi(2,4),mb(3,1)+1,1:nrvars))
!
             qf(    mx(1),fi(2,1),      3:mx(3)-1:2,1:nrvars) = &
      a1*    qf(    mx(1),fi(2,2),      1:mx(3)-3:2,1:nrvars) &
    + a2*    qf(    mx(1),fi(2,3),      2:mx(3)-2:2,1:nrvars) &
    + a3*(b1*qc(mb(1,2)+1,fi(2,4),mb(3,1)+1:mb(3,2),1:nrvars) &
    +     b2*qc(mb(1,2)  ,fi(2,4),mb(3,1)+1:mb(3,2),1:nrvars) &
    +     b3*qc(mb(1,2)-1,fi(2,4),mb(3,1)+1:mb(3,2),1:nrvars))
!
! Red:
         qf(      2:mx(1)-2:2,fi(2,1),      3:mx(3)-1:2,1:nrvars) = &
      a1*qf(      4:mx(1)  :2,fi(2,2),      1:mx(3)-3:2,1:nrvars) &
    + a2*qf(      3:mx(1)-1:2,fi(2,3),      2:mx(3)-2:2,1:nrvars) &
    + a3*qc(mb(1,1):mb(1,2)-1,fi(2,4),mb(3,1)+1:mb(3,2),1:nrvars)
!
             qf(        1,fi(2,1),      2:mx(3)-2:2,1:nrvars) = &
      a1*    qf(        1,fi(2,2),      4:mx(3)  :2,1:nrvars) &
    + a2*    qf(        1,fi(2,3),      3:mx(3)-1:2,1:nrvars) &
    + a3*(b1*qc(mb(1,1)-1,fi(2,4),mb(3,1):mb(3,2)-1,1:nrvars) &
    +     b2*qc(mb(1,1)  ,fi(2,4),mb(3,1):mb(3,2)-1,1:nrvars) &
    +     b3*qc(mb(1,1)+1,fi(2,4),mb(3,1):mb(3,2)-1,1:nrvars))
!
             qf(      3:mx(1)-1:2,fi(2,1),    mx(3),1:nrvars) = &
      a1*    qf(      1:mx(1)-3:2,fi(2,2),    mx(3),1:nrvars) &
    + a2*    qf(      2:mx(1)-2:2,fi(2,3),    mx(3),1:nrvars) &
    + a3*(b1*qc(mb(1,1)+1:mb(1,2),fi(2,4),mb(3,2)+1,1:nrvars) &
    +     b2*qc(mb(1,1)+1:mb(1,2),fi(2,4),mb(3,2)  ,1:nrvars) &
    +     b3*qc(mb(1,1)+1:mb(1,2),fi(2,4),mb(3,2)-1,1:nrvars))
  END IF
!
! Corner grey cells:
  DO parity_1 = 1, 2
    CALL GetFaceIndex(mx,mb,1,parity_1,fi(1,1:5))
    DO parity_3 = 1, 2
      CALL GetFaceIndex(mx,mb,3,parity_3,fi(3,1:5))
               qf(fi(1,3),fi(2,1),fi(3,3),1:nrvars) = &
        a1*    qf(fi(1,3),fi(2,2),fi(3,3),1:nrvars) &
      + a2*    qf(fi(1,3),fi(2,3),fi(3,3),1:nrvars) &
      + a3*(c1*qc(fi(1,5),fi(2,4),fi(3,5),1:nrvars) &
      +     c2*qc(fi(1,4),fi(2,4),fi(3,5),1:nrvars) &
      +     c3*qc(fi(1,5),fi(2,4),fi(3,4),1:nrvars) &
      +     c4*qc(fi(1,4),fi(2,4),fi(3,4),1:nrvars))
    END DO
  END DO
END DO parity_2_loop
!
! Faces 5 and 6 (coordinate 3 fixed):
parity_3_loop: DO parity_3 = 1, 2
  IF(cycle_face(4+parity_3)) CYCLE parity_3_loop
  CALL GetFaceIndex(mx,mb,3,parity_3,fi(3,1:5))
!
  IF(ipass==0 .OR. MODULO(parity_3+ipass,2)==0) THEN
!
! Yellow:
         qf(      2:mx(1)-2:2,      2:mx(2)-2:2,fi(3,1),1:nrvars) = &
      a1*qf(      4:mx(1)  :2,      4:mx(2)  :2,fi(3,2),1:nrvars) &
    + a2*qf(      3:mx(1)-1:2,      3:mx(2)-1:2,fi(3,3),1:nrvars) &
    + a3*qc(mb(1,1):mb(1,2)-1,mb(2,1):mb(2,2)-1,fi(3,4),1:nrvars)
!
             qf(    mx(1),      2:mx(2)-2:2,fi(3,1),1:nrvars) = &
      a1*    qf(    mx(1),      4:mx(2)  :2,fi(3,2),1:nrvars) &
    + a2*    qf(    mx(1),      3:mx(2)-1:2,fi(3,3),1:nrvars) &
    + a3*(b1*qc(mb(1,2)+1,mb(2,1):mb(2,2)-1,fi(3,4),1:nrvars) &
    +     b2*qc(mb(1,2)  ,mb(2,1):mb(2,2)-1,fi(3,4),1:nrvars) &
    +     b3*qc(mb(1,2)-1,mb(2,1):mb(2,2)-1,fi(3,4),1:nrvars))
!
             qf(      2:mx(1)-2:2,    mx(2),fi(3,1),1:nrvars) = &
      a1*    qf(      4:mx(1)  :2,    mx(2),fi(3,2),1:nrvars) &
    + a2*    qf(      3:mx(1)-1:2,    mx(2),fi(3,3),1:nrvars) &
    + a3*(b1*qc(mb(1,1):mb(1,2)-1,mb(2,2)+1,fi(3,4),1:nrvars) &
    +     b2*qc(mb(1,1):mb(1,2)-1,mb(2,2)  ,fi(3,4),1:nrvars) &
    +     b3*qc(mb(1,1):mb(1,2)-1,mb(2,2)-1,fi(3,4),1:nrvars))
!
! Blue:
         qf(      3:mx(1)-1:2,      3:mx(2)-1:2,fi(3,1),1:nrvars) = &
      a1*qf(      1:mx(1)-3:2,      1:mx(2)-3:2,fi(3,2),1:nrvars) &
    + a2*qf(      2:mx(1)-2:2,      2:mx(2)-2:2,fi(3,3),1:nrvars) &
    + a3*qc(mb(1,1)+1:mb(1,2),mb(2,1)+1:mb(2,2),fi(3,4),1:nrvars)
!
             qf(      3:mx(1)-1:2,        1,fi(3,1),1:nrvars) = &
      a1*    qf(      1:mx(1)-3:2,        1,fi(3,2),1:nrvars) &
    + a2*    qf(      2:mx(1)-2:2,        1,fi(3,3),1:nrvars) &
    + a3*(b1*qc(mb(1,1)+1:mb(1,2),mb(2,1)-1,fi(3,4),1:nrvars) &
    +     b2*qc(mb(1,1)+1:mb(1,2),mb(2,1)  ,fi(3,4),1:nrvars) &
    +     b3*qc(mb(1,1)+1:mb(1,2),mb(2,1)+1,fi(3,4),1:nrvars))
!
             qf(        1,      3:mx(2)-1:2,fi(3,1),1:nrvars) = &
      a1*    qf(        1,      1:mx(2)-3:2,fi(3,2),1:nrvars) &
    + a2*    qf(        1,      2:mx(2)-2:2,fi(3,3),1:nrvars) &
    + a3*(b1*qc(mb(1,1)-1,mb(2,1)+1:mb(2,2),fi(3,4),1:nrvars) &
    +     b2*qc(mb(1,1)  ,mb(2,1)+1:mb(2,2),fi(3,4),1:nrvars) &
    +     b3*qc(mb(1,1)+1,mb(2,1)+1:mb(2,2),fi(3,4),1:nrvars))
  END IF
!
  IF(ipass==0 .OR. MODULO(parity_3+ipass,2)==1) THEN
!
! Green:
         qf(      3:mx(1)-1:2,      2:mx(2)-2:2,fi(3,1),1:nrvars) = &
      a1*qf(      1:mx(1)-3:2,      4:mx(2)  :2,fi(3,2),1:nrvars) &
    + a2*qf(      2:mx(1)-2:2,      3:mx(2)-1:2,fi(3,3),1:nrvars) &
    + a3*qc(mb(1,1)+1:mb(1,2),mb(2,1):mb(2,2)-1,fi(3,4),1:nrvars)
!
             qf(        1,      2:mx(2)-2:2,fi(3,1),1:nrvars) = &
      a1*    qf(        1,      4:mx(2)  :2,fi(3,2),1:nrvars) &
    + a2*    qf(        1,      3:mx(2)-1:2,fi(3,3),1:nrvars) &
    + a3*(b1*qc(mb(1,1)-1,mb(2,1):mb(2,2)-1,fi(3,4),1:nrvars) &
    +     b2*qc(mb(1,1)  ,mb(2,1):mb(2,2)-1,fi(3,4),1:nrvars) &
    +     b3*qc(mb(1,1)+1,mb(2,1):mb(2,2)-1,fi(3,4),1:nrvars))
!
             qf(      3:mx(1)-1:2,    mx(2),fi(3,1),1:nrvars) = &
      a1*    qf(      1:mx(1)-3:2,    mx(2),fi(3,2),1:nrvars) &
    + a2*    qf(      2:mx(1)-2:2,    mx(2),fi(3,3),1:nrvars) &
    + a3*(b1*qc(mb(1,1)+1:mb(1,2),mb(2,2)+1,fi(3,4),1:nrvars) &
    +     b2*qc(mb(1,1)+1:mb(1,2),mb(2,2)  ,fi(3,4),1:nrvars) &
    +     b3*qc(mb(1,1)+1:mb(1,2),mb(2,2)-1,fi(3,4),1:nrvars))
!
! Red:
         qf(      2:mx(1)-2:2,      3:mx(2)-1:2,fi(3,1),1:nrvars) = &
      a1*qf(      4:mx(1)  :2,      1:mx(2)-3:2,fi(3,2),1:nrvars) &
    + a2*qf(      3:mx(1)-1:2,      2:mx(2)-2:2,fi(3,3),1:nrvars) &
    + a3*qc(mb(1,1):mb(1,2)-1,mb(2,1)+1:mb(2,2),fi(3,4),1:nrvars)
!
             qf(      2:mx(1)-2:2,        1,fi(3,1),1:nrvars) = &
      a1*    qf(      4:mx(1)  :2,        1,fi(3,2),1:nrvars) &
    + a2*    qf(      3:mx(1)-1:2,        1,fi(3,3),1:nrvars) &
    + a3*(b1*qc(mb(1,1):mb(1,2)-1,mb(2,1)-1,fi(3,4),1:nrvars) &
    +     b2*qc(mb(1,1):mb(1,2)-1,mb(2,1)  ,fi(3,4),1:nrvars) &
    +     b3*qc(mb(1,1):mb(1,2)-1,mb(2,1)+1,fi(3,4),1:nrvars))
!
             qf(    mx(1),      3:mx(2)-1:2,fi(3,1),1:nrvars) = &
      a1*    qf(    mx(1),      1:mx(2)-3:2,fi(3,2),1:nrvars) &
    + a2*    qf(    mx(1),      2:mx(2)-2:2,fi(3,3),1:nrvars) &
    + a3*(b1*qc(mb(1,2)+1,mb(2,1)+1:mb(2,2),fi(3,4),1:nrvars) &
    +     b2*qc(mb(1,2)  ,mb(2,1)+1:mb(2,2),fi(3,4),1:nrvars) &
    +     b3*qc(mb(1,2)-1,mb(2,1)+1:mb(2,2),fi(3,4),1:nrvars))
  END IF
!
! Corner grey cells:
  DO parity_1 = 1, 2
    CALL GetFaceIndex(mx,mb,1,parity_1,fi(1,1:5))
    DO parity_2 = 1, 2
      CALL GetFaceIndex(mx,mb,2,parity_2,fi(2,1:5))
               qf(fi(1,3),fi(2,3),fi(3,1),1:nrvars) = &
        a1*    qf(fi(1,3),fi(2,3),fi(3,2),1:nrvars) &
      + a2*    qf(fi(1,3),fi(2,3),fi(3,3),1:nrvars) &
      + a3*(c1*qc(fi(1,5),fi(2,5),fi(3,4),1:nrvars) &
      +     c2*qc(fi(1,4),fi(2,5),fi(3,4),1:nrvars) &
      +     c3*qc(fi(1,5),fi(2,4),fi(3,4),1:nrvars) &
      +     c4*qc(fi(1,4),fi(2,4),fi(3,4),1:nrvars))
    END DO
  END DO
END DO parity_3_loop
!
END SUBROUTINE Interpolate3D
END MODULE Boundary


