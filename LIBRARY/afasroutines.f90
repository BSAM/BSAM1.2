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
! File:             afasroutines.f90
! Purpose:          Adaptive FAS Multigrid module
! Contains:
! Revision History: Ver. 1.0 Oct. 2006 Steven Wise
! Revision History: Ver. 1.1 May. 2007 Steven Wise
! Revision History: Ver. 1.2 Jul. 2015 Steven Wise
! -----------------------------------------------------------------------
MODULE AFASRoutines
!
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE MultigridIterations
USE NodeInfoDef
USE TreeOps, ONLY: ApplyOnLevel
IMPLICIT NONE
!
TYPE(funcparam):: dummy
INTEGER:: itmg, level
REAL(KIND=r8):: residual
!
CALL FillDown(solutionfield)
!
! If the auxiliary fields are updated in the vcycle, perform here:
DO level = finestlevel, minlevel, -1
  CALL ApplyOnLevel(level,UpdateAuxInVcycle,dummy)
END DO
!
! Maybe implement this later, but it needs more work:
!CALL FillDown(auxiliaryfield)
!
DO level = finestlevel, rootlevel, -1
  CALL ApplyOnLevel(level,GetSourceFunction,dummy)
END DO
!
! If the source function is updated after each vcycle, perform here first:
DO level = finestlevel, rootlevel, -1
  CALL ApplyOnLevel(level,UpdateSourceFunction,dummy)
END DO
!
CALL FillDown(sourcefield)
!
! Perform Adaptive Full Aproximation Scheme Vcycles on the grid hierarchy:
vcycleloop: DO itmg = 1, maxvcycles
!
  CALL AFASVcycle(finestlevel)
!
  CALL FillDown(solutionfield)
!
! If the auxiliary fields are updated in the vcycle, perform here:
  IF(MODULO(itmg-1,updateauxfreq)==0) THEN
    DO level = finestlevel, minlevel, -1
      CALL ApplyOnLevel(level,UpdateAuxInVcycle,dummy)
    END DO
  END IF
!
! Maybe implement this later, but it needs more work:
!  CALL FillDown(auxiliaryfield)
!
! If the source function is updated after each vcycle, perform here:
  DO level = finestlevel, rootlevel, -1
    CALL ApplyOnLevel(level,UpdateSourceFunction,dummy)
  END DO
!
  CALL FillDown(sourcefield)
!
! Note that it is necessary to have an up-to-date source function for 
! calculating the residual:
  residual = ErrorAFAS(finestlevel)
  PRINT *, 'it =', itmg, ' residual =', residual
  IF(residual<qerrortol) EXIT vcycleloop
!
END DO vcycleloop
!
END SUBROUTINE MultigridIterations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION ErrorAFAS(level) RESULT(errorresult)
USE NodeInfoDef
USE TreeOps, ONLY: ApplyOnLevel, GetRootInfo
IMPLICIT NONE
!
TYPE(nodeinfo), POINTER:: rootinfo
!
INTEGER, INTENT(IN):: level
REAL(KIND=r8):: errorresult
!
TYPE(funcparam):: dummy
INTEGER:: i, ierror, nrvars
REAL(KIND=r8), DIMENSION(1:maxnrvars)::  componenterror
!
SELECT CASE(errortype)
  CASE(1)
!
! Blended error calculation:
    integralresult(1:2) = 0.0_r8
!
    CALL ApplyOnLevel(level,L2Error,dummy)
!
    errorresult = SQRT(integralresult(1)/integralresult(2))
  CASE(2)
!
! Component errors:
    ierror = GetRootInfo(rootinfo)
    nrvars = rootinfo%nrvars
    integralresult(2) = 0.0_r8
    componentintegral(1:nrvars) = 0.0_r8
!
    CALL ApplyOnLevel(level,L2ComponentErrors,dummy)
!
    componenterror(1:nrvars) = SQRT(componentintegral(1:nrvars) &
                             / integralresult(2))
!
    DO i = 1, nrvars
      PRINT *, i, componenterror(i)
    END DO
!
    errorresult = MAXVAL(componenterror(1:nrvars))
  CASE DEFAULT
    PRINT *, 'ErrorAFAS: Only errortype = 1,2 are supported.'
    STOP
!
END SELECT
!
END FUNCTION ErrorAFAS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION L2Error(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok
IMPLICIT NONE
!
! Adds the L2 error on this grid to the total.
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
INTEGER:: nrvars
INTEGER, DIMENSION(1:maxdims):: mx
!
L2Error = err_ok
!
nrvars = info%nrvars
mx = 1; mx(1:ndims) = info%mx(1:ndims)
!
CALL Residual(info)
!
integralresult(1) = SquareSum(info%rf(1:mx(1),1:mx(2),1:mx(3),1:nrvars)) &
                  + integralresult(1)
!
integralresult(2) = REAL(nrvars*PRODUCT(mx(1:ndims)),KIND=r8) &
                  + integralresult(2)
!
END FUNCTION L2Error
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION L2ComponentErrors(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok
IMPLICIT NONE
!
! Adds the L2 component error on this grid to the component total.
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
INTEGER:: i, nrvars
INTEGER, DIMENSION(1:maxdims):: mx
!
L2ComponentErrors = err_ok
!
nrvars = info%nrvars
mx = 1; mx(1:ndims) = info%mx(1:ndims)
!
CALL Residual(info)
!
DO i = 1, nrvars
  componentintegral(i) = SquareSum(info%rf(1:mx(1),1:mx(2),1:mx(3),i:i)) &
                       + componentintegral(i)
END DO
!
integralresult(2) = REAL(PRODUCT(mx(1:ndims)),KIND=r8)+integralresult(2)
!
END FUNCTION L2ComponentErrors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION SquareSum(rf) RESULT(squaresumresult)
USE NodeInfoDef
IMPLICIT NONE
!
! Dimensionally invariant square integral.
!
REAL(KIND=r8), DIMENSION(:,:,:,:), INTENT(IN):: rf
REAL(KIND=r8):: squaresumresult
!
squaresumresult = SUM(rf*rf)
!
END FUNCTION SquareSum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE FillDown(field)
USE NodeInfoDef
USE TreeOps, ONLY: ApplyOnLevel
USE Boundary, ONLY: SetGhost
IMPLICIT NONE
!
INTEGER, INTENT(IN):: field
!
TYPE(funcparam):: dummy
INTEGER:: level
!
dummy%iswitch = field
!
DO level = finestlevel, minlevel+1, -1
  CALL ApplyOnLevel(level,FillDownLevel,dummy)
  IF(field==solutionfield) CALL SetGhost(level-1,0)
END DO
!
END SUBROUTINE FillDown
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION FillDownLevel(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok, GetParentInfo
USE GridUtilities, ONLY: Restriction2D, Restriction3D
IMPLICIT NONE
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
TYPE(nodeinfo), POINTER:: parent
INTEGER:: field, ierror, maux, mmaux, nrvars
INTEGER, DIMENSION(1:maxdims):: cmx, mx
INTEGER, DIMENSION(1:maxdims,1:2):: mb
!
FillDownLevel = err_ok
!
IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
field = dummy%iswitch
!
nrvars = info%nrvars
mx = 1; mx(1:ndims) = info%mx(1:ndims)
cmx = 1; cmx(1:ndims) = mx(1:ndims)/2
mb = 1; mb(1:ndims,1:2) = info%mbounds(1:ndims,1:2)
!
!maux = info%maux; mmaux = MAX(maux,1)
!
ierror = GetParentInfo(parent)
!
SELECT CASE(field)
  CASE(sourcefield)
    SELECT CASE(ndims)
      CASE(2)
        parent%f(mb(1,1):mb(1,2),mb(2,1):mb(2,2),1,1:nrvars) &
          = Restriction2D(info%f(1:mx(1),1:mx(2),1,1:nrvars))
      CASE(3)
        parent%f(mb(1,1):mb(1,2),mb(2,1):mb(2,2),mb(3,1):mb(3,2),1:nrvars) &
          = Restriction3D(info%f(1:mx(1),1:mx(2),1:mx(3),1:nrvars))
      CASE DEFAULT
        PRINT *, 'FillDownLevel: Only ndims = 2,3 are supported.'
        STOP
    END SELECT
  CASE(solutionfield)
    SELECT CASE(ndims)
      CASE(2)
        parent%q(mb(1,1):mb(1,2),mb(2,1):mb(2,2),1,1:nrvars) &
          = Restriction2D(info%q(1:mx(1),1:mx(2),1,1:nrvars))
      CASE(3)
        parent%q(mb(1,1):mb(1,2),mb(2,1):mb(2,2),mb(3,1):mb(3,2),1:nrvars) &
          = Restriction3D(info%q(1:mx(1),1:mx(2),1:mx(3),1:nrvars))
      CASE DEFAULT
        PRINT *, 'FillDownLevel: Only ndims = 2,3 are supported.'
        STOP
    END SELECT
!  CASE(auxiliaryfield)
!    IF(maux>0) THEN
!      SELECT CASE(ndims)
!        CASE(2)
!          parent%aux(mb(1,1):mb(1,2),mb(2,1):mb(2,2),1,1:mmaux) &
!            = Restriction2D(info%aux(1:mx(1),1:mx(2),1,1:mmaux))
!        CASE(3)
!          parent%aux(mb(1,1):mb(1,2),mb(2,1):mb(2,2),mb(3,1):mb(3,2),1:mmaux) &
!            = Restriction3D(info%aux(1:mx(1),1:mx(2),1:mx(3),1:mmaux))
!        CASE DEFAULT
!          PRINT *, 'FillDownLevel: Only ndims = 2,3 are supported.'
!          STOP
!      END SELECT
!    END IF
END SELECT
!
END FUNCTION FillDownLevel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION GetSourceFunction(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok
USE Problem, ONLY: Source2D, Source3D
IMPLICIT NONE
!
! Dimensionally invariant source (or partial source) function routine.
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
INTEGER:: ll, maux, mbc, mmaux, nrvars
INTEGER, DIMENSION(1:maxdims):: amx, aul, mx, ul
REAL(KIND=r8):: h
REAL(KIND=r8), DIMENSION(1:maxdims):: xlower
!
GetSourceFunction = err_ok
!
IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
nrvars = info%nrvars
mx = 1; mx(1:ndims) = info%mx(1:ndims)
h = info%dx(1)
xlower(1:ndims) = info%xlower(1:ndims)
!
maux = info%maux; mmaux = MAX(maux,1)
amx = 1; IF(maux>0) amx(1:ndims) = mx(1:ndims)
!
mbc = info%mbc
ll = 1-mbc
 ul = 1;  ul(1:ndims) =  mx(1:ndims)+mbc
aul = 1; aul(1:ndims) = amx(1:ndims)+mbc
!
SELECT CASE(ndims)
  CASE(2)
    info%ftmp(1:mx(1),1:mx(2),1,1:nrvars) &
      = Source2D(info%q(   ll: ul(1),ll: ul(2),1,1:nrvars), &
                 info%qold(ll: ul(1),ll: ul(2),1,1:nrvars), &
                 info%aux( ll:aul(1),ll:aul(2),1,1:mmaux ), &
                 ll,mx(1:2),nrvars,maux,mbc,h,xlower(1:2))
  CASE(3)
    info%ftmp(1:mx(1),1:mx(2),1:mx(3),1:nrvars) &
      = Source3D(info%q(   ll: ul(1),ll: ul(2),ll: ul(3),1:nrvars), &
                 info%qold(ll: ul(1),ll: ul(2),ll: ul(3),1:nrvars), &
                 info%aux( ll:aul(1),ll:aul(2),ll:aul(3),1:mmaux ), &
                 ll,mx(1:3),nrvars,maux,mbc,h,xlower(1:3))
  CASE DEFAULT
    PRINT *, 'GetSourceFunction: Only ndims = 2,3 are supported.'
    STOP
END SELECT
!
END FUNCTION GetSourceFunction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION UpdateAuxInVcycle(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok
USE Problem, ONLY: UpdateAuxVcycle2D, UpdateAuxVcycle3D
IMPLICIT NONE
!
! Used if the auxiliary variables are updated after each Vcycle iteration.
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
INTEGER:: ll, maux, mbc, mmaux, nrvars
INTEGER, DIMENSION(1:maxdims):: amx, aul, mx, ul
REAL(KIND=r8):: h
REAL(KIND=r8), DIMENSION(1:maxdims):: xlower
!
UpdateAuxInVcycle = err_ok
!
IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
nrvars = info%nrvars
mx = 1; mx(1:ndims) = info%mx(1:ndims)
h = info%dx(1)
xlower(1:ndims) = info%xlower(1:ndims)
!
maux = info%maux; mmaux = MAX(maux,1)
amx = 1; IF(maux>0) amx(1:ndims) = mx(1:ndims)
!
mbc = info%mbc
ll = 1-mbc
 ul = 1;  ul(1:ndims) =  mx(1:ndims)+mbc
aul = 1; aul(1:ndims) = amx(1:ndims)+mbc
!
IF(maux<=0) RETURN
!
SELECT CASE(ndims)
  CASE(2)
    CALL UpdateAuxVcycle2D(info%q(   ll: ul(1),ll: ul(2),1,1:nrvars), &
                           info%qold(ll: ul(1),ll: ul(2),1,1:nrvars), &
                           info%aux( ll:aul(1),ll:aul(2),1,1:mmaux ), &
                           ll,mx(1:2),nrvars,maux,mbc,h,xlower(1:2))
  CASE(3)
    CALL UpdateAuxVcycle3D(info%q(   ll: ul(1),ll: ul(2),ll: ul(3),1:nrvars), &
                           info%qold(ll: ul(1),ll: ul(2),ll: ul(3),1:nrvars), &
                           info%aux( ll:aul(1),ll:aul(2),ll:aul(3),1:mmaux ), &
                           ll,mx(1:3),nrvars,maux,mbc,h,xlower(1:3))
  CASE DEFAULT
    PRINT *, 'UpdateAuxInVcycle: Only ndims = 2,3 are supported.'
    STOP
END SELECT
!
END FUNCTION UpdateAuxInVcycle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION UpdateSourceFunction(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok
USE Problem, ONLY: SourceUpdate2D, SourceUpdate3D
IMPLICIT NONE
!
! Dimensionally invariant source update routine.  If no update is needed set 
! SourceUpdate# = 0.
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
INTEGER:: ll, maux, mbc, mmaux, nrvars
INTEGER, DIMENSION(1:maxdims):: amx, aul, mx, ul
REAL(KIND=r8):: h
REAL(KIND=r8), DIMENSION(1:maxdims):: xlower
!
UpdateSourceFunction = err_ok
!
IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
nrvars = info%nrvars
mx = 1; mx(1:ndims) = info%mx(1:ndims)
h = info%dx(1)
xlower(1:ndims) = info%xlower(1:ndims)
!
maux = info%maux; mmaux = MAX(maux,1)
amx = 1; IF(maux>0) amx(1:ndims) = mx(1:ndims)
!
mbc = info%mbc
ll = 1-mbc
 ul = 1;  ul(1:ndims) =  mx(1:ndims)+mbc
aul = 1; aul(1:ndims) = amx(1:ndims)+mbc
!
SELECT CASE(ndims)
  CASE(2)
    info%f(1:mx(1),1:mx(2),1,1:nrvars) &
      = info%ftmp(1:mx(1),1:mx(2),1,1:nrvars) &
      + SourceUpdate2D(info%q(   ll: ul(1),ll: ul(2),1,1:nrvars), &
                       info%qold(ll: ul(1),ll: ul(2),1,1:nrvars), &
                       info%aux( ll:aul(1),ll:aul(2),1,1:mmaux ), &
                       ll,mx(1:2),nrvars,maux,mbc,h,xlower(1:2))
  CASE(3)
    info%f(1:mx(1),1:mx(2),1:mx(3),1:nrvars) &
      = info%ftmp(1:mx(1),1:mx(2),1:mx(3),1:nrvars) &
      + SourceUpdate3D(info%q(   ll: ul(1),ll: ul(2),ll: ul(3),1:nrvars), &
                       info%qold(ll: ul(1),ll: ul(2),ll: ul(3),1:nrvars), &
                       info%aux( ll:aul(1),ll:aul(2),ll:aul(3),1:mmaux ), &
                       ll,mx(1:3),nrvars,maux,mbc,h,xlower(1:3))
  CASE DEFAULT
    PRINT *, 'UpdateSourceFunction: Only ndims = 2,3 are supported.'
    STOP
END SELECT
!
END FUNCTION UpdateSourceFunction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
RECURSIVE SUBROUTINE AFASVcycle(level)
USE NodeInfoDef
USE TreeOps, ONLY: ApplyOnLevel
USE Boundary, ONLY: SetGhost, GetCoarseGhostPoints
IMPLICIT NONE
!
! Dimensionally invariant recursive AFAS vcycle routine.
!
INTEGER, INTENT(IN):: level
!
TYPE(funcparam):: dummy
!
CALL LevelRelax(level)
!
IF(level>minlevel) THEN
!
  CALL ApplyOnLevel(level,RestrictSolution,dummy)
!
  CALL SetGhost(level-1,0)
!
  CALL ApplyOnLevel(level,GetCoarseGhostPoints,dummy)
!
  CALL GetCoarseLevelLoading(level)
!
  CALL AFASVcycle(level-1)
!
  CALL ApplyOnLevel(level,CorrectFine,dummy) 
!
  CALL SetGhost(level,0)
!
  CALL LevelRelax(level)
!
END IF
!
END SUBROUTINE AFASVcycle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE LevelRelax(level)
USE NodeInfoDef
USE TreeOps, ONLY: ApplyOnLevel
USE Boundary, ONLY: SetGhost
IMPLICIT NONE
!
! Dimensionally invariant relaxation.
!
INTEGER, INTENT(IN):: level
!
TYPE(funcparam):: dummy
INTEGER:: redblack, smoothingpass
!
! Assume Ghost points are set on entry to this routine.
DO smoothingpass = 1, nsmoothingpasses
  DO redblack = 1, 2
!
  dummy%iswitch = redblack
  CALL ApplyOnLevel(level,RelaxPatch,dummy)
!
  CALL SetGhost(level,redblack)
!
  END DO
END DO
!
END SUBROUTINE LevelRelax
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION RelaxPatch(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok
USE Problem, ONLY: RelaxGrid2d, RelaxGrid3D
IMPLICIT NONE
!
! Dimensionally invariant relaxation of a patch, on RED squares for
! dummy%iswitch = 1, and BLACK squares for dummy%iswitch = 2.
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
INTEGER:: maux, mmaux, nrvars
INTEGER, DIMENSION(1:maxdims):: amx, mx
REAL(KIND=r8):: h
REAL(KIND=r8), DIMENSION(1:maxdims):: xlower
!
RelaxPatch = err_ok
!
IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
nrvars = info%nrvars
mx = 1; mx(1:ndims) = info%mx(1:ndims)
h = info%dx(1)
xlower(1:ndims) = info%xlower(1:ndims)
!
maux = info%maux; mmaux = MAX(maux,1)
amx = 1; IF(maux>0) amx(1:ndims) = mx(1:ndims)
!
SELECT CASE(ndims)
  CASE(2)
    CALL RelaxGrid2D(info%q(   0: mx(1)+1,0: mx(2)+1,1,1:nrvars), &
                     info%qold(0: mx(1)+1,0: mx(2)+1,1,1:nrvars), &
                     info%aux( 0:amx(1)+1,0:amx(2)+1,1,1:mmaux ), &
                     info%f(   1: mx(1),  1: mx(2)  ,1,1:nrvars), &
                     mx(1:2),nrvars,maux,h,xlower(1:2),dummy%iswitch)
  CASE(3)
    CALL RelaxGrid3D(info%q(   0: mx(1)+1,0: mx(2)+1,0: mx(3)+1,1:nrvars), &
                     info%qold(0: mx(1)+1,0: mx(2)+1,0: mx(3)+1,1:nrvars), &
                     info%aux( 0:amx(1)+1,0:amx(2)+1,0:amx(3)+1,1:mmaux ), &
                     info%f(   1: mx(1),  1: mx(2)  ,1: mx(3)  ,1:nrvars), &
                     mx(1:3),nrvars,maux,h,xlower(1:3),dummy%iswitch)
  CASE DEFAULT
    PRINT *, 'RelaxPatch: Only ndims = 2,3 are supported.'
    STOP
END SELECT
!
END FUNCTION RelaxPatch
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION RestrictSolution(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok, GetParentInfo
USE GridUtilities, ONLY: Restriction2D, Restriction3D
IMPLICIT NONE
!
! Dimensionally invariant restriction.
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
TYPE(nodeinfo), POINTER:: parent
INTEGER:: ierror, nrvars
INTEGER, DIMENSION(1:maxdims):: cmx, mx
INTEGER, DIMENSION(1:maxdims,1:2):: mb
!
RestrictSolution = err_ok
!
IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
nrvars = info%nrvars
mx = 1; mx(1:ndims) = info%mx(1:ndims)
cmx = 1; cmx(1:ndims) = mx(1:ndims)/2
mb = 1; mb(1:ndims,1:2) = info%mbounds(1:ndims,1:2)
!
ierror = GetParentInfo(parent)
!
SELECT CASE(ndims)
  CASE(2)
    info%qc(1:cmx(1),1:cmx(2),1,1:nrvars) &
      = Restriction2D(info%q(1:mx(1),1:mx(2),1,1:nrvars))
  CASE(3)
    info%qc(1:cmx(1),1:cmx(2),1:cmx(3),1:nrvars) &
      = Restriction3D(info%q(1:mx(1),1:mx(2),1:mx(3),1:nrvars))
  CASE DEFAULT
    PRINT *, 'BSAM 1.2 RestrictSolution: Only ndims = 2,3 are supported.'
    STOP
END SELECT
!
parent%q(mb(1,1):mb(1,2),mb(2,1):mb(2,2),mb(3,1):mb(3,2),1:nrvars) &
  = info%qc(1:cmx(1),1:cmx(2),1:cmx(3),1:nrvars)
!
END FUNCTION RestrictSolution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE GetCoarseLevelLoading(level)
USE NodeInfoDef
USE TreeOps, ONLY: ApplyOnLevel
IMPLICIT NONE
!
INTEGER, INTENT(IN):: level
!
TYPE(funcparam):: dummy
!
! 1) Calculate the standard coarse load function in AFAS:
CALL ApplyOnLevel(level,CoarseLoadingFunction,dummy)
!
! 2) Correct the coarse load function to preserve mass:
! Will be implemented in BSAM 1.2.
!
END SUBROUTINE GetCoarseLevelLoading
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION CoarseLoadingFunction(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok, GetParentInfo
USE GridUtilities, ONLY: Restriction2D, Restriction3D
USE Problem, ONLY: Operator2D, Operator3D
IMPLICIT NONE
!
! Dimensionally invariant coarse loading function routine.
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
TYPE(nodeinfo), POINTER:: parent
INTEGER:: ierror, maux, mmaux, nrvars
INTEGER, DIMENSION(1:maxdims):: cmx, mx
INTEGER, DIMENSION(1:maxdims,1:2):: amb, mb
REAL(KIND=r8):: ch
REAL(KIND=r8), DIMENSION(1:maxdims):: xlower
!
CoarseLoadingFunction = err_ok
!
IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
nrvars = info%nrvars
mx = 1; mx(1:ndims) = info%mx(1:ndims)
cmx = 1; cmx(1:ndims) = mx(1:ndims)/2
mb = 1; mb(1:ndims,1:2) = info%mbounds(1:ndims,1:2)
!
ierror = GetParentInfo(parent)
ch = parent%dx(1)
xlower(1:ndims) = info%xlower(1:ndims)
!
maux = info%maux; mmaux = MAX(maux,1)
amb = 1; IF(maux>0) amb(1:ndims,1:2) = mb(1:ndims,1:2)
!
CALL Residual(info)
!
SELECT CASE(ndims)
  CASE(2)
    parent%f(mb(1,1):mb(1,2),mb(2,1):mb(2,2),1,1:nrvars) &
      = Restriction2D(info%rf(1:mx(1),1:mx(2),1,1:nrvars)) &
      + Operator2D(info%qc(0:cmx(1)+1,0:cmx(2)+1,1,1:nrvars), &
                   parent%qold(mb(1,1)-1: mb(1,2)+1, &
                               mb(2,1)-1: mb(2,2)+1,1,1:nrvars), &
                   parent%aux(amb(1,1)-1:amb(1,2)+1, &
                              amb(2,1)-1:amb(2,2)+1,1,1:mmaux ), &
                   cmx(1:2),nrvars,maux,ch,xlower(1:2))
  CASE(3)
    parent%f(mb(1,1):mb(1,2),mb(2,1):mb(2,2),mb(3,1):mb(3,2),1:nrvars) &
      = Restriction3D(info%rf(1:mx(1),1:mx(2),1:mx(3),1:nrvars)) &
      + Operator3D(info%qc(0:cmx(1)+1,0:cmx(2)+1,0:cmx(3)+1,1:nrvars), &
                   parent%qold(mb(1,1)-1: mb(1,2)+1, &
                               mb(2,1)-1: mb(2,2)+1, &
                               mb(3,1)-1: mb(3,2)+1,1:nrvars), &
                   parent%aux(amb(1,1)-1:amb(1,2)+1, &
                              amb(2,1)-1:amb(2,2)+1, &
                              amb(3,1)-1:amb(3,2)+1,1:mmaux ), &
                   cmx(1:3),nrvars,maux,ch,xlower(1:3))
  CASE DEFAULT
    PRINT *, 'CoarseLoadingFunction: Only ndims = 2,3 are supported.'
    STOP
END SELECT
!
END FUNCTION CoarseLoadingFunction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION CoarseLoadingCorrection(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok, GetParentInfo
USE GridUtilities, ONLY: Restriction2D, Restriction3D
USE ProblemDef
IMPLICIT NONE
!
! Dimensionally invariant coarse loading correction function routine.
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
TYPE(nodeinfo), POINTER:: parent
LOGICAL:: masscorrneeded
INTEGER:: i, ierror, j, level, maux, mmaux, nrvars
INTEGER, DIMENSION(1:maxdims):: cmx, mx, pmx
INTEGER, DIMENSION(1:maxdims,1:2):: amb, mb, mg
REAL(KIND=r8):: ch, ch2
REAL(KIND=r8), DIMENSION(1:maxnrvars):: diffcons
REAL(KIND=r8), DIMENSION(1:maxdims):: xlower
REAL(KIND=r8), DIMENSION(:,:,:,:), ALLOCATABLE:: zombiediff1, zombiediff2, zombiediff3
!:
CoarseLoadingCorrection = err_ok
!
IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
nrvars = info%nrvars
mx = 1; mx(1:ndims) = info%mx(1:ndims)
cmx = 1; cmx(1:ndims) = mx(1:ndims)/2
mb = 1; mb(1:ndims,1:2) = info%mbounds(1:ndims,1:2)
!
ierror = GetParentInfo(parent)
pmx = 1; pmx(1:ndims) = parent%mx(1:ndims)
ch = parent%dx(1); ch2 = ch*ch
level = info%level
xlower(1:ndims) = info%xlower(1:ndims)
!
maux = info%maux; mmaux = MAX(maux,1)
amb = 1; IF(maux>0) amb(1:ndims,1:2) = mb(1:ndims,1:2)
!
! These need to come from the user. This is only for testing purposes:
diffcons(1) = dt
diffcons(2) = -eps2
!
SELECT CASE(ndims)
CASE(2)
  ALLOCATE(zombiediff1(1:2     ,1:cmx(2),1:1,1:nrvars), &
           zombiediff2(1:cmx(1),1:2     ,1:1,1:nrvars))
  zombiediff1 = 0.0_r8; zombiediff2 = 0.0_r8
!
! Bottom correction:
  masscorrneeded = .FALSE.
  DO i = 1, cmx(1)
    IF(info%levellandscape(i,0,1) == level-1) THEN
      zombiediff2(i,1,1,1:nrvars) &
        = (info%q(2*i  ,1,1,1:nrvars)-info%q(2*i  ,0,1,1:nrvars) &
        +  info%q(2*i-1,1,1,1:nrvars)-info%q(2*i-1,0,1,1:nrvars) &
        +  parent%q(mb(1,1)+i-1,mb(2,1)-1,1,1:nrvars) &
        -  info%qc(i,1,1,1:nrvars))/ch2
!
      zombiediff2(i,1,1,1) = zombiediff2(i,1,1,2)*diffcons(1)
      zombiediff2(i,1,1,2) = zombiediff2(i,1,1,1)*diffcons(2)
!
      masscorrneeded = .TRUE.
    END IF
  END DO
  IF(masscorrneeded) THEN
    IF(mb(2,1)-1 >= 1) THEN
        parent%f(mb(1,1):mb(1,2),mb(2,1)-1,1,1:nrvars) &
      = parent%f(mb(1,1):mb(1,2),mb(2,1)-1,1,1:nrvars) &
      + zombiediff2(1:cmx(1),1,1,1:nrvars)
    ELSE
      mg = 1
      mg(1,1:2) = parent%mglobal(1,1)+mb(1,1:2)-1
      mg(2,1:2) = parent%mglobal(2,1)-1
      CALL AddStripToList(mg,nrvars,zombiediff2(1:cmx(1),1:1,1:1,1:nrvars))
    END IF
  END IF
!
! Top correction:
  masscorrneeded = .FALSE.
  DO i = 1, cmx(1)
    IF(info%levellandscape(i,cmx(2)+1,1) == level-1) THEN
      zombiediff2(i,2,1,1:nrvars) &
        = (info%q(2*i  ,mx(2),1,1:nrvars)-info%q( 2*i  ,mx(2)+1,1,1:nrvars) &
        +  info%q(2*i-1,mx(2),1,1:nrvars)-info%q( 2*i-1,mx(2)+1,1,1:nrvars) &
        +  parent%q(mb(1,1)+i-1,mb(2,2)+1,1,1:nrvars) &
        -  info%qc(i,cmx(2),1,1:nrvars))/ch2
!
      zombiediff2(i,2,1,1) = zombiediff2(i,2,1,2)*diffcons(1)
      zombiediff2(i,2,1,2) = zombiediff2(i,2,1,1)*diffcons(2)
!
      masscorrneeded = .TRUE.
    END IF
  END DO
  IF(masscorrneeded) THEN
    IF(mb(2,2)+1 <= pmx(2)) THEN
        parent%f(mb(1,1):mb(1,2),mb(2,2)+1,1,1:nrvars) &
      = parent%f(mb(1,1):mb(1,2),mb(2,2)+1,1,1:nrvars) &
      + zombiediff2(1:cmx(1),2,1,1:nrvars)
    ELSE
      mg = 1
      mg(1,1:2) = parent%mglobal(1,1)+mb(1,1:2)-1
      mg(2,1:2) = parent%mglobal(2,1)+pmx(2)
      CALL AddStripToList(mg,nrvars,zombiediff2(1:cmx(1),2:2,1:1,1:nrvars))
    END IF
  END IF
!
! Left Correction:
  masscorrneeded = .FALSE.
  DO j = 1, cmx(2)
    IF(info%levellandscape(0,j,1) == level-1) THEN
      zombiediff1(1,j,1,1:nrvars) &
        = (info%q(1,2*j  ,1,1:nrvars)-info%q(0,2*j  ,1,1:nrvars) &
        +  info%q(1,2*j-1,1,1:nrvars)-info%q(0,2*j-1,1,1:nrvars) &
        +  parent%q(mb(1,1)-1,mb(2,1)+j-1,1,1:nrvars) &
        -  info%qc(1,j,1,1:nrvars))/ch2
!
      zombiediff1(1,j,1,1) = zombiediff1(1,j,1,2)*diffcons(1)
      zombiediff1(1,j,1,2) = zombiediff1(1,j,1,1)*diffcons(2)
!
      masscorrneeded = .TRUE.
    END IF
  END DO
  IF(masscorrneeded) THEN
    IF(mb(1,1)-1 >= 1) THEN
        parent%f(mb(1,1)-1,mb(2,1):mb(2,2),1,1:nrvars) &
      = parent%f(mb(1,1)-1,mb(2,1):mb(2,2),1,1:nrvars) &
      + zombiediff1(1,1:cmx(2),1,1:nrvars)
    ELSE
      mg = 1
      mg(1,1:2) = parent%mglobal(1,1)-1
      mg(2,1:2) = parent%mglobal(2,1)+mb(2,1:2)-1
      CALL AddStripToList(mg,nrvars,zombiediff1(1:1,1:cmx(2),1:1,1:nrvars))
    END IF
  END IF
!
! Right Correction:
  masscorrneeded = .FALSE.
  DO j = 1, cmx(2)
    IF(info%levellandscape(cmx(1)+1,j,1) == level-1) THEN
      zombiediff1(2,j,1,1:nrvars) &
        = (info%q(mx(1),2*j  ,1,1:nrvars)-info%q(mx(1)+1,2*j  ,1,1:nrvars) &
        +  info%q(mx(1),2*j-1,1,1:nrvars)-info%q(mx(1)+1,2*j-1,1,1:nrvars) &
        +  parent%q(mb(1,2)+1,mb(2,1)+j-1,1,1:nrvars) &
        -  info%qc(cmx(1),j,1,1:nrvars))/ch2
!
      zombiediff1(2,j,1,1) = zombiediff1(2,j,1,2)*diffcons(1)
      zombiediff1(2,j,1,2) = zombiediff1(2,j,1,1)*diffcons(2)
!
      masscorrneeded = .TRUE.
    END IF
  END DO
  IF(masscorrneeded) THEN
    IF(mb(1,2)+1 <= pmx(1)) THEN
        parent%f(mb(1,2)+1,mb(2,1):mb(2,2),1,1:nrvars) &
      = parent%f(mb(1,2)+1,mb(2,1):mb(2,2),1,1:nrvars) &
      + zombiediff1(2,1:cmx(2),1,1:nrvars)
    ELSE
      mg = 1
      mg(1,1:2) = parent%mglobal(1,1)+pmx(1)
      mg(2,1:2) = parent%mglobal(2,1)+mb(2,1:2)-1
      CALL AddStripToList(mg,nrvars,zombiediff1(2:2,1:cmx(2),1:1,1:nrvars))
    END IF
  END IF
  DEALLOCATE(zombiediff1,zombiediff2)
CASE(3)
  CONTINUE
CASE DEFAULT
  PRINT *, 'CoarseLoadingFunction: Only ndims = 2,3 are supported.'
  STOP
END SELECT
!
END FUNCTION CoarseLoadingCorrection
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE AddStripToList(mg,nrvars,masscorr)
USE NodeInfoDef
IMPLICIT NONE
!
INTEGER, DIMENSION(1:maxdims,1:2), INTENT(IN):: mg
INTEGER, INTENT(IN):: nrvars
REAL(KIND=r8), DIMENSION(mg(1,1):mg(1,2), &
                         mg(2,1):mg(2,2), &
                         mg(3,1):mg(3,2),1:nrvars), INTENT(IN):: masscorr
!
! This routine adds the mass correction strip to a global linked list. Each
! strip is indexed with global coordinates.
!
nmasscorrlayers =  nmasscorrlayers+1
ALLOCATE(currentlayer)
currentlayer%id = nmasscorrlayers
currentlayer%mg(1:maxdims,1:2) = mg(1:maxdims,1:2)
!
ALLOCATE(currentlayer%masscorr(mg(1,1):mg(1,2), &
                               mg(2,1):mg(2,2), &
                               mg(3,1):mg(3,2),1:nrvars))
!
currentlayer%masscorr(mg(1,1):mg(1,2),mg(2,1):mg(2,2),mg(3,1):mg(3,2),1:nrvars) &
           = masscorr(mg(1,1):mg(1,2),mg(2,1):mg(2,2),mg(3,1):mg(3,2),1:nrvars)
!
currentlayer%prevlayer => lastlayer
lastlayer => currentlayer
!
END SUBROUTINE AddStripToList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE TransferMassCorrectionLayers(level)
USE NodeInfoDef
USE TreeOps, ONLY: ApplyOnLevel
USE Boundary, ONLY: GetPeriodicOffsets
IMPLICIT NONE
!
INTEGER, INTENT(IN):: level
!
TYPE(funcparam):: dummy
LOGICAL:: periodicbuffer
INTEGER:: n, offset, polarity
INTEGER, DIMENSION(1:maxdims,1:2):: mgsave
!
mgsave = 1
currentlayer => lastlayer
searchloop: DO
  IF(.NOT. ASSOCIATED(currentlayer%prevlayer)) EXIT searchloop
!
  dummy%offset = 0
  CALL ApplyOnLevel(level,TransferLoads,dummy)
!
! Buffering of periodic edge tags. Check to see if the buffer area cuts across 
! a periodic boundary.  If so, add offset and apply buffer:
  mgsave(1:ndims,1:2) = currentlayer%mg(1:ndims,1:2)
  IF(periodicboundaryconditions) THEN
    CALL GetPeriodicOffsets(level)
    DO polarity = -1, 1, 2
      DO offset = 1, nperiodicoffsets
        dummy%offset = 0
        DO n = 1, ndims
          dummy%offset(n) = polarity*poffset(n,offset)
          currentlayer%mg(n,1:2) = mgsave(n,1:2)+dummy%offset(n)
        END DO
        CALL ApplyOnLevel(level,TransferLoads,dummy)
      END DO
    END DO
  END IF
  currentlayer%mg(1:ndims,1:2) = mgsave(1:ndims,1:2)
!
  currentlayer => currentlayer%prevlayer
END DO searchloop
!
END SUBROUTINE TransferMassCorrectionLayers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION TransferLoads(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok
IMPLICIT NONE
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
INTEGER:: n, nrvars
INTEGER, DIMENSION(1:maxdims,1:2):: mglobal, mlf, mls, movr
!
TransferLoads = err_ok
IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
nrvars = info%nrvars
mglobal = 1; mglobal(1:ndims,1:2) = info%mglobal(1:ndims,1:2)
!
! 1. Find overlap region in global index space:
movr = 1
movr(1:ndims,1) = MAX(currentlayer%mg(1:ndims,1),mglobal(1:ndims,1))
movr(1:ndims,2) = MIN(currentlayer%mg(1:ndims,2),mglobal(1:ndims,2))
!
! 2. Check for nonempty intersection:
IF(ANY(movr(1:maxdims,2)-movr(1:maxdims,1)<0)) RETURN
!
! 3. Transform common index space to grid index spaces:
mlf = 1; mls = 1
DO n = 1, ndims
  mlf(n,1:2) = movr(n,1:2)-mglobal(n,1)+1
  mls(n,1:2) = movr(n,1:2)-dummy%offset(n)
END DO
!
                 info%f(mlf(1,1):mlf(1,2),  &
                        mlf(2,1):mlf(2,2),  &
                        mlf(3,1):mlf(3,2),1:nrvars)  &
=                info%f(mlf(1,1):mlf(1,2),  &
                        mlf(2,1):mlf(2,2),  &
                        mlf(3,1):mlf(3,2),1:nrvars)  & 
+ currentlayer%masscorr(mls(1,1):mls(1,2),  &
                        mls(2,1):mls(2,2),  &
                        mls(3,1):mls(3,2),1:nrvars)
!
END FUNCTION TransferLoads
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE DeleteMassCorrectionList
USE NodeInfoDef
IMPLICIT NONE
!
searchloop: DO
  currentlayer => lastlayer%prevlayer
  DEALLOCATE(lastlayer%masscorr)
  DEALLOCATE(lastlayer)
  lastlayer => currentlayer
  IF(.NOT. ASSOCIATED(lastlayer%prevlayer)) EXIT searchloop
END DO searchloop
!
END SUBROUTINE DeleteMassCorrectionList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION RelativeTruncationError(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok, GetParentInfo
USE GridUtilities, ONLY: Restriction2D, Restriction3D
USE Problem, ONLY: Operator2D, Operator3D
IMPLICIT NONE
!
! Dimensionally invariant relative trunctation error:
!
!   L_{2h}(I_h^{2h} q_h)-I_h^{2h}(L_h q_h).
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
TYPE(nodeinfo), POINTER:: parent
INTEGER:: ierror, maux, mmaux, nrvars
INTEGER, DIMENSION(1:maxdims):: amx, cmx, mx
INTEGER, DIMENSION(1:maxdims,1:2):: amb, mb
REAL(KIND=r8):: ch, h
REAL(KIND=r8), DIMENSION(1:maxdims):: xlower
!
RelativeTruncationError = err_ok
!
IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
nrvars = info%nrvars
mx = 1; mx(1:ndims) = info%mx(1:ndims)
cmx = 1; cmx(1:ndims) = mx(1:ndims)/2
mb = 1; mb(1:ndims,1:2) = info%mbounds(1:ndims,1:2)
h = info%dx(1)
xlower(1:ndims) = info%xlower(1:ndims)
!
maux = info%maux; mmaux = MAX(maux,1)
amx = 1; IF(maux>0) amx(1:ndims) = mx(1:ndims)
!
ierror = GetParentInfo(parent)
ch = parent%dx(1)
!
maux = info%maux; mmaux = MAX(maux,1)
amb = 1; IF(maux>0) amb(1:ndims,1:2) = mb(1:ndims,1:2)
!
SELECT CASE(ndims)
  CASE(2)
    info%qrte(1:cmx(1),1:cmx(2),1,1:nrvars) &
      = Restriction2D( &
        Operator2D(info%q(   0: mx(1)+1,0: mx(2)+1,1,1:nrvars), &
                   info%qold(0: mx(1)+1,0: mx(2)+1,1,1:nrvars), &
                   info%aux( 0:amx(1)+1,0:amx(2)+1,1,1:mmaux ), &
                   mx(1:2),nrvars,maux,h,xlower(1:2))) &
      - Operator2D(info%qc(0:cmx(1)+1,0:cmx(2)+1,1,1:nrvars), &
                   parent%qold(mb(1,1)-1: mb(1,2)+1, &
                               mb(2,1)-1: mb(2,2)+1,1,1:nrvars), &
                   parent%aux(amb(1,1)-1:amb(1,2)+1, &
                              amb(2,1)-1:amb(2,2)+1,1,1:mmaux ), &
                   cmx(1:2),nrvars,maux,ch,xlower(1:2))
  CASE(3)
    info%qrte(1:cmx(1),1:cmx(2),1:cmx(3),1:nrvars) &
      = Restriction3D( &
        Operator3D(info%q(   0: mx(1)+1,0: mx(2)+1,0: mx(3)+1,1:nrvars), &
                   info%qold(0: mx(1)+1,0: mx(2)+1,0: mx(3)+1,1:nrvars), &
                   info%aux( 0:amx(1)+1,0:amx(2)+1,0:amx(3)+1,1:mmaux ), &
                   mx(1:3),nrvars,maux,h,xlower(1:3))) &
      - Operator3D(info%qc(0:cmx(1)+1,0:cmx(2)+1,0:cmx(3)+1,1:nrvars), &
                   parent%qold(mb(1,1)-1: mb(1,2)+1, &
                               mb(2,1)-1: mb(2,2)+1, &
                               mb(3,1)-1: mb(3,2)+1,1:nrvars), &
                   parent%aux(amb(1,1)-1:amb(1,2)+1, &
                              amb(2,1)-1:amb(2,2)+1, &
                              amb(3,1)-1:amb(3,2)+1,1:mmaux ), &
                   cmx(1:3),nrvars,maux,ch,xlower(1:3))
  CASE DEFAULT
    PRINT *, 'CoarseLoadingFunction: Only ndims = 2,3 are supported.'
    STOP
END SELECT
!
END FUNCTION RelativeTruncationError
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE Residual(info)
USE NodeInfoDef
USE Problem, ONLY: Operator2D, Operator3D
IMPLICIT NONE
!
TYPE(nodeinfo):: info
!
INTEGER:: maux, mmaux, nrvars
INTEGER, DIMENSION(1:maxdims):: amx, mx
REAL(KIND=r8):: h
REAL(KIND=r8), DIMENSION(1:maxdims):: xlower
!
nrvars = info%nrvars
mx = 1; mx(1:ndims) = info%mx(1:ndims)
h = info%dx(1)
xlower(1:ndims) = info%xlower(1:ndims)
!
maux = info%maux; mmaux = MAX(maux,1)
amx = 1; IF(maux>0) amx(1:ndims) = mx(1:ndims)
!
SELECT CASE(ndims)
  CASE(2)
    info%rf(1:mx(1),1:mx(2),1,1:nrvars) &
      = info%f(1:mx(1),1:mx(2),1,1:nrvars) &
      - Operator2D(info%q(   0: mx(1)+1,0: mx(2)+1,1,1:nrvars), &
                   info%qold(0: mx(1)+1,0: mx(2)+1,1,1:nrvars), &
                   info%aux( 0:amx(1)+1,0:amx(2)+1,1,1:mmaux ), &
                   mx(1:2),nrvars,maux,h,xlower(1:2))
  CASE(3)
    info%rf(1:mx(1),1:mx(2),1:mx(3),1:nrvars) &
      = info%f(1:mx(1),1:mx(2),1:mx(3),1:nrvars) &
      - Operator3D(info%q(   0: mx(1)+1,0: mx(2)+1,0: mx(3)+1,1:nrvars), &
                   info%qold(0: mx(1)+1,0: mx(2)+1,0: mx(3)+1,1:nrvars), &
                   info%aux( 0:amx(1)+1,0:amx(2)+1,0:amx(3)+1,1:mmaux ), &
                   mx(1:3),nrvars,maux,h,xlower(1:3))
  CASE DEFAULT
    PRINT *, 'Residual: Only ndims = 2,3 are supported.'
    STOP
END SELECT
!
END SUBROUTINE Residual
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION CorrectFine(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok, GetParentInfo
USE GridUtilities, ONLY:       Prolongation2D,       Prolongation3D, &
                          BiLinProlongationP1,  BiLinProlongationP2, &
                         TriLinProlongationP1, TriLinProlongationP2
IMPLICIT NONE
!
! Dimensionally invariant coarse-grid-correction routine.
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
TYPE(nodeinfo), POINTER:: parent
INTEGER:: ierror, ll, mbc, nrvars
INTEGER, DIMENSION(1:maxdims):: cmx, mx, ul
INTEGER, DIMENSION(1:maxdims,1:2):: mb
!
CorrectFine = err_ok
!
IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
ierror = GetParentInfo(parent)
!
mbc = info%mbc
nrvars = info%nrvars
mx = 1; mx(1:ndims) = info%mx(1:ndims)
cmx = 1; cmx(1:ndims) = mx(1:ndims)/2
mb = 1; mb(1:ndims,1:2) = info%mbounds(1:ndims,1:2)
!
SELECT CASE(ndims)
  CASE(2)
       info%qc(1-mbc:cmx(1)+mbc,1-mbc:cmx(2)+mbc,1,1:nrvars) &
    = parent%q(mb(1,1)-mbc:mb(1,2)+mbc, &
               mb(2,1)-mbc:mb(2,2)+mbc,1,1:nrvars) &
    -  info%qc(1-mbc:cmx(1)+mbc,1-mbc:cmx(2)+mbc,1,1:nrvars)
!
  SELECT CASE(mbc)
    CASE(1)
!
! Bilinear prolongation:
          info%q( 0: mx(1)+1,0: mx(2)+1,1,1:nrvars) &
        = info%q( 0: mx(1)+1,0: mx(2)+1,1,1:nrvars) &
        + BiLinProlongationP1( &
          info%qc(0:cmx(1)+1,0:cmx(2)+1,1,1:nrvars))
!
! Mass conserving bilinear prolongation:
!          info%q( 0: mx(1)+1,0: mx(2)+1,1,1:nrvars) &
!        = info%q( 0: mx(1)+1,0: mx(2)+1,1,1:nrvars) &
!        + BiLinProlongationP1MC( &
!          info%qc(0:cmx(1)+1,0:cmx(2)+1,1,1:nrvars))
!
! Simple injection:
!          info%q( 1: mx(1),1: mx(2),1,1:nrvars) &
!        = info%q( 1: mx(1),1: mx(2),1,1:nrvars) &
!        + Prolongation2D( &
!          info%qc(1:cmx(1),1:cmx(2),1,1:nrvars))
!
    CASE(2)
!
! Bilinear prolongation:
          info%q( -1: mx(1)+2,-1: mx(2)+2,1,1:nrvars) &
        = info%q( -1: mx(1)+2,-1: mx(2)+2,1,1:nrvars) &
        + BiLinProlongationP2( &
          info%qc(-1:cmx(1)+2,-1:cmx(2)+2,1,1:nrvars))
!
! Mass-conserving bilinear prolongation:
!          info%q( -1: mx(1)+2,-1: mx(2)+2,1,1:nrvars) &
!        = info%q( -1: mx(1)+2,-1: mx(2)+2,1,1:nrvars) &
!        + BiLinProlongationP2MC( &
!          info%qc(-1:cmx(1)+2,-1:cmx(2)+2,1,1:nrvars))
  END SELECT
!
  CASE(3)
      info%qc(1-mbc:cmx(1)+mbc,1-mbc:cmx(2)+mbc,1-mbc:cmx(3)+mbc,1:nrvars) &
    = parent%q(mb(1,1)-mbc:mb(1,2)+mbc, &
               mb(2,1)-mbc:mb(2,2)+mbc, &
               mb(3,1)-mbc:mb(3,2)+mbc,1:nrvars) &
    - info%qc(1-mbc:cmx(1)+mbc,1-mbc:cmx(2)+mbc,1-mbc:cmx(3)+mbc,1:nrvars)
!
  SELECT CASE(mbc)
    CASE(1)
!
! Trilinear prolongation:
        info%q( 0: mx(1)+1,0: mx(2)+1,0: mx(3)+1,1:nrvars) &
      = info%q( 0: mx(1)+1,0: mx(2)+1,0: mx(3)+1,1:nrvars) &
      + TriLinProlongationP1( &
        info%qc(0:cmx(1)+1,0:cmx(2)+1,0:cmx(3)+1,1:nrvars))
!
! Mass conserving trilinear prolongation:
!        info%q( 0: mx(1)+1,0: mx(2)+1,0: mx(3)+1,1:nrvars) &
!      = info%q( 0: mx(1)+1,0: mx(2)+1,0: mx(3)+1,1:nrvars) &
!      + TriLinProlongationP1MC( &
!        info%qc(0:cmx(1)+1,0:cmx(2)+1,0:cmx(3)+1,1:nrvars))
!
! Simple injection:
!        info%q( 1: mx(1),1: mx(2),1: mx(3),1:nrvars) &
!      = info%q( 1: mx(1),1: mx(2),1: mx(3),1:nrvars) &
!      + Prolongation3D( &
!        info%qc(1:cmx(1),1:cmx(2),1:cmx(3),1:nrvars))
    CASE(2)
!
! Trilinear prolongation:
        info%q( -1: mx(1)+2,-1: mx(2)+2,-1: mx(3)+2,1:nrvars) &
      = info%q( -1: mx(1)+2,-1: mx(2)+2,-1: mx(3)+2,1:nrvars) &
      + TriLinProlongationP2( &
        info%qc(-1:cmx(1)+2,-1:cmx(2)+2,-1:cmx(3)+2,1:nrvars))
!
! Mass-conserving trilinear prolongation:
!        info%q( -1: mx(1)+2,-1: mx(2)+2,-1: mx(3)+2,1:nrvars) &
!      = info%q( -1: mx(1)+2,-1: mx(2)+2,-1: mx(3)+2,1:nrvars) &
!      + TriLinProlongationP2( &
!        info%qc(-1:cmx(1)+2,-1:cmx(2)+2,-1:cmx(3)+2,1:nrvars))
  END SELECT
!
  CASE DEFAULT
    PRINT *, 'CorrectFine: Only ndims=2,3 are supported.'
    STOP
END SELECT
!
END FUNCTION CorrectFine
END MODULE AFASRoutines