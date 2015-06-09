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
! File:             bsamroutines.f90
! Purpose:          BSAM control module.
! Contains:
! Revision History: Ver. 1.0 Oct. 2006 Steven Wise
! Revision History: Ver. 1.1 May. 2007 Steven Wise
! Revision History: Ver. 1.2 Jul. 2015 Steven Wise
! -----------------------------------------------------------------------
MODULE BSAMRoutines
IMPLICIT NONE
!
SAVE
PRIVATE
PUBLIC BSAMSolver
!
CONTAINS
!
SUBROUTINE BSAMSolver
USE NodeInfoDef
USE TreeOps, ONLY: AddRootLevelNode, ApplyOnForest, ApplyOnLevel, &
                   CreateBelowSeedLevels, DeleteMarkedNode, InitForest, &
                   KillForest
USE BSAMInputOutput, ONLY: ReadQ
USE BSAMStorage, ONLY: DeallocPeriodicBCStorage
USE Problem, ONLY: SetProb, AfterRun
USE Boundary, ONLY: SetGhost
IMPLICIT NONE
!
TYPE(funcparam):: dummy
CHARACTER(LEN=5):: zone
CHARACTER(LEN=8):: date
CHARACTER(LEN=10):: time
INTEGER:: i, ierror, level, ilevel
INTEGER, DIMENSION(1:8):: values
!
NAMELIST/rundata/dt, errortype, getafterstepstats, maxvcycles, &
                 nsmoothingpasses, omega, outframes, outputuniformmesh, &
                 qerrortol, restart, restartframe, syncelliptic, &
                 timeiterations, updateauxfreq
!
! Initializations:
errortype = 1
dt = 0.0_r8
getafterstepstats = .FALSE.
maxvcycles = 20
nsmoothingpasses = 2
omega = 1.0_r8
outframes = 1
outputuniformmesh = .FALSE.
qerrortol = 1.0E-06_r8
restart = .FALSE.
restartframe = 0
syncelliptic = .FALSE.
timeiterations = 1
updateauxfreq = 1
!
! Read general input data:
OPEN(UNIT=75,FILE='rundata.dat',STATUS='OLD',ACTION='READ',IOSTAT=ierror)
IF(ierror/=0) THEN
  PRINT *,'Error opening input file rundata.dat. Program stop.'
  STOP
END IF
READ(UNIT=75,NML=rundata)
CLOSE(75)
!
CALL DATE_AND_TIME(date,time,zone,values)
OPEN(UNIT=76,FILE='output.dat',STATUS='UNKNOWN',ACTION='WRITE', &
     FORM='FORMATTED',POSITION='APPEND')
WRITE(76,1001) date, time
1001 FORMAT(' '/'New run at date ',A8,' and time ',A10/' ')
WRITE(76,NML=rundata)
CLOSE(76)
!
! By default only one rootlevel grid.  This may change in the future:
nrootgrids = 1
!
! Initialize a forest of trees.  One root level grid generated in this call:
CALL InitForest
!
! Read data file to initialize root level grids:
CALL ApplyOnLevel(rootlevel,RootInit,dummy)
!
! Create the levels below the seed needed for multigrid:
CALL CreateBelowSeedLevels(minlevel)
!
! Initialize the levels below the forest seed:
DO ilevel = rootlevel-1, minlevel, -1
  CALL ApplyOnLevel(ilevel,InitSeed,dummy)
END DO
!
! Set user problem parameters:
CALL SetProb
!
! Initialize the rootlevel data fields.  In the case of a restart, we read the
! fields at all above root levels, saving the data to uniformgrid(level)%q':
IF(restart) THEN
  PRINT *, 'Restart of computation from plot frame ', restartframe
  CALL ReadQ
  outputinitialdata = .FALSE.
  syncelliptic = .FALSE.
ELSE
  CALL ApplyOnLevel(rootlevel,Initialq,dummy)
  outputinitialdata = .TRUE.
END IF
!
CALL SetGhost(rootlevel,0)
CALL ApplyOnLevel(rootlevel,SetAuxFields,dummy)
CALL ApplyOnLevel(rootlevel,SetSrcFields,dummy)
CALL ApplyOnLevel(rootlevel,CopyQToQold,dummy)
!
! Initialization complete. Start run:
PRINT *, 'BSAM 1.2 is running ... '
PRINT *, ' '
!
! Carry out the time steps:
CALL TakeTimeSteps
!
! User-specified actions before program ends:
CALL AfterRun
!
! Delete the below-seed-level grids:
DO level = minlevel, maxlevel
  CALL ApplyOnLevel(level,MarkNodeInactive,dummy)
  CALL ApplyOnLevel(level,ReleaseInactiveFields,dummy)
  PRINT *,' Level ',level,' fields have been released.'
END DO
!
! Delete the forest of trees:
CALL KillForest
!
! Delete forest seed and below-seed levels:
CALL ApplyOnLevel(minlevel,MarkNodeToBeDeleted,dummy)
CALL ApplyOnLevel(minlevel,DeleteMarkedNode,dummy)
!
CALL DeallocPeriodicBCStorage
!
END SUBROUTINE BSAMSolver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION MarkNodeToBeDeleted(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok
IMPLICIT NONE
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
MarkNodeToBeDeleted = err_ok
!
info%tobedeleted = .TRUE.
!
END FUNCTION MarkNodeToBeDeleted
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION MarkNodeInactive(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok
IMPLICIT NONE
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
MarkNodeInactive = err_ok
!
info%activegrid = .FALSE.
!
END FUNCTION MarkNodeInactive
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION MarkNodeNonInitial(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok
IMPLICIT NONE
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
MarkNodeNonInitial = err_ok
!
info%initialgrid = .FALSE.
!
END FUNCTION MarkNodeNonInitial
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION RootInit(rootinfo,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok
USE Boundary, ONLY: PeriodicSetup
USE BSAMStorage, ONLY: AllocFields
IMPLICIT NONE
!
TYPE(nodeinfo):: rootinfo
TYPE(funcparam):: dummy
!
CHARACTER(LEN=12):: filename
INTEGER:: i, ierror, maux, mbc, n, nrvars
INTEGER, DIMENSION(1:maxdims):: mx
INTEGER, DIMENSION(1:2*maxdims):: mthbc
INTEGER, DIMENSION(1:maxdims,1:2):: mglobal
REAL(KIND=r8), DIMENSION(1:maxdims):: dx, xlower, xupper
!
NAMELIST/griddata/desiredfillratios, errflagopt, ibuffer, maxlevel, maux, mbc, &
                  mglobal, minimumgridpoints, minlevel, mthbc, mx, ndims, &
                  nrvars, qpo, qtolerance, xlower, xupper
!
PRINT *, 'Reading grid data for root-level grid.'
!
! Set default values !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
desiredfillratios = 8.0E-01_r8
errflagopt = errflagdefault
ibuffer = 0
maxlevel = 0
maux = 0
mbc = 1
mglobal = 1
minimumgridpoints = 2
minlevel = 0
mthbc = 10
mx = 1
ndims = 2
nrvars = 1
qpo = 0.0_r8
qtolerance = 1.0E-06_r8
xlower = 0.0_r8
xupper = 1.0_r8
!
! Read from namelist !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
OPEN(UNIT=75,FILE='griddata.dat',STATUS='OLD',ACTION='READ',IOSTAT=ierror)
IF(ierror/=0) THEN
  PRINT *,'Error opening input file griddata.dat. Program stop.'
  STOP
END IF
READ(75,NML=griddata)
CLOSE(75)
OPEN(UNIT=76,FILE='output.dat',STATUS='OLD',ACTION='WRITE',FORM='FORMATTED', &
     POSITION='APPEND')
WRITE(76,NML=griddata)
CLOSE(76)
!
! Check input variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
IF(ANY(desiredfillratios<0.0) .OR. ANY(desiredfillratios>1.0)) THEN
  PRINT *, 'BSAM 1.2 Error: Code only supports desiredfilratios between 0 and 1.'
  STOP
END IF
!
IF(ANY(ibuffer<0)) THEN
  PRINT *, 'BSAM 1.2 Error: Code only supports ibuffers>=0.'
  STOP
END IF
!
IF(maxlevel<0) THEN
  PRINT *, 'BSAM 1.2 Error: Code only supports maxlevel>=0.'
  STOP
END IF
!
IF(maux<0) THEN
  PRINT *, 'BSAM 1.2 Error: Code only supports maux>=0.'
  STOP
END IF
!
IF(mbc<1 .OR. mbc>2) THEN
  PRINT *, 'BSAM 1.2 Error: Code only supports mbc=1,2.'
  STOP
END IF
!
IF(minlevel>0) THEN
  PRINT *, 'BSAM 1.2 Error: Code only supports minlevel<=0.'
  STOP
END IF
!
IF(ndims<2 .OR. ndims>3) THEN
  PRINT *, 'BSAM 1.2 Error: Code only supports ndims=2,3.'
  STOP
END IF
!
DO i = 1, 2*ndims-1, 2
  IF(nrootgrids==1 .AND. &
    mthbc(i)==2 .AND. mthbc(i+1)/=2) THEN
    PRINT *, 'BSAM 1.2 Error: Incorrect periodic boundary conditions.'
    STOP 
  END IF
  IF(nrootgrids==1 .AND. &
    mthbc(i+1)==2 .AND. mthbc(i)/=2) THEN
    PRINT *, 'BSAM 1.2 Error: Incorrect periodic boundary conditions.'
    STOP
  END IF
END DO
!
DO i = 1, ndims
  IF(mx(i)<=0) THEN
    PRINT *, 'BSAM 1.2 Error: mx<=0 along dim', i, '.'
    STOP
  END IF
END DO
!
DO i = 1, ndims
  IF(MODULO(mx(i),2)/=0) THEN
    PRINT *, 'BSAM 1.2 Error: Initial grid must have dimensions which are a multiple'
    PRINT *, '       of the refinement ratio 2.'
    STOP
  END IF
END DO
!
DO i = 1, ndims
  IF(mglobal(i,2)-mglobal(i,1)+1/=mx(i)) THEN
    PRINT *, 'BSAM 1.2 Error: mglobal(i,2)-mglobal(i,1)+1/=mx(i), i=', i, '.'
    STOP
  END IF
END DO
!
IF(ANY(minimumgridpoints<1)) THEN
  PRINT *, 'BSAM 1.2 Error: Code only supports minimumgridpoints>=1.'
  STOP
END IF
!
IF(nrvars>maxnrvars .OR. nrvars<1) THEN
  PRINT *, 'BSAM 1.2 Error: Code only supports nrvars<=maxnrvars and nrvars>=1.'
  STOP
END IF
!
IF(ANY(qtolerance<=0.0)) THEN
  PRINT *, 'BSAM 1.2 Error: Code only supports qtolerance>0.0.'
  STOP
END IF
!
DO i = 1, ndims
  IF(xupper(i)<=xlower(i)) THEN
    PRINT *, 'BSAM 1.2 Error: Code only supports xupper > xlower.'
    STOP
  END IF
END DO
!
dx = 0.0_r8
dx(1:ndims) = (xupper(1:ndims)-xlower(1:ndims))/REAL(mx(1:ndims),KIND=r8)
!
DO i = 2, ndims
  IF(ABS(dx(1)-dx(i))>1.0E-10_r8) THEN
    PRINT *, 'BSAM 1.2 Error: dx(1)\=dx(i) along dim', i, '.'
    STOP
  END IF
END DO
!
mxmax = 1; mxmax(0,1:ndims) = mx(1:ndims)
DO i = rootlevel+1, maxlevel
  mxmax(i,1:ndims) = 2*mxmax(i-1,1:ndims)
END DO
!
! Copy data to global storage !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
rootinfo%ngrid = 0
rootinfo%level = rootlevel
rootinfo%maxlevel = maxlevel
!
rootinfo%tobedeleted = .FALSE.
rootinfo%initialgrid = .TRUE.
rootinfo%activegrid = .TRUE.
rootinfo%defective = .FALSE.
!
rootinfo%mx = 1
rootinfo%mx(1:ndims) = mx(1:ndims)
rootinfo%mglobal = 1
rootinfo%mglobal(1:ndims,1:2) = mglobal(1:ndims,1:2)
!
rootinfo%gridtime = 0.0_r8
!
rootinfo%xlower = 0.0_r8
rootinfo%xlower(1:ndims) = xlower(1:ndims)
rootinfo%xupper = 0.0_r8
rootinfo%xupper(1:ndims) = xupper(1:ndims)
rootinfo%dx = 0.0_r8
rootinfo%dx(1:ndims) = dx(1:ndims)
rootinfo%mbc = mbc
rootinfo%mthbc = 10
rootinfo%mthbc(1:2*ndims) = mthbc(1:2*ndims)
!
rootinfo%nrvars = nrvars
rootinfo%maux = maux
rootinfo%nroutvars = rootinfo%nrvars
!
rootinfo%mbounds = 1; rootinfo%mbounds(1:ndims,2) = rootinfo%mx(1:ndims)
!
! Allocate storage !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Set up the periodic offsets used in transferring periodic bc's:
CALL PeriodicSetup(rootinfo)
!
CALL AllocFields(rootinfo)
!
rootinfo%levellandscape = 0
!
! Finished initialization of root node info structure
RootInit = err_ok
!
PRINT *, 'Finished reading root-level grid data.'
!
END FUNCTION RootInit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION InitSeed(seedinfo,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok, GetChildInfo
USE BSAMStorage, ONLY: AllocFields
IMPLICIT NONE
!
TYPE(nodeinfo):: seedinfo
TYPE(funcparam) :: dummy
!
TYPE(nodeinfo), POINTER:: child   
INTEGER:: ierror, level, maux, mbc, nrvars
!
! Only data needed for multigrid need be copied from the child grid.  Parent
! (seedinfo) is born from child:
!
ierror = GetChildInfo(child)
!
nrvars = child%nrvars
mbc = child%mbc
level = child%level-1
!
seedinfo%nrvars = nrvars
seedinfo%mbc = mbc
seedinfo%maxlevel = maxlevel
seedinfo%level = level
!
PRINT *, 'level=', level, 'ndims=', ndims
!
IF(ANY(MODULO(child%mx(1:ndims),2)/=0)) THEN
  PRINT *, 'BSAM 1.2 Error in InitSeed: grid on level', child%level, 'will not coarsen.'
  STOP
END IF
!
seedinfo%mx(1:ndims) = child%mx(1:ndims)/2
seedinfo%dx(1:ndims) = child%dx(1:ndims)*REAL(2,KIND=r8)
seedinfo%mglobal(1:ndims,1) = 1
seedinfo%mglobal(1:ndims,2) = child%mglobal(1:ndims,2)/2
child%mbounds(1:ndims,1) = 1; child%mbounds(1:ndims,2) = seedinfo%mx(1:ndims)
!
seedinfo%mbounds(1:ndims,1) = 1; seedinfo%mbounds(1:ndims,2) = seedinfo%mx(1:ndims)
!
IF(ANY(seedinfo%mx(1:ndims)/=seedinfo%mglobal(1:ndims,2))) THEN
  PRINT *, 'BSAM 1.2 Error in InitSeed: on level-1', seedinfo%level
  PRINT *, 'mx(:)/=mglobal(:,2)'
  STOP
END IF
!
seedinfo%ngrid = 0
!
seedinfo%tobedeleted = .FALSE.; seedinfo%level = child%level-1
seedinfo%initialgrid = .FALSE.
seedinfo%activegrid = .TRUE.
seedinfo%defective = .FALSE.
!
seedinfo%gridtime = child%gridtime
!
seedinfo%xlower = child%xlower; seedinfo%xupper = child%xupper
seedinfo%mthbc = child%mthbc
!
seedinfo%maux = child%maux
seedinfo%nroutvars = seedinfo%nrvars
!
CALL AllocFields(seedinfo)
!
! Finished initialization of seed node info structure:
InitSeed = err_ok
!
END FUNCTION InitSeed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION SetAuxFields(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok
USE Problem, ONLY: SetAux
IMPLICIT NONE
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
SetAuxFields = err_ok
IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
IF(info%maux>0) CALL SetAux(info)
!
END FUNCTION SetAuxFields
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!
INTEGER FUNCTION SetSrcFields(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok
USE Problem, ONLY: SetSrc
IMPLICIT NONE
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
SetSrcFields = err_ok
IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
CALL SetSrc(info)
!
END FUNCTION SetSrcFields
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION Initialq(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok
IMPLICIT NONE
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
Initialq = err_ok
IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
CALL QInit(info)
!
END FUNCTION Initialq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE QInit(info)
USE NodeInfoDef
USE PROBLEM, ONLY: QInit2D, QInit3D
IMPLICIT NONE
!
TYPE(nodeinfo):: info
!
INTEGER:: nrvars
INTEGER, DIMENSION(1:maxdims):: mx
REAL(KIND=r8):: h
REAL(KIND=r8), DIMENSION(1:maxdims):: xlower
!
nrvars = info%nrvars
mx(1:ndims) = info%mx(1:ndims)
h = info%dx(1)
xlower(1:ndims) = info%xlower(1:ndims)
!
SELECT CASE(ndims)
  CASE(2)
    CALL QInit2D(info%q(1:mx(1),1:mx(2),1      ,1:nrvars),mx(1:2), &
                 nrvars,h,xlower(1:2))
  CASE(3)
    CALL QInit3D(info%q(1:mx(1),1:mx(2),1:mx(3),1:nrvars),mx(1:3), &
                 nrvars,h,xlower(1:3))
END SELECT
!
END SUBROUTINE QInit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION CopyQToQold(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok
IMPLICIT NONE
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
CopyQToQold = err_ok
IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
info%qold = info%q
!
END FUNCTION CopyQToQold
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE TakeTimeSteps
USE NodeInfoDef
USE TreeOps, ONLY: ApplyOnLevel, DeleteMarkedNode
USE BSAMInputOutput, ONLY: WriteQ, WriteUniformMeshQ
USE AFASRoutines, ONLY: FillDown, MultigridIterations
IMPLICIT NONE
!
! currenttime is the time of the rootlevel grid.
! restarttime is the time at restart.
! finaltime is the calculated final time.
!
TYPE(funcparam):: dummy
LOGICAL:: firstamr
INTEGER:: firstframe, it, itperprint, lastframe, level, n
REAL(KIND=r8):: starttime, dtsave
!
dtsave = dt
!
IF(restart) THEN
  IF(restartframe>=outframes) THEN
    PRINT *, 'BSAM 1.2 Error in rundata.dat: outframes=', outframes, ' restartframe=', &
             restartframe
    PRINT *, 'outframes must be greater than restartframe.'
    STOP
  END IF
!
  currenttime = restarttime
  firstframe = restartframe+1; lastframe = outframes
!
  PRINT 1002, restartframe, outframes
  1002 FORMAT('================================================'/ &
              ' BSAM 1.2: Restart frame ', i4,' out of ', i4,' requested frames'/ &
              '================================================')
ELSE
  currenttime = 0.0_r8
  firstframe = 1; lastframe = outframes
  IF(syncelliptic) THEN
    dt = 0.0_r8
  END IF
END IF
!
itperprint = timeiterations/outframes
!
starttime = currenttime
finaltime = starttime &
          + REAL((lastframe-firstframe+1)*itperprint,KIND=r8)*dtsave
!
firstamr = .TRUE.
!
printloop: DO n = firstframe, lastframe
!
  it = 1
  timesteploop: DO
!
! For testing purposes only:
!
! Print out grid information for memory checking:
!    DO level = rootlevel+1, maxlevel
!      PRINT *, 'Grid information, level =', level
!      gridnumber = 0
!      CALL ApplyOnLevel(level,PrintGridInfo,dummy)
!      READ(*,*)
!    END DO
!
! Mark all above-root-level grids from the old mesh inactive:
    DO level = rootlevel+1, maxlevel
      CALL ApplyOnLevel(level,MarkNodeInactive,dummy)
    END DO
!
! Make a new mesh and move the solution to the new mesh:
    finestlevel = rootlevel
    defectivegridlevel(0:maxlevel) = .FALSE.
    amrrestarts = 0; meshbuildcomplete = .TRUE.
!
    CALL AMR(rootlevel)
!
! After exit from AMR delete all inactive grids (the previous mesh):
    meshbuildcomplete = .TRUE.
    DO level = maxlevel, rootlevel+1, -1
      CALL ApplyOnLevel(level,ReleaseInactiveFields,dummy)
      CALL ApplyOnLevel(level,DeleteMarkedNode,dummy)
    END DO
!
! This ensures that the below rootlevel values of qold are filled after a
! clean start or restart:
    IF(firstamr) THEN
      CALL FillDown(solutionfield)
      firstamr = .FALSE.
      CALL ApplyOnLevel(rootlevel,MarkNodeNonInitial,dummy)
    END IF
!
! Save a copy of data at the last time step:
    DO level = minlevel, finestlevel
      CALL ApplyOnLevel(level,CopyQToQold,dummy)
    END DO
!
! If initial grid, write data:
    IF(outputinitialdata .AND. (.NOT. syncelliptic)) THEN
!
      outputinitialdata = .FALSE.
      PRINT 799
      799 FORMAT('==============================='/ &
                 ' BSAM 1.2: Writing initial data'/ &
                 '===============================')
      IF(getafterstepstats) THEN
        integralresult = 0.0_r8
        totalmeshsize = 0
        DO level = finestlevel, rootlevel, -1
          CALL ApplyOnLevel(level,AfterStepStatistics,dummy)
          CALL ApplyOnLevel(level,GetMeshSize,dummy)
        END DO
        OPEN(UNIT=65,FILE='OUT/stats.dat',STATUS='UNKNOWN',ACTION='WRITE', &
             FORM='FORMATTED',POSITION='APPEND')
        WRITE(65,'(3(F25.12),1x,I8)') currenttime, integralresult(1:2), &
                                      totalmeshsize
        CLOSE(65)
      END IF
      totalmeshsize = 0
      DO level = finestlevel, rootlevel, -1
        CALL ApplyOnLevel(level,GetMeshSize,dummy)
      END DO
      PRINT *, 'Initial composite mesh size =', totalmeshsize
      CALL WriteQ(0,currenttime)
      IF(outputuniformmesh) CALL WriteUniformMeshQ(0,currenttime)
!
    END IF
!
    IF(syncelliptic) THEN
      PRINT *, ' '
      PRINT *, ' '
      PRINT *, 'BSAM 1.2: Synchronizing elliptic fields at start.' 
      PRINT *, ' '
    ELSE
      PRINT *, ' '
      PRINT *, ' '
      PRINT *, 'BSAM 1.2: Advancing to time = ', currenttime+dt, ' out of ', finaltime 
      PRINT *, ' '
    END IF
!
! For testing purposes only:
!
! Print out grid information for memory checking:
!    DO level = rootlevel+1, maxlevel
!      PRINT *, 'Grid information, level =', level
!      gridnumber = 0
!      CALL ApplyOnLevel(level,PrintGridInfo,dummy)
!      READ(*,*)
!    END DO
!
! Perform Multigrid on the multilevel mesh:
    CALL MultigridIterations
!
    currenttime = starttime+REAL((n-firstframe)*itperprint+it,KIND=r8)*dt
!
! Calculate after step statistics:
    IF(getafterstepstats .AND. (.NOT. syncelliptic)) THEN
!
      integralresult = 0.0_r8
      totalmeshsize = 0
      DO level = finestlevel, rootlevel, -1
        CALL ApplyOnLevel(level,AfterStepStatistics,dummy)
        CALL ApplyOnLevel(level,GetMeshSize,dummy)
      END DO
      OPEN(UNIT=65,FILE='OUT/stats.dat',STATUS='UNKNOWN',ACTION='WRITE', &
           FORM='FORMATTED',POSITION='APPEND')
      WRITE(65,'(3(F25.12),1x,I8)') currenttime, integralresult(1:2), &
                                    totalmeshsize
      CLOSE(65)
    END IF
!
    IF(syncelliptic) THEN
      syncelliptic = .false.
      dt = dtsave
      totalmeshsize = 0
      DO level = finestlevel, rootlevel, -1
        CALL ApplyOnLevel(level,GetMeshSize,dummy)
      END DO
      PRINT *, 'BSAM 1.2: Synchronization composite mesh size =', totalmeshsize
    ELSE
      it = it+1
    END IF
!
! Set current time on the grid:
    DO level = finestlevel, rootlevel, -1
      CALL ApplyOnLevel(level,SetCurrentTime,dummy)
    END DO
!
    IF(it>itperprint) EXIT timesteploop
!
  END DO timesteploop
!
  totalmeshsize = 0
  DO level = finestlevel, rootlevel, -1
    CALL ApplyOnLevel(level,GetMeshSize,dummy)
  END DO
  PRINT *, 'Composite mesh size =', totalmeshsize
!
  PRINT 1001, n, outframes, currenttime
  1001 FORMAT( &
    '===================================================================='/ &
    ' Writing frame ',I4,' out of ',I4,' requested frames at t=',E11.4,   / &
    '====================================================================')
!
  CALL WriteQ(n,currenttime)
  IF(outputuniformmesh) CALL WriteUniformMeshQ(n,currenttime)
!
END DO printloop
!
END SUBROUTINE TakeTimeSteps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION PrintGridInfo(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok
IMPLICIT NONE
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
PrintGridInfo = err_ok
!
gridnumber = gridnumber+1
PRINT *, ' '
PRINT *, 'Grid number =', gridnumber, 'level=', info%level
PRINT *, 'mx =', info%mx
PRINT *, 'mb =', info%mbounds
PRINT *, 'ngrid', info%ngrid
PRINT *, 'activegrid', info%activegrid
PRINT *, 'initialgrid', info%initialgrid
PRINT *, 'tobedeleted', info%tobedeleted
PRINT *, 'fieldsallocated', info%fieldsallocated
PRINT *, ' '
!
END FUNCTION PrintGridInfo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION GetMeshSize(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok
IMPLICIT NONE
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
INTEGER, DIMENSION(1:maxdims):: mx, cmx
!
GetMeshSize = err_ok
IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
mx(1:ndims) = info%mx(1:ndims)
cmx(1:ndims) = mx(1:ndims)/2
!
IF(info%level==0) THEN
  totalmeshsize = totalmeshsize+PRODUCT(mx(1:ndims))
ELSE
  totalmeshsize = totalmeshsize+PRODUCT(mx(1:ndims))-PRODUCT(cmx(1:ndims))
END IF
!
END FUNCTION GetMeshSize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION ReleaseInactiveFields(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok
USE BSAMStorage, ONLY: DeAllocFields
IMPLICIT NONE
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
ReleaseInactiveFields = err_ok
!
IF(.NOT. info%activegrid) CALL DeAllocFields(info)
!
END FUNCTION ReleaseInactiveFields
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION ReleaseActiveFields(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok
USE BSAMStorage, ONLY: DeAllocFields
IMPLICIT NONE
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
ReleaseActiveFields = err_ok
!
IF(info%activegrid) CALL DeAllocFields(info)
!
END FUNCTION ReleaseActiveFields
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION SetCurrentTime(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok
IMPLICIT NONE
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
SetCurrentTime = err_ok
!
info%gridtime = currenttime
!
END FUNCTION SetCurrentTime
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION AfterStepStatistics(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok, GetParentInfo
USE PROBLEM, ONLY: AfterStep2D, AfterStep3D
IMPLICIT NONE
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
TYPE(nodeinfo), POINTER:: parent
INTEGER:: ierror, level, nrvars
INTEGER, DIMENSION(1:maxdims):: mx
INTEGER, DIMENSION(1:maxdims,1:2):: mb
REAL(KIND=r8):: h
REAL(KIND=r8), DIMENSION(1:maxdims):: xlower
!
AfterStepStatistics = err_ok
IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
mx(1:ndims) = info%mx(1:ndims)
nrvars = info%nrvars
h = info%dx(1)
xlower(1:ndims) = info%xlower(1:ndims)
level = info%level
mb = 1; mb(1:ndims,1:2) = info%mbounds(1:ndims,1:2)
!
ierror = GetParentInfo(parent)
!
SELECT CASE(ndims)
  CASE(2)
    CALL AfterStep2D(info%q(0:mx(1)+1,0:mx(2)+1,1,1:nrvars), &
                     parent%q(mb(1,1)-1:mb(1,2)+1, &
                              mb(2,1)-1:mb(2,2)+1,1,1:nrvars), &
                     mx(1:2),nrvars,h,xlower(1:2),level)
  CASE(3)
    CALL AfterStep3D(info%q(0:mx(1)+1,0:mx(2)+1,0:mx(3)+1,1:nrvars), &
                     parent%q(mb(1,1)-1:mb(1,2)+1, &
                              mb(2,1)-1:mb(2,2)+1, &
                              mb(3,1)-1:mb(3,2)+1,1:nrvars), &
                     mx(1:3),nrvars,h,xlower(1:3),level)
  CASE DEFAULT
    PRINT *, 'BSAM 1.2: AfterStepStatistics: Only ndims = 2,3 are supported.'
    STOP
END SELECT
!
END FUNCTION AfterStepStatistics
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
RECURSIVE SUBROUTINE AMR(level,estimateswitch)
USE NodeInfoDef
USE TreeOps, ONLY: ApplyOnLevel, DeleteMarkedNode
USE Boundary, ONLY: SetGhost
IMPLICIT NONE
!
INTEGER, INTENT(IN):: level
INTEGER, OPTIONAL, INTENT(IN):: estimateswitch
!
INTEGER:: lvl, noestimates
TYPE(funcparam):: dummy
!
meshbuildcomplete = .TRUE.
!
! Fill in the ghost points:
CALL SetGhost(level,0)
!
! Tag cells for refinement:
IF(level < maxlevel) CALL EstimateLevelErrors(level)
!
! Create new subgrids if necessary:
CALL GridAdapt(level)
!
IF(level < finestlevel) CALL AMR(level+1)
!
END SUBROUTINE AMR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE EstimateLevelErrors(level)
USE NodeInfoDef
USE TreeOps, ONLY: ApplyOnLevel
USE AFASRoutines, ONLY: RestrictSolution, RelativeTruncationError
USE Boundary, ONLY: SetGhost, GetCoarseGhostPoints
IMPLICIT NONE
!
INTEGER, INTENT(IN):: level
!
TYPE(funcparam):: dummy
!
! 1) Calculate the relative truncation error if needed:
IF(errflagopt(level) == errflagdefault) THEN
  CALL ApplyOnLevel(level,RestrictSolution,dummy)
!
  CALL SetGhost(level-1,0)
!
  CALL ApplyOnLevel(level,GetCoarseGhostPoints,dummy)
!
  CALL ApplyOnLevel(level,CopyQToQold,dummy)
!
  CALL ApplyOnLevel(level-1,CopyQToQold,dummy)
!
  CALL ApplyOnLevel(level,RelativeTruncationError,dummy)
END IF
!
! 2) Refine based on the size of the error estimator.  Make a linked list of
!    tagged cells:
ALLOCATE(zerothtaggedcell)
NULLIFY(zerothtaggedcell%prevcell)
!
lasttaggedcell => zerothtaggedcell
ntaggedcells = 0
!
CALL ApplyOnLevel(level,EstimateError,dummy)
!
! 3) Inflate tagged regions and destroy the list:
IF(ntaggedcells > 0) THEN
  CALL InflateEdgeTags(level)
  CALL DeleteTaggedCellsList
END IF
!
NULLIFY(lasttaggedcell)
NULLIFY(currenttaggedcell)
DEALLOCATE(zerothtaggedcell)
!
END SUBROUTINE EstimateLevelErrors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION EstimateError(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok
IMPLICIT NONE
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
EstimateError = err_ok
IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
CALL ErrFlag(info)
!
END FUNCTION EstimateError
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE ErrFlag(info)
USE NodeInfoDef
USE Problem, ONLY: SetErrFlagsUser2D, SetErrFlagsUser3D
IMPLICIT NONE
!
TYPE(nodeinfo):: info
!
INTEGER:: level, nrvars
INTEGER, DIMENSION(1:maxdims):: mx, cmx
INTEGER, DIMENSION(1:maxdims,1:2):: mglobal
REAL(KIND=r8):: h
REAL(KIND=r8), DIMENSION(1:maxdims):: xlower
!
nrvars = info%nrvars
level = info%level
h = info%dx(1)
xlower(1:ndims) = info%xlower(1:ndims)
mx = 1; mx(1:ndims) = info%mx(1:ndims)
cmx = 1; cmx(1:ndims) = mx(1:ndims)/2
mglobal = 1; mglobal(1:ndims,1:2) = info%mglobal(1:ndims,1:2)
!
SELECT CASE(errflagopt(level))
CASE(errflagdefault)
!
! Default error tagging based on the relative truncation error:
  SELECT CASE(ndims)
  CASE(2)
    CALL SetErrFlags2D(info%qrte(1:cmx(1),1:cmx(2),1,1:nrvars), &
                       info%errorflags(1:mx(1),1:mx(2),1), &
                       mx(1:2),cmx(1:2),nrvars,h,level)
  CASE(3)
  CALL SetErrFlags3D(info%qrte(1:cmx(1),1:cmx(2),1:cmx(3),1:nrvars), &
                     info%errorflags(1:mx(1),1:mx(2),1:mx(3)), &
                     mx(1:3),cmx(1:3),nrvars,h,level)
!
  CASE DEFAULT
    PRINT *, 'ErrFlag: only ndims=2d,3d supported'
    STOP
  END SELECT
!
CASE(errflaguser)
!
! User chooses the error tagging proceedure:
  SELECT CASE(ndims)
  CASE(2)
    CALL SetErrFlagsUser2D(info%qrte(1:cmx(1)  ,1:cmx(2)  ,1,1:nrvars), &
                           info%q   (0: mx(1)+1,0: mx(2)+1,1,1:nrvars), &
                           info%errorflags(1:mx(1),1:mx(2),1), &
                           mx(1:2),cmx(1:2),nrvars,h,xlower(1:2),level)
  CASE(3)
    CALL SetErrFlagsUser3D(info%qrte(1:cmx(1)  ,1:cmx(2)  ,1:cmx(3)  ,1:nrvars), &
                           info%q   (0: mx(1)+1,0: mx(2)+1,0: mx(3)+1,1:nrvars), &
                           info%errorflags(1:mx(1),1:mx(2),1:mx(3)), &
                           mx(1:3),cmx(1:3),nrvars,h,xlower(1:3),level)
!
  CASE DEFAULT
    PRINT *, 'ErrFlag: only ndims=2d,3d supported'
    STOP
  END SELECT
!
CASE DEFAULT
  PRINT *, 'ErrFlag: No error flagging algorithm selected'
  STOP
END SELECT
!
! Add the tagged cells to a linked list. Buffering requires a global approach,
! since buffer layers might go into neighboring grids:
IF(ibuffer(level) > 0) THEN
  CALL BufferAndList(info%errorflags(1:mx(1),1:mx(2),1:mx(3)), &
                     mglobal(1:3,1:2),mx(1:3),level)
END IF
!
END SUBROUTINE ErrFlag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE SetErrFlags2D(qrte,errorflags,mx,cmx,nrvars,h,level)
USE NodeInfoDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(1:,1:,1:), INTENT(IN):: qrte
INTEGER, DIMENSION(1:,1:), INTENT(OUT):: errorflags
INTEGER, DIMENSION(1:2), INTENT(IN):: mx
INTEGER, DIMENSION(1:2), INTENT(IN):: cmx
INTEGER, INTENT(IN):: nrvars
REAL(KIND=r8), INTENT(IN):: h
INTEGER, INTENT(IN):: level
!
INTEGER:: i, j
REAL(KIND=r8):: tol
!
tol = qtolerance(level)/h/h
!
errorflags(1:mx(1),1:mx(2)) = 0
!
DO i = 1, cmx(1)
  DO j = 1, cmx(2)
    IF(MAXVAL(ABS(qrte(i,j,1:nrvars)))>tol) THEN
      errorflags(2*i  ,2*j  ) = 1
      errorflags(2*i-1,2*j  ) = 1
      errorflags(2*i  ,2*j-1) = 1
      errorflags(2*i-1,2*j-1) = 1
    END IF
  END DO
END DO
!
END SUBROUTINE SetErrFlags2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE SetErrFlags3D(qrte,errorflags,mx,cmx,nrvars,h,level)
USE NodeInfoDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(1:,1:,1:,1:), INTENT(IN):: qrte
INTEGER, DIMENSION(1:,1:,1:), INTENT(OUT):: errorflags
INTEGER, DIMENSION(1:3), INTENT(IN):: mx
INTEGER, DIMENSION(1:3), INTENT(IN):: cmx
INTEGER, INTENT(IN):: nrvars
REAL(KIND=r8), INTENT(IN):: h
INTEGER, INTENT(IN):: level
!
INTEGER:: i, j, k
REAL(KIND=r8):: tol
!
tol = qtolerance(level)/h/h/h
!
errorflags(1:mx(1),1:mx(2),1:mx(3)) = 0
!
DO i = 1, cmx(1)
  DO j = 1, cmx(2)
    DO k = 1, cmx(3)
      IF(MAXVAL(ABS(qrte(i,j,k,1:nrvars)))>tol) THEN
        errorflags(2*i  ,2*j  ,2*k  ) = 1
        errorflags(2*i-1,2*j  ,2*k  ) = 1
        errorflags(2*i  ,2*j-1,2*k  ) = 1
        errorflags(2*i-1,2*j-1,2*k  ) = 1
        errorflags(2*i  ,2*j  ,2*k-1) = 1
        errorflags(2*i-1,2*j  ,2*k-1) = 1
        errorflags(2*i  ,2*j-1,2*k-1) = 1
        errorflags(2*i-1,2*j-1,2*k-1) = 1
      END IF
    END DO
  END DO
END DO
!
END SUBROUTINE SetErrFlags3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE BufferAndList(errorflags,mglobal,mx,level)
USE NodeInfoDef
IMPLICIT NONE
!
INTEGER, DIMENSION(1:,1:,1:), INTENT(IN OUT):: errorflags
INTEGER, DIMENSION(1:3,1:2), INTENT(IN):: mglobal
INTEGER, DIMENSION(1:3), INTENT(IN):: mx
INTEGER, INTENT(IN):: level
!
INTEGER:: i, j, k
INTEGER, DIMENSION(1:maxdims):: index
INTEGER, DIMENSION(1:maxdims,1:2):: mtg
INTEGER, DIMENSION(1:mx(1),1:mx(2),1:mx(3)):: errorflagstmp
!
errorflagstmp = errorflags
mtg = 1
!
DO k = 1, mx(3)
  index(3) = k
  DO j = 1, mx(2)
    index(2) = j
    DO i = 1, mx(1)
      index(1) = i
      IF(errorflagstmp(i,j,k)==1) THEN
        mtg(1:ndims,1) = index(1:ndims)-ibuffer(level)
        mtg(1:ndims,2) = index(1:ndims)+ibuffer(level)
        IF(ANY(mtg(1:3,1)<1) .OR. ANY(mtg(1:3,2)>mx(1:3))) THEN
!
! If the buffer area overlaps with any edge of the patch, then record the
! global coordinates of the cell:
          Call AddTaggedCellToList(index(1:maxdims)+mglobal(1:maxdims,1)-1)
        ELSE
!
! If the buffer area lies within the patch, apply the buffer:
          errorflags(mtg(1,1):mtg(1,2),mtg(2,1):mtg(2,2),mtg(3,1):mtg(3,2)) = 1
        END IF
      END IF
    END DO
  END DO
END DO
!
END SUBROUTINE BufferAndList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE AddTaggedCellToList(globalindex)
USE NodeInfoDef
IMPLICIT NONE
!
INTEGER, DIMENSION(1:maxdims), INTENT(IN):: globalindex
!
ntaggedcells = ntaggedcells+1
ALLOCATE(currenttaggedcell)
currenttaggedcell%id = ntaggedcells
currenttaggedcell%coordinate(1:maxdims) = 1
currenttaggedcell%coordinate(1:ndims) = globalindex(1:ndims)
currenttaggedcell%prevcell => lasttaggedcell
lasttaggedcell => currenttaggedcell
!
END SUBROUTINE AddTaggedCellToList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE InflateEdgeTags(level)
USE NodeInfoDef
USE TreeOps, ONLY: ApplyOnLevel
USE Boundary, ONLY: GetPeriodicOffsets
IMPLICIT NONE
!
INTEGER, INTENT(IN):: level
!
TYPE(funcparam):: dummy
LOGICAL:: periodicbuffer
INTEGER:: offset, polarity
INTEGER, DIMENSION(1:maxdims):: coordinatesave
!
coordinatesave = 1
currenttaggedcell => lasttaggedcell
searchloop: DO
  IF(.NOT. ASSOCIATED(currenttaggedcell%prevcell)) EXIT searchloop
!
! Ordinary buffering of an edge tag:
  dummy%iswitch = ibuffer(level)
  CALL ApplyOnLevel(level,BufferTaggedCells,dummy)
!
! Buffering of periodic edge tags. Check to see if the buffer area cuts across 
! a periodic boundary.  If so, add offset and apply buffer:
  coordinatesave(1:ndims) = currenttaggedcell%coordinate(1:ndims)
  IF(periodicboundaryconditions) THEN
    CALL GetPeriodicOffsets(level)
    DO polarity = -1, 1, 2
      DO offset = 1, nperiodicoffsets
        currenttaggedcell%coordinate(1:ndims) &
          = coordinatesave(1:ndims)+polarity*poffset(1:ndims,offset)
        dummy%iswitch = ibuffer(level)
        CALL ApplyOnLevel(level,BufferTaggedCells,dummy)
      END DO
    END DO
  END IF
  currenttaggedcell%coordinate(1:ndims) = coordinatesave(1:ndims)
!
  currenttaggedcell => currenttaggedcell%prevcell
END DO searchloop
!
END SUBROUTINE InflateEdgeTags
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION BufferTaggedCells(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok
IMPLICIT NONE
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
INTEGER:: ibuff, n
INTEGER, DIMENSION(1:maxdims,1:2):: mglobal, mglobaltag, mlocal, moverlap
!
BufferTaggedCells = err_ok
IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
ibuff = dummy%iswitch
mglobal = 1; mglobal(1:ndims,1:2) = info%mglobal(1:ndims,1:2)
mglobaltag = 1 
mglobaltag(1:ndims,1) = currenttaggedcell%coordinate(1:ndims)-ibuff
mglobaltag(1:ndims,2) = currenttaggedcell%coordinate(1:ndims)+ibuff
!
! 1. Find overlap region in global index space:
moverlap = 1
moverlap(1:ndims,1) = MAX(mglobaltag(1:ndims,1),mglobal(1:ndims,1))
moverlap(1:ndims,2) = MIN(mglobaltag(1:ndims,2),mglobal(1:ndims,2))
!
! 2. Check for nonempty intersection:
IF(ANY(moverlap(:,2)-moverlap(:,1)<0)) RETURN
!
! 3. Transform common index space to grid index spaces:
mlocal = 1
DO n = 1, ndims
  mlocal(n,1:2) = moverlap(n,1:2)-mglobal(n,1)+1
END DO
!
info%errorflags(mlocal(1,1):mlocal(1,2), &
                mlocal(2,1):mlocal(2,2), &
                mlocal(3,1):mlocal(3,2)) = 1
!
END FUNCTION BufferTaggedCells
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE DeleteTaggedCellsList
USE NodeInfoDef
IMPLICIT NONE
!
searchloop: DO
  currenttaggedcell => lasttaggedcell%prevcell
  DEALLOCATE(lasttaggedcell)
  lasttaggedcell => currenttaggedcell
  IF(.NOT. ASSOCIATED(lasttaggedcell%prevcell)) EXIT searchloop
END DO searchloop
!
END SUBROUTINE DeleteTaggedCellsList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE GridAdapt(level)
USE NodeInfoDef
USE TreeOps, ONLY: ApplyOnLevel, ApplyOnLevelPairs, DeleteMarkedNode
IMPLICIT NONE
!
INTEGER, INTENT(IN):: level
!
! Generate new subgrids of level:
!
TYPE(funcparam):: dummy
INTEGER:: ilevel
!
! Generate new grids in accordance with error flags:
CALL ApplyOnLevel(level,RefineGrid,dummy)
!
! Transfer field values from previous grids on this level to the newly created
! grids:
CALL ApplyOnLevelPairs(level+1,TransferValues,dummy)
!
END SUBROUTINE GridAdapt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION RefineGrid(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok
IMPLICIT NONE
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
INTEGER, DIMENSION(1:maxdims,1:2):: mbounds
!
RefineGrid = err_ok
IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
! Generate new subgrids of the current grid:
mbounds(1:ndims,1) = 1; mbounds(1:ndims,2) = info%mx(1:ndims)
CALL NewSubGrids(info,mbounds)
!
END FUNCTION RefineGrid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE NewSubGrids(info,mbounds)
USE NodeInfoDef
IMPLICIT NONE
!
! Modified implementation of Berger-Rigoutsos algorithm (IEEE Trans. Systems, 
! Man. & Cyber., 21(5):1278-1286, 1991):
!
TYPE(nodeinfo):: info
INTEGER, DIMENSION(1:maxdims,1:2), INTENT(IN OUT):: mbounds
!
LOGICAL havesplit, cansplitgrid
LOGICAL, DIMENSION(1:maxsubgrids):: cansplit
INTEGER, PARAMETER:: maxsplitpasses = 15
INTEGER:: del, dist0, dist1, i, i1, i2, ierror, igrid, inflect, level, maxm, &
          mgp, minm, n, nn, ngrid, npass
INTEGER, DIMENSION(1:maxdims):: isplit, mx
INTEGER, DIMENSION(1:maxdims,1:2,1:maxsubgrids):: msubbounds
INTEGER, DIMENSION(:,:), ALLOCATABLE:: signature ,ddsignature
REAL(KIND=r8):: fillratio, desfillratio
!
mx = info%mx
level = info%level; info%nsubgrids = 0
desfillratio = desiredfillratios(level)
!
mgp = minimumgridpoints(level)
IF(mgp<4) THEN
  PRINT *,'BSAM 1.2: Error on level', level, 'minimumgridpoints cannot be less than 4.'
  STOP
END IF
IF(MODULO(mgp,2)==1) THEN
  PRINT *,'BSAM 1.2: Error on level', level, 'minimumgridpoints must be even.'
  STOP
END IF
!
! Compute fill ratio for this grid:
fillratio = GridFlagRatio(info%errorflags(1:mx(1),1:mx(2),1:mx(3)), &
                          mbounds)
!
! Don't generate a new grid if we don't have flagged points:
IF(fillratio < 1.0E-10_r8) RETURN
!
! Allocate space for signatures:
maxm = MAXVAL(mx(1:ndims))
ALLOCATE(signature(1:maxm,1:ndims),ddsignature(1:maxm,1:ndims),STAT=ierror)
IF(ierror /= 0) THEN
  PRINT *,'BSAM 1.2: Error allocating signatures arrays in NewSubGrids'
  STOP
END IF
signature=0; ddsignature=0
!
! Initialize list of subgrids:
ngrid = 1; cansplit(:) = .TRUE.
msubbounds(1:ndims,1:2,ngrid) = mbounds(1:ndims,1:2)
!
! Loop until no better grid splitting can be found:
igrid = 1
DO WHILE (ngrid<maxsubgrids .AND. igrid<=ngrid)
  npass = 0
  DO WHILE (cansplit(igrid) .AND. npass<maxsplitpasses)
    npass = npass+1
    signature = GetSignatures(info%errorflags(1:mx(1),1:mx(2),1:mx(3)), &
                              msubbounds(:,:,igrid),maxm)
!
! Trim unflagged points on the edges of this grid:
    DO n = 1, ndims
      i1 = msubbounds(n,1,igrid); i2=msubbounds(n,2,igrid)     
!
      DO WHILE(signature(i1,n)==0 .AND. i1<msubbounds(n,2,igrid) .AND. &
               i2-i1+1>mgp)
        i1 = i1+1
      END DO
      DO WHILE(signature(i2,n)==0 .AND. i2>msubbounds(n,1,igrid) .AND. &
               i2-i1+1>mgp)
        i2 = i2-1
      END DO
!
! We only make grid adjustments in increments of 2 (coarse) grid points:
      msubbounds(n,1,igrid) = i1-MODULO(i1+1,2)
      msubbounds(n,2,igrid) = i2+MODULO(i2  ,2)  
    END DO
!
    fillratio = GridFlagRatio(info%errorflags(1:mx(1),1:mx(2),1:mx(3)), &
                              msubbounds(:,:,igrid))
    minm = MINVAL(msubbounds(1:ndims,2,igrid)-msubbounds(1:ndims,1,igrid))+1
    IF(fillratio<desfillratio) THEN
      signature = GetSignatures(info%errorflags(1:mx(1),1:mx(2),1:mx(3)), &
                                msubbounds(:,:,igrid),maxm)
!
! Look for holes along which to split grid:
      isplit = 0; havesplit = .FALSE.
      DO n = 1, ndims
        i1 = msubbounds(n,1,igrid); i2 = msubbounds(n,2,igrid)
        DO i = i1+mgp-1, i2-mgp, 2
!
! i is the second (terminal) index of the first half of the split grid.
          IF(signature(i,n)==0 .AND. MIN(i-i1+1,i2-i)>=mgp) THEN
            isplit(n) = i
            havesplit = .TRUE.
            EXIT
          END IF
        END DO
        IF(havesplit) EXIT
      END DO
!
      IF(.NOT. havesplit) THEN
!
! No split along a hole. Try split along inflection point:
        DO n = 1, ndims
          i1 = msubbounds(n,1,igrid); i2 = msubbounds(n,2,igrid)
          DO i = i1+1, i2-1
            ddsignature(i,n) = signature(i-1,n)-2*signature(i,n) &
                             + signature(i+1,n)
          END DO
        END DO
!
        inflect = 0; dist0 = 0
        DO n = 1, ndims
          i1 = msubbounds(n,1,igrid); i2 = msubbounds(n,2,igrid)
          DO i = i1+mgp, i2-mgp+1, 2
!
! Here i is the first index of the second half of the split grid.
            del = ABS(ddsignature(i,n)-ddsignature(i-1,n))
            IF(del>inflect) THEN
              inflect = del; isplit = 0; isplit(n) = i-1; havesplit = .TRUE.
              dist0 = MIN(i-i1,i2-i+1)
            ELSE IF(del==inflect .AND. inflect>0) THEN
              dist1 = MIN(i-i1,i2-i+1)
              IF(dist1>dist0) THEN
                isplit = 0; isplit(n) = i-1; havesplit = .TRUE.
                dist0 = dist1
              END IF
            END IF
          END DO
        END DO
      END IF
!
      IF(havesplit) THEN
!
! Split the grid along a determined line:
        DO n = 1, ndims
          IF(isplit(n)>0 .AND. &
             MIN(msubbounds(n,2,igrid)-isplit(n), &
                 isplit(n)-msubbounds(n,1,igrid)+1)>=mgp) THEN
!
! Add a new subgrid to the end of the grid list:
            ngrid = ngrid+1
            cansplit(ngrid) = .TRUE.
            msubbounds(1:ndims,1:2,ngrid) = msubbounds(1:ndims,1:2,igrid)
            msubbounds(n,1,ngrid) = isplit(n)+1
!
! Replace current grid with a subgrid:
            msubbounds(n,2,igrid) = isplit(n)
            EXIT
          END IF
        END DO
      ELSE
!
! Mark grid if no split is possible:
        cansplit(igrid) = .FALSE.
      END IF
    ELSE
      cansplit(igrid) = .FALSE.
    END IF
  END DO
  igrid = igrid+1
END DO
!
! Generate the newly determined subgrids:
info%nsubgrids = ngrid
DO i = 1, ngrid
  CALL MakeNewGrid(info,msubbounds(:,:,i))
  IF(MINVAL(msubbounds(1:ndims,2,i)-msubbounds(1:ndims,1,i)+1)<mgp) THEN
    PRINT *, 'BSAM 1.2: Error in NewSubGrids, grid smaller than minimumgridpoints'
    STOP
  END IF
  DO n = 1, ndims
    IF(MODULO(msubbounds(n,2,i)-msubbounds(n,1,i)+1,2)==1) THEN
      PRINT *, 'BSAM 1.2: Error in NewSubGrids, igrid=', i, ', grid length not even number'
      STOP
    END IF
    IF(MODULO(msubbounds(n,1,i),2)==0) THEN
      PRINT *, 'BSAM 1.2: Error in NewSubGrids, igrid=', i, ', index starts on even number'
      STOP
    END IF
    IF(MODULO(msubbounds(n,2,i),2)==1) THEN
      PRINT *, 'BSAM 1.2: Error in NewSubGrids, igrid=', i, ', index end on odd number'
      STOP
    END IF
  END DO
END DO
!
DEALLOCATE(signature,ddsignature,STAT=ierror)
IF(ierror/=0) THEN
  PRINT *,'BSAM 1.2: Error deallocating signatures arrays in NewSubGrids'
  STOP
END IF
!
END SUBROUTINE NewSubGrids
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION GetSignatures(errorflags,msubbounds,maxm) RESULT(gsresult)
USE NodeInfoDef
IMPLICIT NONE
!
INTEGER, DIMENSION(1:,1:,1:), INTENT(IN):: errorflags
INTEGER, DIMENSION(1:maxdims,1:2), INTENT(IN):: msubbounds
INTEGER, INTENT(IN):: maxm
INTEGER, DIMENSION(1:maxm,1:ndims):: gsresult
!
INTEGER:: i, n
INTEGER, DIMENSION(1:maxdims):: i1, i2
!
i1 = 1; i2 = 1
DO n = 1, ndims
  i1(1:ndims) = msubbounds(1:ndims,1)
  i2(1:ndims) = msubbounds(1:ndims,2)
  DO i = msubbounds(n,1), msubbounds(n,2)
    i1(n) = i; i2(n) = i
!
    gsresult(i,n) = SUM(errorflags(i1(1):i2(1),i1(2):i2(2),i1(3):i2(3)))
!
  END DO
END DO
!
END FUNCTION GetSignatures
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION GridFlagRatio(errorflags,mbounds) RESULT(gfrresult)
USE NodeInfoDef
IMPLICIT NONE
!
INTEGER, DIMENSION(1:,1:,1:), INTENT(IN):: errorflags
INTEGER, DIMENSION(1:maxdims,1:2), INTENT(IN OUT):: mbounds
REAL(KIND=r8):: gfrresult
!
REAL(KIND=r8):: flagged, total
!
total = REAL(PRODUCT(mbounds(1:ndims,2)-mbounds(1:ndims,1)+1),KIND=r8)
!
mbounds(ndims+1:maxdims,1:2) = 1
!
flagged = REAL(SUM(errorflags(mbounds(1,1):mbounds(1,2), &
                              mbounds(2,1):mbounds(2,2), &
                              mbounds(3,1):mbounds(3,2))),KIND=r8)
!
IF(flagged<-1.0E-08_r8) THEN
  PRINT *, 'BSAM 1.2: Error in GridFlagRatio: flagged < 0.'
  STOP
END IF
!
gfrresult = flagged/total
!
END FUNCTION GridFlagRatio
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE MakeNewGrid(parent,mbounds)
USE NodeInfoDef
USE TreeOps, ONLY: CreateChild, GetChildInfo
USE BSAMStorage, ONLY: AllocFields
USE Problem, ONLY: SetAux, SetSrc
IMPLICIT NONE
!
TYPE(nodeinfo):: parent
INTEGER, DIMENSION(1:maxdims,1:2), INTENT(IN):: mbounds
!
! Generate a new, finer grid within mbounds of the grid info:
!
TYPE(nodeinfo), POINTER:: child
INTEGER:: i, ierror, n, nb
INTEGER, DIMENSION(1:maxdims):: cmx
INTEGER, DIMENSION(1:maxdims,1:2):: mglobalbounds
REAL(KIND=r8):: rand
!
! Create a child of the currentnode:
CALL CreateChild
!
! Initialize the grid information for the new node:
! Get pointer to youngest child's info:
ierror = GetChildInfo(child)
!
! Start filling in the fields of child's info:
child%tobedeleted = .FALSE.
child%activegrid = .TRUE.
child%defective = .FALSE.
child%fieldsallocated = .FALSE.
child%initialgrid = parent%initialgrid
!
child%maxlevel = parent%maxlevel; child%nsubgrids = 0
child%level = parent%level+1
IF(child%level>finestlevel) finestlevel = child%level
!
child%nrvars = parent%nrvars
!
child%mbc = parent%mbc
!
child%nout = parent%nout; child%nframe = parent%nframe
child%outstyle = parent%outstyle; child%nroutvars = parent%nroutvars
!
child%mx = 1; child%mx(1:ndims) = (mbounds(1:ndims,2)-mbounds(1:ndims,1)+1)*2
child%maux = parent%maux
child%mbounds = 1; child%mbounds(1:ndims,:) = mbounds(1:ndims,:)
child%mglobal = 1
mglobalbounds(1:ndims,1) = parent%mglobal(1:ndims,1)+mbounds(1:ndims,1)-1
mglobalbounds(1:ndims,2) = mglobalbounds(1:ndims,1)+mbounds(1:ndims,2) &
                         - mbounds(1:ndims,1)
child%mglobal(1:ndims,1) = (mglobalbounds(1:ndims,1)-1)*2+1
child%mglobal(1:ndims,2) = mglobalbounds(1:ndims,2)*2
!
! First assume all boundaries are internal:
child%mthbc = internalbc
!
! Now check if we have any physical boundaries:
DO n = 1, ndims
  nb = 2*n-1
!
! If parent boundary condition is physical and child left is same as parent's 
! left:
  IF((parent%mthbc(nb)<internalbc) .AND. &
     (child%mbounds(n,1)==1)) THEN
!
! Then child left boundary is physical and inherited from parent:
    child%mthbc(nb)=parent%mthbc(nb)
  END IF
  nb = nb+1
!
! If parent boundary condition is physical and child right is same as parent's 
! right:
  IF((parent%mthbc(nb)<internalbc) .AND. &
     (child%mbounds(n,2)==parent%mx(n))) THEN
!
! Then child right boundary is physical and inherited from parent
    child%mthbc(nb) = parent%mthbc(nb)
  END IF
END DO
!
! ID the grids randomly:
CALL RANDOM_NUMBER(rand)
child%ngrid = NINT(10000000*rand)
!
child%gridtime = parent%gridtime
!
child%xlower = 0.0_r8; child%xupper = 0.0_r8
child%xlower(1:ndims) = REAL(mbounds(1:ndims,1)-1,KIND=r8)*parent%dx(1:ndims) &
                      + parent%xlower(1:ndims)
child%xupper(1:ndims) = REAL(mbounds(1:ndims,2)  ,KIND=r8)*parent%dx(1:ndims) &
                      + parent%xlower(1:ndims)
child%dx = 0.0_r8; child%dx(1:ndims) = parent%dx(1:ndims)/REAL(2,KIND=r8)
!
! Allocate dynamic space:
CALL AllocFields(child,parent)
!
cmx = 1
cmx(1:ndims) = child%mx(1:ndims)/2
child%levellandscape = parent%level
SELECT CASE(ndims)
CASE(2)
  child%levellandscape(1:cmx(1),1:cmx(2),1       ) = child%level
CASE(3)
  child%levellandscape(1:cmx(1),1:cmx(2),1:cmx(3)) = child%level
END SELECT
!
! Initialize field variables with values from parent or from initial data:
CALL InitFields(parent,child)
!
! Remove this later, if possible!
child%rf = 0.0_r8
parent%rf = 0.0_r8
child%f = 0.0_r8
parent%f = 0.0_r8
!
! Call user routine to initialize auxiliary and source-term array values:
IF(child%maux>0) CALL SetAux(child)
CALL SetSrc(child)
!
child%qold = child%q
!
END SUBROUTINE MakeNewGrid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE InitFields(parent,child)
USE NodeInfoDef
USE GridUtilities, ONLY:  BiLinProlongationP1MC,  BiLinProlongationP2MC, &
                         TriLinProlongationP1MC, TriLinProlongationP2MC
IMPLICIT NONE
!
! Redone to support only bilinear interpolation and 1st layer updating.
!
TYPE(nodeinfo):: parent, child
!
INTEGER:: ierror, interpopt, mbc, nrvars
INTEGER, DIMENSION(1:maxdims):: cmx, mx
INTEGER, DIMENSION(1:maxdims,1:2):: mb
!
nrvars = parent%nrvars
mbc = parent%mbc
!
! Check whether we can apply user-provided QInit for initialization.
IF(child%initialgrid) THEN
  CALL QInit(child)
  RETURN
END IF
!
mx = 1; mx(1:ndims) = child%mx(1:ndims)
cmx = 1; cmx(1:ndims) = mx(1:ndims)/2
mb = 1; mb(1:ndims,1:2) = child%mbounds(1:ndims,1:2)
!
SELECT CASE(ndims)
  CASE(2)
    child%qc(1-mbc:cmx(1)+mbc,1-mbc:cmx(2)+mbc,1,1:nrvars) &
      = parent%q(mb(1,1)-mbc:mb(1,2)+mbc, &
                 mb(2,1)-mbc:mb(2,2)+mbc,1,1:nrvars)
!
    SELECT CASE(mbc)
      CASE(1)
            child%q(  0: mx(1)+1, 0: mx(2)+1,1,1:nrvars) &
          = BiLinProlongationP1MC( &
            child%qc( 0:cmx(1)+1, 0:cmx(2)+1,1,1:nrvars))
      CASE(2)
            child%q( -1: mx(1)+2,-1: mx(2)+2,1,1:nrvars) &
          = BiLinProlongationP2MC( &
            child%qc(-1:cmx(1)+2,-1:cmx(2)+2,1,1:nrvars))
      CASE DEFAULT
        PRINT *, 'InitFields: only supports mbc=1,2.'
    END SELECT
!
  CASE(3)
    child%qc(1-mbc:cmx(1)+mbc,1-mbc:cmx(2)+mbc,1-mbc:cmx(3)+mbc,1:nrvars) &
      = parent%q(mb(1,1)-mbc:mb(1,2)+mbc, &
                 mb(2,1)-mbc:mb(2,2)+mbc, &
                 mb(3,1)-mbc:mb(3,2)+mbc,1:nrvars)
!
    SELECT CASE(mbc)
      CASE(1)
            child%q(  0: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1:nrvars) &
          = TriLinProlongationP1MC( &
            child%qc( 0:cmx(1)+1, 0:cmx(2)+1, 0:cmx(3)+1,1:nrvars))
      CASE(2)
            child%q( -1: mx(1)+2,-1: mx(2)+2,-1: mx(3)+2,1:nrvars) &
          = TriLinProlongationP2MC( &
            child%qc(-1:cmx(1)+2,-1:cmx(2)+2,-1:cmx(3)+2,1:nrvars))
      CASE DEFAULT
        PRINT *, 'InitFields: only supports mbc=1,2.'
    END SELECT
!
  CASE DEFAULT
    PRINT *, 'InitFields: only supports ndims=2,3.'
END SELECT
!
END SUBROUTINE InitFields
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION TransferValues(grid1,grid2,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok
IMPLICIT NONE
!
TYPE(nodeinfo):: grid1, grid2
TYPE(funcparam):: dummy
!
! Transfer values from previous grids on this level to newly created grids:
!    
TransferValues = err_ok
!
IF(grid1%activegrid .AND. (.NOT. grid2%activegrid)) THEN
!
! Look for overlap and transfer grid values from Grid2 to Grid1
  CALL Transferq(grid2,grid1)
END IF
!
IF(grid2%activegrid .AND. (.NOT. grid1%activegrid)) THEN
!
! Look for overlap and transfer grid values from Grid1 to Grid2
  CALL Transferq(grid1,grid2)
END IF
!
! We either have two new grids or two old grids, nothing to be done:
RETURN
END FUNCTION TransferValues
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE Transferq(sourceinfo,targetinfo)
USE NodeInfoDef
IMPLICIT NONE
!
TYPE(nodeinfo):: sourceinfo
TYPE(nodeinfo):: targetinfo
!
! Look for overlap and transfer grid values from source to target
!    
INTEGER, DIMENSION(1:maxdims,1:2):: mtarget, msource
!
IF(.NOT. sourceinfo%fieldsallocated) THEN
  PRINT *, 'BSAM 1.2 Error: Trying to transfer values from unallocated'
  PRINT *, '                sourceinfo in Transferq.'
  STOP
END IF
!
IF(.NOT. targetinfo%fieldsallocated) THEN
  PRINT *, 'BSAM 1.2 Error: Trying to transfer values to unallocated'
  PRINT *, '                targetinfo in Transferq.'
  STOP
END IF
!
! Look for overlap (Later, maybe copy over ghost cells as well):
msource = sourceinfo%mglobal; mtarget = targetinfo%mglobal
!
CALL TransferOverlap(msource,mtarget,sourceinfo,targetinfo)
!
END SUBROUTINE Transferq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE TransferOverlap(msource,mtarget,sourceinfo,targetinfo)
USE NodeInfoDef
IMPLICIT NONE
!
INTEGER, DIMENSION(1:maxdims,1:2), INTENT(IN):: msource
INTEGER, DIMENSION(1:maxdims,1:2), INTENT(IN):: mtarget
TYPE(nodeinfo):: sourceinfo
TYPE(nodeinfo):: targetinfo
!
INTEGER:: n
INTEGER, DIMENSION(1:maxdims,1:2):: moverlap, ms, mt
!
! Transfer values from source to target in overlap region:
! 1. Find overlap region in global index space:
moverlap = 1
moverlap(1:ndims,1) = MAX(msource(1:ndims,1),mtarget(1:ndims,1))
moverlap(1:ndims,2) = MIN(msource(1:ndims,2),mtarget(1:ndims,2))
!
! 2. Check for nonempty intersection:
IF(ANY(moverlap(:,2)-moverlap(:,1)<0)) RETURN
!
! 3. Transform common index space to grid index spaces:
ms = 1; mt = 1
DO n = 1, ndims
  ms(n,:) = moverlap(n,:)-sourceinfo%mglobal(n,1)+1
  mt(n,:) = moverlap(n,:)-targetinfo%mglobal(n,1)+1
END DO
!
! 4. Carry out the transfer:
    targetinfo%q(mt(1,1):mt(1,2),mt(2,1):mt(2,2),mt(3,1):mt(3,2),:) &
  = sourceinfo%q(ms(1,1):ms(1,2),ms(2,1):ms(2,2),ms(3,1):ms(3,2),:)
!
END SUBROUTINE TransferOverlap
END MODULE BSAMRoutines
