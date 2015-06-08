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
! File:             bsaminputoutput.f90
! Purpose:          BSAM I/O module.
! Contains:
! Revision History: Ver. 1.0 Oct. 2006 Steven Wise
! Revision History: Ver. 1.1 May. 2007 Steven Wise
! Revision History: Ver. 1.2 Jul. 2015 Steven Wise
! -----------------------------------------------------------------------
MODULE BSAMInputOutput
USE NodeInfoDef
IMPLICIT NONE
!
SAVE
PRIVATE
PUBLIC WriteQ, ReadQ, WriteUniformMeshQ
!
CONTAINS
!
SUBROUTINE WriteQ(nframe,time)
USE NodeinfoDef
USE TreeOps, ONLY: ApplyOnForest, ApplyOnLevel
IMPLICIT NONE
!
INTEGER, INTENT(IN):: nframe
REAL(KIND=r8), INTENT(IN):: time
!
TYPE(funcparam):: dummy
CHARACTER(LEN=16):: filename
INTEGER:: level
!
! Matlab format:
WRITE(filename,'(A7,I5.5,A4)') './OUT/m', nframe, '.dat'
OPEN(UNIT=54,FILE=filename,ACTION='READWRITE',STATUS='REPLACE',FORM='FORMATTED')
WRITE(54,'(F25.12)') time
WRITE(54,'(I3)') finestlevel
DO level = rootlevel, finestlevel
  CALL ApplyOnLevel(level,OutputQ,dummy)
END DO
CLOSE(54)
!
END SUBROUTINE WriteQ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION OutputQ(info,dummy)
USE NodeinfoDef
USE TreeOps, ONLY: err_ok
IMPLICIT NONE
!
TYPE(nodeinfo) :: info
TYPE(funcparam) :: dummy
!
INTEGER:: nrvars
INTEGER, DIMENSION(1:maxdims):: mx
!
Outputq = err_ok
!
IF(.NOT. info%fieldsallocated) RETURN
!
mx = info%mx; nrvars = info%nrvars
!
WRITE(54,*) ' '
WRITE(54,'(I3,3(1x,I3))') info%level, ndims, 2, nrvars
!
SELECT CASE(ndims)
  CASE(2)
    WRITE(54,2001) info%dx(1), info%dx(2)
    WRITE(54,2001) info%xlower(1), info%xlower(2)
    WRITE(54,2001) info%xupper(1), info%xupper(2)
    2001 FORMAT(F25.12,1x,F25.12)
    WRITE(54,'(I8,1X,I8)') info%mx(1), info%mx(2)
    WRITE(54,'(I8,3(1X,I8))') info%mglobal(1,1), info%mglobal(1,2), &
                              info%mglobal(2,1), info%mglobal(2,2)
    CALL WriteQ2D(info%q(0:mx(1)+1,0:mx(2)+1,1,1:nrvars),mx(1:2), &
                  nrvars)
  CASE(3)
    WRITE(54,2002) info%dx(1), info%dx(2), info%dx(3)
    WRITE(54,2002) info%xlower(1), info%xlower(2), info%xlower(3)
    WRITE(54,2002) info%xupper(1), info%xupper(2), info%xupper(3)
    2002 FORMAT(F25.12,2(1x,F25.12))
    WRITE(54,'(I8,2(1X,I8))') info%mx(1), info%mx(2), info%mx(3)
    WRITE(54,'(I8,5(1X,I8))') info%mglobal(1,1), info%mglobal(1,2), &
                              info%mglobal(2,1), info%mglobal(2,2), &
                              info%mglobal(3,1), info%mglobal(3,2)
    CALL WriteQ3D(info%q(0:mx(1)+1,0:mx(2)+1,0:mx(3)+1,1:nrvars),mx(1:3), &
                  nrvars)
END SELECT
!
END FUNCTION OutputQ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE WriteQ2D(q,mx,nrvars)
USE NodeInfoDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: q
INTEGER, DIMENSION(1:2), INTENT(IN):: mx
INTEGER, INTENT(IN):: nrvars
!
INTEGER:: i, j, k
REAL(KIND=r8):: rdummy
!
DO j = 0, mx(2)+1
  DO i = 0, mx(1)+1
!
    WRITE(54,3001) (q(i,j,k),k=1,nrvars)
    3001 FORMAT(10(F25.12,2X))
    BACKSPACE(54)
    READ(54,'(F25.12)') rdummy
!
  END DO
END DO 
!
END SUBROUTINE WriteQ2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE WriteQ3D(q,mx,nrvars)
USE NodeInfoDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:,0:,1:), INTENT(IN):: q
INTEGER, DIMENSION(1:3), INTENT(IN):: mx
INTEGER, INTENT(IN):: nrvars
!
INTEGER:: i, j, k, l
REAL(KIND=r8):: rdummy
!
DO k = 0, mx(3)+1
  DO j = 0, mx(2)+1
    DO i = 0, mx(1)+1
!
      WRITE(54,3001) (q(i,j,k,l),l=1,nrvars)
      3001 FORMAT(10(F25.12,2X))
      BACKSPACE(54)
      READ(54,'(F25.12)') rdummy
!
    END DO
  END DO
END DO
!
END SUBROUTINE WriteQ3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE ReadQ
USE NodeInfoDef
USE TreeOps, ONLY: CreateChild, CurrentNodeToYoungest, GetChildInfo, &
                   GetRootInfo
USE BSAMStorage, ONLY: AllocFields
IMPLICIT NONE
!
TYPE(nodeinfo), POINTER:: rootinfo, info
LOGICAL:: fileexist
CHARACTER(LEN=1):: string
CHARACTER(LEN=16):: filename
INTEGER:: ierror, ioerror, level, lvl, mbc, restartndims, npatch, nrvars, &
          restartfinestlevel, r
INTEGER, DIMENSION(1:maxdims):: mx
INTEGER, DIMENSION(1:maxdims,1:2):: mg
REAL(KIND=r8), PARAMETER:: small = 1.0E-08_r8
REAL(KIND=r8):: rdummy
REAL(KIND=r8), DIMENSION(:,:,:,:), ALLOCATABLE:: patch
REAL(KIND=r8), DIMENSION(1:maxdims):: dx, x, xlower, xupper
!
! On read-in the usual parent-child relationship is broken.  For simplicity,
! the youngest grid on level=l-1 is the parent of all level=l grids.  In
! particular, coarse-fine boundary conditions can not be properly enforced.
! The first ghost layer values are intact on read-in.  However, even these are
! not needed at restart.
!
! Initialized values for the root grid:
ierror = GetRootInfo(rootinfo)
info => rootinfo
!
mx = 1; mx(1:ndims) = rootinfo%mx(1:ndims)
mbc = rootinfo%mbc
!
WRITE(filename,'(A7,I5.5,A4)') './OUT/m', restartframe, '.dat'
INQUIRE(FILE=filename,EXIST=fileexist)
IF(.NOT. fileexist) THEN
  PRINT *, 'Readq: Error; input file does not exist: ', filename
  STOP
END IF
OPEN(UNIT=54,FILE=filename,STATUS='OLD',FORM='FORMATTED',ACTION='READ', &
     IOSTAT=ioerror)
IF(ioerror/=0) THEN
  PRINT *,'Readq: Error opening restart file ', filename
  STOP
END IF
READ(54,'(F25.12)') restarttime
!
info%gridtime = restarttime
!
READ(54,'(I3)') restartfinestlevel
!
finestlevel = restartfinestlevel
!
IF(restartfinestlevel>maxlevel) THEN
  PRINT *, 'Readq: Error; restartfinestlevel>maxlevel.'
  STOP
END IF
!
npatch = 0
!
level_loop: DO lvl = rootlevel, restartfinestlevel
!
  IF(lvl>rootlevel) CALL CurrentNodeToYoungest(lvl-1)
!
  patch_loop: DO
!  
    READ(54,'(A1)',IOSTAT=ioerror) string; IF(ioerror<0) EXIT level_loop
    READ(54,'(I3,3(1x,I3))') level, restartndims, r, nrvars
!
! Check to see whether we've read all the patches on this level:
    IF(level/=lvl) THEN
      level = lvl
      BACKSPACE(54); BACKSPACE(54)
      EXIT patch_loop
    END IF
!
    IF(level>rootlevel) THEN
      CALL CreateChild
      NULLIFY(info)
      ierror = GetChildInfo(info)
    END IF
!
    IF(restartndims/=ndims) THEN
      PRINT *, 'Readq: Error reading data file ', filename
      PRINT *, 'Spatial dimension (ndims) inconsistancy on restart.'
      STOP
    END IF
!
    IF(r/=2) THEN
      PRINT *, 'Readq: Error reading data file ', filename
      PRINT *, 'Refinement ratio (r) inconsistancy on restart.'
      PRINT *, 'Currently only r = 2 is supported for restarting.'
      STOP
    END IF
!
    SELECT CASE(ndims)
      CASE(2)
        READ(54,2001) dx(1), dx(2)
        READ(54,2001) xlower(1), xlower(2)
        READ(54,2001) xupper(1), xupper(2)
        2001 FORMAT(F25.12,1x,F25.12)
        READ(54,'(I8,1X,I8)') mx(1), mx(2)
        READ(54,'(I8,3(1X,I8))') mg(1,1), mg(1,2), mg(2,1), mg(2,2)
      CASE(3)
        READ(54,2002) dx(1), dx(2), dx(3)
        READ(54,2002) xlower(1), xlower(2), xlower(3)
        READ(54,2002) xupper(1), xupper(2), xupper(3)
        2002 FORMAT(F25.12,2(1x,F25.12))
        READ(54,'(I8,2(1X,I8))') mx(1), mx(2), mx(3)
        READ(54,'(I8,5(1X,I8))') mg(1,1), mg(1,2), mg(2,1), mg(2,2), &
                                 mg(3,1), mg(3,2)
    END SELECT
!
! Root-level grid constructs have already been set.  Check for errors:
    IF(level==rootlevel) THEN
!
      IF(nrvars/=rootinfo%nrvars) THEN
        PRINT *, 'Readq: Error reading data file ', filename
        PRINT *, 'Number of variables (nrvars) inconsistancy on restart.'
        STOP
      END IF
!
      IF(ANY(dx(1:ndims)-rootinfo%dx(1:ndims)>small)) THEN
        PRINT *, 'Readq: Error reading data file ', filename
        PRINT *, 'Root-level grid spacing (dx) inconsistancy on restart.'
        STOP
      END IF
!
      IF(ANY(ABS(xlower(1:ndims)-rootinfo%xlower(1:ndims))>small)) THEN
        PRINT *, 'Readq: Error reading data file ', filename
        PRINT *, 'Root-level grid location (xlower) inconsistancy on restart.'
        STOP
      END IF
!
      IF(ANY(ABS(xupper(1:ndims)-rootinfo%xupper(1:ndims))>small)) THEN
        PRINT *, 'Readq: Error reading data file ', filename
        PRINT *, 'Root-level grid location (xupper) inconsistancy on restart.'
        STOP
      END IF
!
      IF(ANY(mx(1:ndims)/=rootinfo%mx(1:ndims))) THEN
        PRINT *, 'Readq: Error reading data file ', filename
        PRINT *, 'Root-level grid size (mx) inconsistancy on restart.'
        STOP
      END IF
!
      IF(ANY(mg(1:ndims,1)/=1) .OR. &
         ANY(mg(1:ndims,2)/=rootinfo%mx(1:ndims))) THEN
        PRINT *, 'Readq: Error reading data file ', filename
        PRINT *, 'Root-level grid size (mglobal) inconsistancy on restart.'
        STOP
      END IF
    END IF
!
    info%tobedeleted = .FALSE.; info%activegrid = .FALSE.
    info%initialgrid = .FALSE.
!
! These rootlevel constructs should not be changed:
    IF(level>rootlevel) THEN
      info%fieldsallocated = .FALSE.
      info%maxlevel = maxlevel; info%nsubgrids = 0
      info%level = level
!
      info%nrvars = nrvars
!
      info%mbc = rootinfo%mbc
!
      info%nout = rootinfo%nout; info%nframe = rootinfo%nframe
      info%outstyle = rootinfo%outstyle; info%nroutvars = rootinfo%nroutvars
!
      info%mx = 1; info%mx(1:ndims) = mx(1:ndims)
      info%maux = rootinfo%maux
      info%mglobal = 1
      info%mglobal(1:ndims,1) = mg(1:ndims,1)
      info%mglobal(1:ndims,2) = mg(1:ndims,2)
!
! These two items are broken on restart, but are not needed:
      info%mthbc = internalbc
      info%mbounds = 1
!
      npatch = npatch+1
!      info%ngrid = npatch
      info%ngrid = -13
!
      info%gridtime = restarttime
!
      info%xlower = 0.0_r8; info%xupper = 0.0_r8
      info%xlower(1:ndims) = xlower(1:ndims)
      info%xupper(1:ndims) = xupper(1:ndims)
      info%dx = 0.0_r8; info%dx(1:ndims) = dx(1:ndims)
!
      CALL AllocFields(info)
    END IF
!
! Read-in field data:
    SELECT CASE(ndims)
      CASE(2)
        CALL ReadQ2D(info%q(0:mx(1)+1,0:mx(2)+1,1,1:nrvars),mx(1:2), &
                     nrvars)
      CASE(3)
        CALL ReadQ3D(info%q(0:mx(1)+1,0:mx(2)+1,0:mx(3)+1,1:nrvars),mx(1:3), &
                     nrvars)
    END SELECT
!
  END DO patch_loop
END DO level_loop
!
CLOSE(54)
!
END SUBROUTINE ReadQ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE ReadQ2D(q,mx,nrvars)
USE NodeInfoDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(OUT):: q
INTEGER, DIMENSION(1:2), INTENT(IN):: mx
INTEGER, INTENT(IN):: nrvars
!
INTEGER:: i, j, k
!
DO j = 0, mx(2)+1
  DO i = 0, mx(1)+1
!
    READ(54,3001) (q(i,j,k),k=1,nrvars)
    3001 FORMAT(10(F25.12,2X))
!
  END DO
END DO 
!
END SUBROUTINE ReadQ2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE ReadQ3D(q,mx,nrvars)
USE NodeInfoDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:,0:,1:), INTENT(OUT):: q
INTEGER, DIMENSION(1:3), INTENT(IN):: mx
INTEGER, INTENT(IN):: nrvars
!
INTEGER:: i, j, k, l
!
DO k = 0, mx(3)+1
  DO j = 0, mx(2)+1
    DO i = 0, mx(1)+1
!
      READ(54,3001) (q(i,j,k,l),l=1,nrvars)
      3001 FORMAT(10(F25.12,2X))
!
    END DO
  END DO
END DO
!
END SUBROUTINE ReadQ3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE WriteUniformMeshQ(nframe,time)
USE NodeInfoDef
USE TreeOps, ONLY: ApplyOnLevel, GetRootInfo
USE GridUtilities, ONLY: BiLinProlongationP1MC, TriLinProlongationP1MC
USE BSAMStorage, ONLY: AllocUniformGrids, DeallocUniformGrids
IMPLICIT NONE
!
INTEGER, INTENT(IN):: nframe
REAL(KIND=r8), INTENT(IN):: time
!
TYPE(nodeinfo), POINTER:: rootinfo
TYPE(funcparam):: dummy
CHARACTER(LEN=16):: filename
INTEGER:: ierror, i, j, k, l, level, mbc, nrvars
INTEGER, DIMENSION(1:maxdims):: high, low, mx, mxuc, mxuf
REAL(KIND=r8):: rdummy
REAL(KIND=r8), DIMENSION(1:maxdims):: dx, xlower, xupper
!
ierror = GetRootInfo(rootinfo)
!
nrvars = rootinfo%nrvars
mbc = rootinfo%mbc
dx = 0.0_r8; dx(1:ndims) = rootinfo%dx(1:ndims)
xlower = 0.0_r8; xlower(1:ndims) = rootinfo%xlower(1:ndims)
xupper = 0.0_r8; xupper(1:ndims) = rootinfo%xupper(1:ndims)
!
mx = 1; mx(1:ndims) = rootinfo%mx(1:ndims)
!
CALL AllocUniformGrids(mx,mbc,nrvars)
!
low = 1; low(1:ndims) = 0
high = 1; high(1:ndims) = mx(1:ndims)+1
!
uniformgrid(rootlevel)%q(low(1):high(1),low(2):high(2),low(3):high(3),1:nrvars) &
            = rootinfo%q(low(1):high(1),low(2):high(2),low(3):high(3),1:nrvars)
!
mxuf = mx
!
DO level = 1, finestlevel
!
  mxuc = mxuf
  mxuf = 1; mxuf(1:ndims) = 2*mxuc(1:ndims)
  dx(1:ndims) = dx(1:ndims)/2
!
  SELECT CASE(ndims)
    CASE(2)
          uniformgrid(level  )%q(0:mxuf(1)+1,0:mxuf(2)+1,1,1:nrvars) &
        = BiLinProlongationP1MC( &
          uniformgrid(level-1)%q(0:mxuc(1)+1,0:mxuc(2)+1,1,1:nrvars))
    CASE(3)
          uniformgrid(level  )%q(0:mxuf(1)+1,0:mxuf(2)+1,0:mxuf(3)+1,1:nrvars) &
        = TriLinProlongationP1MC( &
          uniformgrid(level-1)%q(0:mxuc(1)+1,0:mxuc(2)+1,0:mxuc(3)+1,1:nrvars))
  END SELECT
!
  CALL ApplyOnLevel(level,CopyPatchToUniformGrid,dummy)
END DO
!
WRITE(filename,'(A7,I5.5,A4)') './OUT/u', nframe, '.dat'
OPEN(UNIT=54,FILE=filename,STATUS='REPLACE',FORM='FORMATTED')
!
WRITE(54,'(F25.12)') time
WRITE(54,'(I3)') 0
WRITE(54,*) ' '
WRITE(54,'(I3,3(1x,I3))') 0, ndims, 2, nrvars
!
SELECT CASE(ndims)
  CASE(2)
    WRITE(54,2001) dx(1), dx(2)
    WRITE(54,2001) xlower(1), xlower(2)
    WRITE(54,2001) xupper(1), xupper(2)
    2001 FORMAT(F25.12,1x,F25.12)
    WRITE(54,'(I8,1X,I8)') mxuf(1), mxuf(2)
    WRITE(54,'(I8,3(1X,I8))') 1, mxuf(1), 1, mxuf(2)
    CALL WriteQ2D(uniformgrid(finestlevel)%q(0:mxuf(1)+1, &
                                             0:mxuf(2)+1, &
                                             1          , &
                                             1:nrvars    ),mxuf(1:2),nrvars)
  CASE(3)
    WRITE(54,2002) dx(1), dx(2), dx(3)
    WRITE(54,2002) xlower(1), xlower(2), xlower(3)
    WRITE(54,2002) xupper(1), xupper(2), xupper(3)
    2002 FORMAT(F25.12,2(1x,F25.12))
    WRITE(54,'(I8,2(1X,I8))') mxuf(1), mxuf(2), mxuf(3)
    WRITE(54,'(I8,5(1X,I8))') 1, mxuf(1), 1, mxuf(2), 1, mxuf(3)
    CALL WriteQ3D(uniformgrid(finestlevel)%q(0:mxuf(1)+1, &
                                             0:mxuf(2)+1, &
                                             0:mxuf(3)+1, &
                                             1:nrvars    ),mxuf(1:3),nrvars)
END SELECT
!
CALL DeallocUniformGrids
!
END SUBROUTINE WriteUniformMeshQ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION CopyPatchToUniformGrid(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok
IMPLICIT NONE
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
INTEGER:: level, nrvars
INTEGER, DIMENSION(1:maxdims):: h1, h2, l1, l2, mx
INTEGER, DIMENSION(1:maxdims,1:2):: mg
!
CopyPatchToUniformGrid = err_ok
!
IF(info%tobedeleted) RETURN
!
nrvars = info%nrvars
level = info%level
!
mg = 1; mg(1:ndims,1:2) = info%mglobal(1:ndims,1:2)
mx = 1; mx(1:ndims) = info%mx(1:ndims)
!
l1 = 1; l1(1:ndims) = mg(1:ndims,1)-1
h1 = 1; h1(1:ndims) = mg(1:ndims,2)+1
l2 = 1; l2(1:ndims) = 0
h2 = 1; h2(1:ndims) = mx(1:ndims)+1
!
uniformgrid(level)%q(l1(1):h1(1),l1(2):h1(2),l1(3):h1(3),1:nrvars) &
  =           info%q(l2(1):h2(1),l2(2):h2(2),l2(3):h2(3),1:nrvars)
!
END FUNCTION CopyPatchToUniformGrid
!
END MODULE BSAMInputOutput