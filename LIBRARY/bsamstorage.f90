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
! File:             bsamstorage.f90
! Purpose:          BSAM memory allocation module.
! Contains:
! Revision History: Ver. 1.0 Oct. 2006 Steven Wise
! Revision History: Ver. 1.1 May. 2007 Steven Wise
! Revision History: Ver. 1.2 Jul. 2015 Steven Wise
! -----------------------------------------------------------------------
MODULE BSAMStorage
!
CONTAINS
!
SUBROUTINE AllocFields(info,parent)
USE NodeInfoDef
IMPLICIT NONE
!
TYPE(nodeinfo):: info
TYPE(nodeinfo), INTENT(IN), OPTIONAL:: parent
!
INTEGER:: ierror, maux, mmaux, mbc, n, nrvars
INTEGER, DIMENSION(1:maxdims):: amx, cmx, mx
!
mx = info%mx; mbc = info%mbc
nrvars = info%nrvars; maux = info%maux
cmx = 1
DO n = 1, ndims
  IF(MODULO(mx(n),2)==1) THEN
    cmx(n) = 1
  ELSE
    cmx(n) = mx(n)/2
  END IF
END DO
!
mmaux = MAX(maux,1)
amx = 1; IF(maux>0) amx(1:ndims) = mx(1:ndims)
!
! Allocate space for all POINTER components from nodeinfodef:
SELECT CASE(ndims)
  CASE(2)
    ALLOCATE( &
      info%errorflags(1    : mx(1)    ,1    : mx(2)    ,1:1         ), &
      info%q(         1-mbc: mx(1)+mbc,1-mbc: mx(2)+mbc,1:1,1:nrvars), &
      info%qold(      1-mbc: mx(1)+mbc,1-mbc: mx(2)+mbc,1:1,1:nrvars), &
      info%qc(        1-mbc:cmx(1)+mbc,1-mbc:cmx(2)+mbc,1:1,1:nrvars), &
      info%qrte(      1    :cmx(1)    ,1    :cmx(2)    ,1:1,1:nrvars), &
      info%aux(       1-mbc:amx(1)+mbc,1-mbc:amx(2)+mbc,1:1,1:mmaux ), &
      info%f(         1    : mx(1)    ,1    : mx(2)    ,1:1,1:nrvars), &
      info%rf(        1    : mx(1)    ,1    : mx(2)    ,1:1,1:nrvars), &
      info%ftmp(      1    : mx(1)    ,1    : mx(2)    ,1:1,1:nrvars), &
        STAT=ierror)
!
    IF(ierror/=0) THEN
      PRINT *, 'BSAM 1.2: Error in allocation of nodeinfodef components in AllocFields'
      STOP
    END IF
!
	ALLOCATE(info%levellandscape(0:cmx(1)+1,0:cmx(2)+1,1:1),STAT=ierror)
!
    IF(ierror/=0) THEN
      PRINT *, 'BSAM 1.2: Error in allocation of nodeinfodef components in AllocFields'
      STOP
    END IF
!
  CASE(3)
    ALLOCATE( &
      info%errorflags(1    : mx(1)    ,1    : mx(2)    ,1    : mx(3)             ), &
      info%q(         1-mbc: mx(1)+mbc,1-mbc: mx(2)+mbc,1-mbc: mx(3)+mbc,1:nrvars), &
      info%qold(      1-mbc: mx(1)+mbc,1-mbc: mx(2)+mbc,1-mbc: mx(3)+mbc,1:nrvars), &
      info%qc(        1-mbc:cmx(1)+mbc,1-mbc:cmx(2)+mbc,1-mbc:cmx(3)+mbc,1:nrvars), &
      info%qrte(      1    :cmx(1)    ,1    :cmx(2)    ,1    :cmx(3)    ,1:nrvars), &
      info%aux(       1-mbc:amx(1)+mbc,1-mbc:amx(2)+mbc,1-mbc:amx(3)+mbc,1:mmaux ), &
      info%f(         1    : mx(1)    ,1    : mx(2)    ,1    : mx(3)    ,1:nrvars), &
      info%rf(        1    : mx(1)    ,1    : mx(2)    ,1    : mx(3)    ,1:nrvars), &
      info%ftmp(      1    : mx(1)    ,1    : mx(2)    ,1    : mx(3)    ,1:nrvars), &
        STAT=ierror)
!
    IF(ierror/=0) THEN
      PRINT *, 'BSAM 1.2: Error in allocation of nodeinfodef components in AllocFields'
      STOP
    END IF
!
	ALLOCATE(info%levellandscape(0:cmx(1)+1,0:cmx(2)+1,0:cmx(3)+1),STAT=ierror)
!
    IF(ierror/=0) THEN
      PRINT *, 'BSAM 1.2: Error in allocation of nodeinfodef components in AllocFields'
      STOP
    END IF
!
  CASE DEFAULT
    PRINT *, 'AllocFields: Only ndims = 2,3 are supported.'
    STOP
END SELECT
!
! Initialize floating point arrays:
info%q = 0.0_r8; info%qold = 0.0_r8; info%qc = 0.0_r8; info%qrte = 0.0_r8
info%aux = 0.0_r8; info%f = 0.0_r8; info%rf = 0.0_r8; info%ftmp = 0.0_r8
!
! Initialize error flags:
info%errorflags = 0
!
info%fieldsallocated = .TRUE.
!
END SUBROUTINE AllocFields
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE DeAllocFields(info)
USE NodeInfoDef
IMPLICIT NONE
!
TYPE(nodeinfo):: info
!
IF(.NOT. info%fieldsallocated) RETURN
!
IF(ASSOCIATED(info%errorflags))     DEALLOCATE(info%errorflags)
IF(ASSOCIATED(info%q))              DEALLOCATE(info%q)
IF(ASSOCIATED(info%qold))           DEALLOCATE(info%qold)
IF(ASSOCIATED(info%qc))             DEALLOCATE(info%qc)
IF(ASSOCIATED(info%qrte))           DEALLOCATE(info%qrte)
IF(ASSOCIATED(info%aux))            DEALLOCATE(info%aux)
IF(ASSOCIATED(info%f))              DEALLOCATE(info%f)
IF(ASSOCIATED(info%rf))             DEALLOCATE(info%rf)
IF(ASSOCIATED(info%ftmp))           DEALLOCATE(info%ftmp)
IF(ASSOCIATED(info%levellandscape)) DEALLOCATE(info%levellandscape)
!
! Nullify the pointers:
NULLIFY(info%errorflags,info%q,info%qold,info%qc,info%qrte,info%aux,info%f, &
        info%rf,info%ftmp,info%levellandscape)
!
info%fieldsallocated = .FALSE.  ! Don't compute on this node:
info%tobedeleted     = .TRUE.   ! Mark for garbage collection:
!
END SUBROUTINE DeAllocFields
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE AllocPeriodicBCStorage(np,nrvars)
USE NodeInfoDef
IMPLICIT NONE
!
INTEGER, INTENT(IN):: np
INTEGER, INTENT(IN):: nrvars
!
INTEGER:: ierror
!
ALLOCATE(periodicoffsetindex(1:np),poffset(1:maxdims,1:np), &
         qoffset(1:np,1:nrvars),STAT=ierror)
!
IF(ierror/=0) THEN
  PRINT *,'BSAM 1.2 Error in AllocPeriodicBCStorage.'
  STOP
END IF
!
END SUBROUTINE AllocPeriodicBCStorage
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE DeallocPeriodicBCStorage
USE NodeInfoDef
IMPLICIT NONE
!
IF(ALLOCATED(periodicoffsetindex)) DEALLOCATE(periodicoffsetindex)
IF(ALLOCATED(poffset)) DEALLOCATE(poffset)
IF(ALLOCATED(qoffset)) DEALLOCATE(qoffset)
!
END SUBROUTINE DeallocPeriodicBCStorage
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE AllocUniformGrids(mxroot,mbc,nrvars)
USE NodeInfoDef
IMPLICIT NONE
!
INTEGER, DIMENSION(1:maxdims), INTENT(IN):: mxroot ! root-level grid size.
INTEGER, INTENT(IN):: mbc
INTEGER, INTENT(IN):: nrvars
!
INTEGER:: ierror, level, n, r
INTEGER, DIMENSION(1:maxdims):: high, low, mx
!
mx = 1; mx(1:ndims) = mxroot(1:ndims)
low = 1; low(1:ndims) = 1-mbc
high = 1
!
DO level = 0, maxlevel
!
  IF(level>0) mx(1:ndims) = mx(1:ndims)*2
  uniformgrid(level)%mx = 1; uniformgrid(level)%mx(1:ndims) = mx(1:ndims)
!
  high(1:ndims) = mx(1:ndims)+mbc
!
  ALLOCATE( &
    uniformgrid(level)%q(low(1):high(1),low(2):high(2),low(3):high(3),1:nrvars), &
    STAT=ierror)
!
  IF(ierror/=0) THEN
    PRINT *, 'AllocUniformGrids: Error allocating uniformgrid on level', &
             level
    STOP
  END IF
!  
END DO
!
END SUBROUTINE AllocUniformGrids
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE DeallocUniformGrids
USE NodeInfoDef
IMPLICIT NONE
!
INTEGER:: level
!
DO level = rootlevel, maxlevel
  uniformgrid(level)%mx = 0
  IF(ASSOCIATED(uniformgrid(level)%q)) DEALLOCATE(uniformgrid(level)%q)
END DO
!
END SUBROUTINE DeallocUniformGrids
!
END MODULE BSAMStorage