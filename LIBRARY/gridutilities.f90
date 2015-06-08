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
! associated text is reproduced on all copies.  For all other uses,
! including distribution of modified versions, please contact the authors.
!
! Commercial use is strictly forbidden without permission.
!
! This software is made available "as is" without any assurance that it
! will work for your purposes. The software may in fact have defects,
! so use the software at your own risk.
!
! -----------------------------------------------------------------------
! File:             gridutilities.f90
! Purpose:          BSAM grid utilities module.
! Contains:
! Revision History: Ver. 1.0 Oct. 2006 Steven Wise
! Revision History: Ver. 1.1 May. 2007 Steven Wise
! Revision History: Ver. 1.2 Jul. 2015 Steven Wise
! -----------------------------------------------------------------------
MODULE GridUtilities
!
! Contains the commonly-used 2 and 3D grid utilities.  These only operate on 
! uniform, rectangular grid patches.
!
CONTAINS
!
FUNCTION ULap2D(a) RESULT(ulapresult)
USE NodeInfoDef
IMPLICIT NONE
!
! 2D UNDIVIDED laplacian of a(0:mx(1)+1,0:mx(2)+1). The result is stored in
! ulapresult(1:mx(1),1:mx(2)).
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: a
REAL(KIND=r8), DIMENSION(1:SIZE(a,1)-2, &
                         1:SIZE(a,2)-2):: ulapresult
!
INTEGER, DIMENSION(1:2):: mx
!
mx(1) = SIZE(a,1)-2
mx(2) = SIZE(a,2)-2
!
ulapresult(1:mx(1),1:mx(2)) =        a(2:mx(1)+1,1:mx(2)  ) &
                            +        a(0:mx(1)-1,1:mx(2)  ) &
                            +        a(1:mx(1)  ,2:mx(2)+1) &
                            +        a(1:mx(1)  ,0:mx(2)-1) &
                            - 4.0_r8*a(1:mx(1)  ,1:mx(2)  )
!
END FUNCTION ULap2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION UDiv2D(f1,f2) RESULT(udivresult)
USE NodeInfoDef
IMPLICIT NONE
!
! UNDIVIDED divergence of the 2D flux function 
!
!  f = (f1(0:mx(1),1:mx(2)),f2(1:mx(1),0:mx(2))).  
!
! The result is stored in udivresult(1:mx(1),1:mx(2)).
!
REAL(KIND=r8), DIMENSION(0:,1:), INTENT(IN):: f1
REAL(KIND=r8), DIMENSION(1:,0:), INTENT(IN):: f2
REAL(KIND=r8), DIMENSION(1:SIZE(f2,1), &
                         1:SIZE(f1,2)):: udivresult
!
INTEGER, DIMENSION(1:2):: mx
!
mx(1) = SIZE(f2,1)
mx(2) = SIZE(f1,2)
!
udivresult(1:mx(1),1:mx(2)) = f1(1:mx(1)  ,1:mx(2)  ) &
                            - f1(0:mx(1)-1,1:mx(2)  ) &
                            + f2(1:mx(1)  ,1:mx(2)  ) &
                            - f2(1:mx(1)  ,0:mx(2)-1)
!
END FUNCTION UDiv2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION Restriction2D(a) RESULT(restrictionresult)
USE NodeInfoDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(:,:,:), INTENT(IN):: a
REAL(KIND=r8), DIMENSION(1:SIZE(a,1)/2, &
                         1:SIZE(a,2)/2, &
                         1:SIZE(a,3)):: restrictionresult
!
INTEGER:: nrvars
INTEGER, DIMENSION(1:2):: cmx, mx
!
mx(1)  = SIZE(a,1)
mx(2)  = SIZE(a,2)
nrvars = SIZE(a,3)
cmx(1:2) = mx(1:2)/2
!
restrictionresult(1:cmx(1),1:cmx(2),1:nrvars) &
  = 0.25_r8*(a(2:mx(1)  :2,2:mx(2)  :2,1:nrvars) &
  +          a(1:mx(1)-1:2,2:mx(2)  :2,1:nrvars) &
  +          a(2:mx(1)  :2,1:mx(2)-1:2,1:nrvars) &
  +          a(1:mx(1)-1:2,1:mx(2)-1:2,1:nrvars))
!
END FUNCTION Restriction2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION Prolongation2D(a) RESULT(presult)
USE NodeInfoDef
IMPLICIT NONE
!
! This prolongation algorithm assumes no ghost layers are present.
!
REAL(KIND=r8), DIMENSION(:,:,:), INTENT(IN):: a
REAL(KIND=r8), DIMENSION(1:SIZE(a,1)*2, &
                         1:SIZE(a,2)*2, &
                         1:SIZE(a,3)):: presult
!
INTEGER:: nrvars
INTEGER, DIMENSION(1:2):: cmx, mx
!
cmx(1) = SIZE(a,1)
cmx(2) = SIZE(a,2)
nrvars = SIZE(a,3)
mx(1:2) = cmx(1:2)*2
!
presult(2:mx(1)  :2,2:mx(2)  :2,1:nrvars) = a(1:cmx(1),1:cmx(2),1:nrvars)
presult(1:mx(1)-1:2,2:mx(2)  :2,1:nrvars) = a(1:cmx(1),1:cmx(2),1:nrvars)
presult(2:mx(1)  :2,1:mx(2)-1:2,1:nrvars) = a(1:cmx(1),1:cmx(2),1:nrvars)
presult(1:mx(1)-1:2,1:mx(2)-1:2,1:nrvars) = a(1:cmx(1),1:cmx(2),1:nrvars)
!
END FUNCTION Prolongation2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION BiLinProlongationP1(a) RESULT(presult)
USE NodeInfoDef
IMPLICIT NONE
!
! This bilinear prolongation algorithm assumes exactly one ghost layer, i.e., 
! mbc = 1.
!
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: a
REAL(KIND=r8), DIMENSION(0:(SIZE(a,1)-2)*2+1, &
                         0:(SIZE(a,2)-2)*2+1, &
                         1: SIZE(a,3)       ):: presult
!
INTEGER:: nrvars
INTEGER, DIMENSION(1:2):: cmx, mx
REAL(KIND=r8), DIMENSION(0:(SIZE(a,1)-2)*2+1, &
                         0: SIZE(a,2)-2   +1, &
                         1: SIZE(a,3)       ):: b
!
cmx(1) = SIZE(a,1)-2
cmx(2) = SIZE(a,2)-2
nrvars = SIZE(a,3)
mx(1:2) = cmx(1:2)*2
!
! Linear interpolation in the x-direction first:
            b(0: mx(1)  :2,0:cmx(2)+1  ,1:nrvars) &
  =  3.0_r8*a(0:cmx(1)    ,0:cmx(2)+1  ,1:nrvars) &
  +         a(1:cmx(1)+1  ,0:cmx(2)+1  ,1:nrvars)
            b(1: mx(1)+1:2,0:cmx(2)+1  ,1:nrvars) &
  =         a(0:cmx(1)    ,0:cmx(2)+1  ,1:nrvars) &
  +  3.0_r8*a(1:cmx(1)+1  ,0:cmx(2)+1  ,1:nrvars)
!
! Linear interpolation in the y-direction second:
      presult(0: mx(1)+1  ,0: mx(2)  :2,1:nrvars) &
  = (3.0_r8*b(0: mx(1)+1  ,0:cmx(2)    ,1:nrvars) &
  +         b(0: mx(1)+1  ,1:cmx(2)+1  ,1:nrvars))/16.0_r8
      presult(0: mx(1)+1  ,1: mx(2)+1:2,1:nrvars) &
  = (       b(0: mx(1)+1  ,0:cmx(2)    ,1:nrvars) &
  +  3.0_r8*b(0: mx(1)+1  ,1:cmx(2)+1  ,1:nrvars))/16.0_r8
!
END FUNCTION BiLinProlongationP1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION BiLinProlongationP2(a) RESULT(presult)
USE NodeInfoDef
IMPLICIT NONE
!
! This bilinear prolongation algorithm assumes exactly two ghost layers, i.e., 
! mbc = 2.
!
REAL(KIND=r8), DIMENSION(-1:,-1:,1:), INTENT(IN):: a
REAL(KIND=r8), DIMENSION(-1:(SIZE(a,1)-4)*2+2, &
                         -1:(SIZE(a,2)-4)*2+2, &
                          1: SIZE(a,3)       ):: presult
!
INTEGER:: nrvars
INTEGER, DIMENSION(1:2):: cmx, mx
REAL(KIND=r8), DIMENSION(-1:(SIZE(a,1)-4)*2+2, &
                         -1: SIZE(a,2)-4   +2, &
                          1: SIZE(a,3)       ):: b
!
cmx(1) = SIZE(a,1)-4
cmx(2) = SIZE(a,2)-4
nrvars = SIZE(a,3)
mx(1:2) = cmx(1:2)*2
!
! Linear interpolation in the x-direction first:
            b(-1: mx(1)+1:2,-1:cmx(2)+2  ,1:nrvars) &
  =  3.0_r8*a( 0:cmx(1)+1  ,-1:cmx(2)+2  ,1:nrvars) &
  +         a(-1:cmx(1)    ,-1:cmx(2)+2  ,1:nrvars)
!
            b( 0: mx(1)+2:2,-1:cmx(2)+2  ,1:nrvars) &
  =  3.0_r8*a( 0:cmx(1)+1  ,-1:cmx(2)+2  ,1:nrvars) &
  +         a( 1:cmx(1)+2  ,-1:cmx(2)+2  ,1:nrvars)
!
! Linear interpolation in the y-direction second:
      presult(-1: mx(1)+2  ,-1: mx(2)+1:2,1:nrvars) &
  = (3.0_r8*b(-1: mx(1)+2  , 0:cmx(2)+1  ,1:nrvars) &
  +         b(-1: mx(1)+2  ,-1:cmx(2)    ,1:nrvars))/16.0_r8
!  
      presult(-1: mx(1)+2  , 0: mx(2)+2:2,1:nrvars) &
  = (3.0_r8*b(-1: mx(1)+2  , 0:cmx(2)+1  ,1:nrvars) &
  +         b(-1: mx(1)+2  , 1:cmx(2)+2  ,1:nrvars))/16.0_r8
!
END FUNCTION BiLinProlongationP2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION BiLinProlongationP1MC(a) RESULT(presult)
USE NodeInfoDef
IMPLICIT NONE
!
! This mass-corrected bilinear prolongation algorithm assumes exactly one ghost
! layer, i.e., mbc = 1.
!
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: a
REAL(KIND=r8), DIMENSION(0:(SIZE(a,1)-2)*2+1, &
                         0:(SIZE(a,2)-2)*2+1, &
                         1: SIZE(a,3)       ):: presult
!
INTEGER:: nrvars
INTEGER, DIMENSION(1:2):: cmx, mx
REAL(KIND=r8), DIMENSION(0:(SIZE(a,1)-2)*2+1, &
                         0: SIZE(a,2)-2   +1, &
                         1: SIZE(a,3)       ):: b
REAL(KIND=r8), DIMENSION(1: SIZE(a,1)-2     , &
                         1: SIZE(a,2)-2     , &
                         1: SIZE(a,3)       ):: cor
!
cmx(1) = SIZE(a,1)-2
cmx(2) = SIZE(a,2)-2
nrvars = SIZE(a,3)
mx(1:2) = cmx(1:2)*2
!
! Linear interpolation in the x-direction first:
            b(0: mx(1)  :2,0:cmx(2)+1  ,1:nrvars) &
  =  3.0_r8*a(0:cmx(1)    ,0:cmx(2)+1  ,1:nrvars) &
  +         a(1:cmx(1)+1  ,0:cmx(2)+1  ,1:nrvars)
            b(1: mx(1)+1:2,0:cmx(2)+1  ,1:nrvars) &
  =         a(0:cmx(1)    ,0:cmx(2)+1  ,1:nrvars) &
  +  3.0_r8*a(1:cmx(1)+1  ,0:cmx(2)+1  ,1:nrvars)
!
! Linear interpolation in the y-direction second:
      presult(0: mx(1)+1  ,0: mx(2)  :2,1:nrvars) &
  = (3.0_r8*b(0: mx(1)+1  ,0:cmx(2)    ,1:nrvars) &
  +         b(0: mx(1)+1  ,1:cmx(2)+1  ,1:nrvars))/16.0_r8
      presult(0: mx(1)+1  ,1: mx(2)+1:2,1:nrvars) &
  = (       b(0: mx(1)+1  ,0:cmx(2)    ,1:nrvars) &
  +  3.0_r8*b(0: mx(1)+1  ,1:cmx(2)+1  ,1:nrvars))/16.0_r8
!
! The mass correction is computed only for regular cells; not for ghost cells.
! Compute mass correction:
  cor(1:cmx(1),1:cmx(2),1:nrvars) &
  = a(1:cmx(1),1:cmx(2),1:nrvars) &
  - 0.25_r8*(presult(2:mx(1)  :2,2:mx(2)  :2,1:nrvars) &
  +          presult(1:mx(1)-1:2,2:mx(2)  :2,1:nrvars) &
  +          presult(2:mx(1)  :2,1:mx(2)-1:2,1:nrvars) &
  +          presult(1:mx(1)-1:2,1:mx(2)-1:2,1:nrvars))
!
! Add the mass correction:
    presult(2:mx(1)  :2,2:mx(2)  :2,1:nrvars) &
  = presult(2:mx(1)  :2,2:mx(2)  :2,1:nrvars)+cor(1:cmx(1),1:cmx(2),1:nrvars)
    presult(1:mx(1)-1:2,2:mx(2)  :2,1:nrvars) &
  = presult(1:mx(1)-1:2,2:mx(2)  :2,1:nrvars)+cor(1:cmx(1),1:cmx(2),1:nrvars)
    presult(2:mx(1)  :2,1:mx(2)-1:2,1:nrvars) &
  = presult(2:mx(1)  :2,1:mx(2)-1:2,1:nrvars)+cor(1:cmx(1),1:cmx(2),1:nrvars)
    presult(1:mx(1)-1:2,1:mx(2)-1:2,1:nrvars) &
  = presult(1:mx(1)-1:2,1:mx(2)-1:2,1:nrvars)+cor(1:cmx(1),1:cmx(2),1:nrvars)
!
END FUNCTION BiLinProlongationP1MC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION BiLinProlongationP2MC(a) RESULT(presult)
USE NodeInfoDef
IMPLICIT NONE
!
! This mass-corrected bilinear prolongation algorithm assumes exactly two ghost
! layers, i.e., mbc = 2.
!
REAL(KIND=r8), DIMENSION(-1:,-1:,1:), INTENT(IN):: a
REAL(KIND=r8), DIMENSION(-1:(SIZE(a,1)-4)*2+2, &
                         -1:(SIZE(a,2)-4)*2+2, &
                          1: SIZE(a,3)       ):: presult
!
INTEGER:: nrvars
INTEGER, DIMENSION(1:2):: cmx, mx
REAL(KIND=r8), DIMENSION(-1:(SIZE(a,1)-4)*2+2, &
                         -1: SIZE(a,2)-4   +2, &
                          1: SIZE(a,3)       ):: b
REAL(KIND=r8), DIMENSION( 1: SIZE(a,1)-4     , &
                          1: SIZE(a,2)-4     , &
                          1: SIZE(a,3)       ):: cor
!
cmx(1) = SIZE(a,1)-4
cmx(2) = SIZE(a,2)-4
nrvars = SIZE(a,3)
mx(1:2) = cmx(1:2)*2
!
! Linear interpolation in the x-direction first:
            b(-1: mx(1)+1:2,-1:cmx(2)+2  ,1:nrvars) &
  =  3.0_r8*a( 0:cmx(1)+1  ,-1:cmx(2)+2  ,1:nrvars) &
  +         a(-1:cmx(1)    ,-1:cmx(2)+2  ,1:nrvars)
!
            b( 0: mx(1)+2:2,-1:cmx(2)+2  ,1:nrvars) &
  =  3.0_r8*a( 0:cmx(1)+1  ,-1:cmx(2)+2  ,1:nrvars) &
  +         a( 1:cmx(1)+2  ,-1:cmx(2)+2  ,1:nrvars)
!
! Linear interpolation in the y-direction second:
      presult(-1: mx(1)+2  ,-1: mx(2)+1:2,1:nrvars) &
  = (3.0_r8*b(-1: mx(1)+2  , 0:cmx(2)+1  ,1:nrvars) &
  +         b(-1: mx(1)+2  ,-1:cmx(2)    ,1:nrvars))/16.0_r8
!  
      presult(-1: mx(1)+2  , 0: mx(2)+2:2,1:nrvars) &
  = (3.0_r8*b(-1: mx(1)+2  , 0:cmx(2)+1  ,1:nrvars) &
  +         b(-1: mx(1)+2  , 1:cmx(2)+2  ,1:nrvars))/16.0_r8
!
! The mass correction is computed only for regular cells; not for ghost cells.
! Compute mass correction:
  cor(1:cmx(1),1:cmx(2),1:nrvars) &
  = a(1:cmx(1),1:cmx(2),1:nrvars) &
  - 0.25_r8*(presult(2:mx(1)  :2,2:mx(2)  :2,1:nrvars) &
  +          presult(1:mx(1)-1:2,2:mx(2)  :2,1:nrvars) &
  +          presult(2:mx(1)  :2,1:mx(2)-1:2,1:nrvars) &
  +          presult(1:mx(1)-1:2,1:mx(2)-1:2,1:nrvars))
!
! Add the mass correction:
    presult(2:mx(1)  :2,2:mx(2)  :2,1:nrvars) &
  = presult(2:mx(1)  :2,2:mx(2)  :2,1:nrvars)+cor(1:cmx(1),1:cmx(2),1:nrvars)
    presult(1:mx(1)-1:2,2:mx(2)  :2,1:nrvars) &
  = presult(1:mx(1)-1:2,2:mx(2)  :2,1:nrvars)+cor(1:cmx(1),1:cmx(2),1:nrvars)
    presult(2:mx(1)  :2,1:mx(2)-1:2,1:nrvars) &
  = presult(2:mx(1)  :2,1:mx(2)-1:2,1:nrvars)+cor(1:cmx(1),1:cmx(2),1:nrvars)
    presult(1:mx(1)-1:2,1:mx(2)-1:2,1:nrvars) &
  = presult(1:mx(1)-1:2,1:mx(2)-1:2,1:nrvars)+cor(1:cmx(1),1:cmx(2),1:nrvars)
!
END FUNCTION BiLinProlongationP2MC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION ULap3D(a) RESULT(ulapresult)
USE NodeInfoDef
IMPLICIT NONE
!
! 3D UNDIVIDED laplacian of a(0:mx(1)+1,0:mx(2)+1,0:mx(3)+1). The result is 
! stored in ulapresult(1:mx(1),1:mx(2)).
!
REAL(KIND=r8), DIMENSION(0:,0:,0:), INTENT(IN):: a
REAL(KIND=r8), DIMENSION(1:SIZE(a,1)-2, &
                         1:SIZE(a,2)-2, &
                         1:SIZE(a,3)-2):: ulapresult
!
INTEGER, DIMENSION(1:3):: mx
!
mx(1) = SIZE(a,1)-2
mx(2) = SIZE(a,2)-2
mx(3) = SIZE(a,3)-2
!
ulapresult(1:mx(1),1:mx(2),1:mx(3)) &
  =        a(2:mx(1)+1,1:mx(2)  ,1:mx(3)  ) &
  +        a(0:mx(1)-1,1:mx(2)  ,1:mx(3)  ) &
  +        a(1:mx(1)  ,2:mx(2)+1,1:mx(3)  ) &
  +        a(1:mx(1)  ,0:mx(2)-1,1:mx(3)  ) &
  +        a(1:mx(1)  ,1:mx(2)  ,2:mx(3)+1) &
  +        a(1:mx(1)  ,1:mx(2)  ,0:mx(3)-1) &
  - 6.0_r8*a(1:mx(1)  ,1:mx(2)  ,1:mx(3)  )
!
END FUNCTION ULap3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION UDiv3D(f1,f2,f3) RESULT(udivresult)
USE NodeInfoDef
IMPLICIT NONE
!
! UNDIVIDED divergence of the 3D flux function 
!
!  (f1(0:mx(1),1:mx(2),1:mx(3)),
!   f2(1:mx(1),0:mx(2),1:mx(3)),
!   f3(1:mx(1),1:mx(2),0:mx(3))).  
!
! The result is stored in udivresult(1:mx(1),1:mx(2),1:mx(3)).
!
REAL(KIND=r8), DIMENSION(0:,1:,1:), INTENT(IN):: f1
REAL(KIND=r8), DIMENSION(1:,0:,1:), INTENT(IN):: f2
REAL(KIND=r8), DIMENSION(1:,1:,0:), INTENT(IN):: f3
REAL(KIND=r8), DIMENSION(1:SIZE(f2,1), &
                         1:SIZE(f3,2), &
                         1:SIZE(f1,3)):: udivresult
!
INTEGER, DIMENSION(1:3):: mx
!
mx(1) = SIZE(f2,1)
mx(2) = SIZE(f3,2)
mx(3) = SIZE(f1,3)
!
udivresult(1:mx(1),1:mx(2),:mx(3)) = f1(1:mx(1)  ,1:mx(2)  ,1:mx(3)  ) &
                                   - f1(0:mx(1)-1,1:mx(2)  ,1:mx(3)  ) &
                                   + f2(1:mx(1)  ,1:mx(2)  ,1:mx(3)  ) &
                                   - f2(1:mx(1)  ,0:mx(2)-1,1:mx(3)  ) &
                                   + f3(1:mx(1)  ,1:mx(2)  ,1:mx(3)  ) &
                                   - f3(1:mx(1)  ,1:mx(2)  ,0:mx(3)-1)
!
END FUNCTION UDiv3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION Restriction3D(a) RESULT(restrictionresult)
USE NodeInfoDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(:,:,:,:), INTENT(IN):: a
REAL(KIND=r8), DIMENSION(1:SIZE(a,1)/2, &
                         1:SIZE(a,2)/2, &
                         1:SIZE(a,3)/2, &
                         1:SIZE(a,4)):: restrictionresult
!
INTEGER:: nrvars
INTEGER, DIMENSION(1:3):: cmx, mx
!
mx(1)  = SIZE(a,1)
mx(2)  = SIZE(a,2)
mx(3)  = SIZE(a,3)
nrvars = SIZE(a,4)
cmx(1:3) = mx(1:3)/2
!
restrictionresult(1:cmx(1),1:cmx(2),1:cmx(3),1:nrvars) &
  = 0.125_r8*(a(2:mx(1)  :2,2:mx(2)  :2,2:mx(3)  :2,1:nrvars) &
  +           a(1:mx(1)-1:2,2:mx(2)  :2,2:mx(3)  :2,1:nrvars) &
  +           a(2:mx(1)  :2,1:mx(2)-1:2,2:mx(3)  :2,1:nrvars) &
  +           a(1:mx(1)-1:2,1:mx(2)-1:2,2:mx(3)  :2,1:nrvars) &
  +           a(2:mx(1)  :2,2:mx(2)  :2,1:mx(3)-1:2,1:nrvars) &
  +           a(1:mx(1)-1:2,2:mx(2)  :2,1:mx(3)-1:2,1:nrvars) &
  +           a(2:mx(1)  :2,1:mx(2)-1:2,1:mx(3)-1:2,1:nrvars) &
  +           a(1:mx(1)-1:2,1:mx(2)-1:2,1:mx(3)-1:2,1:nrvars))
!
END FUNCTION Restriction3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION Prolongation3D(a) RESULT(presult)
USE NodeInfoDef
IMPLICIT NONE
!
! This prolongation algorithm assumes no ghost layers are present.
!
REAL(KIND=r8), DIMENSION(:,:,:,:), INTENT(IN):: a
REAL(KIND=r8), DIMENSION(1:SIZE(a,1)*2, &
                         1:SIZE(a,2)*2, &
                         1:SIZE(a,3)*2, &
                         1:SIZE(a,4)):: presult
!
INTEGER:: nrvars
INTEGER, DIMENSION(1:3):: cmx, mx
!
cmx(1) = SIZE(a,1)
cmx(2) = SIZE(a,2)
cmx(3) = SIZE(a,3)
nrvars = SIZE(a,4)
mx(1:3) = cmx(1:3)*2
!
presult(2: mx(1)  :2,2: mx(2)  :2,2: mx(3)  :2,1:nrvars) &
  =   a(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
presult(1: mx(1)-1:2,2: mx(2)  :2,2: mx(3)  :2,1:nrvars) &
  =   a(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
presult(2: mx(1)  :2,1: mx(2)-1:2,2: mx(3)  :2,1:nrvars) &
  =   a(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
presult(1: mx(1)-1:2,1: mx(2)-1:2,2: mx(3)  :2,1:nrvars) &
  =   a(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
presult(2: mx(1)  :2,2: mx(2)  :2,1: mx(3)-1:2,1:nrvars) &
  =   a(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
presult(1: mx(1)-1:2,2: mx(2)  :2,1: mx(3)-1:2,1:nrvars) &
  =   a(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
presult(2: mx(1)  :2,1: mx(2)-1:2,1: mx(3)-1:2,1:nrvars) &
  =   a(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
presult(1: mx(1)-1:2,1: mx(2)-1:2,1: mx(3)-1:2,1:nrvars) &
  =   a(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
!
END FUNCTION Prolongation3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION TriLinProlongationP1(a) RESULT(presult)
USE NodeInfoDef
IMPLICIT NONE
!
! This trilinear prolongation algorithm assumes exactly one ghost layer, i.e., 
! mbc = 1.
!
REAL(KIND=r8), DIMENSION(0:,0:,0:,1:), INTENT(IN):: a
REAL(KIND=r8), DIMENSION(0:(SIZE(a,1)-2)*2+1, &
                         0:(SIZE(a,2)-2)*2+1, &
                         0:(SIZE(a,3)-2)*2+1, &
                         1: SIZE(a,4)):: presult
!
INTEGER:: nrvars
INTEGER, DIMENSION(1:3):: cmx, mx
REAL(KIND=r8), DIMENSION(0:(SIZE(a,1)-2)*2+1, &
                         0: SIZE(a,2)-2   +1, &
                         0: SIZE(a,3)-2   +1, &
                         1: SIZE(a,4)       ):: b
REAL(KIND=r8), DIMENSION(0:(SIZE(a,1)-2)*2+1, &
                         0:(SIZE(a,2)-2)*2+1, &
                         0: SIZE(a,3)-2   +1, &
                         1: SIZE(a,4)       ):: c
!
cmx(1) = SIZE(a,1)-2
cmx(2) = SIZE(a,2)-2
cmx(3) = SIZE(a,3)-2
nrvars = SIZE(a,4)
mx(1:3) = cmx(1:3)*2
!
! Later, use presult in place of b to save storage space:
!
! Linear interpolation in the x-direction first:
b(0:mx(1)  :2,0:cmx(2)+1,0:cmx(3)+1,1:nrvars) &
  = 3.0_r8*a(0:cmx(1)  ,0:cmx(2)+1,0:cmx(3)+1,1:nrvars) &
  +        a(1:cmx(1)+1,0:cmx(2)+1,0:cmx(3)+1,1:nrvars)
b(1:mx(1)+1:2,0:cmx(2)+1,0:cmx(3)+1,1:nrvars) &
  =        a(0:cmx(1)  ,0:cmx(2)+1,0:cmx(3)+1,1:nrvars) &
  + 3.0_r8*a(1:cmx(1)+1,0:cmx(2)+1,0:cmx(3)+1,1:nrvars)
!
! Linear interpolation in the y-direction second:
c(0:mx(1)+1,0:mx(2)  :2,0:cmx(3)+1,1:nrvars) &
  = 3.0_r8*b(0:mx(1)+1,0:cmx(2)  ,0:cmx(3)+1,1:nrvars) &
  +        b(0:mx(1)+1,1:cmx(2)+1,0:cmx(3)+1,1:nrvars)
c(0:mx(1)+1,1:mx(2)+1:2,0:cmx(3)+1,1:nrvars) &
  =        b(0:mx(1)+1,0:cmx(2)  ,0:cmx(3)+1,1:nrvars) &
  + 3.0_r8*b(0:mx(1)+1,1:cmx(2)+1,0:cmx(3)+1,1:nrvars)
!
! Linear interpolation in the z-direction third:
presult(0:mx(1)+1,0:mx(2)+1,0:mx(3)  :2,1:nrvars) &
  = (3.0_r8*c(0:mx(1)+1,0:mx(2)+1,0:cmx(3)  ,1:nrvars) &
  +         c(0:mx(1)+1,0:mx(2)+1,1:cmx(3)+1,1:nrvars))/64.0_r8
presult(0:mx(1)+1,0:mx(2)+1,1:mx(3)+1:2,1:nrvars) &
  = (       c(0:mx(1)+1,0:mx(2)+1,0:cmx(3)  ,1:nrvars) &
  +  3.0_r8*c(0:mx(1)+1,0:mx(2)+1,1:cmx(3)+1,1:nrvars))/64.0_r8
!
END FUNCTION TriLinProlongationP1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION TriLinProlongationP2(a) RESULT(presult)
USE NodeInfoDef
IMPLICIT NONE
!
! This trilinear prolongation algorithm assumes exactly two ghost layers, i.e., 
! mbc = 2.
!
REAL(KIND=r8), DIMENSION(-1:,-1:,-1:,1:), INTENT(IN):: a
REAL(KIND=r8), DIMENSION(-1:(SIZE(a,1)-4)*2+2, &
                         -1:(SIZE(a,2)-4)*2+2, &
                         -1:(SIZE(a,3)-4)*2+2, &
                          1: SIZE(a,4)       ):: presult
!
INTEGER:: nrvars
INTEGER, DIMENSION(1:3):: cmx, mx
REAL(KIND=r8), DIMENSION(-1:(SIZE(a,1)-4)*2+2, &
                         -1: SIZE(a,2)-4   +2, &
                         -1: SIZE(a,3)-4   +2, &
                          1: SIZE(a,4)       ):: b
REAL(KIND=r8), DIMENSION(-1:(SIZE(a,1)-4)*2+2, &
                         -1:(SIZE(a,2)-4)*2+2, &
                         -1: SIZE(a,3)-4   +2, &
                          1: SIZE(a,4)       ):: c
!
cmx(1) = SIZE(a,1)-4
cmx(2) = SIZE(a,2)-4
cmx(3) = SIZE(a,3)-4
nrvars = SIZE(a,4)
mx(1:3) = cmx(1:3)*2
!
! Linear interpolation in the x-direction first:
           b(-1: mx(1)+1:2,-1:cmx(2)+2  ,-1:cmx(3)+2,1:nrvars) &
  = 3.0_r8*a( 0:cmx(1)+1  ,-1:cmx(2)+2  ,-1:cmx(3)+2,1:nrvars) &
  +        a(-1:cmx(1)    ,-1:cmx(2)+2  ,-1:cmx(3)+2,1:nrvars)
!
           b( 0: mx(1)+2:2,-1:cmx(2)+2  ,-1:cmx(3)+2,1:nrvars) &
  = 3.0_r8*a( 0:cmx(1)+1  ,-1:cmx(2)+2  ,-1:cmx(3)+2,1:nrvars) &
  +        a( 1:cmx(1)+2  ,-1:cmx(2)+2  ,-1:cmx(3)+2,1:nrvars)
!
! Linear interpolation in the y-direction second:
           c(-1: mx(1)+2  ,-1: mx(2)+1:2,-1:cmx(3)+2,1:nrvars) &
  = 3.0_r8*b(-1: mx(1)+2  , 0:cmx(2)+1  ,-1:cmx(3)+2,1:nrvars) &
  +        b(-1: mx(1)+2  ,-1:cmx(2)    ,-1:cmx(3)+2,1:nrvars)
!  
           c(-1: mx(1)+2  , 0: mx(2)+2:2,-1:cmx(3)+2,1:nrvars) &
  = 3.0_r8*b(-1: mx(1)+2  , 0:cmx(2)+1  ,-1:cmx(3)+2,1:nrvars) &
  +        b(-1: mx(1)+2  , 1:cmx(2)+2  ,-1:cmx(3)+2,1:nrvars)
!
! Linear interpolation in the z-direction second:
      presult(-1:mx(1)+2,-1:mx(2)+2,-1: mx(3)+1:2,1:nrvars) &
  = (3.0_r8*c(-1:mx(1)+2,-1:mx(2)+2, 0:cmx(3)+1  ,1:nrvars) &
  +         c(-1:mx(1)+2,-1:mx(2)+2,-1:cmx(3)    ,1:nrvars))/64.0_r8
      presult(-1:mx(1)+2,-1:mx(2)+2, 0: mx(3)+2:2,1:nrvars) &
  = (3.0_r8*c(-1:mx(1)+2,-1:mx(2)+2, 0:cmx(3)+1  ,1:nrvars) &
  +         c(-1:mx(1)+2,-1:mx(2)+2, 1:cmx(3)+2  ,1:nrvars))/64.0_r8
!
END FUNCTION TriLinProlongationP2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION TriLinProlongationP1MC(a) RESULT(presult)
USE NodeInfoDef
IMPLICIT NONE
!
! This trilinear prolongation algorithm assumes exactly one ghost layer, i.e., 
! mbc = 1.
!
REAL(KIND=r8), DIMENSION(0:,0:,0:,1:), INTENT(IN):: a
REAL(KIND=r8), DIMENSION(0:(SIZE(a,1)-2)*2+1, &
                         0:(SIZE(a,2)-2)*2+1, &
                         0:(SIZE(a,3)-2)*2+1, &
                         1: SIZE(a,4)):: presult
!
INTEGER:: nrvars
INTEGER, DIMENSION(1:3):: cmx, mx
REAL(KIND=r8), DIMENSION(0:(SIZE(a,1)-2)*2+1, &
                         0: SIZE(a,2)-2   +1, &
                         0: SIZE(a,3)-2   +1, &
                         1: SIZE(a,4)       ):: b
REAL(KIND=r8), DIMENSION(0:(SIZE(a,1)-2)*2+1, &
                         0:(SIZE(a,2)-2)*2+1, &
                         0: SIZE(a,3)-2   +1, &
                         1: SIZE(a,4)       ):: c
REAL(KIND=r8), DIMENSION(1: SIZE(a,1)-2     , &
                         1: SIZE(a,2)-2     , &
                         1: SIZE(a,3)-2     , &
                         1: SIZE(a,4)       ):: cor
!
cmx(1) = SIZE(a,1)-2
cmx(2) = SIZE(a,2)-2
cmx(3) = SIZE(a,3)-2
nrvars = SIZE(a,4)
mx(1:3) = cmx(1:3)*2
!
! Later, use presult in place of b to save storage space:
!
! Linear interpolation in the x-direction first:
b(0:mx(1)  :2,0:cmx(2)+1,0:cmx(3)+1,1:nrvars) &
  = 3.0_r8*a(0:cmx(1)  ,0:cmx(2)+1,0:cmx(3)+1,1:nrvars) &
  +        a(1:cmx(1)+1,0:cmx(2)+1,0:cmx(3)+1,1:nrvars)
b(1:mx(1)+1:2,0:cmx(2)+1,0:cmx(3)+1,1:nrvars) &
  =        a(0:cmx(1)  ,0:cmx(2)+1,0:cmx(3)+1,1:nrvars) &
  + 3.0_r8*a(1:cmx(1)+1,0:cmx(2)+1,0:cmx(3)+1,1:nrvars)
!
! Linear interpolation in the y-direction second:
c(0:mx(1)+1,0:mx(2)  :2,0:cmx(3)+1,1:nrvars) &
  = 3.0_r8*b(0:mx(1)+1,0:cmx(2)  ,0:cmx(3)+1,1:nrvars) &
  +        b(0:mx(1)+1,1:cmx(2)+1,0:cmx(3)+1,1:nrvars)
c(0:mx(1)+1,1:mx(2)+1:2,0:cmx(3)+1,1:nrvars) &
  =        b(0:mx(1)+1,0:cmx(2)  ,0:cmx(3)+1,1:nrvars) &
  + 3.0_r8*b(0:mx(1)+1,1:cmx(2)+1,0:cmx(3)+1,1:nrvars)
!
! Linear interpolation in the z-direction third:
presult(0:mx(1)+1,0:mx(2)+1,0:mx(3)  :2,1:nrvars) &
  = (3.0_r8*c(0:mx(1)+1,0:mx(2)+1,0:cmx(3)  ,1:nrvars) &
  +         c(0:mx(1)+1,0:mx(2)+1,1:cmx(3)+1,1:nrvars))/64.0_r8
presult(0:mx(1)+1,0:mx(2)+1,1:mx(3)+1:2,1:nrvars) &
  = (       c(0:mx(1)+1,0:mx(2)+1,0:cmx(3)  ,1:nrvars) &
  +  3.0_r8*c(0:mx(1)+1,0:mx(2)+1,1:cmx(3)+1,1:nrvars))/64.0_r8
!
! The mass correction is computed only for regular cells; not for ghost cells.
! Compute mass correction:
  cor(1:cmx(1),1:cmx(2),1:cmx(3),1:nrvars) &
  = a(1:cmx(1),1:cmx(2),1:cmx(3),1:nrvars) &
  - 0.125_r8*(presult(2:mx(1)  :2,2:mx(2)  :2,2:mx(3)  :2,1:nrvars) &
  +           presult(1:mx(1)-1:2,2:mx(2)  :2,2:mx(3)  :2,1:nrvars) &
  +           presult(2:mx(1)  :2,1:mx(2)-1:2,2:mx(3)  :2,1:nrvars) &
  +           presult(1:mx(1)-1:2,1:mx(2)-1:2,2:mx(3)  :2,1:nrvars) &
  +           presult(2:mx(1)  :2,2:mx(2)  :2,1:mx(3)-1:2,1:nrvars) &
  +           presult(1:mx(1)-1:2,2:mx(2)  :2,1:mx(3)-1:2,1:nrvars) &
  +           presult(2:mx(1)  :2,1:mx(2)-1:2,1:mx(3)-1:2,1:nrvars) &
  +           presult(1:mx(1)-1:2,1:mx(2)-1:2,1:mx(3)-1:2,1:nrvars))
!
! Add the mass correction:
    presult(2: mx(1)  :2,2: mx(2)  :2,2: mx(3)  :2,1:nrvars) &
  = presult(2: mx(1)  :2,2: mx(2)  :2,2: mx(3)  :2,1:nrvars) &
  +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
    presult(1: mx(1)-1:2,2: mx(2)  :2,2: mx(3)  :2,1:nrvars) &
  = presult(1: mx(1)-1:2,2: mx(2)  :2,2: mx(3)  :2,1:nrvars) &
  +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
    presult(2: mx(1)  :2,1: mx(2)-1:2,2: mx(3)  :2,1:nrvars) &
  = presult(2: mx(1)  :2,1: mx(2)-1:2,2: mx(3)  :2,1:nrvars) &
  +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
    presult(1: mx(1)-1:2,1: mx(2)-1:2,2: mx(3)  :2,1:nrvars) &
  = presult(1: mx(1)-1:2,1: mx(2)-1:2,2: mx(3)  :2,1:nrvars) &
  +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
    presult(2: mx(1)  :2,2: mx(2)  :2,1: mx(3)-1:2,1:nrvars) &
  = presult(2: mx(1)  :2,2: mx(2)  :2,1: mx(3)-1:2,1:nrvars) &
  +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
    presult(1: mx(1)-1:2,2: mx(2)  :2,1: mx(3)-1:2,1:nrvars) &
  = presult(1: mx(1)-1:2,2: mx(2)  :2,1: mx(3)-1:2,1:nrvars) &
  +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
    presult(2: mx(1)  :2,1: mx(2)-1:2,1: mx(3)-1:2,1:nrvars) &
  = presult(2: mx(1)  :2,1: mx(2)-1:2,1: mx(3)-1:2,1:nrvars) &
  +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
    presult(1: mx(1)-1:2,1: mx(2)-1:2,1: mx(3)-1:2,1:nrvars) &
  = presult(1: mx(1)-1:2,1: mx(2)-1:2,1: mx(3)-1:2,1:nrvars) &
  +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
!
END FUNCTION TriLinProlongationP1MC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION TriLinProlongationP2MC(a) RESULT(presult)
USE NodeInfoDef
IMPLICIT NONE
!
! This trilinear prolongation algorithm assumes exactly two ghost layers, i.e., 
! mbc = 2.
!
REAL(KIND=r8), DIMENSION(-1:,-1:,-1:,1:), INTENT(IN):: a
REAL(KIND=r8), DIMENSION(-1:(SIZE(a,1)-4)*2+2, &
                         -1:(SIZE(a,2)-4)*2+2, &
                         -1:(SIZE(a,3)-4)*2+2, &
                          1: SIZE(a,4)       ):: presult
!
INTEGER:: nrvars
INTEGER, DIMENSION(1:3):: cmx, mx
REAL(KIND=r8), DIMENSION(-1:(SIZE(a,1)-4)*2+2, &
                         -1: SIZE(a,2)-4   +2, &
                         -1: SIZE(a,3)-4   +2, &
                          1: SIZE(a,4)       ):: b
REAL(KIND=r8), DIMENSION(-1:(SIZE(a,1)-4)*2+2, &
                         -1:(SIZE(a,2)-4)*2+2, &
                         -1: SIZE(a,3)-4   +2, &
                          1: SIZE(a,4)       ):: c
REAL(KIND=r8), DIMENSION( 1: SIZE(a,1)-4     , &
                          1: SIZE(a,2)-4     , &
                          1: SIZE(a,3)-4     , &
                          1: SIZE(a,4)       ):: cor
!
cmx(1) = SIZE(a,1)-4
cmx(2) = SIZE(a,2)-4
cmx(3) = SIZE(a,3)-4
nrvars = SIZE(a,4)
mx(1:3) = cmx(1:3)*2
!
! Linear interpolation in the x-direction first:
           b(-1: mx(1)+1:2,-1:cmx(2)+2  ,-1:cmx(3)+2,1:nrvars) &
  = 3.0_r8*a( 0:cmx(1)+1  ,-1:cmx(2)+2  ,-1:cmx(3)+2,1:nrvars) &
  +        a(-1:cmx(1)    ,-1:cmx(2)+2  ,-1:cmx(3)+2,1:nrvars)
!
           b( 0: mx(1)+2:2,-1:cmx(2)+2  ,-1:cmx(3)+2,1:nrvars) &
  = 3.0_r8*a( 0:cmx(1)+1  ,-1:cmx(2)+2  ,-1:cmx(3)+2,1:nrvars) &
  +        a( 1:cmx(1)+2  ,-1:cmx(2)+2  ,-1:cmx(3)+2,1:nrvars)
!
! Linear interpolation in the y-direction second:
           c(-1: mx(1)+2  ,-1: mx(2)+1:2,-1:cmx(3)+2,1:nrvars) &
  = 3.0_r8*b(-1: mx(1)+2  , 0:cmx(2)+1  ,-1:cmx(3)+2,1:nrvars) &
  +        b(-1: mx(1)+2  ,-1:cmx(2)    ,-1:cmx(3)+2,1:nrvars)
!  
           c(-1: mx(1)+2  , 0: mx(2)+2:2,-1:cmx(3)+2,1:nrvars) &
  = 3.0_r8*b(-1: mx(1)+2  , 0:cmx(2)+1  ,-1:cmx(3)+2,1:nrvars) &
  +        b(-1: mx(1)+2  , 1:cmx(2)+2  ,-1:cmx(3)+2,1:nrvars)
!
! Linear interpolation in the z-direction third:
      presult(-1:mx(1)+2,-1:mx(2)+2,-1: mx(3)+1:2,1:nrvars) &
  = (3.0_r8*c(-1:mx(1)+2,-1:mx(2)+2, 0:cmx(3)+1  ,1:nrvars) &
  +         c(-1:mx(1)+2,-1:mx(2)+2,-1:cmx(3)    ,1:nrvars))/64.0_r8
      presult(-1:mx(1)+2,-1:mx(2)+2, 0: mx(3)+2:2,1:nrvars) &
  = (3.0_r8*c(-1:mx(1)+2,-1:mx(2)+2, 0:cmx(3)+1  ,1:nrvars) &
  +         c(-1:mx(1)+2,-1:mx(2)+2, 1:cmx(3)+2  ,1:nrvars))/64.0_r8
!
! The mass correction is computed only for regular cells; not for ghost cells.
! Compute mass correction:
  cor(1:cmx(1),1:cmx(2),1:cmx(3),1:nrvars) &
  = a(1:cmx(1),1:cmx(2),1:cmx(3),1:nrvars) &
  - 0.125_r8*(presult(2:mx(1)  :2,2:mx(2)  :2,2:mx(3)  :2,1:nrvars) &
  +           presult(1:mx(1)-1:2,2:mx(2)  :2,2:mx(3)  :2,1:nrvars) &
  +           presult(2:mx(1)  :2,1:mx(2)-1:2,2:mx(3)  :2,1:nrvars) &
  +           presult(1:mx(1)-1:2,1:mx(2)-1:2,2:mx(3)  :2,1:nrvars) &
  +           presult(2:mx(1)  :2,2:mx(2)  :2,1:mx(3)-1:2,1:nrvars) &
  +           presult(1:mx(1)-1:2,2:mx(2)  :2,1:mx(3)-1:2,1:nrvars) &
  +           presult(2:mx(1)  :2,1:mx(2)-1:2,1:mx(3)-1:2,1:nrvars) &
  +           presult(1:mx(1)-1:2,1:mx(2)-1:2,1:mx(3)-1:2,1:nrvars))
!
! Add the mass correction:
    presult(2: mx(1)  :2,2: mx(2)  :2,2: mx(3)  :2,1:nrvars) &
  = presult(2: mx(1)  :2,2: mx(2)  :2,2: mx(3)  :2,1:nrvars) &
  +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
    presult(1: mx(1)-1:2,2: mx(2)  :2,2: mx(3)  :2,1:nrvars) &
  = presult(1: mx(1)-1:2,2: mx(2)  :2,2: mx(3)  :2,1:nrvars) &
  +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
    presult(2: mx(1)  :2,1: mx(2)-1:2,2: mx(3)  :2,1:nrvars) &
  = presult(2: mx(1)  :2,1: mx(2)-1:2,2: mx(3)  :2,1:nrvars) &
  +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
    presult(1: mx(1)-1:2,1: mx(2)-1:2,2: mx(3)  :2,1:nrvars) &
  = presult(1: mx(1)-1:2,1: mx(2)-1:2,2: mx(3)  :2,1:nrvars) &
  +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
    presult(2: mx(1)  :2,2: mx(2)  :2,1: mx(3)-1:2,1:nrvars) &
  = presult(2: mx(1)  :2,2: mx(2)  :2,1: mx(3)-1:2,1:nrvars) &
  +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
    presult(1: mx(1)-1:2,2: mx(2)  :2,1: mx(3)-1:2,1:nrvars) &
  = presult(1: mx(1)-1:2,2: mx(2)  :2,1: mx(3)-1:2,1:nrvars) &
  +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
    presult(2: mx(1)  :2,1: mx(2)-1:2,1: mx(3)-1:2,1:nrvars) &
  = presult(2: mx(1)  :2,1: mx(2)-1:2,1: mx(3)-1:2,1:nrvars) &
  +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
    presult(1: mx(1)-1:2,1: mx(2)-1:2,1: mx(3)-1:2,1:nrvars) &
  = presult(1: mx(1)-1:2,1: mx(2)-1:2,1: mx(3)-1:2,1:nrvars) &
  +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
!
END FUNCTION TriLinProlongationP2MC
END MODULE GridUtilities