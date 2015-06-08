!==========================================================================
! BSAM 1.2: Block-Structured Adaptive Multigrid Solver
!==========================================================================
!
! WiseSoft: Innovators, Brothers.
!
! (c) Copyright Sorin Mitran, 2000
! Department of Applied Mathematics
! University of Washington
! mitran@amath.washington.edu
!
! (c) Copyright Sorin Mitran, 2002
! Department of Mathematics
! University of North Carolina at Chapel Hill
! mitran@amath.unc.edu
!
! Portions of the code
!
! (c) Copyright Steven M. Wise, 2006
! Department of Mathematics
! University of California at Irvine
! swise@math.uci.edu
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
! This software is made available "as is" without any assurance that it
! will work for your purposes. The software may in fact have defects,
! so use the software at your own risk.
!
! -----------------------------------------------------------------------
! File:             treeops.f90
! Purpose:          Tree operations module.
! Contains:
! Revision History: Ver. 1.0 Oct. 2000 Sorin Mitran
! Revision History: Ver. 1.1 May. 2007 Steven Wise
! Revision History: Ver. 1.2 Jul. 2015 Steven Wise
! -----------------------------------------------------------------------

!*************
MODULE TreeOps
!*************
  ! User definitions of data structures associated with each node
  USE NodeInfoDef
  IMPLICIT NONE
  PRIVATE  ! Everything is implicitly PRIVATE
  PUBLIC   InitForest,KillForest,AddRootLevelNode,          &
           CreateChild,KillNode,DeleteMarkedNode,           &
           ApplyOnLeaves,ApplyOnLevels,ApplyOnForest,       &
           ApplyOnLevel,ApplyOnLevelPairs,ApplyOnChildren,  &
           GetTreeOpsErrCode,GetNodeInfo,SetNodeInfo,       &
           GetCurrentNode,GetNodeNo,GetChildInfo,           &
           GetParentInfo,GetParent,GetSibling,GetChild,     &
           GetLevel,GetSiblingIndex,GetNrOfChildren,        &
           GetRootInfo,SetRootInfo,ExistLevel,              &
           PushForest,PopForest,CreateBelowSeedLevels,      &
           CurrentNodeToYoungest
!
! Error handling
  INTEGER  ErrCode      ! The most recent error code
  ! Error codes
  PUBLIC   err_OK,                                          &
           err_UndefinedNode,                               &
           err_NoParent,err_NoSibling,err_NoChild
  INTEGER, PARAMETER :: ErrPrint          = 1000
  INTEGER, PARAMETER :: err_OK            =    0
  INTEGER, PARAMETER :: err_UndefinedNode =  100
  INTEGER, PARAMETER :: err_NoParent      =  201
  INTEGER, PARAMETER :: err_NoSibling     =  202
  INTEGER, PARAMETER :: err_NoChild       =  203
  INTEGER, PARAMETER :: err_BadLevel      =  301

  INTEGER, PARAMETER :: IntTrue=0, IntFalse=-1

  !

  TYPE, PUBLIC :: Node
    PRIVATE
    TYPE (NodeInfo), POINTER :: Info
    TYPE (Node),     POINTER :: Parent
    TYPE (Node),     POINTER :: Sibling  ! Next elder sibling
    TYPE (Node),     POINTER :: Child
    TYPE (Node),     POINTER :: Neighbor ! Nodes not of this parent but
                                         ! on the same level
    INTEGER level                        ! This node's level
    INTEGER ChildNo                      ! Child's ordinal. Eldest=1, Next=2, ...
    INTEGER NrOfChildren
    INTEGER LeafDist                     ! Distance from leaf
                                         ! This is *NOT* maintained dynamically as
                                         ! nodes are created/destroyed
    INTEGER NodeNo                       ! Node number. Useful as ID in debugging
  END TYPE Node

  TYPE :: NodePointer
    TYPE (Node), POINTER :: NodePtr
  END TYPE NodePointer

  TYPE, PUBLIC :: InfoPointer
    TYPE (NodeInfo), POINTER :: InfoPtr
  END TYPE InfoPointer

  INTEGER, PARAMETER, PUBLIC :: FOREST_SEED = -999
  INTEGER, PARAMETER :: NO_CHILDREN = 0
  INTEGER, PARAMETER :: FIRST_CHILD = 1
  INTEGER, PARAMETER :: NOT_A_CHILD = 0
  INTEGER, PARAMETER :: BAD_CHILD_NO = 0
  INTEGER, PARAMETER :: ZERO_LEAF_NODE_DIST = 0

  LOGICAL, PARAMETER :: NoInfoInit = .FALSE.
  LOGICAL, PARAMETER :: PreEvalNext = .TRUE.
  LOGICAL, PARAMETER :: PostEvalNext = .FALSE.

  ! Global variables showing current state of tree traversal
  TYPE (Node), POINTER, SAVE :: ForestSeed
  TYPE (Node), POINTER, SAVE :: Root
  TYPE (Node), POINTER, SAVE :: CurrentNode
  TYPE (NodePointer), DIMENSION(-MaxDepth:MaxDepth) :: Stack
  TYPE (NodePointer), DIMENSION(-MaxDepth:MaxDepth) :: YoungestOnLevel
  TYPE (NodePointer), DIMENSION(-MaxDepth:MaxDepth) :: EldestOnLevel
  INTEGER, DIMENSION(-MaxDepth:MaxDepth) :: NrOfNodes
  INTEGER CurrentLevel
  INTEGER GlobalNodeNo

  ! Stack mechanism to allow nested tree traversals
  INTEGER, PARAMETER :: MaxForestStacks=32
  INTEGER, SAVE :: StackLevel
  TYPE (NodePointer), SAVE, DIMENSION(0:MaxForestStacks) :: ForestStack

  ! Global I/O variables
  INTEGER InputUnit,OutputUnit

  ! Global variables showing memory usage
  INTEGER, PARAMETER :: iPrec = SELECTED_INT_KIND(9)
  !INTEGER (KIND=iPrec) :: MemFree, MemAllocated

  !-------
  CONTAINS
  !-------
    ! Interface routines are defined first.
    SUBROUTINE InitForest(InfoInit)
      LOGICAL, INTENT(IN), OPTIONAL :: InfoInit
      INTEGER :: iError
      INTEGER L
      !
      ALLOCATE (Root,ForestSeed,STAT=iError)
      IF (iError /= 0) THEN
        PRINT *,"Error allocating tree/forest roots in InitForest."
        STOP
      END IF
      StackLevel=0
      ! Seed of all trees
      NULLIFY(ForestSeed%Parent); NULLIFY(ForestSeed%Sibling)
      NULLIFY(ForestSeed%Neighbor)
      ForestSeed%Child=>Root; ForestSeed%level=FOREST_SEED
      ForestSeed%ChildNo=NOT_A_CHILD; ForestSeed%NrOfChildren=1
      ForestSeed%LeafDist=ZERO_LEAF_NODE_DIST; ForestSeed%NodeNo=FOREST_SEED
      ! First tree root
      NULLIFY(Root%Sibling); NULLIFY(Root%Child); NULLIFY(Root%Neighbor)
      Root%Parent=>ForestSeed
      Root%level=rootlevel; Root%ChildNo=FIRST_CHILD
      Root%NrOfChildren=NO_CHILDREN; Root%LeafDist=ZERO_LEAF_NODE_DIST
      GlobalNodeNo=1; Root%NodeNo=GlobalNodeNo
      DO L=-MaxDepth,MaxDepth
        NULLIFY(YoungestOnLevel(L)%NodePtr)
        NULLIFY(EldestOnLevel(L)%NodePtr)
        NrOfNodes(L)=0
      END DO
      YoungestOnLevel(rootlevel)%NodePtr=>Root
      EldestOnLevel(rootlevel)%NodePtr=>Root
      NrOfNodes(rootlevel)=1
      ALLOCATE(Root%Info,STAT=iError)  ! Allocate space for NodeInfo
      IF (iError /= 0) THEN
        PRINT *,"Error allocating tree root Info in InitForest."
        STOP
      END IF
      ! Set current node pointer
      CurrentNode => Root; CurrentLevel=rootlevel
      !
      IF (.NOT. PRESENT(InfoInit)) RETURN
      IF (.NOT. InfoInit) RETURN
    END SUBROUTINE InitForest
    
    SUBROUTINE CreateBelowSeedLevels(MinLevel,InfoInit)
      INTEGER, INTENT(IN):: MinLevel
      LOGICAL, INTENT(IN), OPTIONAL :: InfoInit
      INTEGER :: iError
      !
      IF(MinLevel<-MaxDepth) THEN
        PRINT *,"Error MinLevel less than MaxDepth in CreateBelowSeedLevels."
        STOP
      END IF
      !
      ALLOCATE(ForestSeed%Info,STAT=iError)  ! Allocate space for NodeInfo
      IF (iError /= 0) THEN
        PRINT *,"Error allocating forest seed Info in CreateBelowSeedLevels."
        STOP
      END IF
      !
      ForestSeed%Level = -1
      !
      YoungestOnLevel(-1)%NodePtr => ForestSeed
      EldestOnLevel(-1)%NodePtr => ForestSeed
      NrOfNodes(-1) = 1
      
      IF(-1>MinLevel) THEN      
        ForestSeed%ChildNo = 1
        CurrentNode => ForestSeed; CurrentLevel=-1
        CALL CreateBelowSeedLevelNode(MinLevel,InfoInit)
        ! Set current node pointer back to root.
        CurrentNode => Root; CurrentLevel=rootlevel
      END IF
      !
      IF (.NOT. PRESENT(InfoInit)) RETURN
      IF (.NOT. InfoInit) RETURN
    END SUBROUTINE CreateBelowSeedLevels
    
    RECURSIVE SUBROUTINE CreateBelowSeedLevelNode(MinLevel,InfoInit)
      ! Interface declarations
      INTEGER, INTENT(IN):: MinLevel
      LOGICAL, INTENT(IN), OPTIONAL:: InfoInit
      ! Internal declarations
      TYPE (Node), POINTER:: BelowSeedNode
      INTEGER:: iError, BelowSeedLevel
      !
      ! There are only one of these grids at each level below the root.
      ! Current node is assumed to be root node or lower.
      !
      ALLOCATE(BelowSeedNode,STAT=iError)
      IF (iError /= err_OK) THEN
        PRINT *,"Error allocating BelowSeedNode in treeops::CreateBelowSeedLevelNode."
        STOP
      END IF
      !
      NULLIFY(BelowSeedNode%Parent); NULLIFY(BelowSeedNode%Sibling)
      NULLIFY(BelowSeedNode%Child);  NULLIFY(BelowSeedNode%Neighbor)
      !
      BelowSeedNode%NodeNo = FOREST_SEED
      BelowSeedLevel = CurrentNode%level-1       ! BelowSeedNode is one level up.
      !
      NrOfNodes(BelowSeedLevel) = 1
      BelowSeedNode%Child => CurrentNode         ! Parent is created for CurrentNode
      BelowSeedNode%NrOfChildren = 1             ! BelowSeedNodes have only one child.
      CurrentNode%Parent => BelowSeedNode
      !
      BelowSeedNode%ChildNo = NOT_A_CHILD
      !
      YoungestOnLevel(BelowSeedLevel)%NodePtr => BelowSeedNode  ! Start of level stack
      EldestOnLevel(BelowSeedLevel)%NodePtr => BelowSeedNode    ! End of level stack
      !
      BelowSeedNode%level = BelowSeedLevel       ! Set the below-root level counter
      BelowSeedNode%LeafDist = FOREST_SEED
      !
      ! Finished with internal tree structure maintenance
      ! Allocate space for NodeInfo
      ALLOCATE(BelowSeedNode%Info,STAT=iError)
      IF (iError /= err_OK) THEN
        PRINT *,"Error allocating BelowSeedNode Info."
        STOP
      END IF
      !
      IF(BelowSeedLevel>MinLevel) THEN
        BelowSeedNode%ChildNo = 1
        CurrentNode => BelowSeedNode; CurrentLevel = BelowSeedLevel
        CALL CreateBelowSeedLevelNode(MinLevel,InfoInit)
      END IF
      !
      IF (.NOT. PRESENT(InfoInit)) RETURN
      IF (.NOT. InfoInit) RETURN
    END SUBROUTINE CreateBelowSeedLevelNode

    SUBROUTINE KillForest
      CALL KillNode(Root)
    END SUBROUTINE KillForest

    SUBROUTINE AddRootLevelNode(InfoInit)
      LOGICAL, INTENT(IN), OPTIONAL :: InfoInit
      ! Internal declarations
      TYPE (Node), POINTER :: NewRootLevelNode
      INTEGER :: iError
      !
      ALLOCATE (NewRootLevelNode,STAT=iError)
      IF (iError /= 0) THEN
        PRINT *,"Error allocating node in AddRootLevelNode."
        STOP
      END IF
      NewRootLevelNode%level=rootlevel
      NrOfNodes(rootlevel)=NrOfNodes(rootlevel)+1
      GlobalNodeNo=GlobalNodeNo+1; NewRootLevelNode%NodeNo=GlobalNodeNo
      NULLIFY(NewRootLevelNode%Parent); NULLIFY(NewRootLevelNode%Child)
      NewRootLevelNode%NrOfChildren=NO_CHILDREN
      NewRootLevelNode%ChildNo=Root%ChildNo+1
      NewRootLevelNode%LeafDist=ZERO_LEAF_NODE_DIST
      ! The newly created node becomes the first node on this level and
      ! the start of the tree
      !  1. Set the Sibling and Neighbor of the newly created node to
      !     point to the old Root node
      NewRootLevelNode%Sibling => Root; NewRootLevelNode%Neighbor => Root
      !  2. Set the global Root to point to the newly created root level node
      Root => NewRootLevelNode; ForestSeed%Child=>Root; CurrentNode=>Root
      !  3. The newly created node starts this level's neighbor list
      YoungestOnLevel(rootlevel)%NodePtr=>NewRootLevelNode
      ALLOCATE(NewRootLevelNode%Info,STAT=iError)  ! Allocate space for NodeInfo
      IF (iError /= 0) THEN
        PRINT *,"Error allocating Info in AddRootLevelNode."
        STOP
      END IF
      IF (.NOT. PRESENT(InfoInit)) RETURN
      IF (.NOT. InfoInit) RETURN
    END SUBROUTINE AddRootLevelNode

    SUBROUTINE CreateChild(InfoInit,ReadMode)
      ! Interface declarations
      LOGICAL, INTENT(IN), OPTIONAL :: InfoInit
      LOGICAL, INTENT(IN), OPTIONAL :: ReadMode
      ! Internal declarations
      LOGICAL UpdateYoungest
      TYPE (Node), POINTER :: NewNode
      INTEGER :: iError,ThisLevel
      !
      ALLOCATE(NewNode,STAT=iError)
      IF (iError /= err_OK) THEN
        PRINT *,"Error allocating NewNode in treeops::CreateChild."
        STOP
      END IF
      NULLIFY(NewNode%Parent); NULLIFY(NewNode%Sibling)
      NULLIFY(NewNode%Child);  NULLIFY(NewNode%Neighbor)
      GlobalNodeNo=GlobalNodeNo+1; NewNode%NodeNo=GlobalNodeNo
      ThisLevel = CurrentNode%level+1   ! Child is one level down
      IF (CurrentLevel /= CurrentNode%level) THEN
        WRITE(1,*)'Internal inconsistency in level counters in treeops::CreateChild'
        WRITE(1,1000)CurrentLevel,CurrentNode%level
        1000 FORMAT('Global level=',i2,' CurrentNode%level=',i2)
        STOP
      END IF
      IF (ThisLevel > MaxDepth) THEN
        PRINT *,'Error in treeops::CreateChild. Maximum tree depth exceeded'
        STOP
      END IF
      NrOfNodes(ThisLevel)=NrOfNodes(ThisLevel)+1
      NewNode%Parent => CurrentNode         ! Child is born of Current node
      CurrentNode%NrOfChildren=CurrentNode%NrOfChildren+1
      NewNode%Sibling => CurrentNode%Child  ! Point to next elder sibling
      ! First node of this parent?
      IF (.NOT. ASSOCIATED(NewNode%Sibling)) THEN
        NewNode%ChildNo=FIRST_CHILD
      ELSE
        NewNode%ChildNo=NewNode%Sibling%ChildNo+1
      END IF
      ! Determine if we're reading nodes from a file, in which case
      ! we'll be adding elements to the end of the level list
      IF (.NOT. PRESENT(ReadMode)) THEN
        UpdateYoungest = .TRUE.
      ELSE
        IF (ReadMode) THEN
          UpdateYoungest = .FALSE.
        ELSE
          UpdateYoungest = .TRUE.
        END IF
      END IF
      ! Is this the first child created on this level?
      IF (NrOfNodes(ThisLevel)==1) THEN
        ! Yes. Start stack spanning this level
        YoungestOnLevel(ThisLevel)%NodePtr => NewNode  ! Start of level stack
        EldestOnLevel(ThisLevel)%NodePtr => NewNode    ! End of level stack
        NULLIFY(NewNode%Neighbor)   ! ApplyOnLevel will end on this node
      ELSE
        ! No. Update the Neighbor list spanning this level
        IF (UpdateYoungest) THEN
          ! We're creating a new node during program execution
          ! Set the new node's Neighbor pointer to the previous YoungestOnLevel
          NewNode%Neighbor => YoungestOnLevel(thisLevel)%NodePtr
        ELSE
          ! We're creating a new node while reading from file
          EldestOnLevel(ThisLevel)%NodePtr%Neighbor => NewNode
          NULLIFY(NewNode%Neighbor)
          EldestOnLevel(ThisLevel)%NodePtr => NewNode
        END IF
      END IF
      CurrentNode%Child => NewNode           ! Parent points to youngest child
      NULLIFY(NewNode%Child)                 ! Just born doesn't have children
      NewNode%NrOfChildren=0
      NewNode%level = ThisLevel              ! Set the child's level counter
      NewNode%LeafDist = ZERO_LEAF_NODE_DIST ! New node is a leaf
      ! Save the most recently created child on this level. Used in maintaining
      ! the Neighbor list spanning a level
      IF (UpdateYoungest) YoungestOnLevel(ThisLevel)%NodePtr => NewNode
      ! Finished with internal tree structure maintenance
      ! Allocate space for NodeInfo
      ALLOCATE(NewNode%Info,STAT=iError)
      IF (iError /= err_OK) THEN
        PRINT *,"Error allocating NewNode Info."
        STOP
      END IF
      IF (.NOT. PRESENT(InfoInit)) RETURN
      IF (.NOT. InfoInit) RETURN
    END SUBROUTINE CreateChild

    ! Kill aNode and all of its children
    RECURSIVE SUBROUTINE KillNode(aNode)
      ! Interface declarations
      TYPE (Node), POINTER :: aNode
      ! Internal declarations
      TYPE (Node), POINTER :: Child,Sibling,Prev
      LOGICAL NeighborUpdated
      INTEGER :: iError,ThisLevel
      ! If aNode has children kill those first
      Child => aNode%Child
      DO
        IF (.NOT. ASSOCIATED(Child)) EXIT
        Sibling => Child%Sibling
        Call KillNode(Child)
        Child => Sibling
      END DO
      ! Remove leaf node
      IF (ASSOCIATED(aNode%Parent)) THEN
        aNode%Parent%NrOfChildren = aNode%Parent%NrOfChildren - 1
        IF (ASSOCIATED(aNode%Parent%Child,aNode)) THEN
          ! Update parent child pointer to next eldest
          aNode%Parent%Child => aNode%Sibling
        ELSE
          ! Search for position of this node withing sibling list
          Sibling => aNode%Parent%Child
          DO
            IF (ASSOCIATED(Sibling%Sibling,aNode)) EXIT
            Sibling => Sibling%Sibling
          END DO
          ! Remove this node from the sibling list
          Sibling%Sibling => aNode%Sibling
        END IF
      END IF
      ThisLevel = aNode%level
      NrOfNodes(ThisLevel)=NrOfNodes(ThisLevel)-1
      ! Update the Neighbor list spanning this level
      NeighborUpdated=.FALSE.
      ! Was this node the last node on this level?
      IF (NrOfNodes(ThisLevel)==0) THEN
        NULLIFY(YoungestOnLevel(ThisLevel)%NodePtr)
        NULLIFY(EldestOnLevel(ThisLevel)%NodePtr)
        NeighborUpdated=.TRUE.
      END IF
      ! Was it the start of the Neighbor list?
      IF ((.NOT. NeighborUpdated) .AND. &
          (ASSOCIATED(aNode,YoungestOnLevel(ThisLevel)%NodePtr))) THEN
         ! Set the start of the list to the next node on this level
         YoungestOnLevel(ThisLevel)%NodePtr => &
            YoungestOnLevel(ThisLevel)%NodePtr%Neighbor
         NeighborUpdated=.TRUE.
      END IF
      IF (.NOT. NeighborUpdated) THEN
        ! aNode definitely has a previous element. Find it.
        Prev => YoungestOnLevel(ThisLevel)%NodePtr
        DO WHILE (.NOT. ASSOCIATED(Prev%Neighbor,aNode))
          Prev => Prev%Neighbor
        END DO
      END IF
      ! Was it the end of the Neighbor stack?
      IF ((.NOT. NeighborUpdated) .AND. &
          (ASSOCIATED(aNode,EldestOnLevel(ThisLevel)%NodePtr))) THEN
        NULLIFY(Prev%Neighbor)
        EldestOnLevel(ThisLevel)%NodePtr => Prev
        NeighborUpdated = .TRUE.
      END IF
      IF (.NOT. NeighborUpdated) THEN
        ! If we're here, aNode is neither the last nor the first in Neighbor list
        Prev%Neighbor => aNode%Neighbor ! Skip aNode in Neighbor list
      END IF
      DEALLOCATE(aNode%Info,STAT=iError)
      IF (iError /= 0) THEN
        PRINT *,"Error deallocating aNode Info."
        STOP
      END IF
      DEALLOCATE(aNode,STAT=iError)
      IF (iError /= 0) THEN
        PRINT *,"Error deallocating aNode."
        STOP
      END IF
    END SUBROUTINE KillNode
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION DeleteMarkedNode(info,dummy)
IMPLICIT NONE
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
DeleteMarkedNode = err_ok
!
IF(.NOT. info%tobedeleted) RETURN
!
CALL KillNode(currentnode)
!
END FUNCTION DeleteMarkedNode
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE CurrentNodeToYoungest(level)
IMPLICIT NONE
!
INTEGER, INTENT(IN):: level
!
currentnode => youngestonlevel(level)%nodeptr
currentlevel = level
!
END SUBROUTINE CurrentNodeToYoungest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE ApplyOnLevel(level,f,fparam)
IMPLICIT NONE
!
INTEGER:: level
INTERFACE
  INTEGER FUNCTION f(info,param)
    USE NodeInfoDef
    TYPE(nodeinfo):: info
    TYPE(funcparam):: param
  END FUNCTION f
END INTERFACE
TYPE(funcparam):: fparam
!
TYPE(node), POINTER:: nextnode
INTEGER:: ferrcode
!
IF(ABS(level)>maxdepth) THEN
  PRINT *,'Error in ApplyOnLevel. Maximum tree depth exceeded.'
  STOP
END IF
!
currentnode => youngestonlevel(level)%nodeptr
currentlevel = level
!
DO
  IF(.NOT. ASSOCIATED(currentnode)) EXIT
  nextnode => currentnode%neighbor
  ferrcode = f(currentnode%info,fparam)
  IF(ferrcode/=err_ok) THEN
    errcode = ferrcode
    IF(ferrcode<errprint) PRINT *, 'Error in ApplyOnLevel. ferrcode=', &
                                    ferrcode
    RETURN
  END IF
  currentnode => nextnode
END DO
!
END SUBROUTINE ApplyOnLevel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE ApplyOnLevelPairs(L,f,fparam)
      ! Interface declarations
      INTEGER L
      INTERFACE
        INTEGER FUNCTION f(Info1,Info2,Param)
          USE NodeInfoDef
          ! Interface declarations
          TYPE (NodeInfo) :: Info1,Info2
          TYPE (FuncParam) :: Param
        END FUNCTION f
      END INTERFACE
      TYPE (FuncParam) :: fparam
      ! Internal declarations
      TYPE (Node), POINTER :: Node1,Node2
      INTEGER fErrCode
      IF (ABS(L) > MaxDepth) THEN
        PRINT *,'Error in ApplyOnLevelPairs. Maximum tree depth exceeded.'
        STOP
      END IF
      Node1 => YoungestOnLevel(L)%NodePtr
      DO
        IF (.NOT. ASSOCIATED(Node1)) EXIT
        Node2 => Node1%Neighbor
        DO
          IF (.NOT. ASSOCIATED(Node2)) EXIT
          fErrCode=f(Node1%Info,Node2%Info,fparam)
          IF (fErrCode /= err_OK) THEN
            ErrCode = fErrCode
            IF (fErrCode<ErrPrint ) PRINT *,'Error in ApplyOnLevelPairs. fErrCode=',fErrCode
            RETURN
          END IF
          Node2=>Node2%Neighbor
        END DO
        Node1=>Node1%Neighbor
      END DO
    END SUBROUTINE ApplyOnLevelPairs

    SUBROUTINE ApplyOnChildren(f,fparam)
      INTERFACE
        INTEGER FUNCTION f(Info,Param)
          USE NodeInfoDef
          ! Interface declarations
          TYPE (NodeInfo) :: Info
          TYPE (FuncParam) :: Param
        END FUNCTION f
      END INTERFACE
      TYPE (FuncParam) :: fparam
      INTEGER fErrCode
      TYPE (Node), POINTER :: SaveCurrentNode
      ! Apply f on all children of the current node
      SaveCurrentNode=>CurrentNode
      CurrentNode=>CurrentNode%Child
      DO
        IF (.NOT. ASSOCIATED(CurrentNode)) EXIT
        fErrCode=f(CurrentNode%Info,fparam)
        IF (fErrCode /= err_OK) THEN
          ErrCode = fErrCode
          IF (fErrCode<ErrPrint) PRINT *,'Error in ApplyOnChildren. fErrCode=',fErrCode
          RETURN
        END IF
        CurrentNode=>CurrentNode%Sibling
      END DO
      CurrentNode=>SaveCurrentNode
    END SUBROUTINE ApplyOnChildren

    SUBROUTINE ApplyOnForest(f,fparam,PreEval)
    ! Apply func on all nodes within tree
      INTERFACE
        INTEGER FUNCTION f(Info,Param)
          USE NodeInfoDef
          ! Interface declarations
          TYPE (NodeInfo) :: Info
          TYPE (FuncParam) :: Param
        END FUNCTION f
      END INTERFACE
      TYPE (FuncParam) :: fparam
      LOGICAL, OPTIONAL :: PreEval
      LOGICAL EvalNextBefore_f
      ! First executable statement
      CurrentNode => Root; CurrentLevel=rootlevel
      IF (PRESENT(PreEval)) THEN
        EvalNextBefore_f=PreEval
      ELSE
        EvalNextBefore_f=.TRUE.
      END IF
      CALL ForestTraversal(Root,f,TrueCond,fparam,EvalNextBefore_f)
    END SUBROUTINE ApplyOnForest

    SUBROUTINE ApplyOnLevels(f,fparam,PreEval)
    ! Level  0: Root
    ! Level  1: 1st generation children
    ! ....
    ! Level  n: n-th generation children
    ! Level -1: Parents of leaves
    ! Level -2: Grandparents of leaves
      INTERFACE
        INTEGER FUNCTION f(Info,Param)
          USE NodeInfoDef
          ! Interface declarations
          TYPE (NodeInfo) :: Info
          TYPE (FuncParam) :: Param
        END FUNCTION f
      END INTERFACE
      TYPE (FuncParam) :: fparam
      LOGICAL, OPTIONAL :: PreEval
      LOGICAL EvalNextBefore_f
      CurrentNode => Root; CurrentLevel=rootlevel
      IF (PRESENT(PreEval)) THEN
        EvalNextBefore_f=PreEval
      ELSE
        EvalNextBefore_f=.TRUE.
      END IF
      CALL ForestTraversal(Root,f,LevelCond,fparam,EvalNextBefore_f)
    END SUBROUTINE ApplyOnLevels

    SUBROUTINE ApplyOnLeaves(f,fparam,PreEval)
    ! Leaves = youngest generation in existence
      INTERFACE
        INTEGER FUNCTION f(Info,Param)
          USE NodeInfoDef
          ! Interface declarations
          TYPE (NodeInfo) :: Info
          TYPE (FuncParam) :: Param
        END FUNCTION f
      END INTERFACE
      TYPE (FuncParam) :: fparam
      LOGICAL, OPTIONAL :: PreEval
      LOGICAL EvalNextBefore_f
      ! First executable statement
      CurrentNode => Root; CurrentLevel=rootlevel
      IF (PRESENT(PreEval)) THEN
        EvalNextBefore_f=PreEval
      ELSE
        EvalNextBefore_f=.TRUE.
      END IF
      CALL ForestTraversal(Root,f,LeafCond,fparam,EvalNextBefore_f)
    END SUBROUTINE ApplyOnLeaves

    ! Routines interior to the module

    SUBROUTINE ForestTraversal(aNode,f,cond,fparam,EvalNextBefore_f)
      ! Traverse the root level nodes
      TYPE (Node), POINTER :: aNode
      INTERFACE
        INTEGER FUNCTION f(Info,Param)
          USE NodeInfoDef
          ! Interface declarations
          TYPE (NodeInfo) :: Info
          TYPE (FuncParam) :: Param
        END FUNCTION f
        LOGICAL FUNCTION cond()
        END FUNCTION cond
      END INTERFACE
      TYPE (FuncParam) :: fparam
      LOGICAL EvalNextBefore_f
      ! Internal declarations
      TYPE (Node), POINTER :: Next
      Next => aNode
      DO
        IF (.NOT. ASSOCIATED(Next)) EXIT
        CALL TreeTraversal(Next,f,cond,fparam,EvalNextBefore_f)
        Next=>Next%Sibling
      END DO
    END SUBROUTINE ForestTraversal

    RECURSIVE SUBROUTINE TreeTraversal(aNode,f,cond,fparam,EvalNextBefore_f)
      ! The core tree function in terms of which all others are defined:
      !   1) Traverse tree
      !   2) If condition *Cond* is satisfied, apply function *Func* to node
      !      *aNode*
      !   3) Update *Level* counter
      TYPE (Node), TARGET, INTENT(IN) :: aNode
      INTERFACE
        INTEGER FUNCTION f(Info,Param)
          USE NodeInfoDef
          ! Interface declarations
          TYPE (NodeInfo) :: Info
          TYPE (FuncParam) :: Param
        END FUNCTION f
        LOGICAL FUNCTION cond()
        END FUNCTION cond
      END INTERFACE
      TYPE (FuncParam) :: fparam
      LOGICAL :: EvalNextBefore_f
      ! Internal declarations
      TYPE (Node), POINTER :: Next
      TYPE (Node), POINTER, SAVE :: aNodeSave
      INTEGER fErrCode,Applyf,level
      ! First executable statement
      IF (ErrCode /= err_OK) RETURN   ! Go up recursion upon error
      CurrentNode => aNode            ! Set current node and maintain stack
      Stack(CurrentLevel)%NodePtr => CurrentNode ! so other functions know what to work on
      IF (EvalNextBefore_f) Next => CurrentNode%Child
      ! Do work on this node if condition cond() is satisfied
      IF (cond()) THEN
        fErrCode=f(CurrentNode%Info,fparam)
        IF (fErrCode /= err_OK) THEN
          ErrCode = fErrCode
          IF (fErrCode<ErrPrint ) PRINT *,'Error on applying f. fErrCode=',fErrCode
          RETURN
        END IF
      END IF
      IF (.NOT. EvalNextBefore_f) Next => CurrentNode%Child
      ! Find next node to work on
      DO
        IF (.NOT. ASSOCIATED(Next)) THEN
          ! Reached a leaf
          CurrentLevel=MAX(CurrentLevel-1,rootlevel); level=CurrentLevel
          ! Go up tree, exit clause from DO loop
          EXIT
        ELSE
          CurrentLevel=CurrentLevel+1; level=CurrentLevel  ! Node has children; go down
          CALL TreeTraversal(Next,f,cond,fparam,EvalNextBefore_f)
          ! Reset global context from local context when returning from recursion
          CurrentNode => aNode
          Next => Next%Sibling       ! After children have been exhausted
                                     ! work on next sibling
        END IF
      END DO
    END SUBROUTINE TreeTraversal

    LOGICAL FUNCTION TrueCond()
      TrueCond = .TRUE.
    END FUNCTION TrueCond

    LOGICAL FUNCTION LevelCond()
      LevelCond = .TRUE.
    END FUNCTION LevelCond

    LOGICAL FUNCTION LeafCond()
      IF (.NOT. ASSOCIATED(CurrentNode%Child)) THEN
        LeafCond = .TRUE.
      ELSE
        LeafCond = .FALSE.
      END IF
    END FUNCTION LeafCond

    INTEGER FUNCTION SetChildNo(Info,Param)
      TYPE (NodeInfo) :: Info
      TYPE (FuncParam) :: Param
      !
      TYPE (Node), POINTER :: Youngest
      INTEGER ListNo,NrOfChildren
      IF (CurrentNode%Level==rootlevel) THEN
        CurrentNode%ChildNo=FIRST_CHILD
        SetChildNo=err_OK
        RETURN
      END IF
      IF (CurrentNode%ChildNo==BAD_CHILD_NO) THEN
        ! Number all children within this family
        Youngest => CurrentNode%Parent%Child
        ListNo=1
        NrOfChildren=CurrentNode%Parent%NrOfChildren
        Youngest%ChildNo=NrOfChildren
        DO
          IF (.NOT. ASSOCIATED(Youngest%Sibling)) EXIT
          Youngest => Youngest%Sibling
          ListNo=ListNo+1
          Youngest%ChildNo=NrOfChildren-ListNo+1
        END DO
      END IF
      SetChildNo=err_OK
    END FUNCTION SetChildNo

    INTEGER FUNCTION GetTreeOpsErrCode()
      GetTreeOpsErrCode = ErrCode
    END FUNCTION GetTreeOpsErrCode

    INTEGER FUNCTION GetNodeInfo(aNode,aNodeInfo)
      TYPE (Node), POINTER :: aNode
      TYPE (NodeInfo) :: aNodeInfo
      IF (ASSOCIATED(aNode)) THEN
        aNodeInfo = aNode%Info
        GetNodeInfo = err_OK
      ELSE
        GetNodeInfo = err_UndefinedNode
      END IF
    END FUNCTION GetNodeInfo

    INTEGER FUNCTION SetNodeInfo(aNode,aNodeInfo)
      TYPE (Node), POINTER :: aNode
      TYPE (NodeInfo) :: aNodeInfo
      IF (ASSOCIATED(aNode)) THEN
        aNode%Info = aNodeInfo
        SetNodeInfo = err_OK
      ELSE
        SetNodeInfo = err_UndefinedNode
      END IF
    END FUNCTION SetNodeInfo

    INTEGER FUNCTION GetCurrentNodeInfo(aNodeInfo)
      TYPE (NodeInfo) :: aNodeInfo
      IF (ASSOCIATED(CurrentNode)) THEN
        aNodeInfo = CurrentNode%Info
        GetCurrentNodeInfo = err_OK
      ELSE
        GetCurrentNodeInfo = err_UndefinedNode
      END IF
    END FUNCTION GetCurrentNodeInfo

    INTEGER FUNCTION GetNodeNo()
      GetNodeNo=CurrentNode%NodeNo
    END FUNCTION GetNodeNo

    INTEGER FUNCTION SetCurrentNodeInfo(aNodeInfo)
      TYPE (NodeInfo) :: aNodeInfo
      IF (ASSOCIATED(CurrentNode)) THEN
        CurrentNode%Info = aNodeInfo
        SetCurrentNodeInfo = err_OK
      ELSE
        SetCurrentNodeInfo = err_UndefinedNode
      END IF
    END FUNCTION SetCurrentNodeInfo

    INTEGER FUNCTION GetRootInfo(aNodeInfo)
      TYPE (NodeInfo), POINTER :: aNodeInfo
      IF (ASSOCIATED(Root)) THEN
        aNodeInfo => Root%Info
        GetRootInfo = err_OK
      ELSE
        GetRootInfo = err_UndefinedNode
      END IF
    END FUNCTION GetRootInfo

    INTEGER FUNCTION SetRootInfo(aNodeInfo)
      TYPE (NodeInfo) :: aNodeInfo
      IF (ASSOCIATED(Root)) THEN
        Root%Info = aNodeInfo
        SetRootInfo = err_OK
      ELSE
        SetRootInfo = err_UndefinedNode
      END IF
    END FUNCTION SetRootInfo

    INTEGER FUNCTION GetParentInfo(aNodeInfo)
      TYPE (NodeInfo), POINTER :: aNodeInfo
      IF (ASSOCIATED(CurrentNode) .AND. ASSOCIATED(CurrentNode%Parent) .AND.  &
          CurrentNode%Parent%level > FOREST_SEED) THEN
        aNodeInfo => CurrentNode%Parent%Info
        GetParentInfo = err_OK
      ELSE
        GetParentInfo = err_NoParent
      END IF
    END FUNCTION GetParentInfo

    INTEGER FUNCTION GetCurrentNode(aNode)
      TYPE (Node), POINTER :: aNode
      GetCurrentNode = err_OK
      aNode => CurrentNode
    END FUNCTION GetCurrentNode

    INTEGER FUNCTION GetParent(aNode,Parent)
      TYPE (Node), POINTER :: aNode
      TYPE (Node), POINTER :: Parent
      IF (ASSOCIATED(aNode)) THEN
        IF (ASSOCIATED(aNode%Parent)) THEN
          Parent => aNode%Parent
          GetParent = err_OK
        ELSE
          NULLIFY(Parent)
          GetParent = err_NoParent
        END IF
      ELSE
        GetParent = err_UndefinedNode
      END IF
    END FUNCTION GetParent

    INTEGER FUNCTION GetSibling(aNode,Sibling)
      TYPE (Node), POINTER :: aNode
      TYPE (Node), POINTER :: Sibling
      IF (ASSOCIATED(aNode)) THEN
        IF (ASSOCIATED(aNode%Sibling)) THEN
          Sibling => aNode%Sibling
          GetSibling = err_OK
        ELSE
          NULLIFY(Sibling)
          GetSibling = err_NoSibling
        END IF
      ELSE
        GetSibling = err_UndefinedNode
      END IF
    END FUNCTION GetSibling

    INTEGER FUNCTION GetChild(aNode,Child)
      TYPE (Node), POINTER :: aNode
      TYPE (Node), POINTER :: Child
      IF (ASSOCIATED(aNode)) THEN
        IF (ASSOCIATED(aNode%Child)) THEN
          Child => aNode%Child
          GetChild = err_OK
        ELSE
          NULLIFY(Child)
          GetChild = err_NoChild
        END IF
      ELSE
        GetChild = err_UndefinedNode
      END IF
    END FUNCTION GetChild

    INTEGER FUNCTION GetChildInfo(aNodeInfo)
      TYPE (NodeInfo), POINTER :: aNodeInfo
      IF (ASSOCIATED(CurrentNode) .AND. ASSOCIATED(CurrentNode%Child)) THEN
        aNodeInfo => CurrentNode%Child%Info
        GetChildInfo = err_OK
      ELSE
        GetChildInfo = err_NoChild
      END IF
    END FUNCTION GetChildInfo

    INTEGER FUNCTION GetLevel(ThisLevel)
      INTEGER ThisLevel
      ThisLevel = CurrentLevel
      IF (ThisLevel /= CurrentNode%level) THEN
        PRINT *,'Internal inconsistency in level counter'
        STOP
      END IF
      GetLevel = err_OK
    END FUNCTION GetLevel

    INTEGER FUNCTION GetSiblingIndex(SiblingIndex)
      INTEGER SiblingIndex
      INTEGER NrOfSiblings
      TYPE (Node), POINTER :: NextSibling,ThisSibling
      ThisSibling => CurrentNode
      SiblingIndex=1; NrOfSiblings=1
      IF (.NOT. ASSOCIATED(CurrentNode%Parent)) THEN
        ! We are on the forest seed level
        GetSiblingIndex = err_OK
        RETURN
      ELSE
        ! We are below the forest seed level
        NextSibling => CurrentNode%Parent%Child
        DO
          IF (ASSOCIATED(NextSibling,ThisSibling)) EXIT
          SiblingIndex=SiblingIndex+1
          NextSibling => NextSibling%Sibling
        END DO
        GetSiblingIndex = err_OK
        NextSibling => CurrentNode%Parent%Child
        DO
          IF (.NOT. ASSOCIATED(nextSibling%Sibling)) EXIT
          NrOfSiblings=NrOfSiblings+1
          NextSibling => NextSibling%Sibling
        END DO
      END IF
      IF ((NrOfSiblings+1-SiblingIndex) /= CurrentNode%ChildNo) THEN
        PRINT *,'Internal inconsistency in ChildNo counter'
        PRINT *,'Level=',CurrentNode%level
        PRINT *,'ChildNo=',CurrentNode%ChildNo
        PRINT *,'SiblingIndex=',SiblingIndex
        PRINT *,'NrOfSiblings=',NrOfSiblings
        STOP
      END IF
    END FUNCTION GetSiblingIndex

    INTEGER FUNCTION GetNrOfChildren(NrOfChildren)
      INTEGER NrOfChildren
      TYPE (Node), POINTER :: NextChild
      NrOfChildren=0
      NextChild => CurrentNode%Child
      DO
        IF (.NOT. ASSOCIATED(NextChild)) EXIT
        NrOfChildren=NrOfChildren+1
        NextChild => NextChild%Sibling
      END DO
      IF (NrOfChildren /= CurrentNode%NrOfChildren) THEN
        WRITE(1,1001)
        1001 FORMAT('GetNrOfChildren: Internal inconsistency in NrOfChildren')
        WRITE(1,1002)NrOfChildren,CurrentNode%NrOfChildren
        1002 FORMAT('Computed from links:',i3,' stored in Node:',i3)
      END IF
      GetNrOfChildren = err_OK
    END FUNCTION GetNrOfChildren

    LOGICAL FUNCTION ExistLevel(level)
      INTEGER level
      IF (ABS(level)>MaxDepth) THEN
        ErrCode=err_BadLevel
        PRINT *,'Bad level in call to ExistLevel'
        STOP
      END IF
      IF (ASSOCIATED(EldestOnLevel(level)%NodePtr)) THEN
        ExistLevel=.TRUE.
      ELSE
        ExistLevel=.FALSE.
      END IF
    END FUNCTION ExistLevel

    SUBROUTINE PushForest
      ! Save current tree traversal state to allow a subsidiary tree traversal
      ! to take place
      ForestStack(StackLevel)%NodePtr => CurrentNode
      StackLevel=StackLevel+1
      IF (StackLevel>MaxForestStacks) THEN
        PRINT *,'Too many forest traversal recursions. Increase MaxForestStacks in treeops.f90'
        STOP
      END IF
      CurrentNode=>Root
    END SUBROUTINE PushForest

    SUBROUTINE PopForest
      ! Restore tree traversal state on return from a subsidiary tree traversal
      StackLevel=StackLevel-1
      IF (StackLevel<0) THEN
        PRINT *,'PopForest requested without prior PushForest'
        STOP
      END IF
      CurrentNode => ForestStack(StackLevel)%NodePtr
    END SUBROUTINE PopForest

!*****************
END MODULE TreeOps
!*****************
