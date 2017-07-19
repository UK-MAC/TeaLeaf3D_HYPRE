#define EXTERNAL_FACE -1

MODULE PETScTeaLeaf

  USE definitions_module

  IMPLICIT NONE

#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"
#include "finclude/petscdm.h"
#include "finclude/petscdmda.h"
#include "finclude/petscdmda.h90"

  INTEGER :: perr
  INTEGER :: mpisize

  KSP :: kspObj
  Vec :: Sol
  Vec :: X
  Vec :: B
  Mat :: A
  Vec :: XLoc
  Vec :: RHSLoc
  DM  :: petscDA

CONTAINS

SUBROUTINE setup_petsc(eps,max_iters)

  USE data_module

#include "finclude/petscsys.h"

  INTEGER :: c,cx,cy
  REAL(kind=8) :: eps
  INTEGER :: max_iters

  CALL PetscInitialize(PETSC_NULL_CHARACTER,perr)

  ! px, py, pz set in tea_decompose to be chunks_x, chunks_y, and chunks_z
  ! tea_decompose MUST be called first

  CALL DMDACreate3D(PETSC_COMM_WORLD,                                   &
                    DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE, &
                    DMDA_STENCIL_STAR,                                  &
                    grid%x_cells,grid%y_cells,grid%z_cells,             &
                    px,py,pz,                                           &
                    1,1,                                              &
                    lx(1:px),ly(1:py),lz(1:pz),                         &
                    petscDA,perr)

  ! Setup the KSP Solver
  CALL KSPCreate(MPI_COMM_WORLD,kspObj,perr)
  CALL KSPSetTolerances(kspObj,eps,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,max_iters,perr)

  CALL KSPSetCheckNormIteration(kspObj,-1,perr)
  CALL KSPSetFromOptions(kspObj,perr)

  IF(parallel%boss) THEN
    WRITE(g_out,*)
    WRITE(g_out,*)'PETSc Setup:'
    WRITE(g_out,*)'Absolute Tolerance set to', eps
    WRITE(g_out,*)'max_iters set to', max_iters
  ENDIF

  CALL MPI_Comm_Size(MPI_COMM_WORLD,mpisize,perr)

  IF(mpisize .EQ. 1) THEN
    CALL DMSetMatType(petscDA,'seqaij',perr)
  ELSE
    CALL DMSetMatType(petscDA,'mpiaij',perr)
  ENDIF

  CALL DMCreateMatrix(petscDA,A,perr)

  ! Setup the initial solution vector
  CALL DMCreateGlobalVector(petscDA,X,perr)

  ! Duplicate Vector to setup B
  CALL DMCreateGlobalVector(petscDA,B,perr)

  ! Local Vector for RHS Vector
  CALL DMCreateLocalVector(petscDA,RHSLoc,perr)

  ! Local Vector for X Vector
  CALL DMCreateLocalVector(petscDA,XLoc,perr)

  total_petsc_iter = 0

END SUBROUTINE setup_petsc

SUBROUTINE cleanup_petsc()

#include "finclude/petscsys.h"

  CALL PetscFinalize(perr)

END SUBROUTINE cleanup_petsc

SUBROUTINE setupSol_petsc(c,rx,ry,rz)

    PetscErrorCode :: ierr
    INTEGER       :: c                                ! What chunk are we solving
    INTEGER       :: x_min,x_max,y_min,y_max,z_min,z_max
    INTEGER       :: g_xmin, g_xmax, g_ymin, g_ymax, g_zmin, g_zmax
    INTEGER       :: left,right,top,bottom,back,front
    INTEGER :: i,j,k,ilen,count
    REAL(KIND=8),pointer,dimension(:) :: rowdata
    INTEGER,pointer,dimension(:) :: rowloc
    REAL(KIND=8) :: val
    REAL(KIND=8),INTENT(IN)   :: rx,ry,rz
    PetscScalar,pointer :: xv(:,:,:)

    count = 0
    x_min = chunks(c)%field%x_min
    x_max = chunks(c)%field%x_max
    y_min = chunks(c)%field%y_min
    y_max = chunks(c)%field%y_max
    z_min = chunks(c)%field%z_min
    z_max = chunks(c)%field%z_max
    g_xmin = 1
    g_xmax = grid%x_cells
    g_ymin = 1
    g_ymax = grid%y_cells
    g_zmin = 1
    g_zmax = grid%z_cells
    left = chunks(c)%field%left
    right = chunks(c)%field%right
    top = chunks(c)%field%top
    bottom = chunks(c)%field%bottom
    back = chunks(c)%field%back
    front = chunks(c)%field%front

    ! u0 = chunks(c)%field%u

    CALL VecZeroEntries(X,perr)

    CALL DMDAVecGetArrayF90(petscDA,X,xv,perr)
    DO k = back, front
      DO j = bottom, top
        DO i = left, right
            xv(i-1,j-1,k-1) = chunks(c)%field%u((i-left)+1,(j-bottom)+1,(k-back)+1)
        ENDDO
      ENDDO
    ENDDO
    CALL DMDAVecRestoreArrayF90(petscDA,X,xv,perr)

END SUBROUTINE setupSol_petsc

SUBROUTINE setupRHS_petsc(c,rx,ry,rz)

    PetscErrorCode :: ierr
    INTEGER       :: c                                ! What chunk are we solving
    INTEGER       :: left,right,top,bottom,back,front
    INTEGER       :: x_min,x_max,y_min,y_max,z_min,z_max
    INTEGER       :: g_xmin, g_xmax, g_ymin, g_ymax, g_zmin, g_zmax
    INTEGER :: i,j,k,ilen,count,i_local,j_local,k_local
    REAL(KIND=8),pointer,dimension(:) :: rowdata
    INTEGER,pointer,dimension(:) :: rowloc
    REAL(KIND=8) :: val
    REAL(KIND=8),INTENT(IN)   :: rx,ry,rz
    PetscScalar,pointer :: bv(:,:,:)

    count = 0
    
    left = chunks(c)%field%left
    right = chunks(c)%field%right
    top = chunks(c)%field%top
    bottom = chunks(c)%field%bottom
    back = chunks(c)%field%back
    front = chunks(c)%field%front
    x_min = chunks(c)%field%x_min
    x_max = chunks(c)%field%x_max
    y_min = chunks(c)%field%y_min
    y_max = chunks(c)%field%y_max
    z_min = chunks(c)%field%z_min
    z_max = chunks(c)%field%z_max
    g_xmin = 1
    g_xmax = grid%x_cells
    g_ymin = 1
    g_ymax = grid%y_cells
    g_zmin = 1
    g_zmax = grid%z_cells

    CALL VecZeroEntries(B,perr)

    CALL DMDAVecGetArrayF90(petscDA,B,bv,perr)
    
    DO k = back, front
      DO j = bottom, top
        DO i = left, right

          i_local = (i - left)+1
          j_local = (j - bottom)+1
          k_local = (k - back)+1

          val = chunks(c)%field%u(i_local,j_local,k_local)

          bv(i-1,j-1,k-1) = val

        ENDDO
      ENDDO
    ENDDO
    
    CALL DMDAVecRestoreArrayF90(petscDA,B,bv,perr)

END SUBROUTINE setupRHS_petsc

SUBROUTINE getSolution_petsc(c)

    INTEGER :: i,j,k,ilen,rowloc
    PetscScalar :: sol_v(1)
    PetscOffset iss
    INTEGER :: count
    INTEGER       :: x_min,x_max,y_min,y_max,z_min,z_max
    INTEGER       :: c                                ! What chunk are we solving
    PetscScalar,pointer :: xv(:,:,:)
    INTEGER       :: left,right,top,bottom,back,front

    left = chunks(c)%field%left
    right = chunks(c)%field%right
    top = chunks(c)%field%top
    bottom = chunks(c)%field%bottom
    back = chunks(c)%field%back
    front = chunks(c)%field%front
    
    CALL DMDAVecGetArrayF90(petscDA,X,xv,perr)

    DO k = back, front
      DO j = bottom, top
        DO i = left, right
          chunks(c)%field%u((i-left)+1,(j-bottom)+1,(k-back+1)) = xv(i-1,j-1,k-1)
        ENDDO
      ENDDO
    ENDDO
    
    CALL DMDAVecRestoreArrayF90(petscDA,X,xv,perr)

END SUBROUTINE getSolution_petsc

!! STILL TO FIX setupMatA_petsc

SUBROUTINE setupMatA_petsc(c,rx,ry,rz)

  INTEGER       :: c                                ! What chunk are we solving
  INTEGER       :: left,right,top,bottom,back,front
  INTEGER       :: x_min,x_max,y_min,y_max,z_min,z_max
  INTEGER       :: g_xmin, g_xmax, g_ymin, g_ymax, g_zmin, g_zmax
  REAL(KIND=8),INTENT(IN)   :: rx,ry,rz
  INTEGER       :: i,j,k,i_g,j_g,k_g
  INTEGER       :: count
  REAL(KIND=8)  :: c1,c2,c3,c4,c5,c6,c7  ! 7 point stencil

  !! THIS ALL NEEDS TO CHANGE FOR 7 POINT STENCIL
  MatStencil  :: row(4,1)       ! 6 x Number of stencils in entry (adding 1 at a time for now, so 1)
  MatStencil  :: column(4,7)    ! 6 x Stencil Size: (i,j,k,m,n,l) * stencil Size (7 point stencil)
  PetscScalar :: stencil(7)     ! 6 x Stencil Size: (i,j,k,m,n,l) * stencil Size (7 point stencil)

  left = chunks(c)%field%left
  right = chunks(c)%field%right
  top = chunks(c)%field%top
  bottom = chunks(c)%field%bottom
  back = chunks(c)%field%back
  front = chunks(c)%field%front
  
  x_min = chunks(c)%field%x_min
  x_max = chunks(c)%field%x_max
  y_min = chunks(c)%field%y_min
  y_max = chunks(c)%field%y_max
  z_min = chunks(c)%field%z_min
  z_max = chunks(c)%field%z_max  

  g_xmin = 1
  g_xmax = grid%x_cells
  g_ymin = 1
  g_ymax = grid%y_cells
  g_zmin = 1
  g_zmax = grid%z_cells
  
  ! Zero Matrix
  CALL MatZeroEntries(A,perr)
  
  ! Fixed the loop number and local 
  
  DO k_g = back, front
    DO j_g = bottom,top
      DO i_g = left,right

          ! -1 is applied to all indexes to shift them to zero-index, presuming the fortran code maintains
          ! an index from 1 to max

          i = (i_g - left) + 1
          j = (j_g - bottom) + 1
          k = (k_g - back) + 1

          count = 1
          row(MatStencil_i,1) = i_g - 1
          row(MatStencil_j,1) = j_g - 1
          row(MatStencil_k,1) = k_g - 1

          ! Check these are correct!!!
          c2 = (-1*rx) * chunks(c)%field%Vector_Kx(i,j,k)
          c3 = (-1*rx) * chunks(c)%field%Vector_Kx(i+1,j,k)
          c4 = (-1*ry) * chunks(c)%field%Vector_Ky(i,j,k)
          c5 = (-1*ry) * chunks(c)%field%Vector_Ky(i,j+1,k)
          c6 = (-1*rz) * chunks(c)%field%Vector_Kz(i,j,k)
          c7 = (-1*rz) * chunks(c)%field%Vector_Kz(i,j,k+1)

          ! Dirichlet boundary
          IF(i_g == 1) THEN ! Global X Low Boundary
            c2=0
          ENDIF

          IF(i_g == grid%x_cells) THEN
            c3=0
          ENDIF

          IF(j_g == 1) THEN
            c4=0
          ENDIF
        
          IF(j_g == grid%y_cells) THEN
            c5=0
          ENDIF
          IF(k_g == 1) THEN
            c6=0
          ENDIF
        
          IF(k_g == grid%y_cells) THEN
            c7=0
          ENDIF

          c1 = (1.0_8-c2-c3-c4-c5-c6-c7)

          ! If you are not on the Dirichlet Boundary
          IF(.NOT.(k_g == 1)) THEN     ! Not back External Boundary
            stencil(count) = c6
            column(MatStencil_i,count) = (i_g-1)
            column(MatStencil_j,count) = (j_g-1) 
            column(MatStencil_k,count) = (k_g-1) - 1
            count = count + 1
          ENDIF
          IF(.NOT.(j_g == 1)) THEN     ! Not Bottom External Boundary
            stencil(count) = c4
            column(MatStencil_i,count) = (i_g-1)
            column(MatStencil_j,count) = (j_g-1) - 1
            column(MatStencil_k,count) = (k_g-1) 
            count = count + 1
          ENDIF

          IF(.NOT.(i_g == 1)) THEN     ! Not Left External Boundary
            stencil(count) = c2
            column(MatStencil_i,count) = (i_g-1) - 1
            column(MatStencil_j,count) = j_g - 1
            column(MatStencil_k,count) = (k_g-1)
            count = count + 1
          ENDIF

          ! Cell Centered Value
          stencil(count) = c1
          column(MatStencil_i,count) = i_g - 1
          column(MatStencil_j,count) = j_g - 1
          column(MatStencil_k,count) = k_g - 1
          count = count + 1

          ! If you are not on the Dirichlet Boundary
          IF(.NOT.(i_g == grid%x_cells)) THEN     ! Not Right External Boundary
            stencil(count) = c3
            column(MatStencil_i,count) = (i_g+1) - 1
            column(MatStencil_j,count) = j_g - 1
            column(MatStencil_k,count) = k_g - 1
            count = count + 1
          ENDIF

          IF(.NOT.(j_g == grid%y_cells)) THEN     ! Not Top External Boundary
            stencil(count) = c5
            column(MatStencil_i,count) = i_g - 1
            column(MatStencil_j,count) = (j_g+1) - 1
            column(MatStencil_k,count) = k_g - 1
            count = count + 1
          ENDIF
          
          IF(.NOT.(k_g == grid%z_cells)) THEN     ! Not Front External Boundary
            stencil(count) = c7
            column(MatStencil_i,count) = i_g - 1
            column(MatStencil_j,count) = j_g - 1
            column(MatStencil_k,count) = (k_g+1)  - 1
            count = count + 1
          ENDIF          

          CALL MatSetValuesStencil(A,1,row,count-1,column,stencil,INSERT_VALUES,perr)

      ENDDO
    ENDDO
  ENDDO
  
  CALL MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,perr)
  CALL MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,perr)

END SUBROUTINE setupMatA_petsc

!! STILL TO FIX setupMatA_petsc

SUBROUTINE solve_petsc(numit,error)

    USE MPI
    USE data_module

    INTEGER,INTENT(INOUT) :: numit
    REAL(KIND=8),INTENT(INOUT) :: error
    INTEGER :: errcode, mpierr
    KSPConvergedReason :: reason
    PetscReal :: residual  ! Final residual

    CALL KSPSetOperators(kspObj,A,A,perr)
    CALL KSPSolve(kspObj,B,X,perr)
    CALL KSPGetIterationNumber(kspObj, numit, perr)

    CALL KSPGetConvergedReason(kspObj,reason,perr)
    CALL KSPGetResidualNorm(kspObj,residual,perr)
    error=residual

    IF(reason < 0) THEN
      WRITE(g_out,*) ' Error: Did not Converge. Calling MPI_Abort'
      WRITE(g_out,*) ' Divergence Reason:'

      IF(reason == -2) THEN
        WRITE(g_out,*) ' Diverged Null'
      ELSE IF(reason == -3) THEN
        WRITE(g_out,*) ' Diverged ITS, took ', numit, ' iterations'
      ELSE IF(reason == -4) THEN
        WRITE(g_out,*) ' Diverged DTOL'
      ELSE IF(reason == -5) THEN
        WRITE(g_out,*) ' Diverged Breakdown'
      ELSE IF(reason == -6) THEN
        WRITE(g_out,*) ' Diverged Breakdown BICG'
      ELSE IF(reason == -7) THEN
        WRITE(g_out,*) ' Diverged NonSymmetric'
      ELSE IF(reason == -8) THEN
        WRITE(g_out,*) ' Diverged Indefinite PC'
      ELSE IF(reason == -9) THEN
        WRITE(g_out,*) ' Diverged nanorinf'
      ELSE IF(reason == -10) THEN
        WRITE(g_out,*) ' Diverged indefinite mat'
      ENDIF

      CALL MPI_Abort(MPI_COMM_WORLD,errcode,mpierr)
    ENDIF 

    total_petsc_iter = total_petsc_iter + numit

END SUBROUTINE solve_petsc



SUBROUTINE printXVec(fileName)

    USE MPI

    IMPLICIT NONE

#   include "finclude/petscviewer.h"

    CHARACTER(LEN=*) :: fileName
    CHARACTER(LEN=8) :: id
    PetscViewer :: viewer

    ! Write Global Solution Vector
    CALL PetscViewerCreate(MPI_COMM_WORLD,viewer,perr)
    CALL PetscViewerSetType(viewer,"ascii",perr)
    CALL PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_INDEX,perr)
    CALL PetscViewerFileSetName(viewer,fileName,perr)
    CALL VecView(X,viewer,perr)
    CALL PetscViewerDestroy(viewer,perr)

END SUBROUTINE printXVec

SUBROUTINE printBVec(fileName)

    USE MPI

    IMPLICIT NONE

#   include "finclude/petscviewer.h"

    CHARACTER(LEN=*) :: fileName
    CHARACTER(LEN=8) :: id
    PetscViewer :: viewer

    CALL PetscViewerCreate(MPI_COMM_WORLD,viewer,perr)
    CALL PetscViewerSetType(viewer,"ascii",perr)
    CALL PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_INDEX,perr)
    CALL PetscViewerFileSetName(viewer,fileName,perr)
    CALL VecView(B,viewer,perr)
    CALL PetscViewerDestroy(viewer,perr)

END SUBROUTINE printBVec

SUBROUTINE printMatA(fileName)

    USE MPI

    IMPLICIT NONE

#   include "finclude/petscviewer.h"

    CHARACTER(LEN=*) :: fileName
    CHARACTER(LEN=8) :: id
    PetscViewer :: viewer

    CALL PetscViewerCreate(MPI_COMM_WORLD,viewer,perr)
    CALL PetscViewerSetType(viewer,"ascii",perr)
    CALL PetscViewerSetFormat(viewer,PETSC_VIEWER_DEFAULT,perr)
    CALL PetscViewerFileSetName(viewer,fileName,perr)
    CALL MatView(A,viewer,perr)
    CALL PetscViewerDestroy(viewer,perr)

END SUBROUTINE printMatA

END MODULE PETScTeaLeaf
