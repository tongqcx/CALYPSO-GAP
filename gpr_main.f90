module GPR_MAIN
use constants
use math
use my_line

real(dp), dimension(:,:), allocatable :: c_subYY_sqrtInverseLambda
real(dp), dimension(:,:), allocatable :: factor_c_subYsubY
real(dp), dimension(:,:), allocatable :: a
real(dp), dimension(:),   allocatable :: globalY
type(LA_Matrix)                       :: LA_c_subYsubY, LA_q_subYsubY
integer                               :: n_globalSparseX , n_globalY
integer                               :: error


REAL(DP),DIMENSION(:),ALLOCATABLE       :: lamda
REAL(DP),DIMENSION(:),ALLOCATABLE       :: lamdaobe
REAL(DP),DIMENSION(:,:),ALLOCATABLE     :: cmm
REAL(DP),DIMENSION(:,:,:),ALLOCATABLE   :: cmo
REAL(DP),DIMENSION(:),ALLOCATABLE       :: sparsecut
REAL(DP),DIMENSION(:),ALLOCATABLE       :: sparseX

contains
SUBROUTINE GPR(cmm, cmo, lamdaobe, alpha)
implicit none
real(dp),intent(in),dimension(:,:)      :: cmm
real(dp),intent(in),dimension(:,:)      :: cmo
real(dp),intent(in),dimension(:)        :: lamdaobe
real(dp),intent(out),dimension(:)       :: alpha
!--- local ---
integer                                 :: i,j
n_globalSparseX = size(cmo, 1)
n_globalY = size(cmo, 2)
!==========================================================================
allocate(factor_c_subYsubY(n_globalSparseX,n_globalSparseX))
allocate(a(n_globalY+n_globalSparseX,n_globalSparseX))
allocate(globalY(n_globalY+n_globalSparseX))

globalY = 0.d0
call LA_Matrix_initialise(LA_c_subYsubY,cmm)
call LA_Matrix_Factorise(LA_c_subYsubY,factor_c_subYsubY,error=error)
call LA_Matrix_finalise(LA_c_subYsubY)


do i = 1, n_globalSparseX-1
   do j = i+1, n_globalSparseX
      factor_c_subYsubY(j,i) = 0.0_qp
   end do
end do

a(1:n_globalY,:) = transpose(cmo(:,:))
a(n_globalY+1:,:) = factor_c_subYsubY
do i = 1, n_globalY
    globalY(i) = lamdaobe(i)
enddo
call LA_matrix_initialise(LA_q_subYsubY,a)
call LA_Matrix_QR_Solve_Vector(LA_q_subYsubY,globalY,alpha)
call LA_matrix_finalise(LA_q_subYsubY)
deallocate(factor_c_subYsubY)
deallocate(a)
deallocate(globalY)
END subroutine

SUBROUTINE INI_GAP(nsparse, nobf)
integer,intent(in)      :: nsparse
integer,intent(in)      :: nobf
!--local--
integer                 :: ninteraction
ninteraction = nspecies * (nspecies + 1)/2.d0
allocate(cmm(nsparse, nsparse))
allocate(cmo(nsparse, nobf, ninteraction))
allocate(sparseX(nsparse))
END SUBROUTINE

SUBROUTINE INI_GAP_2B()
dr3 = (rcut - rmin)/(n_sparse - 1)
do i = 1, n_sparse
    sparseX(i) = rmin + (i - 1)*dr3
enddo
!$OMP parallel do schedule(dynamic) default(shared) private(i,j)
    do i  = 1, n_sparse
        cmm(i,i) = delta**2 
        fc_i = fcutij(sparseX(i))
        cmm(i,i) = cmm(i,i)*fc_i*fc_i + sigma_jitter
        do j = i + 1, n_sparse
            fc_j = fcutij(sparseX(j))
            cmm(i,j) = covariance(sparseX(i),sparseX(j))
            cmm(i,j) = cmm(i,j) * fc_i * fc_j
            cmm(j,i) = cmm(i,j)
        enddo
    enddo
END SUBROUTINE

END module
