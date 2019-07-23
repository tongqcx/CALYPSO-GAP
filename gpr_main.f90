module GPR_MAIN
use constants
use io
use math
use linearalgebra
use struct

interface covariance
    module procedure covariance_2B, covariance_MB
end interface covariance
interface INI_GAP
    module procedure INI_GAP_2B, INI_GAP_MB
end interface INI_GAP

real(dp), dimension(:,:), allocatable :: c_subYY_sqrtInverseLambda
real(dp), dimension(:,:), allocatable :: factor_c_subYsubY
real(dp), dimension(:,:), allocatable :: a
real(dp), dimension(:),   allocatable :: globalY
type(LA_Matrix)                       :: LA_c_subYsubY, LA_q_subYsubY
integer                               :: n_globalSparseX , n_globalY
integer                               :: error



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

!SUBROUTINE INI_GAP(GAP, nobf)
!type(GAP_type),intent(in) :: GAP
!integer,intent(in)        :: nobf
!!--local--
!!integer                 :: ninteraction
!
!allocate(GAP%cmm(GAP%nsparse, GAP%nsparse))
!allocate(GAP%cmo(GAP%nsparse, nobf, ninteraction))
!allocate(GAP%sparseX(GAP%nsparse. GAP%dd))
!allocate(GAP%obe(nobf))
!allocate(GAP%coeff(GAP%nsparse,ninteraction))
!allocate(GAP%lamda(nobf))
!allocate(GAP%lamdaobe(nobf))
!allocate(GAP%sparsecut(GAP%nsparse))
!END SUBROUTINE

SUBROUTINE INI_GAP_2B(GAP, nsparse, nobf)
type(GAP_type),intent(inout) :: GAP
integer,intent(in)           :: nsparse, nobf

!local
real(DP)                     :: dr3
GAP%dd = 1
GAP%nsparse = nsparse

allocate(GAP%cmm(GAP%nsparse, GAP%nsparse))
allocate(GAP%cmo(GAP%nsparse, nobf, ninteraction))
allocate(GAP%sparseX(GAP%nsparse, GAP%dd))
allocate(GAP%obe(nobf))
allocate(GAP%coeff(GAP%nsparse,ninteraction))
allocate(GAP%lamda(nobf))
allocate(GAP%lamdaobe(nobf))
allocate(GAP%sparsecut(GAP%nsparse))

dr3 = (rcut - rmin)/(GAP%nsparse - 1)
do i = 1, GAP%nsparse
    sparseX(i,1) = rmin + (i - 1)*dr3
enddo
!$OMP parallel do schedule(dynamic) default(shared) private(i,j,fc_i,fc_j)
    do i  = 1, nsparse
        GAP%cmm(i,i) = delta**2 
        fc_i = fcutij(sparseX(i,1))
        GAP%cmm(i,i) = GAP%cmm(i,i)*fc_i*fc_i + sigma_jitter
        do j = i + 1, nsparse
            !fc_j = fcutij(sparseX(j))
            GAP%cmm(i,j) = covariance(sparseX(i, :),sparseX(j, :))
            !cmm(i,j) = cmm(i,j) * fc_i * fc_j
            GAP%cmm(j,i) = GAP%cmm(i,j)
        enddo
    enddo
call write_array(GAP%cmm,'cmm.dat')
END SUBROUTINE

SUBROUTINE INI_GAP_MB(GAP, nsparse, dd, nobf)
type(GAP_type),intent(inout) :: GAP
integer,intent(in)           :: nsparse, dd, nobf

!local
GAP%nsparse = nsparse
GAP%dd = dd

allocate(GAP%cmm(GAP%nsparse, GAP%nsparse))
allocate(GAP%cmo(GAP%nsparse, nobf, ninteraction))
allocate(GAP%sparseX(GAP%nsparse, GAP%dd))
allocate(GAP%obe(nobf))
allocate(GAP%coeff(GAP%nsparse,ninteraction))
allocate(GAP%lamda(nobf))
allocate(GAP%lamdaobe(nobf))
allocate(GAP%sparsecut(GAP%nsparse))

END SUBROUTINE INI_GAP_MB

FUNCTION  covariance_2B(x,y)
implicit none
real(8),intent(in)  :: x
real(8),intent(in)  :: y
real(8)             :: covariance_2B

!integer  i 
REAL(DP)            :: fc_i, fc_j
fc_i = fcutij(x)
fc_j = fcutij(y)
covariance = 0.d0
covariance = covariance + ((x-y)/theta)**2
covariance = delta**2*exp(-0.5d0*covariance) * fc_i * fc_j
END FUNCTION covariance_2B

FUNCTION  covariance_MB(x,y, theta)
implicit none
real(DP),intent(in)  :: x(:)
real(DP),intent(in)  :: y(:)
real(DP),intent(in)  :: theta(:)
real(DP)             :: covariance_MB

integer  i 
!REAL(DP)            :: fc_i, fc_j
covariance = 0.d0
do i = 1, size(x)
    covariance = covariance + ((x(i)-y(i))/theta(i))**2
enddo
covariance = delta**2*exp(-0.5d0*covariance) 
END FUNCTION covariance_MB

FUNCTION  DcovarianceDx(x,y)
implicit none
real(8),intent(in)  :: x
real(8),intent(in)  :: y
real(8)             :: DcovarianceDx
!integer  i 
REAL(DP)            :: fc_i, fc_j, dfc_i, exp_part

DcovarianceDx = 0.d0
exp_part = 0.d0
fc_i = fcutij(x)
fc_j = fcutij(y)
dfc_i = dfcutij(x)
exp_part = exp_part + ((x-y)/theta)**2
exp_part = delta**2 * exp(-0.5d0 * exp_part)
DcovarianceDx = exp_part * -1.d0 * (x-y)/theta**2
DcovarianceDx = DcovarianceDx * fc_i + exp_part * dfc_i
DcovarianceDx = DcovarianceDx * fc_j
END FUNCTION DcovarianceDx

subroutine matmuldiag(x,y)
real(8),intent(in)   :: x(:)
real(8),intent(inout) :: y(:,:)
do i = 1,size(x)
    y(i,:) = x(i)*y(i,:)
enddo
end subroutine matmuldiag

subroutine matmuldiag_T(y,x)
real(8),intent(in)   :: x(:)
real(8),intent(inout) :: y(:,:)
do i = 1,size(x)
    y(:,i) = x(i)*y(:,i)
enddo
end subroutine matmuldiag_T

SUBROUTINE gp_predict(at)
implicit none
type(Structure),intent(inout)          :: at


! local 
integer                               :: i,j,k, k1, k2, k3, k4
integer                               :: interaction_index
REAL(DP)                              :: rij, fcut_ij , dfcut_ij
REAL(DP)                              :: ene
REAL(DP),allocatable,dimension(:,:,:) :: stress_i

at%energy_cal = 0.d0
at%force_cal = 0.d0
at%stress = 0.d0
allocate(stress_i(3,3,at%natoms))
stress_i = 0.d0

!$OMP parallel do schedule(dynamic) default(shared) private(i,j,k,rij, fcut_ij, interaction_index, k1, k2, k3, k4, dfcut_ij, ene)
do i = 1,at%natoms
    do j = 1, nspecies
        do k = 1, at%atom(i)%count(j)
            rij = at%atom(i)%neighbor(j,k,4)
            fcut_ij = fcutij(rij)
            dfcut_ij = dfcutij(rij)
            interaction_index = at%interaction_mat(at%index(i),j)
            ene = 0.d0
            do k1 = 1, nsparse
!***********  get total energy                   
                ene = ene + covariance(rij, sparseX(k1)) * fcut_ij * coeff(k1,interaction_index) 

!***********  get atomic force                    
                do k2 = 1,3
                    at%force_cal(i,k2) = at%force_cal(i,k2) + &
                    (dfcut_ij * covariance(rij, sparseX(k1)) + DcovarianceDx(rij, sparseX(k1)) * fcut_ij) * &
                    (at%atom(i)%pos(k2) - at%atom(i)%neighbor(j,k,k2))/rij * coeff(k1,interaction_index) 

!***********  get atomic cell stress
                    do k3 = 1,3
                        stress_i(k3,k2,i) = stress_i(k3,k2,i) + (at%atom(i)%pos(k3) - at%atom(i)%neighbor(j,k,k3)) * &
                        (dfcut_ij * covariance(rij, sparseX(k1)) + DcovarianceDx(rij, sparseX(k1)) * fcut_ij) * &
                        (at%atom(i)%pos(k2) - at%atom(i)%neighbor(j,k,k2))/rij * coeff(k1,interaction_index) 

                    enddo ! k3
                enddo ! k2
            enddo ! k1
            at%atomic_energy(i) = at%atomic_energy(i) + ene
        enddo ! k
    enddo ! j
!    print*, at%atomic_energy(i)
enddo
!print*, at%atomic_energy
at%energy_cal = sum(at%atomic_energy) * 0.5d0 + at%natoms * ene_cons
!print*, at%energy_cal
at%force_cal = -1.0 * at%force_cal
at%volume = volume(at%lat)
at%stress = sum(stress_i, dim=3) * -0.5d0
at%stress = at%stress / at%volume * (1.0/GPa2eVPang) 
at%stress_cal(1) = at%stress(1,1)
at%stress_cal(2) = at%stress(1,2)
at%stress_cal(3) = at%stress(1,3)
at%stress_cal(4) = at%stress(2,2)
at%stress_cal(5) = at%stress(2,3)
at%stress_cal(6) = at%stress(3,3)
deallocate(stress_i)
END SUBROUTINE

SUBROUTINE get_cmo(GAP)
type(GAP_type),intent(in)    ::   GAP
integer   ::   i, j, k1, k2, k3, interaction_index
REAL(DP)  ::   rij

!$OMP parallel do schedule(dynamic) default(shared) private(i, j, k1, k2, k3, interaction_index, rij)
do i = 1, GAP%nsparse
    do j = 1, nconfig
!***********************************************
        do k1 = 1, at(j)%natoms
            do k2 = 1, nspecies
                do k3 = 1, at(j)%atom(k1)%count(k2)
                    interaction_index = at(j)%interaction_mat(at(j)%index(k1),k2)
                    rij = at(j)%atom(k1)%neighbor(k2,k3,4)
                    GAP%cmo(i,j,interaction_index) = GAP%cmo(i,j,interaction_index) + covariance(sparseX(i,1), rij) * 0.5d0
                enddo
            enddo
        enddo
    enddo
enddo
!cmo(:,:,2) = cmo(:,:,2)/2.d0
END SUBROUTINE

END module
