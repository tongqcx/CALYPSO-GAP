! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! X  the module of gaussian process regression 
! X  for two-body interaction
! X  Qunchao Tong 2019.07.24 20:55
! X
! X
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
module GPR_2B
use constants
use io
use math
use linearalgebra
use struct
use gpr_base

CONTAINS
SUBROUTINE GAP_INI_2B(GAP, AT, DATA_C)
type(GAP_type),intent(inout)             :: GAP
type(Structure),intent(in),dimension(:)  :: at
type(DATA_type),intent(in)               :: DATA_C

!local
real(DP)                     :: dr3
GAP%dd = 1
GAP%nsparse = DATA_C%nsparse_2B
GAP%nglobalY = DATA_C%ne
GAP%delta = DATA_C%delta_2B
GAP%ninteraction = data_c%ninteraction
GAP%sigma_e = DATA_C%sigma_e_2B
print*, '********************************************************'
print*, 'Parameters for TWO_BODY interaction'
print*, 'The size of sparse set:', GAP%nsparse
print*, 'The dimension of ACSF descriptors', GAP%dd
print*, 'The size of globalY:', GAP%nglobalY
print*, 'The value of delta:', GAP%delta
print*, 'The number of interaction:', GAP%ninteraction
print*, '********************************************************'

allocate(GAP%theta(1))
allocate(GAP%cmm(GAP%nsparse, GAP%nsparse))
allocate(GAP%cmo(GAP%nsparse, GAP%nglobalY, GAP%ninteraction))
allocate(GAP%obe(GAP%nglobalY))
allocate(GAP%coeff(GAP%nsparse, GAP%ninteraction))
allocate(GAP%lamda(GAP%nglobalY))
allocate(GAP%lamdaobe(GAP%nglobalY, 1))
allocate(GAP%sparsecut(GAP%nsparse))
allocate(GAP%sparseX(GAP%nsparse, GAP%dd))

GAP%theta(1)  = data_c%theta_2b

dr3 = (data_c%rcut - data_c%rmin)/(GAP%nsparse - 1)
do i = 1, GAP%nsparse
    GAP%sparseX(i,1) = data_c%rmin + (i - 1)*dr3
    GAP%sparsecut(i) = fcut_ij(GAP%sparseX(i,1))
enddo
!$OMP parallel do schedule(dynamic) default(shared) private(i,j,fc_i,fc_j)
    do i  = 1, GAP%nsparse
        GAP%cmm(i,i) = GAP%delta**2 
        fc_i = fcut_ij(GAP%sparseX(i,1))
        GAP%cmm(i,i) = GAP%cmm(i,i)*fc_i*fc_i + DATA_C%sigma_jitter
        do j = i + 1, GAP%nsparse

            GAP%cmm(i,j) = covariance_2B(GAP%delta, GAP%sparseX(i,1), GAP%sparseX(j,1),GAP%theta(1))
            GAP%cmm(j,i) = GAP%cmm(i,j)
        enddo
    enddo
call write_array(GAP%cmm,'cmm.dat')
do i = 1, DATA_C%ne
    GAP%lamda(i) = (GAP%sigma_e * (sqrt(1.d0 * at(i)%natoms)))**2
enddo
END SUBROUTINE

SUBROUTINE GAP_COEFF_2B(GAP, DATA_C)
type(GAP_type),intent(inout)             :: GAP
type(DATA_type),intent(in)               :: DATA_C

do i = 1, DATA_C%ne
    GAP%obe(i) = DATA_C%obe(i) 
    GAP%lamdaobe(i, 1) = DATA_C%obe(i) * sqrt(1.d0/GAP%lamda(i))
enddo
call write_array(GAP%lamdaobe(:, 1), 'lamdaobe.dat')
do k = 1, DATA_C%ninteraction
    call matmuldiag_T(GAP%cmo(:,:,k),sqrt(1.0/GAP%lamda))
    call gpr(GAP%cmm, GAP%cmo(:,:,k), GAP%lamdaobe(:,1), GAP%coeff(:,k))
enddo
END SUBROUTINE GAP_COEFF_2B

SUBROUTINE GAP_WRITE_PARAS_2B(GAP)
type(GAP_type),intent(in)             :: GAP

open(2234, file='gap_paras_2b.dat')
write(2234,*) GAP%nsparse, GAP%ninteraction
write(2234,*)
write(2234,*)
write(2234,*)
write(2234,*) GAP%delta
write(2234,*) GAP%theta
do i = 1,GAP%nsparse
    write(2234,'(I3,F30.15,$)') i, GAP%sparseX(i,1)
    write(2234,'(F30.15, $)') GAP%sparsecut(i)
    do k = 1,DATA_C%ninteraction
        write(2234,'(F30.15,$)') GAP%coeff(i,k)
    enddo
    write(2234,*)
enddo
close(2234)
ENDSUBROUTINE GAP_WRITE_PARAS_2B

SUBROUTINE GAP_READ_PARAS_2B(GAP)
type(GAP_type),intent(inout)             :: GAP
open(2234, file='gap_paras_2b.dat')
read(2234,*) GAP%nsparse,  GAP%ninteraction
read(2234,*)
read(2234,*)
read(2234,*)
allocate(GAP%sparsecut(GAP%nsparse))
allocate(GAP%sparseX(GAP%nsparse, 1))
allocate(GAP%coeff(GAP%nsparse, GAP%ninteraction))
allocate(GAP%theta(1))
read(2234,*) GAP%delta
read(2234,*) GAP%theta
do i_sparsex = 1,GAP%nsparse
    read(2234,*) i, GAP%sparseX(i,1), GAP%sparsecut(i), GAP%coeff(i,:)
enddo
END SUBROUTINE GAP_READ_PARAS_2B


SUBROUTINE GAP_predict_2B(GAP, at, DATA_C)
implicit none
type(GAP_type),intent(in)              :: GAP
type(Structure),intent(inout)          :: at
type(DATA_type),intent(in)             :: DATA_C


! local 
integer                               :: i,j,k, k1, k2, k3
integer                               :: interaction_index, n
REAL(DP)                              :: rij, fc_ij , dfc_ij
REAL(DP)                              :: ene
REAL(DP)                              :: cov, dcov
REAL(DP),allocatable,dimension(:,:,:) :: stress_i

at%energy_cal = 0.d0
at%force_cal = 0.d0
at%stress = 0.d0
at%atomic_energy = 0.d0
allocate(stress_i(3,3,at%natoms))
stress_i = 0.d0

!$OMP parallel do schedule(dynamic) default(shared) private(i,j,k,rij, fc_ij, interaction_index, k1, k2, k3, dfc_ij, ene, cov, dcov)
do i = 1,at%natoms
    do j = 1, at%nspecies
        do k = 1, at%atom(i)%count(j)
            rij = at%atom(i)%neighbor(j, k, 4)
            n = int(at%atom(i)%neighbor(j,k, 5))
            fc_ij = fcut_ij(rij)
            dfc_ij = dfcut_ij(rij)
            interaction_index = DATA_C%interaction_mat(at%index(i),at%index(n))
            ene = 0.d0
            do k1 = 1, GAP%nsparse
!***********  get total energy                   
!print*, GAP%delta, rij, GAP%sparseX(k1, 1), GAP%theta(1)
                cov = covariance_2B(GAP%delta, rij, GAP%sparseX(k1, 1), GAP%theta(1))
                dcov = dcovdx_2B(GAP%delta, rij, GAP%sparseX(k1, 1), GAP%theta(1))
                ene = ene + cov * fc_ij * GAP%coeff(k1,interaction_index) 
!                print*, cov , fc_ij , GAP%coeff(k1,interaction_index)
!                pause

!***********  get atomic force                    
                do k2 = 1,3
                    at%force_cal(i,k2) = at%force_cal(i,k2) + (dfc_ij * cov + dcov * fc_ij) * &
                    (at%atom(i)%pos(k2) - at%atom(i)%neighbor(j,k,k2))/rij * GAP%coeff(k1,interaction_index) 

!***********  get atomic cell stress
                    do k3 = 1,3
                        stress_i(k3,k2,i) = stress_i(k3,k2,i) + (at%atom(i)%pos(k3) - at%atom(i)%neighbor(j,k,k3)) * &
                        (dfc_ij * cov + dcov * fc_ij) * &
                        (at%atom(i)%pos(k2) - at%atom(i)%neighbor(j,k,k2))/rij * GAP%coeff(k1,interaction_index) 

                    enddo ! k3
                enddo ! k2
            enddo ! k1
            at%atomic_energy(i) = at%atomic_energy(i) + ene
        enddo ! k
    enddo ! j
!    print*, at%atomic_energy(i)
enddo
at%energy_cal = sum(at%atomic_energy) * 0.5d0 + at%natoms * ene_cons
at%force_cal = -1.0 * at%force_cal
at%volume = volume(at%lat)
at%stress = sum(stress_i, dim=3) * -0.5d0
!at%stress = at%stress / at%volume * (1.0/GPa2eVPang) * 10.d0
at%stress_cal(1) = at%stress(1,1)
at%stress_cal(2) = at%stress(1,2)
at%stress_cal(3) = at%stress(1,3)
at%stress_cal(4) = at%stress(2,2)
at%stress_cal(5) = at%stress(2,3)
at%stress_cal(6) = at%stress(3,3)
deallocate(stress_i)
END SUBROUTINE

SUBROUTINE GAP_cmo_2B(GAP, AT, DATA_C)
implicit none
type(GAP_type),intent(inout)             :: GAP
type(Structure),intent(in),dimension(:)  :: at
type(DATA_type),intent(in)               :: DATA_C
integer                                  :: i, j, k1, k2, k3
integer                                  :: n, interaction_index
REAL(DP)                                 :: rij

GAP%cmo = 0.d0
!$OMP parallel do schedule(dynamic) default(shared) private(i, j, k1, k2, k3, interaction_index, rij, n)
do i = 1, GAP%nsparse
    do j = 1, GAP%nglobalY
!***********************************************
        do k1 = 1, at(j)%natoms
            do k2 = 1, at(j)%nspecies
                do k3 = 1, at(j)%atom(k1)%count(k2)
                    rij = at(j)%atom(k1)%neighbor(k2, k3, 4)
                    n = int(at(j)%atom(k1)%neighbor(k2, k3, 5))
                    interaction_index = data_c%interaction_mat(at(j)%index(k1),at(j)%index(n))
                    GAP%cmo(i,j,interaction_index) = GAP%cmo(i,j,interaction_index) + &
                    covariance_2B(GAP%delta, GAP%sparseX(i,1), rij, GAP%theta(1)) * 0.5d0
                enddo
            enddo
        enddo
    enddo
enddo
!cmo(:,:,2) = cmo(:,:,2)/2.d0
END SUBROUTINE

END module
