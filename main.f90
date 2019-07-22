Program  test_nGAP
use constants
use io
use struct
use gpr_main
implicit none
integer            :: i,j,k,ii,jj,kk,k1, k2,k3
integer            :: na
integer            :: nconfig
real(dp)           :: fc_i, fc_j, rij
logical            :: alive
integer            :: interaction_index

call read_input()

if (ltrain) then
print*, 'nsparse', nsparse
open(2211,file='config')
read(2211,*) nconfig
allocate(at(nconfig))
close(2211)
call read_structure('config',at)

call ini_gap(nconfig)
call ini_gap_2b()

! Build matrix cmo
print*, 'nspecies',nspecies
cmo(:,:,:) = 0.d0
k = 0

!!!$OMP parallel do schedule(dynamic) default(shared) private(i,j,k1,k2,k,k3,fc_i,fc_j, ii, jj)
!do i = 1, nsparse
!    do j = 1, nconfig
!        k = 0
!        do k1 = 1, nspecies
!            do k2 = k1 , nspecies
!                k = k + 1
!                do ii = at(j)%pos_index(k1,1) , at(j)%pos_index(k1,2)
!                    do jj = 1, at(j)%atom(ii)%count(k2)
!                        fc_i = fcutij(sparseX(i))
!                        fc_j = fcutij(at(j)%atom(ii)%neighbor(k2,jj,4))
!                        cmo(i,j,k) = cmo(i,j,k) + covariance(sparseX(i), at(j)%atom(ii)%neighbor(k2,jj,4)) * 0.5d0
!        !                cmo(i,j,k) = cmo(i,j,k) + covariance(sparseX(i), at(j)%atom(ii)%neighbor(k2,jj,4)) * fc_j * 0.5d0
!                    enddo
!                enddo
!            enddo
!        enddo
!    enddo
!enddo

!!$OMP parallel do schedule(dynamic) default(shared) private(i,j,k1,k2,k,fc_i,fc_j, ii, jj, rij)
do i = 1, nsparse
    do j = 1, nconfig
!***********************************************
        do k1 = 1, at(j)%natoms
            do k2 = 1, nspecies
                do k3 = 1, at(j)%atom(k1)%count(k2)
                    interaction_index = at(j)%interaction_mat(at(j)%index(k1),k2)
                    rij = at(j)%atom(k1)%neighbor(k2,k3,4)
                    cmo(i,j,interaction_index) = cmo(i,j,interaction_index) + covariance(sparseX(i), rij) * 0.5d0
                enddo
            enddo
        enddo
    enddo
enddo
cmo(:,:,2) = cmo(:,:,2)/2.d0

call write_array(cmo(:,:,1),'cmo1.dat')
call write_array(cmo(:,:,2),'cmo2.dat')
call write_array(cmo(:,:,3),'cmo3.dat')

do i = 1,nconfig
    obe(i) = at(i)%energy_ref - at(i)%natoms * ene_cons
    at(i)%sigma_e = sigma_e * sqrt(1.d0 * at(i)%natoms)
    lamda(i) = at(i)%sigma_e**2
    lamdaobe(i) = obe(i) * sqrt((1.0/lamda(i)))
enddo

do k = 1, ninteraction
    call matmuldiag_T(cmo(:,:,k),sqrt(1.0/lamda))
!call write_array(lamdaobe,'lamdaobe.dat')
    call gpr(cmm, cmo(:,:,k), lamdaobe, coeff(:,k))
enddo

open(2234,file='coeffx.dat')
do i = 1,nsparse
    sparsecut(i) = fcutij(sparseX(i))
    write(2234,'(I3,F25.8,$)') i, sparseX(i)
    write(2234,'(F25.8, $)') sparsecut(i)
    do k = 1,ninteraction
        write(2234,'(F25.8,$)') coeff(i,k)
    enddo
    write(2234,*)
enddo
close(2234)
deallocate(at)
ENDIF  ! ltrain

if (ltest) then
print*, '************************************************'
print*, '*************** BEGIN PREDICTING ***************'
print*, '************************************************'
if (.not.ltrain) then
    inquire(file="coeffx.dat",exist=alive)
    if(.not.alive) then
        print*, "coeffx.dat file does not exist!"
        stop
    endif
    allocate(sparseX(nsparse))
    allocate(sparsecut(nsparse))
    allocate(coeff(nsparse, ninteraction))
    open(111,file='coeffx.dat')
    do i = 1, nsparse
        read(111,*) j, sparseX(i), sparsecut(i), coeff(i,:)
    enddo
    close(111)
endif
  
open(2211,file='test')
read(2211,*) nconfig
allocate(at(nconfig))
close(2211)
CALL  SYSTEM_CLOCK(it1)
call read_structure('test',at)
CALL  SYSTEM_CLOCK(it2)
print*, "CPU time used (sec) For converting coord: ",(it2 - it1)/10000.0

CALL  SYSTEM_CLOCK(it1)
!$OMP parallel do schedule(dynamic) default(shared) private(i)
do i = 1, nconfig
    call gp_predict(at(i))
enddo
CALL  SYSTEM_CLOCK(it2)
print*, "CPU time used (sec) For GP Predict: ",(it2 - it1)/10000.0

rmse_energy = 0.d0
rmse_force = 0.d0
rmse_stress = 0.d0
nforce = 0
do i = 1, nconfig
    rmse_energy = rmse_energy + (at(i)%energy_cal/at(i)%natoms - at(i)%energy_ref/at(i)%natoms)**2
    do j = 1, at(i)%natoms
        do k = 1, 3
            nforce = nforce + 1
            rmse_force = rmse_force + (at(i)%force_cal(j,k) - at(i)%force_ref(j,k))**2
        enddo
    enddo
    do j = 1,6
        rmse_stress = rmse_stress + (at(i)%stress_cal(j) - at(i)%stress_ref(j))**2
    enddo
enddo
print *, 'RMSE ENERGY:', sqrt(rmse_energy/nconfig)
print *, 'RMSE FORCE:', sqrt(rmse_force/nforce)
print *, 'RMSE STRESS Units GPa:', sqrt(rmse_stress/nconfig/6.d0)
open(181,file="predited.datf")
do ii = 1, nconfig
        write(181,*) "----------------------------------------------------"
        write(181,'(A9,X,I5,X,A30)') "Structure",ii
        write(181,*) "----------------------------------------------------"
        do i  =1,at(ii)%natoms
            do j = 1,3
            if (j==1) then
            write(181,'(I5,I5,X,A,A,X,3F15.6)') ii,i, "F","X",at(ii)%force_cal(i,j),at(ii)%force_ref(i,j),&
                                    abs(at(ii)%force_cal(i,j) - at(ii)%force_ref(i,j))
            elseif(j ==2 ) then
            write(181,'(I5,I5,X,A,A,X,3F15.6)') ii,i, "F","Y",at(ii)%force_cal(i,j),at(ii)%force_ref(i,j),&
                                    abs(at(ii)%force_cal(i,j) - at(ii)%force_ref(i,j))
            else
            write(181,'(I5,I5,X,A,A,X,3F15.6)') ii,i, "F","Z",at(ii)%force_cal(i,j),at(ii)%force_ref(i,j),&
                                    abs(at(ii)%force_cal(i,j) - at(ii)%force_ref(i,j))
            endif
            enddo
        enddo
        do j = 1,6
            if (j==1) then
            write(181,'(A2,X,3F20.5)') "XX",at(ii)%stress_cal(j),at(ii)%stress_ref(j),&
                                     abs(at(ii)%stress_cal(j) - at(ii)%stress_ref(j))
            elseif(j==2) then
            write(181,'(A2,X,3F20.5)') "XY",at(ii)%stress_cal(j),at(ii)%stress_ref(j),&
                                     abs(at(ii)%stress_cal(j) - at(ii)%stress_ref(j))
            elseif(j==3) then
            write(181,'(A2,X,3F20.5)') "XZ",at(ii)%stress_cal(j),at(ii)%stress_ref(j),&
                                     abs(at(ii)%stress_cal(j) - at(ii)%stress_ref(j))
            elseif(j==4) then
            write(181,'(A2,X,3F20.5)') "YY",at(ii)%stress_cal(j),at(ii)%stress_ref(j),&
                                     abs(at(ii)%stress_cal(j) - at(ii)%stress_ref(j))
            elseif(j==5) then
            write(181,'(A2,X,3F20.5)') "YZ",at(ii)%stress_cal(j),at(ii)%stress_ref(j),&
                                     abs(at(ii)%stress_cal(j) - at(ii)%stress_ref(j))
            else
            write(181,'(A2,X,3F20.5)') "ZZ",at(ii)%stress_cal(j),at(ii)%stress_ref(j),&
                                     abs(at(ii)%stress_cal(j) - at(ii)%stress_ref(j))
            endif
        enddo
        write(181,'(A6,X,I5,X,3F15.6)') "ENERGY",ii,at(ii)%energy_cal/at(ii)%natoms, at(ii)%energy_ref/at(ii)%natoms,&
               abs(at(ii)%energy_cal-at(ii)%energy_ref)/at(ii)%natoms
        !write(181,'(A10X3F15.8)') 'E-V:', at(i)%volume/at(i)%na, (at(i)%stress_cal(1) + at(i)%stress_cal(4) + at(i)%stress_cal(6))/30.0, &
        !at(i)%energy_cal/at(i)%na
enddo
ENDIF ! ltest
end program

