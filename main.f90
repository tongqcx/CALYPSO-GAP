Program  test_nGAP
use constants
use io
use struct
use gpr_base
use gpr_mb
use gpr_2b
implicit none
integer            :: i,j,k,ii,jj,kk,k1, k2,k3
integer            :: na
integer            :: nconfig
real(dp)           :: fc_i, fc_j, rij
logical            :: alive
integer            :: interaction_index
integer            :: ran_seed(1)
ran_seed=1
call random_seed(put=ran_seed)

call read_input()

if (ltrain) then
open(2211,file='config')
read(2211,*) nconfig
allocate(at(nconfig))
close(2211)

call read_structure('config',at, data_c)

!*****************************************************
!
!*****************************************************
call gap_ini_mb(GAP_MB, AT, ACSF, DATA_C)
CALL  SYSTEM_CLOCK(it1)
call gap_cmo_mb(GAP_MB, AT, DATA_C)
CALL  SYSTEM_CLOCK(it2)
print*, "CPU time used (sec) For matrix CMO: ",(it2 - it1)/10000.0
call gap_coeff_mb(GAP_MB,AT,DATA_C)
call gap_write_paras_mb(GAP_MB)


!*****************************************************
!
!call gap_ini_2b(GAP_2B, AT, DATA_C)
!call gap_cmo_2b(GAP_2B, AT, DATA_C)
!call gap_coeff_2b(GAP_2B, DATA_C)
!print*, GAP_2B%coeff(:,1)
!stop





!open(2234,file='coeffx.dat')
!do i = 1,GAP_2B%nsparse
!    GAP_2B%sparsecut(i) = fcut_ij(GAP_2B%sparseX(i,1))
!    write(2234,'(I3,F25.8,$)') i, GAP_2B%sparseX(i,1)
!    write(2234,'(F25.8, $)') GAP_2B%sparsecut(i)
!    do k = 1,DATA_C%ninteraction
!        write(2234,'(F25.8,$)') GAP_2B%coeff(i,k)
!    enddo
!    write(2234,*)
!enddo
!close(2234)
deallocate(at)
call destropy_data_type(DATA_C)
ENDIF  ! ltrain

if (ltest) then
print*, '************************************************'
print*, '*************** BEGIN PREDICTING ***************'
print*, '************************************************'
if (.not.ltrain) then
    call gap_read_paras_mb(GAP_MB, ACSF)
endif
  
open(2211,file='test')
read(2211,*) nconfig
allocate(at(nconfig))
close(2211)
CALL  SYSTEM_CLOCK(it1)
call read_structure('test',at, data_c)
CALL  SYSTEM_CLOCK(it2)
print*, "CPU time used (sec) For converting coord: ",(it2 - it1)/10000.0

CALL  SYSTEM_CLOCK(it1)
!$OMP parallel do schedule(dynamic) default(shared) private(i)
do i = 1, nconfig
!    call gap_predict_2B(GAP_2B, at(i), DATA_C)
    call gap_predict_MB(GAP_MB, at(i))
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
    at(i)%stress_ref = at(i)%stress_ref * (1.0/GPa2eVPang) * 10.0 / at(i)%volume
    do j = 1,6
        rmse_stress = rmse_stress + (at(i)%stress_cal(j)/10.0 - at(i)%stress_ref(j)/10.0)**2
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

