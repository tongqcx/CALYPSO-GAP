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
logical            :: alive, T_MB, T_2B
integer            :: interaction_index
integer            :: ran_seed(1)
T_MB = .true.
T_2B = .false.
ran_seed=1
call random_seed(put=ran_seed)

call read_input()

if (ltrain) then
open(2211,file='config')
read(2211,*) nconfig
allocate(at(nconfig))
close(2211)

call read_structure('config',at, data_c)
call set_gpr_ob(at, data_c)

!*****************************************************
!
!*****************************************************
if (T_MB) then
    call gap_ini_mb(GAP_MB, AT, ACSF, DATA_C)
    call gap_cmo_mb(GAP_MB, AT, DATA_C)
    call gap_coeff_mb(GAP_MB,AT,DATA_C)
    call gap_write_paras_mb(GAP_MB)
endif


!*****************************************************
!
if (T_2B) then
    call gap_ini_2b(GAP_2B, AT, DATA_C)
    call gap_cmo_2b(GAP_2B, AT, DATA_C)
    call write_array(GAP_2B%cmo(:,:,1),'cmo.dat')
    call gap_coeff_2b(GAP_2B, DATA_C)
    call gap_write_paras_2b(GAP_2B)
endif
!print*, GAP_2B%coeff(:,1)
!stop
deallocate(at)
ENDIF  ! ltrain

if (ltest) then
print*, '*************** BEGIN PREDICTING ***************'
if (.not.ltrain) then
    if (T_MB) call gap_read_paras_mb(GAP_MB, ACSF)
    if (T_2B) call gap_read_paras_2b(GAP_2B)
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
    if (T_2B) call gap_predict_2B(GAP_2B, at(i), DATA_C)
    if (T_MB) call gap_predict_MB(GAP_MB, at(i))
enddo
CALL  SYSTEM_CLOCK(it2)
print*, "CPU time used (sec) For GP Predict: ",(it2 - it1)/10000.0

call get_rmse(at)

ENDIF ! ltest
end program

