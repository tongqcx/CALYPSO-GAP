!###################################
! 2019.10.27
!
!
!###################################
SUBROUTINE MLP_MAIN(imode)
! for gap fitting
use gap_constants
use io
use struct
use gpr_base
use gpr_mb
! for geometry optimization
use read_module
use gap_module
use relax_module
implicit none
integer,intent(in) :: imode 
! local
integer            :: i,j,k,ii,jj,kk,k1, k2,k3
integer            :: nconfig, itype
real(dp)           :: fc_i, fc_j, rij
logical            :: alive, T_MB
integer            :: interaction_index
integer            :: ran_seed(1)

!imode = 2
if (imode==0) then
    print*,  ''
    print*, 'Geometry Optimization Using GAP potentials'
    print*,  ''
    call read_gulp()
    call relax_main_conj(NA, SPECIES, LAT, POS, STRESS, maxcycle, ftol, gtol)
    stop
else
    print*, ''
    print*, 'Fitting GAP potentials'
    print*, ''
    ran_seed=1
    call random_seed(put=ran_seed)
    call read_input()
    T_MB = DATA_C%ltrain_mb
    
    if (ltrain) then
    open(2211,file='config')
    read(2211,*) nconfig
    allocate(at(nconfig))
    close(2211)
    
    call read_structure('config',at, data_c)
    call ini_at_calc(AT, data_c)
    
    !*****************************************************
    !
    !*****************************************************
    call set_gpr_ob(at, data_c)
    call gap_ini_mb(GAP_MB, AT, ACSF, DATA_C)
    call gap_cmo_mb(GAP_MB, AT, DATA_C)
    deallocate(at)
    call gap_coeff_mb(GAP_MB, DATA_C)
    call gap_write_paras_mb(GAP_MB, ACSF)
    itype = 2
    !*****************************************************
    ENDIF  ! ltrain
    
    if (ltest) then
    print*, '*************** BEGIN PREDICTING ***************'
    if (.not.ltrain) then
        if (T_MB) call gap_read_paras_mb(GAP_MB, ACSF)
    endif
      
    open(2211,file='test')
    read(2211,*) nconfig
    allocate(at(nconfig))
    close(2211)
    call read_structure('test',at, data_c)
    
    CALL  SYSTEM_CLOCK(it1)
    !$OMP parallel do schedule(dynamic) default(shared) private(i)
    do i = 1, nconfig
        if (T_MB) call gap_predict_MB(GAP_MB, at(i), DATA_C, .true.)
    enddo
    CALL  SYSTEM_CLOCK(it2)
    print*, "CPU time used (sec) For GP Predict: ",(it2 - it1)/10000.0
    call get_rmse(at,itype)
    ENDIF ! ltest
    stop
endif
end SUBROUTINE

