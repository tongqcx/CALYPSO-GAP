! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! X  the module of gaussian process regression 
! X  for many-body interaction
! X  Qunchao Tong 2019.07.24 20:55
! X
! X
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
MODULE GPR_MB
use constants
use io
use math
use linearalgebra
use struct
use gpr_base

CONTAINS
SUBROUTINE GAP_INI_MB(GAP, AT, ACSF, DATA_C)
type(GAP_type),intent(inout)                :: GAP
type(Structure),intent(inout),dimension(:)  :: at
type(ACSF_type),intent(inout)               :: ACSF
type(DATA_type),intent(in)                  :: DATA_C

!local
n_config = size(at)

call READ_ACSF('neural.in', ACSF)
GAP%nsparse = DATA_C%nsparse_mb
GAP%dd = ACSF%NSF * 2     ! D_tot = D_topology + D_species
GAP%nglobalY = DATA_C%nob
GAP%sigma_e = DATA_C%sigma_e_mb
GAP%sigma_f = DATA_C%sigma_f_mb
GAP%sigma_s = DATA_C%sigma_s_mb
GAP%delta = DATA_C%delta_mb
GAP%sparse_method = DATA_C%sparse_method
GAP%sigma_atom = DATA_C%sigma_atom
GAP%sparse_dis_len = DATA_C%sparse_dis_len
print*, '==============================================='
print*, 'Parameters for MANY_BODY interaction'
print*, 'The size of sparse set:', GAP%nsparse
print*, 'The dimension of ACSF descriptors', GAP%dd
write(*,'(A30X3F8.5)'), 'The value of sigma E/F/S:', GAP%sigma_e, GAP%sigma_f, GAP%sigma_s
print*, 'The size of globalY:', GAP%nglobalY, DATA_C%nob
print*, 'sparse_dis_len:', GAP%sparse_dis_len
print*, 'sparse_method:',GAP%sparse_method
print*, 'sigma_atom:',GAP%sigma_atom
if ( .not. data_c%lstress) print*, 'CAUTIONS: Stress is not used in fitting GAP potentials'

CALL  SYSTEM_CLOCK(it1)
!$OMP parallel do schedule(dynamic) default(shared) private(i)
do i = 1, n_config
    call car2acsf(at(i), GAP, ACSF, DATA_C)
enddo
CALL  SYSTEM_CLOCK(it2)
print*, "CPU time used (sec) For convert COORD: ",(it2 - it1)/10000.0

allocate(GAP%theta(GAP%dd))
allocate(GAP%obe(GAP%nglobalY))
allocate(GAP%lamda(GAP%nglobalY))
allocate(GAP%lamdaobe(GAP%nglobalY, 1))

if (GAP%sparse_method == 1) then
    allocate(GAP%DescriptorX(DATA_C%natoms, GAP%dd))
    allocate(GAP%cmm(GAP%nsparse, GAP%nsparse))
    allocate(GAP%cmo(GAP%nsparse, GAP%nglobalY, 1))
    allocate(GAP%sparseX(GAP%nsparse, GAP%dd))
    allocate(GAP%sparsex_index(GAP%nsparse))
    allocate(GAP%coeff(GAP%nsparse, 1))
    spaese_index = 0
    k = 0
    do i_struc = 1, n_config
        do i_atom = 1, at(i_struc)%natoms
            k = k + 1
            GAP%descriptorx(k,:) = at(i_struc)%xx(:,i_atom)
        enddo
    enddo
    !call write_array(GAP%descriptorx,'des.dat')
    call cur_decomposition(transpose(GAP%descriptorx), GAP%sparseX_index)
    do i = 1, GAP%nsparse
        GAP%sparseX(i,:) = GAP%descriptorx(GAP%sparsex_index(i),:)
    enddo
    !call write_array(GAP%sparseX,'sparseX.dat')
    call GAP_SET_THETA(GAP%sparseX, GAP%theta)
    
    
    GAP%cmm = 0.0
    do i = 1, GAP%nsparse
        GAP%cmm(i,i) = 1.d0 + DATA_C%sigma_jitter
        do  j = i+1 , GAP%nsparse
            GAP%cmm(i,j) = covariance_MB(GAP%delta, GAP%sparseX(i,:),GAP%sparseX(j,:), GAP%theta)
            GAP%cmm(j,i) = GAP%cmm(i,j)
        enddo
    enddo
elseif (GAP%sparse_method == 2 .or. GAP%sparse_method == 3) then
    call GAP_SPARSE(GAP,AT)
else
    print*, 'Unknow sparse_method'
endif
!call write_array(GAP%cmm,'cmm.dat')



!!!!!!!!!!!!!!!!!!!!!!
! initial GAP%lamda
!!!!!!!!!!!!!!!!!!!!!!
kf = 0
do i = 1, n_config
    kf = kf + 1
    GAP%lamda(kf) = (GAP%sigma_e * (sqrt(1.d0 * at(i)%natoms)))**2
    GAP%lamda(kf) = sqrt(1.d0/GAP%lamda(kf))
    do j = 1,at(i)%natoms
        do k1 = 1,3
            kf = kf + 1
            GAP%lamda(kf) = GAP%sigma_f**2
            GAP%lamda(kf) = sqrt(1.d0/GAP%lamda(kf))
        enddo
    enddo
    if (data_c%lstress) then
        do j = 1,6
            kf = kf + 1
            GAP%lamda(kf) = GAP%sigma_s**2
            GAP%lamda(kf) = sqrt(1.d0/GAP%lamda(kf))
        enddo
    endif
enddo
print*, 'GAP INI FINISHED'
END SUBROUTINE GAP_INI_MB

SUBROUTINE GAP_SPARSE(GAP, AT)
type(GAP_type),intent(inout)                :: GAP
type(Structure),intent(inout),dimension(:)  :: at

! local
REAL(DP),dimension(:),allocatable           :: theta_temp
REAL(DP),dimension(:,:),allocatable         :: MM
REAL(DP),dimension(:,:),allocatable         :: calc_det
logical                                     :: lexit, has_upper_bound
REAL(DP)                                    :: my_length

!spaese_dis_len = 1.d0
!open(3310, file='ddd.dat')
!do i1 = 1, size(AT)
!    write(3310,*), i1
!    do i2 = 1, at(i1)%natoms
!        do j = 1, GAP%dd
!            write(3310,'(F20.10,$)') at(i1)%xx(j,i2)
!        enddo
!        write(3310,*)
!    enddo
!enddo
!close(3310)
!stop


n_config = size(AT)
allocate(theta_temp(GAP%dd))
allocate(MM(max_mm_len,GAP%dd))
lexit = .false.
MM = 0.d0
k = 0
do i1 = 1,n_config
    do i2 = 1,at(i1)%natoms
        call random_number(rx)
        if (rx < 0.1d0) then
            k = k + 1
            mm(k,:) = at(i1)%xx(:,i2)
        endif
        if (k >= max_mm_len) then
            lexit = .true.
            exit
        endif
    enddo
    if (lexit) exit
enddo
!call GAP_SET_THETA(mm, GAP%theta)
deallocate(MM)
allocate(MM(max_mm_len,GAP%dd))
theta_temp = 1.d0
has_upper_bound = .FALSE.
CALL  SYSTEM_CLOCK(it1)
do while (.true.)
    lexit = .false.
    mm = 0.d0
    k = 1
    mm(k,:) = at(1)%xx(:,1)
    do i1 = 1, n_config
        do i2 = 1, at(i1)%natoms
            ladd = 0
            do j = 1, k
                my_length = length(at(i1)%xx(:,i2) - mm(j,:),theta_temp)
                if (my_length >= GAP%sparse_dis_len) ladd = ladd + 1
            enddo
            if (ladd == k) then
                k = k + 1
                mm(k,:) = at(i1)%xx(:,i2)
            endif
            if (k == max_mm_len) then
                lexit = .true.
                !print*, i1,'nconfig',i2,'atom'
                !print*, k,"is large than max_mm_len",max_mm_len
                exit
            endif
        enddo
        if (lexit) exit
    enddo
    call GAP_SET_THETA(mm(1:k,:), GAP%theta)
!call write_array(GAP%theta, 'theta,dat')
!call write_array(mm(1:k,:),'MM.dat')
    GAP%nsparse = k
    GAP%theta = GAP%theta * GAP%sigma_atom
    if (.not. allocated(GAP%cmm))  allocate(GAP%cmm(GAP%nsparse, GAP%nsparse))
    if (.not. allocated(calc_det))  allocate(calc_det(GAP%nsparse, GAP%nsparse))
    GAP%cmm = 0.0
    calc_det  = 0.0
    do i = 1, GAP%nsparse
        GAP%cmm(i,i) = 1.d0 + DATA_C%sigma_jitter
        do  j = i+1 , GAP%nsparse
            GAP%cmm(i,j) = covariance_MB(GAP%delta, mm(i,:),mm(j,:), GAP%theta)
            GAP%cmm(j,i) = GAP%cmm(i,j)
        enddo
    enddo
!    call write_array(GAP%cmm, 'cmm.datx')
    calc_det = GAP%cmm
    det_cmm = my_det(calc_det)**(1.d0/GAP%nsparse)
    if (GAP%sparse_method == 2) then
        !if (det_cmm > Inverse_error_min .and. det_cmm < Inverse_error_max) then
        if (GAP%nsparse > 120 .and. GAP%nsparse < 200) then
            print*, "The number of atomic environment in Sparse set:",GAP%nsparse
            print*, 'sparse_dis_len:', GAP%sparse_dis_len
            print*, 'Det of CMM:',det_cmm
            exit
        !elseif (det_cmm <= Inverse_error_min) then
        elseif (GAP%nsparse >= 200) then

            if (has_upper_bound) then
                GAP%sparse_dis_len = (GAP%sparse_dis_len + upper_bound)/2.d0
            else
                GAP%sparse_dis_len = GAP%sparse_dis_len * 2.d0
            endif
            print *,'Increasing sparse_dis_cut to:',GAP%sparse_dis_len

            if (allocated(GAP%cmm))  deallocate(GAP%cmm)
            if (allocated(calc_det))  deallocate(calc_det)
        !elseif (det_cmm >= Inverse_error_max )then
        elseif (GAP%nsparse <= 120 )then

            has_upper_bound = .TRUE.
            upper_bound = GAP%sparse_dis_len
            GAP%sparse_dis_len = GAP%sparse_dis_len * 0.75d0
            print *,'*Decreasing sparse_dis_cut to:',GAP%sparse_dis_len

            if (allocated(GAP%cmm))  deallocate(GAP%cmm)
            if (allocated(calc_det))  deallocate(calc_det)
        endif
    else
        if (det_cmm > Inverse_error_min) then
            print*, "The number of atomic environment in Sparse set:",GAP%nsparse
            print*, 'sparse_dis_len:', GAP%sparse_dis_len
            print*, 'Det of CMM:',det_cmm
            exit
        else
            GAP%sparse_dis_len = GAP%sparse_dis_len + 0.2d0
            print *,'Increasing sparse_dis_cut',GAP%sparse_dis_len, GAP%nsparse
            if (allocated(GAP%cmm))  deallocate(GAP%cmm)
            if (allocated(calc_det))  deallocate(calc_det)
        endif
    endif
enddo
CALL  SYSTEM_CLOCK(it2)
print*, 'Sparsing FINISHED',(it2 - it1)/10000.0,'Seconds'
allocate(GAP%cmo(GAP%nsparse, GAP%nglobalY, 1))
allocate(GAP%sparseX(GAP%nsparse, GAP%dd))
allocate(GAP%coeff(GAP%nsparse, 1))
do i = 1, GAP%nsparse
    GAP%sparseX(i,:) = MM(i,:)
enddo
deallocate(MM)
deallocate(theta_temp)
deallocate(calc_det)
END SUBROUTINE GAP_SPARSE

SUBROUTINE GAP_WRITE_PARAS_MB(GAP)
type(GAP_type),intent(in)             :: GAP

open(2234,file='gap_paras_mb.dat')
write(2234,*) GAP%nsparse, GAP%dd
write(2234,*)
write(2234,*)
write(2234,*)
do i = 1, GAP%dd
    write(2234,'(F25.10,$)') GAP%theta(i)
enddo
write(2234,*)
do i = 1, GAP%nsparse
    do j = 1, GAP%dd
        write(2234,'(F25.10,$)') GAP%sparseX(i,j)
    enddo
    write(2234,*)
enddo
do i = 1, GAP%nsparse
    write(2234,'(F25.10,$)') GAP%coeff(i,1)
enddo
write(2234,*)
close(2234)
END SUBROUTINE GAP_WRITE_PARAS_MB

SUBROUTINE GAP_READ_PARAS_MB(GAP, ACSF)
type(GAP_type),intent(inout)             :: GAP
type(ACSF_type),intent(inout)            :: ACSF

! local
logical                                  :: alive

inquire(file="gap_paras_mb.dat",exist=alive)
if(.not.alive) then
    print*, "FILES-> gap_paras_mb.dat does not exist!"
    stop
endif
open(2234,file='gap_paras_mb.dat')
read(2234,*) GAP%nsparse, GAP%dd
allocate(GAP%sparseX(GAP%nsparse, GAP%dd))
allocate(GAP%theta(GAP%dd))
allocate(GAP%coeff(GAP%nsparse, 1))
read(2234,*)
read(2234,*)
read(2234,*)
read(2234,*) GAP%theta(:)
do i = 1, GAP%nsparse
    read(2234,*)  GAP%sparseX(i,:)
enddo
read(2234,*) GAP%coeff(:,1)
close(2234)
call READ_ACSF('neural.in', ACSF)
GAP%delta = 1.d0
END SUBROUTINE GAP_READ_PARAS_MB

SUBROUTINE GAP_COEFF_MB(GAP, DATA_C)
type(GAP_type),intent(inout)             :: GAP
type(DATA_type),intent(in)               :: DATA_C

do i = 1, DATA_C%nob
    GAP%lamdaobe(i,1) = GAP%lamda(i) * DATA_C%ob(i)
enddo
call gpr(GAP%cmm, GAP%cmo(:,:,1), GAP%lamdaobe(:,1), GAP%coeff(:,1))
END SUBROUTINE GAP_COEFF_MB

SUBROUTINE GAP_CMO_MB(GAP, at, DATA_C)
type(GAP_type),intent(inout)             :: GAP
type(Structure),intent(in),dimension(:)  :: at
type(DATA_type),intent(in)               :: DATA_C


!local
!REAL(DP),allocatable,dimension(:)        :: cov
INTEGER                                  :: i,j,k1,kf
REAL(DP)                                 :: cov(3 * max_atoms + 7) ! the maximum number of atoms is 500  modified 2019.08.15
CALL  SYSTEM_CLOCK(it1)
!allocate(cov(DATA_C%nob))
!!$OMP parallel do schedule(dynamic) default(shared) private(i_sparse, i_struc ,i_ob, kf, cov)
!$OMP parallel private(i_sparse, i_struc ,i_ob, kf, cov)
!$OMP DO
do i_sparse = 1, GAP%nsparse
    kf = 1
!    print*, i_sparse
    do i_struc = 1, DATA_C%ne
        call new_COV(GAP%delta, GAP%sparseX(i_sparse,:), GAP%theta, AT(i_struc)%xx, AT(i_struc)%dxdy, AT(i_struc)%strs, cov(:))
        if (data_c%lstress) then
            do i_ob = 1, 3*at(i_struc)%natoms + 7
                GAP%cmo(i_sparse, kf, 1) = cov( i_ob)
                kf = kf + 1
            enddo
        else
            do i_ob = 1, 3*at(i_struc)%natoms + 1
                GAP%cmo(i_sparse, kf, 1) = cov( i_ob)
                kf = kf + 1
            enddo
        endif
    enddo
enddo
!$OMP END PARALLEL 

!deallocate(cov)
call matmuldiag_T(GAP%cmo(:,:,1), GAP%lamda)
CALL  SYSTEM_CLOCK(it2)
print*, 'GAP_MB CMO FINISHED',(it2 - it1)/10000.0,'Seconds'
END SUBROUTINE GAP_CMO_MB

SUBROUTINE GAP_PREDICT_MB(GAP,AT, DATA_C, lcar2acsf)
type(GAP_type),intent(in)                   :: GAP
type(Structure),intent(inout)               :: at
type(DATA_type),intent(in)                  :: DATA_C
logical,intent(in)                          :: lcar2acsf

if (lcar2acsf)  call car2acsf(at, GAP, ACSF, DATA_C)

if (.not. allocated(at%kk))   allocate(at%kk(at%natoms, GAP%dd))
if (.not. allocated(at%ckm))  allocate(at%ckm(at%natoms, GAP%nsparse))
if (.not. allocated(at%dedg)) allocate(at%dedg(at%natoms, GAP%dd))

at%kk = 0.d0
at%ckm = 0.d0
at%dedg = 0.d0

do i = 1,at%natoms
    at%kk(i,:) = at%xx(:,i)
enddo
call get_cov(GAP%delta, at%kk, GAP%sparseX, GAP%theta, at%ckm)

at%atomic_energy = matmul(at%ckm, GAP%coeff(:,1))
at%energy_cal_mb = sum(at%atomic_energy)
do i = 1, at%natoms
    do j = 1, GAP%dd
        do k = 1, GAP%nsparse
            at%dedg(i,j) = at%dedg(i,j) - 1.d0 * (at%kk(i,j) - GAP%sparsex(k,j))/GAP%theta(j)**2 &
            * at%ckm(i,k) * GAP%coeff(k,1)
        enddo
    enddo
enddo

!!for force calculation
at%force_cal_mb = 0.d0
do i = 1, at%natoms
    do n = 1, at%natoms
        do j = 1,3
            do k = 1, GAP%dd
                at%force_cal_mb(i,j) = at%force_cal_mb(i,j) - at%dedg(n,k) * at%dxdy(k,n,i,j)
            enddo
        enddo
    enddo
enddo
!! for stress calculation
at%stress_cal_mb = 0.d0
at%volume = volume(at%lat)
k1 = 1
do i = 1, 3
    do j = i, 3
        do n = 1, at%natoms
            do k = 1, GAP%dd
                at%stress_cal_mb(k1) = at%stress_cal_mb(k1) - at%dedg(n,k) * at%strs(i,j,k,n)
            enddo
        enddo
        k1 = k1 + 1
    enddo
enddo
!at%stress_cal = at%stress_cal * (1.0/GPa2eVPang) * 10.0 / at%volume
END SUBROUTINE GAP_PREDICT_MB

SUBROUTINE   new_cov(delta, x, theta, xx, dxdy, strs, covf)
implicit none
REAL(DP),intent(in)                      :: delta   
real(DP),intent(in),dimension(:)         :: x, theta
real(DP),intent(in),dimension(:,:)       :: xx
real(DP),intent(in),dimension(:,:,:,:)   :: dxdy, strs
real(DP),intent(out),dimension(:)        :: covf
real(DP),allocatable,dimension(:,:)      :: for
real(DP)                                 :: ene
real(DP)                                 :: factor
real(DP)                                 :: stress(6)
integer                                  :: i,j,k1,k2,k,na,nf,kf

na = size(xx,2)
nf = size(xx,1)
allocate(for(na,3))
ene = 0.d0
for = 0.d0
stress = 0.d0
factor = 0.d0
covf = 0.d0

do i = 1,na  ! number of atoms
    ene = ene + covariance_mb(delta, x, xx(:,i), theta)
    do j = 1,nf  ! number of symmetry function
        factor = (x(j) - xx(j,i))/theta(j)**2 * covariance_mb(delta, x,xx(:,i), theta)
        do k1 = 1,na
            do k2 = 1,3
                !write(111,'(4I4X3F10.6)') i,j,k1,k2,(x(j) - xx(j,1,i))/theta(j)**2,covariance(x,xx(:,1,i)),dxdy(j,i,k1,k2)
                for(k1,k2) = for(k1,k2) - factor*dxdy(j,i,k1,k2)
            enddo
        enddo
        stress(1) = stress(1) - factor * strs(1,1,j,i)
        stress(2) = stress(2) - factor * strs(1,2,j,i)
        stress(3) = stress(3) - factor * strs(1,3,j,i)
        stress(4) = stress(4) - factor * strs(2,2,j,i)
        stress(5) = stress(5) - factor * strs(2,3,j,i)
        stress(6) = stress(6) - factor * strs(3,3,j,i)
    enddo
enddo
kf = 2
covf(1) = ene
do k1 = 1,na
    do k2 = 1,3
        covf(kf) = for(k1,k2)
        kf = kf + 1
    enddo
enddo
do i = 1,6
    covf(kf) = stress(i)
    kf = kf + 1
enddo
deallocate(for)
end subroutine new_cov

SUBROUTINE CAR2ACSF(at, GAP, ACSF, DATA_C)

implicit real(DP) (a-h,o-z)

type(Structure),intent(inout)          :: at
type(GAP_type),intent(in)              :: GAP
type(ACSF_type),intent(in)             :: acsf
type(DATA_type),intent(in)             :: DATA_C

!local
REAL(DP),dimension(3)                  :: xyz, xyz_j, xyz_k

allocate(at%xx(GAP%dd, at%natoms))
allocate(at%dxdy(GAP%dd, at%natoms, at%natoms, 3))
allocate(at%strs(3, 3, GAP%dd, at%natoms))
! @@@@ this three array must be initial
at%xx = 0.d0
at%dxdy = 0.d0
at%strs = 0.d0
! @@@@

rmin = 0.5d0
nnn = ACSF%nsf
do ii = 1, nnn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!G1 = SUM_j{exp(-alpha*rij**2)*fc(rij)}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (ACSF%sf(ii)%ntype.eq.1) then
        cutoff = ACSF%sf(ii)%cutoff
        alpha = ACSF%sf(ii)%alpha
        do i = 1, at%natoms
            ! ******************
            ! Be careful!!!
            ! i_type loop from 1 to data_c%nspecies not at%nspecies
            ! 2019.09.04 STUPID!!!
            ! ******************
            do i_type = 1, data_c%nspecies
                do i_neighbor = 1, at%atom(i)%count(i_type)
                    rij = at%atom(i)%neighbor(i_type,i_neighbor,4)
                    if (rij.gt.cutoff) cycle
                    xyz = at%atom(i)%neighbor(i_type,i_neighbor,1:3)
                    n = int(at%atom(i)%neighbor(i_type,i_neighbor,5))
                    weights = at%mlp_weights(n)
                    deltaxj = -1.d0*(at%atom(i)%pos(1) - xyz(1))
                    deltayj = -1.d0*(at%atom(i)%pos(2) - xyz(2))
                    deltazj = -1.d0*(at%atom(i)%pos(3) - xyz(3))
                    drijdxi = -1.d0*deltaxj/rij
                    drijdyi = -1.d0*deltayj/rij
                    drijdzi = -1.d0*deltazj/rij
                    drijdxj = -1.d0*drijdxi
                    drijdyj = -1.d0*drijdyi
                    drijdzj = -1.d0*drijdzi
                    fcutij = 0.5d0 * (dcos(pi*rij/cutoff) + 1.d0)
                    temp1=0.5d0*(-dsin(pi*rij/cutoff))*(pi/cutoff)
                    dfcutijdxi=temp1*drijdxi
                    dfcutijdyi=temp1*drijdyi
                    dfcutijdzi=temp1*drijdzi
                    dfcutijdxj=-1.d0*dfcutijdxi
                    dfcutijdyj=-1.d0*dfcutijdyi
                    dfcutijdzj=-1.d0*dfcutijdzi
        !            if (i==2 .and. ii==2) print*, dexp(-1.d0*alpha*rij**2)*fcutij
                    at%xx(ii,i) = at%xx(ii,i) + dexp(-1.d0*alpha*rij**2)*fcutij
                    at%xx(ii + nnn, i) = at%xx(ii + nnn, i) + dexp(-1.d0*alpha*rij**2)*fcutij * weights !!!!!!! 
                    !write(*,'(I4X4F25.10)') i, rij, at%xx(ii,i),at%xx(ii + nnn,i)  , weights

                    !dxx/dx
                    temp1=-2.d0*alpha*rij*dexp(-1.d0*alpha*rij**2)*fcutij
                    temp2= dexp(-1.d0*alpha*rij**2)

                    at%dxdy(ii,i,i,1)=at%dxdy(ii,i,i,1)+(drijdxi*temp1 + temp2*dfcutijdxi)
                    at%dxdy(ii+nnn,i,i,1)=at%dxdy(ii+nnn,i,i,1) + &
                    (drijdxi*temp1+ temp2*dfcutijdxi) * weights
                
                    temp3=drijdxj*temp1 + temp2*dfcutijdxj
                    at%dxdy(ii,i,n,1)=at%dxdy(ii,i,n,1)+temp3
                
                    temp4=temp3*weights
                    at%dxdy(ii + nnn,i,n,1)=at%dxdy(ii + nnn,i,n,1)+temp4
                
                    at%strs(1,1,ii,i)=at%strs(1,1,ii,i)+deltaxj*temp3
                    at%strs(2,1,ii,i)=at%strs(2,1,ii,i)+deltayj*temp3
                    at%strs(3,1,ii,i)=at%strs(3,1,ii,i)+deltazj*temp3
                
                    at%strs(1,1,ii+nnn,i)=at%strs(1,1,ii+nnn,i)+deltaxj*temp4
                    at%strs(2,1,ii+nnn,i)=at%strs(2,1,ii+nnn,i)+deltayj*temp4
                    at%strs(3,1,ii+nnn,i)=at%strs(3,1,ii+nnn,i)+deltazj*temp4
                    !dxx/dy
                    at%dxdy(ii,i,i,2)=at%dxdy(ii,i,i,2)+(drijdyi*temp1+temp2*dfcutijdyi)
                    at%dxdy(ii+nnn,i,i,2)=at%dxdy(ii+nnn,i,i,2)+(drijdyi*temp1+temp2*dfcutijdyi)*weights
                    temp3= drijdyj*temp1 + temp2*dfcutijdyj
                    at%dxdy(ii,i,n,2)=at%dxdy(ii,i,n,2)+temp3
                    temp4=temp3 * weights
                    at%dxdy(ii + nnn,i,n,2)=at%dxdy(ii + nnn ,i,n,2)+temp4
                
                    at%strs(1,2,ii,i)=at%strs(1,2,ii,i)+deltaxj*temp3
                    at%strs(2,2,ii,i)=at%strs(2,2,ii,i)+deltayj*temp3
                    at%strs(3,2,ii,i)=at%strs(3,2,ii,i)+deltazj*temp3
                
                    at%strs(1,2,ii + nnn,i)=at%strs(1,2,ii + nnn,i)+deltaxj*temp4
                    at%strs(2,2,ii + nnn,i)=at%strs(2,2,ii + nnn,i)+deltayj*temp4
                    at%strs(3,2,ii + nnn,i)=at%strs(3,2,ii + nnn,i)+deltazj*temp4
                    !dxx/dz
                    at%dxdy(ii,i,i,3)=at%dxdy(ii,i,i,3)+&
                           (drijdzi*temp1&
                          + temp2*dfcutijdzi)
                    at%dxdy(ii + nnn,i,i,3)=at%dxdy(ii + nnn,i,i,3)+&
                           (drijdzi*temp1&
                          + temp2*dfcutijdzi)*weights
                    temp3=drijdzj*temp1 + temp2*dfcutijdzj
                    at%dxdy(ii,i,n,3)=at%dxdy(ii,i,n,3)+temp3
                    temp4=temp3*weights
                    at%dxdy(ii + nnn,i,n,3)=at%dxdy(ii + nnn,i,n,3)+temp4
                
                    at%strs(1,3,ii,i)=at%strs(1,3,ii,i)+deltaxj*temp3
                    at%strs(2,3,ii,i)=at%strs(2,3,ii,i)+deltayj*temp3
                    at%strs(3,3,ii,i)=at%strs(3,3,ii,i)+deltazj*temp3
                
                    at%strs(1,3,ii + nnn,i)=at%strs(1,3,ii + nnn,i)+deltaxj*temp4
                    at%strs(2,3,ii + nnn,i)=at%strs(2,3,ii + nnn,i)+deltayj*temp4
                    at%strs(3,3,ii + nnn,i)=at%strs(3,3,ii + nnn,i)+deltazj*temp4
                enddo ! i_neighbor
            enddo ! i_type
        enddo ! i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! lamda = 1
! eta = 1
! G2 = SUM_jk{(1+lamda*costheta_ijk)^eta*
! exp(-alpha*(rij**2+rik**2+rjk**2))*fc(rij)*fc(rik)*fc(rjk)}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    elseif (ACSF%sf(ii)%ntype.eq.2) then
        cutoff = ACSF%sf(ii)%cutoff
        alpha = ACSF%sf(ii)%alpha
!        print*, 'cutoff',cutoff,'alpha',alpha
        do i = 1, at%natoms
!            lllll = 0
            do j_type = 1, data_c%nspecies
                do j_neighbor = 1, at%atom(i)%count(j_type)
                    rij = at%atom(i)%neighbor(j_type,j_neighbor,4)
                    if (rij.gt.cutoff) cycle
                    xyz_j = at%atom(i)%neighbor(j_type,j_neighbor,1:3)
                    n = int(at%atom(i)%neighbor(j_type,j_neighbor,5))
                    weights_j = at%mlp_weights(n)
                    !print*,  xyz_j,'j'
                    deltaxj = -1.d0*(at%atom(i)%pos(1) - xyz_j(1))
                    deltayj = -1.d0*(at%atom(i)%pos(2) - xyz_j(2))
                    deltazj = -1.d0*(at%atom(i)%pos(3) - xyz_j(3))
                    drijdxi = -1.d0*deltaxj/rij
                    drijdyi = -1.d0*deltayj/rij
                    drijdzi = -1.d0*deltazj/rij
                    drijdxj = -1.d0*drijdxi
                    drijdyj = -1.d0*drijdyi
                    drijdzj = -1.d0*drijdzi
                    drijdxk = 0.d0
                    drijdyk = 0.d0
                    drijdzk = 0.d0
                    fcutij=0.5d0*(dcos(pi*rij/cutoff)+1.d0)
                    temp1=0.5d0*(-dsin(pi*rij/cutoff))*(pi/cutoff)
                    dfcutijdxi=temp1*drijdxi
                    dfcutijdyi=temp1*drijdyi
                    dfcutijdzi=temp1*drijdzi
                    dfcutijdxj=-1.d0*dfcutijdxi
                    dfcutijdyj=-1.d0*dfcutijdyi
                    dfcutijdzj=-1.d0*dfcutijdzi
                    dfcutijdxk=0.0d0
                    dfcutijdyk=0.0d0
                    dfcutijdzk=0.0d0
                    do k_type = 1, data_c%nspecies
                        do k_neighbor = 1, at%atom(i)%count(k_type)
                            !if ((k_type <= j_type) .and. (k_neighbor <= j_neighbor)) cycle
                            ! ******************
                            ! Be careful
                            ! ******************
                            if (k_type < j_type) cycle
                            if ((k_type==j_type) .and. (k_neighbor <= j_neighbor)) cycle
                            rik = at%atom(i)%neighbor(k_type,k_neighbor,4)
                            if (rik.gt.cutoff) cycle
                   !         lllll = lllll + 1
                            xyz_k = at%atom(i)%neighbor(k_type,k_neighbor,1:3)
                    !        print*, xyz_k,'k'
                            m = int(at%atom(i)%neighbor(k_type,k_neighbor,5))
                            weights_k = at%mlp_weights(m)

                            deltaxk = -1.d0*(at%atom(i)%pos(1) - xyz_k(1))
                            deltayk = -1.d0*(at%atom(i)%pos(2) - xyz_k(2))
                            deltazk = -1.d0*(at%atom(i)%pos(3) - xyz_k(3))
                            drikdxi = -deltaxk/rik
                            drikdyi = -deltayk/rik
                            drikdzi = -deltazk/rik
                            drikdxk = -1.d0*drikdxi
                            drikdyk = -1.d0*drikdyi
                            drikdzk = -1.d0*drikdzi
                            drikdxj = 0.d0
                            drikdyj = 0.d0
                            drikdzj = 0.d0
                            fcutik=0.5d0*(dcos(pi*rik/cutoff)+1.d0)
                            temp1=0.5d0*(-dsin(pi*rik/cutoff))*(pi/cutoff)
                            dfcutikdxi=temp1*drikdxi
                            dfcutikdyi=temp1*drikdyi
                            dfcutikdzi=temp1*drikdzi
                            dfcutikdxj=0.0d0
                            dfcutikdyj=0.0d0
                            dfcutikdzj=0.0d0
                            dfcutikdxk=-1.d0*dfcutikdxi
                            dfcutikdyk=-1.d0*dfcutikdyi
                            dfcutikdzk=-1.d0*dfcutikdzi
                            rjk = (xyz_j(1) - xyz_k(1))**2 + (xyz_j(2) - xyz_k(2))**2 + (xyz_j(3) - xyz_k(3))**2
                            rjk = dsqrt(rjk)

                            if (rjk.gt.cutoff) cycle  ! CAUTAINS STUPID!!!
                            if (rjk < Rmin) then
                                print*, 'Rjk', rjk,' smaller than Rmin'
                                stop
                            endif
                            drjkdxj = (xyz_j(1) - xyz_k(1))/rjk
                            drjkdyj = (xyz_j(2) - xyz_k(2))/rjk
                            drjkdzj = (xyz_j(3) - xyz_k(3))/rjk
                            drjkdxk = -1.d0*drjkdxj
                            drjkdyk = -1.d0*drjkdyj
                            drjkdzk = -1.d0*drjkdzj
                            drjkdxi = 0.d0
                            drjkdyi = 0.d0
                            drjkdzi = 0.d0
                            fcutjk=0.5d0*(dcos(pi*rjk/cutoff)+1.d0)
                            temp1=0.5d0*(-dsin(pi*rjk/cutoff))*(pi/cutoff)
                            dfcutjkdxj=temp1*drjkdxj
                            dfcutjkdyj=temp1*drjkdyj
                            dfcutjkdzj=temp1*drjkdzj
                            dfcutjkdxk=-1.d0*dfcutjkdxj
                            dfcutjkdyk=-1.d0*dfcutjkdyj
                            dfcutjkdzk=-1.d0*dfcutjkdzj
                            dfcutjkdxi=0.0d0
                            dfcutjkdyi=0.0d0
                            dfcutjkdzi=0.0d0

                            f=rjk**2 - rij**2 -rik**2
                            g=-2.d0*rij*rik
                            costheta=f/g
                            !!!!  2^1-eta (1+lamda coseta_ijk)^eta 
                            !!!!  eta = 1 lamda = +1.d0
                            costheta=1.d0 + costheta
                            dfdxi=-2.d0*rij*drijdxi - 2.d0*rik*drikdxi
                            dfdyi=-2.d0*rij*drijdyi - 2.d0*rik*drikdyi
                            dfdzi=-2.d0*rij*drijdzi - 2.d0*rik*drikdzi

                            dfdxj=2.d0*rjk*drjkdxj - 2.d0*rij*drijdxj
                            dfdyj=2.d0*rjk*drjkdyj - 2.d0*rij*drijdyj
                            dfdzj=2.d0*rjk*drjkdzj - 2.d0*rij*drijdzj

                            dfdxk=2.d0*rjk*drjkdxk - 2.d0*rik*drikdxk
                            dfdyk=2.d0*rjk*drjkdyk - 2.d0*rik*drikdyk
                            dfdzk=2.d0*rjk*drjkdzk - 2.d0*rik*drikdzk

                            dgdxi=-2.d0*(drijdxi*rik + rij*drikdxi)
                            dgdyi=-2.d0*(drijdyi*rik + rij*drikdyi)
                            dgdzi=-2.d0*(drijdzi*rik + rij*drikdzi)

                            dgdxj=-2.d0*drijdxj*rik
                            dgdyj=-2.d0*drijdyj*rik
                            dgdzj=-2.d0*drijdzj*rik

                            dgdxk=-2.d0*rij*drikdxk
                            dgdyk=-2.d0*rij*drikdyk
                            dgdzk=-2.d0*rij*drikdzk

                            temp1=1.d0/g**2
                            dcosthetadxi=(dfdxi*g - f*dgdxi)*temp1
                            dcosthetadyi=(dfdyi*g - f*dgdyi)*temp1
                            dcosthetadzi=(dfdzi*g - f*dgdzi)*temp1
                            dcosthetadxj=(dfdxj*g - f*dgdxj)*temp1
                            dcosthetadyj=(dfdyj*g - f*dgdyj)*temp1
                            dcosthetadzj=(dfdzj*g - f*dgdzj)*temp1
                            dcosthetadxk=(dfdxk*g - f*dgdxk)*temp1
                            dcosthetadyk=(dfdyk*g - f*dgdyk)*temp1
                            dcosthetadzk=(dfdzk*g - f*dgdzk)*temp1
                            expxyz=dexp(-alpha*(rij**2+rik**2+rjk**2))
                            temp1=-alpha*2.0d0*expxyz
                            dexpxyzdxi=(rij*drijdxi+rik*drikdxi+rjk*drjkdxi)*temp1
                            dexpxyzdyi=(rij*drijdyi+rik*drikdyi+rjk*drjkdyi)*temp1
                            dexpxyzdzi=(rij*drijdzi+rik*drikdzi+rjk*drjkdzi)*temp1
                            dexpxyzdxj=(rij*drijdxj+rik*drikdxj+rjk*drjkdxj)*temp1
                            dexpxyzdyj=(rij*drijdyj+rik*drikdyj+rjk*drjkdyj)*temp1
                            dexpxyzdzj=(rij*drijdzj+rik*drikdzj+rjk*drjkdzj)*temp1
                            dexpxyzdxk=(rij*drijdxk+rik*drikdxk+rjk*drjkdxk)*temp1
                            dexpxyzdyk=(rij*drijdyk+rik*drikdyk+rjk*drjkdyk)*temp1
                            dexpxyzdzk=(rij*drijdzk+rik*drikdzk+rjk*drjkdzk)*temp1
                            at%xx(ii,i)=at%xx(ii,i)+costheta*expxyz*fcutij*fcutik*fcutjk
                            at%xx(ii + nnn,i)=at%xx(ii + nnn,i)+&
                            costheta*expxyz*fcutij*fcutik*fcutjk*weights_j*weights_k

                            temp1=(dcosthetadxi*expxyz*fcutij*fcutik*fcutjk&
                                  +costheta*dexpxyzdxi*fcutij*fcutik*fcutjk&
                                  +costheta*expxyz*dfcutijdxi*fcutik*fcutjk&
                                  +costheta*expxyz*fcutij*dfcutikdxi*fcutjk&
                                  +costheta*expxyz*fcutij*fcutik*dfcutjkdxi)
                            temp2=(dcosthetadxj*expxyz*fcutij*fcutik*fcutjk&
                                  +costheta*dexpxyzdxj*fcutij*fcutik*fcutjk&
                                  +costheta*expxyz*dfcutijdxj*fcutik*fcutjk&
                                  +costheta*expxyz*fcutij*dfcutikdxj*fcutjk&
                                  +costheta*expxyz*fcutij*fcutik*dfcutjkdxj)
                            temp3=(dcosthetadxk*expxyz*fcutij*fcutik*fcutjk&
                                  +costheta*dexpxyzdxk*fcutij*fcutik*fcutjk&
                                  +costheta*expxyz*dfcutijdxk*fcutik*fcutjk&
                                  +costheta*expxyz*fcutij*dfcutikdxk*fcutjk&
                                  +costheta*expxyz*fcutij*fcutik*dfcutjkdxk)
                            temp4 = temp1 * weights_j * weights_k
                            temp5 = temp2 * weights_j * weights_k
                            temp6 = temp3 * weights_j * weights_k
                            at%dxdy(ii,i,i,1)=at%dxdy(ii,i,i,1)+temp1
                            at%dxdy(ii,i,n,1)=at%dxdy(ii,i,n,1)+temp2
                            at%dxdy(ii,i,m,1)=at%dxdy(ii,i,m,1)+temp3
                            at%dxdy(ii + nnn,i,i,1)=at%dxdy(ii + nnn,i,i,1)+temp4
                            at%dxdy(ii + nnn,i,n,1)=at%dxdy(ii + nnn,i,n,1)+temp5
                            at%dxdy(ii + nnn,i,m,1)=at%dxdy(ii + nnn,i,m,1)+temp6

                            at%strs(1,1,ii,i)=at%strs(1,1,ii,i)+deltaxj*temp2+deltaxk*temp3
                            at%strs(2,1,ii,i)=at%strs(2,1,ii,i)+deltayj*temp2+deltayk*temp3
                            at%strs(3,1,ii,i)=at%strs(3,1,ii,i)+deltazj*temp2+deltazk*temp3
                            at%strs(1,1,ii + nnn,i)=at%strs(1,1,ii + nnn,i)+deltaxj*temp5+deltaxk*temp6
                            at%strs(2,1,ii + nnn,i)=at%strs(2,1,ii + nnn,i)+deltayj*temp5+deltayk*temp6
                            at%strs(3,1,ii + nnn,i)=at%strs(3,1,ii + nnn,i)+deltazj*temp5+deltazk*temp6
                            ! dxxii/dy_i
                            temp1=(dcosthetadyi*expxyz*fcutij*fcutik*fcutjk&
                                  +costheta*dexpxyzdyi*fcutij*fcutik*fcutjk&
                                  +costheta*expxyz*dfcutijdyi*fcutik*fcutjk&
                                  +costheta*expxyz*fcutij*dfcutikdyi*fcutjk&
                                  +costheta*expxyz*fcutij*fcutik*dfcutjkdyi)
                            temp2=(dcosthetadyj*expxyz*fcutij*fcutik*fcutjk&
                                  +costheta*dexpxyzdyj*fcutij*fcutik*fcutjk&
                                  +costheta*expxyz*dfcutijdyj*fcutik*fcutjk&
                                  +costheta*expxyz*fcutij*dfcutikdyj*fcutjk&
                                  +costheta*expxyz*fcutij*fcutik*dfcutjkdyj)
                            temp3=(dcosthetadyk*expxyz*fcutij*fcutik*fcutjk&
                                  +costheta*dexpxyzdyk*fcutij*fcutik*fcutjk&
                                  +costheta*expxyz*dfcutijdyk*fcutik*fcutjk&
                                  +costheta*expxyz*fcutij*dfcutikdyk*fcutjk&
                                  +costheta*expxyz*fcutij*fcutik*dfcutjkdyk)
                            temp4 = temp1 * weights_j * weights_k
                            temp5 = temp2 * weights_j * weights_k
                            temp6 = temp3 * weights_j * weights_k
                            at%dxdy(ii,i,i,2)=at%dxdy(ii,i,i,2)+temp1
                            at%dxdy(ii,i,n,2)=at%dxdy(ii,i,n,2)+temp2
                            at%dxdy(ii,i,m,2)=at%dxdy(ii,i,m,2)+temp3

                            at%dxdy(ii + nnn,i,i,2)=at%dxdy(ii + nnn,i,i,2)+temp4
                            at%dxdy(ii + nnn,i,n,2)=at%dxdy(ii + nnn,i,n,2)+temp5
                            at%dxdy(ii + nnn,i,m,2)=at%dxdy(ii + nnn,i,m,2)+temp6
                            at%strs(1,2,ii,i)=at%strs(1,2,ii,i)+deltaxj*temp2+deltaxk*temp3
                            at%strs(2,2,ii,i)=at%strs(2,2,ii,i)+deltayj*temp2+deltayk*temp3
                            at%strs(3,2,ii,i)=at%strs(3,2,ii,i)+deltazj*temp2+deltazk*temp3

                            at%strs(1,2,ii + nnn,i)=at%strs(1,2,ii + nnn,i)+deltaxj*temp5+deltaxk*temp6
                            at%strs(2,2,ii + nnn,i)=at%strs(2,2,ii + nnn,i)+deltayj*temp5+deltayk*temp6
                            at%strs(3,2,ii + nnn,i)=at%strs(3,2,ii + nnn,i)+deltazj*temp5+deltazk*temp6
                            ! dxxii/dz_i
                            temp1=(dcosthetadzi*expxyz*fcutij*fcutik*fcutjk&
                                  +costheta*dexpxyzdzi*fcutij*fcutik*fcutjk&
                                  +costheta*expxyz*dfcutijdzi*fcutik*fcutjk&
                                  +costheta*expxyz*fcutij*dfcutikdzi*fcutjk&
                                  +costheta*expxyz*fcutij*fcutik*dfcutjkdzi)
                            temp2=(dcosthetadzj*expxyz*fcutij*fcutik*fcutjk&
                                  +costheta*dexpxyzdzj*fcutij*fcutik*fcutjk&
                                  +costheta*expxyz*dfcutijdzj*fcutik*fcutjk&
                                  +costheta*expxyz*fcutij*dfcutikdzj*fcutjk&
                                  +costheta*expxyz*fcutij*fcutik*dfcutjkdzj)
                            temp3=(dcosthetadzk*expxyz*fcutij*fcutik*fcutjk&
                                  +costheta*dexpxyzdzk*fcutij*fcutik*fcutjk&
                                  +costheta*expxyz*dfcutijdzk*fcutik*fcutjk&
                                  +costheta*expxyz*fcutij*dfcutikdzk*fcutjk&
                                  +costheta*expxyz*fcutij*fcutik*dfcutjkdzk)
                            temp4 = temp1 * weights_j * weights_k
                            temp5 = temp2 * weights_j * weights_k
                            temp6 = temp3 * weights_j * weights_k
                            at%dxdy(ii,i,i,3)=at%dxdy(ii,i,i,3)+temp1
                            at%dxdy(ii,i,n,3)=at%dxdy(ii,i,n,3)+temp2
                            at%dxdy(ii,i,m,3)=at%dxdy(ii,i,m,3)+temp3

                            at%dxdy(ii + nnn,i,i,3)=at%dxdy(ii + nnn,i,i,3)+temp4
                            at%dxdy(ii + nnn,i,n,3)=at%dxdy(ii + nnn,i,n,3)+temp5
                            at%dxdy(ii + nnn,i,m,3)=at%dxdy(ii + nnn,i,m,3)+temp6
                            at%strs(1,3,ii,i)=at%strs(1,3,ii,i)+deltaxj*temp2+deltaxk*temp3
                            at%strs(2,3,ii,i)=at%strs(2,3,ii,i)+deltayj*temp2+deltayk*temp3
                            at%strs(3,3,ii,i)=at%strs(3,3,ii,i)+deltazj*temp2+deltazk*temp3

                            at%strs(1,3,ii + nnn,i)=at%strs(1,3,ii + nnn,i)+deltaxj*temp5+deltaxk*temp6
                            at%strs(2,3,ii + nnn,i)=at%strs(2,3,ii + nnn,i)+deltayj*temp5+deltayk*temp6
                            at%strs(3,3,ii + nnn,i)=at%strs(3,3,ii + nnn,i)+deltazj*temp5+deltazk*temp6
                        enddo ! k_neighbor
                    enddo ! k_type       
                enddo ! j_neighbor
            enddo ! j_type
!            print*, 'lllll',lllll
        enddo ! i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!G3 = SUM_j{exp(-alpha*(rij-rshift)**2)*fc(rij)}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    elseif (ACSF%sf(ii)%ntype.eq.3) then
        cutoff = ACSF%sf(ii)%cutoff
        rshift = ACSF%sf(ii)%alpha
        alpha = 4.d0
        do i = 1, at%natoms
            do i_type = 1, data_c%nspecies
                do i_neighbor = 1, at%atom(i)%count(i_type)
                    rij = at%atom(i)%neighbor(i_type,i_neighbor,4)

                    if (rij.gt.cutoff) cycle
                    xyz = at%atom(i)%neighbor(i_type,i_neighbor,1:3)
                    n = int(at%atom(i)%neighbor(i_type,i_neighbor,5))
                    weights = at%mlp_weights(n)
                    deltaxj = -1.d0*(at%atom(i)%pos(1) - xyz(1))
                    deltayj = -1.d0*(at%atom(i)%pos(2) - xyz(2))
                    deltazj = -1.d0*(at%atom(i)%pos(3) - xyz(3))
                    drijdxi = -1.d0*deltaxj/rij
                    drijdyi = -1.d0*deltayj/rij
                    drijdzi = -1.d0*deltazj/rij
                    drijdxj = -1.d0*drijdxi
                    drijdyj = -1.d0*drijdyi
                    drijdzj = -1.d0*drijdzi
                    fcutij = 0.5d0 * (dcos(pi*rij/cutoff) + 1.d0)
                    temp1=0.5d0*(-dsin(pi*rij/cutoff))*(pi/cutoff)
                    dfcutijdxi=temp1*drijdxi
                    dfcutijdyi=temp1*drijdyi
                    dfcutijdzi=temp1*drijdzi
                    dfcutijdxj=-1.d0*dfcutijdxi
                    dfcutijdyj=-1.d0*dfcutijdyi
                    dfcutijdzj=-1.d0*dfcutijdzi
                    at%xx(ii,i)=at%xx(ii,i)+dexp(-1.d0*alpha*(rij-rshift)**2)*fcutij
                    at%xx(ii + nnn,i)=at%xx(ii + nnn,i)+dexp(-1.d0*alpha*(rij-rshift)**2)*fcutij*weights
                    temp1=-2.d0*alpha*(rij-rshift)
                    temp2=dexp(-1.d0*alpha*(rij-rshift)**2)
                    ! dxx/dx
                    at%dxdy(ii,i,i,1)=at%dxdy(ii,i,i,1)+&
                           (temp1*drijdxi*temp2*fcutij&
                          + temp2*dfcutijdxi)

                    at%dxdy(ii + nnn,i,i,1)=at%dxdy(ii + nnn,i,i,1)+&
                           (temp1*drijdxi*temp2*fcutij&
                          + temp2*dfcutijdxi)*weights
                    temp3=temp1*drijdxj*temp2*fcutij + temp2*dfcutijdxj
                    at%dxdy(ii,i,n,1)=at%dxdy(ii,i,n,1)+temp3
                    temp4 = temp3 * weights
                    at%dxdy(ii + nnn,i,n,1)=at%dxdy(ii + nnn,i,n,1)+temp4
                    at%strs(1,1,ii,i)=at%strs(1,1,ii,i)+deltaxj*temp3
                    at%strs(2,1,ii,i)=at%strs(2,1,ii,i)+deltayj*temp3
                    at%strs(3,1,ii,i)=at%strs(3,1,ii,i)+deltazj*temp3
                    at%strs(1,1,ii + nnn,i)=at%strs(1,1,ii + nnn,i)+deltaxj*temp4
                    at%strs(2,1,ii + nnn,i)=at%strs(2,1,ii + nnn,i)+deltayj*temp4
                    at%strs(3,1,ii + nnn,i)=at%strs(3,1,ii + nnn,i)+deltazj*temp4
                    ! dxx/dy
                    at%dxdy(ii,i,i,2)=at%dxdy(ii,i,i,2)+&
                           (temp1*drijdyi*temp2*fcutij&
                          + temp2*dfcutijdyi)
                    at%dxdy(ii + nnn,i,i,2)=at%dxdy(ii + nnn,i,i,2)+&
                           (temp1*drijdyi*temp2*fcutij&
                          + temp2*dfcutijdyi)*weights
                    temp3= temp1*drijdyj*temp2*fcutij + temp2*dfcutijdyj
                    at%dxdy(ii,i,n,2)=at%dxdy(ii,i,n,2)+temp3
                    temp4 = temp3 * weights
                    at%dxdy(ii + nnn,i,n,2)=at%dxdy(ii + nnn,i,n,2)+temp4
                    at%strs(1,2,ii,i)=at%strs(1,2,ii,i)+deltaxj*temp3
                    at%strs(2,2,ii,i)=at%strs(2,2,ii,i)+deltayj*temp3
                    at%strs(3,2,ii,i)=at%strs(3,2,ii,i)+deltazj*temp3

                    at%strs(1,2,ii + nnn,i)=at%strs(1,2,ii + nnn,i)+deltaxj*temp4
                    at%strs(2,2,ii + nnn,i)=at%strs(2,2,ii + nnn,i)+deltayj*temp4
                    at%strs(3,2,ii + nnn,i)=at%strs(3,2,ii + nnn,i)+deltazj*temp4
                    ! dxx/dz
                    at%dxdy(ii,i,i,3)=at%dxdy(ii,i,i,3)+&
                           (temp1*drijdzi*temp2*fcutij&
                          + temp2*dfcutijdzi)
                    at%dxdy(ii + nnn,i,i,3)=at%dxdy(ii + nnn,i,i,3)+&
                           (temp1*drijdzi*temp2*fcutij&
                          + temp2*dfcutijdzi)*weights
                    temp3=temp1*drijdzj*temp2*fcutij + temp2*dfcutijdzj
                    at%dxdy(ii,i,n,3)=at%dxdy(ii,i,n,3)+temp3
                    temp4 = temp3 * weights
                    at%dxdy(ii + nnn,i,n,3)=at%dxdy(ii + nnn,i,n,3)+temp4
                    at%strs(1,3,ii,i)=at%strs(1,3,ii,i)+deltaxj*temp3
                    at%strs(2,3,ii,i)=at%strs(2,3,ii,i)+deltayj*temp3
                    at%strs(3,3,ii,i)=at%strs(3,3,ii,i)+deltazj*temp3

                    at%strs(1,3,ii + nnn,i)=at%strs(1,3,ii + nnn,i)+deltaxj*temp4
                    at%strs(2,3,ii + nnn,i)=at%strs(2,3,ii + nnn,i)+deltayj*temp4
                    at%strs(3,3,ii + nnn,i)=at%strs(3,3,ii + nnn,i)+deltazj*temp4
                enddo ! i_neighbor
            enddo ! i_type
        enddo ! i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! lamda = -1.d0
! eta = 1
! G2 = SUM_jk{(1+lamda*costheta_ijk)^eta*
! exp(-alpha*(rij**2+rik**2+rjk**2))*fc(rij)*fc(rik)*fc(rjk)}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    elseif (ACSF%sf(ii)%ntype.eq.4) then
        cutoff = ACSF%sf(ii)%cutoff
        alpha = ACSF%sf(ii)%alpha
        do i = 1, at%natoms
            do j_type = 1, data_c%nspecies
                do j_neighbor = 1, at%atom(i)%count(j_type)
                    rij = at%atom(i)%neighbor(j_type,j_neighbor,4)
                    if (rij.gt.cutoff) cycle
                    xyz_j = at%atom(i)%neighbor(j_type,j_neighbor,1:3)
                    n = int(at%atom(i)%neighbor(j_type,j_neighbor,5))
                    weights_j = at%mlp_weights(n)
                    deltaxj = -1.d0*(at%atom(i)%pos(1) - xyz_j(1))
                    deltayj = -1.d0*(at%atom(i)%pos(2) - xyz_j(2))
                    deltazj = -1.d0*(at%atom(i)%pos(3) - xyz_j(3))
                    drijdxi = -1.d0*deltaxj/rij
                    drijdyi = -1.d0*deltayj/rij
                    drijdzi = -1.d0*deltazj/rij
                    drijdxj = -1.d0*drijdxi
                    drijdyj = -1.d0*drijdyi
                    drijdzj = -1.d0*drijdzi
                    drijdxk = 0.d0
                    drijdyk = 0.d0
                    drijdzk = 0.d0

                    fcutij=0.5d0*(dcos(pi*rij/cutoff)+1.d0)
                    temp1=0.5d0*(-dsin(pi*rij/cutoff))*(pi/cutoff)
                    dfcutijdxi=temp1*drijdxi
                    dfcutijdyi=temp1*drijdyi
                    dfcutijdzi=temp1*drijdzi
                    dfcutijdxj=-1.d0*dfcutijdxi
                    dfcutijdyj=-1.d0*dfcutijdyi
                    dfcutijdzj=-1.d0*dfcutijdzi
                    dfcutijdxk=0.0d0
                    dfcutijdyk=0.0d0
                    dfcutijdzk=0.0d0
                    do k_type = 1, data_c%nspecies
                        do k_neighbor = 1, at%atom(i)%count(k_type)
                            !if ((k_type <= j_type) .and. (k_neighbor <= j_neighbor)) cycle
                            if (k_type < j_type) cycle
                            if ((k_type==j_type) .and. (k_neighbor <= j_neighbor)) cycle
                            rik = at%atom(i)%neighbor(k_type,k_neighbor,4)
                            if (rik.gt.cutoff) cycle
                            xyz_k = at%atom(i)%neighbor(k_type,k_neighbor,1:3)
                            m = int(at%atom(i)%neighbor(k_type,k_neighbor,5))
                            weights_k = at%mlp_weights(m)

                            deltaxk = -1.d0*(at%atom(i)%pos(1) - xyz_k(1))
                            deltayk = -1.d0*(at%atom(i)%pos(2) - xyz_k(2))
                            deltazk = -1.d0*(at%atom(i)%pos(3) - xyz_k(3))
                            drikdxi = -deltaxk/rik
                            drikdyi = -deltayk/rik
                            drikdzi = -deltazk/rik
                            drikdxk = -1.d0*drikdxi
                            drikdyk = -1.d0*drikdyi
                            drikdzk = -1.d0*drikdzi
                            drikdxj = 0.d0
                            drikdyj = 0.d0
                            drikdzj = 0.d0
                            fcutik=0.5d0*(dcos(pi*rik/cutoff)+1.d0)
                            temp1=0.5d0*(-dsin(pi*rik/cutoff))*(pi/cutoff)
                            dfcutikdxi=temp1*drikdxi
                            dfcutikdyi=temp1*drikdyi
                            dfcutikdzi=temp1*drikdzi
                            dfcutikdxj=0.0d0
                            dfcutikdyj=0.0d0
                            dfcutikdzj=0.0d0
                            dfcutikdxk=-1.d0*dfcutikdxi
                            dfcutikdyk=-1.d0*dfcutikdyi
                            dfcutikdzk=-1.d0*dfcutikdzi
                            rjk = (xyz_j(1) - xyz_k(1))**2 + (xyz_j(2) - xyz_k(2))**2 + (xyz_j(3) - xyz_k(3))**2
                            rjk = dsqrt(rjk)

                            if (rjk.gt.cutoff) cycle  ! Be careful STUPID!!!
                            if (rjk < Rmin) then
                                print*, 'Rjk', rjk,' smaller than Rmin'
                                stop
                            endif
                            drjkdxj = (xyz_j(1) - xyz_k(1))/rjk
                            drjkdyj = (xyz_j(2) - xyz_k(2))/rjk
                            drjkdzj = (xyz_j(3) - xyz_k(3))/rjk
                            drjkdxk = -1.d0*drjkdxj
                            drjkdyk = -1.d0*drjkdyj
                            drjkdzk = -1.d0*drjkdzj
                            drjkdxi = 0.d0
                            drjkdyi = 0.d0
                            drjkdzi = 0.d0
                            fcutjk=0.5d0*(dcos(pi*rjk/cutoff)+1.d0)
                            temp1=0.5d0*(-dsin(pi*rjk/cutoff))*(pi/cutoff)
                            dfcutjkdxj=temp1*drjkdxj
                            dfcutjkdyj=temp1*drjkdyj
                            dfcutjkdzj=temp1*drjkdzj
                            dfcutjkdxk=-1.d0*dfcutjkdxj
                            dfcutjkdyk=-1.d0*dfcutjkdyj
                            dfcutjkdzk=-1.d0*dfcutjkdzj
                            dfcutjkdxi=0.0d0
                            dfcutjkdyi=0.0d0
                            dfcutjkdzi=0.0d0

                            f=rjk**2 - rij**2 -rik**2
                            g=-2.d0*rij*rik
                            costheta=f/g
                            costheta=1.d0 - costheta  ! avoid negative values
                            dfdxi=-2.d0*rij*drijdxi - 2.d0*rik*drikdxi
                            dfdyi=-2.d0*rij*drijdyi - 2.d0*rik*drikdyi
                            dfdzi=-2.d0*rij*drijdzi - 2.d0*rik*drikdzi

                            dfdxj=2.d0*rjk*drjkdxj - 2.d0*rij*drijdxj
                            dfdyj=2.d0*rjk*drjkdyj - 2.d0*rij*drijdyj
                            dfdzj=2.d0*rjk*drjkdzj - 2.d0*rij*drijdzj

                            dfdxk=2.d0*rjk*drjkdxk - 2.d0*rik*drikdxk
                            dfdyk=2.d0*rjk*drjkdyk - 2.d0*rik*drikdyk
                            dfdzk=2.d0*rjk*drjkdzk - 2.d0*rik*drikdzk

                            dgdxi=-2.d0*(drijdxi*rik + rij*drikdxi)
                            dgdyi=-2.d0*(drijdyi*rik + rij*drikdyi)
                            dgdzi=-2.d0*(drijdzi*rik + rij*drikdzi)

                            dgdxj=-2.d0*drijdxj*rik
                            dgdyj=-2.d0*drijdyj*rik
                            dgdzj=-2.d0*drijdzj*rik

                            dgdxk=-2.d0*rij*drikdxk
                            dgdyk=-2.d0*rij*drikdyk
                            dgdzk=-2.d0*rij*drikdzk

                            temp1=1.d0/g**2
                            !!!! Be careful costheta = 1.d0 - costheta 2019.07.25
                            dcosthetadxi=-1.d0 * (dfdxi*g - f*dgdxi)*temp1  
                            dcosthetadyi=-1.d0 * (dfdyi*g - f*dgdyi)*temp1 
                            dcosthetadzi=-1.d0 * (dfdzi*g - f*dgdzi)*temp1 
                            dcosthetadxj=-1.d0 * (dfdxj*g - f*dgdxj)*temp1 
                            dcosthetadyj=-1.d0 * (dfdyj*g - f*dgdyj)*temp1 
                            dcosthetadzj=-1.d0 * (dfdzj*g - f*dgdzj)*temp1 
                            dcosthetadxk=-1.d0 * (dfdxk*g - f*dgdxk)*temp1 
                            dcosthetadyk=-1.d0 * (dfdyk*g - f*dgdyk)*temp1 
                            dcosthetadzk=-1.d0 * (dfdzk*g - f*dgdzk)*temp1 

                            expxyz=dexp(-alpha*(rij**2+rik**2+rjk**2))
                            temp1=-alpha*2.0d0*expxyz
                            dexpxyzdxi=(rij*drijdxi+rik*drikdxi+rjk*drjkdxi)*temp1
                            dexpxyzdyi=(rij*drijdyi+rik*drikdyi+rjk*drjkdyi)*temp1
                            dexpxyzdzi=(rij*drijdzi+rik*drikdzi+rjk*drjkdzi)*temp1
                            dexpxyzdxj=(rij*drijdxj+rik*drikdxj+rjk*drjkdxj)*temp1
                            dexpxyzdyj=(rij*drijdyj+rik*drikdyj+rjk*drjkdyj)*temp1
                            dexpxyzdzj=(rij*drijdzj+rik*drikdzj+rjk*drjkdzj)*temp1
                            dexpxyzdxk=(rij*drijdxk+rik*drikdxk+rjk*drjkdxk)*temp1
                            dexpxyzdyk=(rij*drijdyk+rik*drikdyk+rjk*drjkdyk)*temp1
                            dexpxyzdzk=(rij*drijdzk+rik*drikdzk+rjk*drjkdzk)*temp1

                            at%xx(ii,i)=at%xx(ii,i)+costheta*expxyz*fcutij*fcutik*fcutjk
                            at%xx(ii + nnn,i)=at%xx(ii + nnn,i)+&
                            costheta*expxyz*fcutij*fcutik*fcutjk*weights_j*weights_k

                            temp1=(dcosthetadxi*expxyz*fcutij*fcutik*fcutjk&
                                  +costheta*dexpxyzdxi*fcutij*fcutik*fcutjk&
                                  +costheta*expxyz*dfcutijdxi*fcutik*fcutjk&
                                  +costheta*expxyz*fcutij*dfcutikdxi*fcutjk&
                                  +costheta*expxyz*fcutij*fcutik*dfcutjkdxi)
                            temp2=(dcosthetadxj*expxyz*fcutij*fcutik*fcutjk&
                                  +costheta*dexpxyzdxj*fcutij*fcutik*fcutjk&
                                  +costheta*expxyz*dfcutijdxj*fcutik*fcutjk&
                                  +costheta*expxyz*fcutij*dfcutikdxj*fcutjk&
                                  +costheta*expxyz*fcutij*fcutik*dfcutjkdxj)
                            temp3=(dcosthetadxk*expxyz*fcutij*fcutik*fcutjk&
                                  +costheta*dexpxyzdxk*fcutij*fcutik*fcutjk&
                                  +costheta*expxyz*dfcutijdxk*fcutik*fcutjk&
                                  +costheta*expxyz*fcutij*dfcutikdxk*fcutjk&
                                  +costheta*expxyz*fcutij*fcutik*dfcutjkdxk)
                            temp4 = temp1 * weights_j * weights_k
                            temp5 = temp2 * weights_j * weights_k
                            temp6 = temp3 * weights_j * weights_k
                            at%dxdy(ii,i,i,1)=at%dxdy(ii,i,i,1)+temp1
                            at%dxdy(ii,i,n,1)=at%dxdy(ii,i,n,1)+temp2
                            at%dxdy(ii,i,m,1)=at%dxdy(ii,i,m,1)+temp3
                            at%dxdy(ii + nnn,i,i,1)=at%dxdy(ii + nnn,i,i,1)+temp4
                            at%dxdy(ii + nnn,i,n,1)=at%dxdy(ii + nnn,i,n,1)+temp5
                            at%dxdy(ii + nnn,i,m,1)=at%dxdy(ii + nnn,i,m,1)+temp6

                            at%strs(1,1,ii,i)=at%strs(1,1,ii,i)+deltaxj*temp2+deltaxk*temp3
                            at%strs(2,1,ii,i)=at%strs(2,1,ii,i)+deltayj*temp2+deltayk*temp3
                            at%strs(3,1,ii,i)=at%strs(3,1,ii,i)+deltazj*temp2+deltazk*temp3
                            at%strs(1,1,ii + nnn,i)=at%strs(1,1,ii + nnn,i)+deltaxj*temp5+deltaxk*temp6
                            at%strs(2,1,ii + nnn,i)=at%strs(2,1,ii + nnn,i)+deltayj*temp5+deltayk*temp6
                            at%strs(3,1,ii + nnn,i)=at%strs(3,1,ii + nnn,i)+deltazj*temp5+deltazk*temp6
                            ! dxxii/dy_i
                            temp1=(dcosthetadyi*expxyz*fcutij*fcutik*fcutjk&
                                  +costheta*dexpxyzdyi*fcutij*fcutik*fcutjk&
                                  +costheta*expxyz*dfcutijdyi*fcutik*fcutjk&
                                  +costheta*expxyz*fcutij*dfcutikdyi*fcutjk&
                                  +costheta*expxyz*fcutij*fcutik*dfcutjkdyi)
                            temp2=(dcosthetadyj*expxyz*fcutij*fcutik*fcutjk&
                                  +costheta*dexpxyzdyj*fcutij*fcutik*fcutjk&
                                  +costheta*expxyz*dfcutijdyj*fcutik*fcutjk&
                                  +costheta*expxyz*fcutij*dfcutikdyj*fcutjk&
                                  +costheta*expxyz*fcutij*fcutik*dfcutjkdyj)
                            temp3=(dcosthetadyk*expxyz*fcutij*fcutik*fcutjk&
                                  +costheta*dexpxyzdyk*fcutij*fcutik*fcutjk&
                                  +costheta*expxyz*dfcutijdyk*fcutik*fcutjk&
                                  +costheta*expxyz*fcutij*dfcutikdyk*fcutjk&
                                  +costheta*expxyz*fcutij*fcutik*dfcutjkdyk)
                            temp4 = temp1 * weights_j * weights_k
                            temp5 = temp2 * weights_j * weights_k
                            temp6 = temp3 * weights_j * weights_k
                            at%dxdy(ii,i,i,2)=at%dxdy(ii,i,i,2)+temp1
                            at%dxdy(ii,i,n,2)=at%dxdy(ii,i,n,2)+temp2
                            at%dxdy(ii,i,m,2)=at%dxdy(ii,i,m,2)+temp3

                            at%dxdy(ii + nnn,i,i,2)=at%dxdy(ii + nnn,i,i,2)+temp4
                            at%dxdy(ii + nnn,i,n,2)=at%dxdy(ii + nnn,i,n,2)+temp5
                            at%dxdy(ii + nnn,i,m,2)=at%dxdy(ii + nnn,i,m,2)+temp6
                            at%strs(1,2,ii,i)=at%strs(1,2,ii,i)+deltaxj*temp2+deltaxk*temp3
                            at%strs(2,2,ii,i)=at%strs(2,2,ii,i)+deltayj*temp2+deltayk*temp3
                            at%strs(3,2,ii,i)=at%strs(3,2,ii,i)+deltazj*temp2+deltazk*temp3

                            at%strs(1,2,ii + nnn,i)=at%strs(1,2,ii + nnn,i)+deltaxj*temp5+deltaxk*temp6
                            at%strs(2,2,ii + nnn,i)=at%strs(2,2,ii + nnn,i)+deltayj*temp5+deltayk*temp6
                            at%strs(3,2,ii + nnn,i)=at%strs(3,2,ii + nnn,i)+deltazj*temp5+deltazk*temp6
                            ! dxxii/dz_i
                            temp1=(dcosthetadzi*expxyz*fcutij*fcutik*fcutjk&
                                  +costheta*dexpxyzdzi*fcutij*fcutik*fcutjk&
                                  +costheta*expxyz*dfcutijdzi*fcutik*fcutjk&
                                  +costheta*expxyz*fcutij*dfcutikdzi*fcutjk&
                                  +costheta*expxyz*fcutij*fcutik*dfcutjkdzi)
                            temp2=(dcosthetadzj*expxyz*fcutij*fcutik*fcutjk&
                                  +costheta*dexpxyzdzj*fcutij*fcutik*fcutjk&
                                  +costheta*expxyz*dfcutijdzj*fcutik*fcutjk&
                                  +costheta*expxyz*fcutij*dfcutikdzj*fcutjk&
                                  +costheta*expxyz*fcutij*fcutik*dfcutjkdzj)
                            temp3=(dcosthetadzk*expxyz*fcutij*fcutik*fcutjk&
                                  +costheta*dexpxyzdzk*fcutij*fcutik*fcutjk&
                                  +costheta*expxyz*dfcutijdzk*fcutik*fcutjk&
                                  +costheta*expxyz*fcutij*dfcutikdzk*fcutjk&
                                  +costheta*expxyz*fcutij*fcutik*dfcutjkdzk)
                            temp4 = temp1 * weights_j * weights_k
                            temp5 = temp2 * weights_j * weights_k
                            temp6 = temp3 * weights_j * weights_k
                            at%dxdy(ii,i,i,3)=at%dxdy(ii,i,i,3)+temp1
                            at%dxdy(ii,i,n,3)=at%dxdy(ii,i,n,3)+temp2
                            at%dxdy(ii,i,m,3)=at%dxdy(ii,i,m,3)+temp3

                            at%dxdy(ii + nnn,i,i,3)=at%dxdy(ii + nnn,i,i,3)+temp4
                            at%dxdy(ii + nnn,i,n,3)=at%dxdy(ii + nnn,i,n,3)+temp5
                            at%dxdy(ii + nnn,i,m,3)=at%dxdy(ii + nnn,i,m,3)+temp6
                            at%strs(1,3,ii,i)=at%strs(1,3,ii,i)+deltaxj*temp2+deltaxk*temp3
                            at%strs(2,3,ii,i)=at%strs(2,3,ii,i)+deltayj*temp2+deltayk*temp3
                            at%strs(3,3,ii,i)=at%strs(3,3,ii,i)+deltazj*temp2+deltazk*temp3

                            at%strs(1,3,ii + nnn,i)=at%strs(1,3,ii + nnn,i)+deltaxj*temp5+deltaxk*temp6
                            at%strs(2,3,ii + nnn,i)=at%strs(2,3,ii + nnn,i)+deltayj*temp5+deltayk*temp6
                            at%strs(3,3,ii + nnn,i)=at%strs(3,3,ii + nnn,i)+deltazj*temp5+deltazk*temp6
                        enddo ! k_neighbor
                    enddo ! k_type       
                enddo ! j_neighbor
            enddo ! j_type
        enddo ! i
    else
        print *, 'Unknown function type',ii, ACSF%sf(ii)%ntype
    endif
enddo  ! types
!print*, 'CCC',at%xx(:,1)
END SUBROUTINE 

END MODULE
