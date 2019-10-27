module io
use gap_constants
use struct
contains

subroutine read_input()

logical :: lex
character(len=150) :: rtp,cerr, ele
integer :: lv1,lv2

!>>> INITIAL PARAMETERS LIST <<< 
    DATA_C%nspecies = 1
    DATA_C%d_width = 1.d0
    DATA_C%rmin = 0.5d0
    DATA_C%rcut = 6.d0
    DATA_C%ncross = 3
    DATA_C%lread_ele_weight = .false.
    sigma_jitter = 1.0d-8
    ltrain = .true.
    ltest = .true.
!{
    DATA_C%nsparse_2B = 100
    DATA_C%theta_2B = 1.d0
    DATA_C%delta_2B = 2.35d0
    DATA_C%sigma_e_2B = 1.d0
    DATA_C%sigma_f_2B = 1.d0
    DATA_C%sigma_s_2B = 1.d0
    DATA_C%ltrain_2B = .false.
!}

!{
    DATA_C%nsparse_MB = 300
    DATA_C%sparse_dis_len = 1.d0
    DATA_C%sigma_e_MB = 0.01d0
    DATA_C%sigma_f_MB = 0.1d0
    DATA_C%sigma_s_MB = 1.d0
    DATA_C%delta_MB = 1.d0
    DATA_C%sparse_method = 2
    DATA_C%sigma_atom = 1.6d0
    DATA_C%ltrain_mb = .true.
    DATA_C%lstress = .true.
!}
    
!>>>                         <<<
open(file='control',unit=60)
do while(.not.eof(60))
    read(60,'(a150)') rtp
    rtp=adjustl(rtp)
    lv1 = index(rtp,'{')
    if (lv1 == 1 ) cycle
    lv1 = index(rtp,'}')
    if (lv1 == 1 ) cycle
    lv1=index(rtp,'#')
    lv2=index(rtp,'!')
    if(lv1==0) lv1=lv2
    if(lv2<lv1.and.lv2/=0) lv1=lv2
    if(lv1>1) rtp=rtp(:lv1-1)
    cerr=rtp
    if(rtp(:1)=='#'.or.rtp(:1)=='!'.or.rtp(:1)==' ') cycle
    lv2=index(rtp,'=')+1
    lv1=lv2-2
    call u2l(rtp(:lv1))
    if(len_trim(rtp(lv2:))==0) cycle
    select case(rtp(:lv1))
    case('numberofspecies')
        read(rtp(lv2:),*) DATA_C%nspecies
        DATA_C%ninteraction = DATA_C%nspecies * (DATA_C%nspecies + 1)/2.d0
    !    print*, 'DATA_C%ninteraction',DATA_C%ninteraction
    !    stop
        if (.not. allocated(DATA_C%elements)) allocate(DATA_C%elements(DATA_C%nspecies))
        if (.not. allocated(DATA_C%elements_weight)) allocate(DATA_C%elements_weight(DATA_C%nspecies))
        if (.not. allocated(DATA_C%elements_count)) allocate(DATA_C%elements_count(DATA_C%nspecies))
!{
    case('nsparse_2b')
        read(rtp(lv2:),*) DATA_C%nsparse_2b
    case('sigma_e_2b')
        read(rtp(lv2:),*) DATA_C%sigma_e_2b
    case('sigma_f_2b')
        read(rtp(lv2:),*) DATA_C%sigma_f_2b
    case('sigma_s_2b')
        read(rtp(lv2:),*) DATA_C%sigma_s_2b
    case('theta_2b')
        read(rtp(lv2:),*) DATA_C%theta_2b
    case('delta_2b')
        read(rtp(lv2:),*) DATA_C%delta_2b
    case('ltrain_2b')
        read(rtp(lv2:),*) DATA_C%ltrain_2b
!}
!{
    case('nsparse_mb')
        read(rtp(lv2:),*) DATA_C%nsparse_mb
    case('sigmae')
        read(rtp(lv2:),*) DATA_C%sigma_e_mb
    case('sigmaf')
        read(rtp(lv2:),*) DATA_C%sigma_f_mb
    case('sigmas')
        read(rtp(lv2:),*) DATA_C%sigma_s_mb
    case('delta_mb')
        read(rtp(lv2:),*) DATA_C%delta_mb
    case('sparsecut')
        read(rtp(lv2:),*) DATA_C%sparse_dis_len
    case('isparse')
        read(rtp(lv2:),*) DATA_C%sparse_method
    case('sigma_atom')
        read(rtp(lv2:),*) DATA_C%sigma_atom
    case('ltrain_mb')
        read(rtp(lv2:),*) DATA_C%ltrain_mb
    case('lstress')
        read(rtp(lv2:),*) DATA_C%lstress
!}
    case('rcut')
        read(rtp(lv2:),*) DATA_C%rcut
    case('rmin')
        read(rtp(lv2:),*) DATA_C%rmin
    case('d_width')
        read(rtp(lv2:),*) DATA_C%d_width
    case('ncross')
        read(rtp(lv2:),*) DATA_C%ncross
    case('sigma_jitter')
        read(rtp(lv2:),*) DATA_C%sigma_jitter
    case('nameofatoms')
        read(rtp(lv2:),*) DATA_C%elements
    case('weightofatoms')
        DATA_C%lread_ele_weight = .TRUE.
        read(rtp(lv2:),*) DATA_C%elements_weight
        !print *, DATA_C%lread_ele_weight
        !print*, DATA_C%elements_weight
    case('ltrain')
        read(rtp(lv2:),*) ltrain
    case('ltest')
        read(rtp(lv2:),*) ltest
    case default
        !write(*,*) "Warning: '"//trim(cerr)//"' ignored."
    end select
end do
close(60)


if (.not. allocated(data_c%interaction_mat)) allocate(data_c%interaction_mat(data_c%nspecies, data_c%nspecies))
iindex = 0
do i = 1, data_c%nspecies
    do j = i,data_c%nspecies
        iindex = iindex + 1
        data_c%interaction_mat(i,j) = iindex
        data_c%interaction_mat(j,i) = data_c%interaction_mat(i,j)
    enddo
enddo

end subroutine read_input

subroutine u2l(ce)

character(len=*) :: ce
integer :: i

do i=1,len_trim(ce)
    if(ce(i:i)>'@'.and.ce(i:i)<'[') ce(i:i)=char(iachar(ce(i:i))+32)
end do
end subroutine u2l

SUBROUTINE read_structure(filename, at, data_c)
character(*),intent(in)                     :: filename
type(Structure),intent(inout),dimension(:)  :: at
type(DATA_TYPE),intent(inout)               :: data_c
 
integer                                     :: n_config
n_config = size(at)

open(2244,file=trim(adjustl(filename)))
read(2244,*)
data_c%nf = 0
data_c%natoms = 0
DATA_C%elements_count = 0
do i = 1, n_config
!    print *, i
    read(2244,*)  at(i)%natoms, at(i)%nspecies
    data_c%nf = data_c%nf + 3*at(i)%natoms
    data_c%natoms = data_c%natoms + at(i)%natoms
    call ini_structure(at(i), data_c)
    do j = 1,3
        read(2244,*) at(i)%lat(j,:)
    enddo
    at(i)%volume = volume(at(i)%lat)
    read(2244,*) at(i)%stress_ref(:)
    at(i)%stress_ref = at(i)%stress_ref/10.d0 ! kB to GPa
    at(i)%stress_ref = at(i)%stress_ref * GPa2eVPang * at(i)%volume ! to virial stress
    do j = 1, at(i)%natoms
        read(2244,*) at(i)%symbols(j), at(i)%pos(j,:), at(i)%force_ref(j,:)
        if (DATA_C%lread_ele_weight) then
            call get_read_weights(DATA_C, at(i)%symbols(j), at(i)%mlp_weights(j))
        else
            call get_default_weights(at(i)%symbols(j),at(i)%mlp_weights(j))
        endif
        call get_elements_count(DATA_C, at(i)%symbols(j))
    enddo
    read(2244,*) at(i)%energy_ref
    call build_neighbor(at(i), data_c)
enddo
close(2244)

print*, '*********'
do i = 1, DATA_C%nspecies
    print*, DATA_C%elements(i), DATA_C%elements_count(i), DATA_C%elements_weight(i)
enddo
print*, '*********'

END SUBROUTINE read_structure

SUBROUTINE READ_ACSF(filename, acsf)
implicit none
character(*),intent(in)                          :: filename
type(ACSF_type),intent(inout)                    :: acsf
! local
INTEGER                                          :: i

open(2244,file=trim(adjustl(filename)))
read(2244,*)  acsf%global_cutoff
read(2244,*)  acsf%nsf
allocate(acsf%sf(acsf%nsf))
do i = 1, acsf%nsf
    read(2244,*) acsf%sf(i)%ntype, acsf%sf(i)%alpha, acsf%sf(i)%cutoff
enddo
close(2244)
END SUBROUTINE READ_ACSF

end module
