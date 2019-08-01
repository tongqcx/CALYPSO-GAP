module io
use constants
use struct
contains

subroutine read_input()

logical :: lex
character(len=150) :: rtp,cerr, ele
integer :: lv1,lv2

inquire(file='contral',exist=lex)
if(.not.lex) then
    open(file='contral',unit=60)
    write(60,fmt=*) '    NSPECIES =  '
    write(60,fmt=*) '  NSPARSE_2B =  '
    write(60,fmt=*) '  NSPARSE_MB =  '
    write(60,fmt=*) '     SIGMA_E =  '
    write(60,fmt=*) '     SIGMA_F =  '
    write(60,fmt=*) '     SIGMA_S =  '
    write(60,fmt=*) '     D_WIDTH =  '
    write(60,fmt=*) '        Rcut =  '
    write(60,fmt=*) '        Rmin =  '
    write(60,fmt=*) '       THETA =  '
    write(60,fmt=*) '    DELTA_2B =  '
    write(60,fmt=*) '    DELTA_MB =  '
    write(60,fmt=*) 'SIGMA_JITTER =  '
    write(60,fmt=*) '    ELEMENTS =  '
    write(60,fmt=*) '    TRAINING =  '
    write(60,fmt=*) '     TESTING =  '
    close(60)
    stop
end if
!>>> INITIAL PARAMETERS LIST <<< 
    DATA_C%nspecies = 1
    DATA_C%nsparse_2B = 100
    DATA_C%nsparse_MB = 1000
    DATA_C%sigma_e = 1.d0
    DATA_C%sigma_f = 1.d0
    DATA_C%sigma_s = 1.d0
    DATA_C%d_width = 1.d0
    DATA_C%rcut = 9.d0
    DATA_C%theta_2B = 1.d0
    DATA_C%delta_2B = 1.d0
    DATA_C%delta_MB = 1.d0
    sigma_jitter = 1.0d-8
    ltrain = .true.
    ltest = .true.
    
!>>>                         <<<
open(file='contral',unit=60)
do while(.not.eof(60))
    read(60,'(a150)') rtp
    rtp=adjustl(rtp)
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
    case('nspecies')
        read(rtp(lv2:),*) DATA_C%nspecies
        DATA_C%ninteraction = DATA_C%nspecies * (DATA_C%nspecies + 1)/2.d0
        if (.not. allocated(DATA_C%elements)) allocate(DATA_C%elements(DATA_C%nspecies))
    case('nsparse_2b')
        read(rtp(lv2:),*) DATA_C%nsparse_2b
    case('nsparse_mb')
        read(rtp(lv2:),*) DATA_C%nsparse_mb
    case('sigma_e')
        read(rtp(lv2:),*) DATA_C%sigma_e
    case('sigma_f')
        read(rtp(lv2:),*) DATA_C%sigma_f
    case('sigma_s')
        read(rtp(lv2:),*) DATA_C%sigma_s
    case('rcut')
        read(rtp(lv2:),*) DATA_C%rcut
    case('rmin')
        read(rtp(lv2:),*) DATA_C%rmin
    case('theta_2b')
        read(rtp(lv2:),*) DATA_C%theta_2b
    case('delta_2b')
        read(rtp(lv2:),*) DATA_C%delta_2b
    case('delta_mb')
        read(rtp(lv2:),*) DATA_C%delta_mb
    case('d_width')
        read(rtp(lv2:),*) DATA_C%d_width
    case('sigma_jitter')
        read(rtp(lv2:),*) DATA_C%sigma_jitter
    case('elements')
        read(rtp(lv2:),*) DATA_C%elements
    case('training')
        read(rtp(lv2:),*) ltrain
    case('testing')
        read(rtp(lv2:),*) ltest
    case default
        write(*,*) "Warning: '"//trim(cerr)//"' ignored."
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
do i = 1, n_config
!    print *, i
    read(2244,*)  at(i)%natoms, at(i)%nspecies
    data_c%nf = data_c%nf + 3*at(i)%natoms
    data_c%natoms = data_c%natoms + at(i)%natoms
    call ini_structure(at(i))
    do j = 1,3
        read(2244,*) at(i)%lat(j,:)
    enddo
    at(i)%volume = volume(at(i)%lat)
    read(2244,*) at(i)%stress_ref(:)
    at(i)%stress_ref = at(i)%stress_ref/10.d0 ! kB to GPa
    at(i)%stress_ref = at(i)%stress_ref * GPa2eVPang * at(i)%volume ! to virial stress
    do j = 1, at(i)%natoms
        read(2244,*) at(i)%symbols(j), at(i)%pos(j,:), at(i)%force_ref(j,:)
        call get_ele_weights(at(i)%symbols(j),at(i)%mlp_weights(j))
    enddo
    read(2244,*) at(i)%energy_ref
    call build_neighbor(at(i), data_c)
enddo
close(2244)
data_c%ns = n_config*6
data_c%ne = n_config
data_c%nob = data_c%ne + data_c%nf + data_c%ns
allocate(data_c%ob(data_c%nob))
allocate(data_c%obe(data_c%ne))
kf = 0
do i = 1, n_config
    kf = kf + 1
    data_c%obe(i) = at(i)%energy_ref
    data_c%ob(kf) = at(i)%energy_ref
    do j = 1, at(i)%natoms
        do k1 = 1, 3
            kf = kf + 1
            data_c%ob(kf) = at(i)%force_ref(j,k1)
        enddo
    enddo
    do j = 1, 6
        kf = kf + 1
        data_c%ob(kf) = at(i)%stress_ref(j)
    enddo
enddo
print*, 'Number of structures/energy/stress:', data_c%ne, data_c%ne, data_c%ns
print*, 'Number of atoms/forces:', data_c%natoms, data_c%nf
print*, 'Number of Observable variables:',data_c%nob
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
