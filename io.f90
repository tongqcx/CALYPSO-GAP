module io
use constants
use struct
contains

subroutine read_input()

logical :: lex
character(len=150) :: rtp,cerr, ele
integer :: lv1,lv2

inquire(file='control',exist=lex)
if(.not.lex) then
    open(file='control',unit=60)
    write(60,fmt=*) '    NSPECIES =  '
    write(60,fmt=*) '     NSPARSE =  '
    write(60,fmt=*) '     SIGMA_E =  '
    write(60,fmt=*) '     D_WIDTH =  '
    write(60,fmt=*) '        Rcut =  '
    write(60,fmt=*) '        Rmin =  '
    write(60,fmt=*) '       THETA =  '
    write(60,fmt=*) '       DELTA =  '
    write(60,fmt=*) 'SIGMA_JITTER =  '
    write(60,fmt=*) '    ELEMENTS =  '
    write(60,fmt=*) '    TRAINING =  '
    write(60,fmt=*) '     TESTING =  '
    close(60)
    stop
end if
!>>> INITIAL PARAMETERS LIST <<< 
    nspecies = 1
    nsparse = 100
    sigma_e = 1.d0
    d_width = 1.d0
    rcut = 9.d0
    theta = 1.d0
    delta = 1.d0
    sigma_jitter = 1.0d-8
    ltrain = .true.
    ltest = .true.
    
!>>>                         <<<
open(file='control',unit=60)
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
        read(rtp(lv2:),*) nspecies
        ninteraction = nspecies * (nspecies + 1)/2.d0
    case('nsparse')
        read(rtp(lv2:),*) nsparse
    case('sigma_e')
        read(rtp(lv2:),*) sigma_e
    case('rcut')
        read(rtp(lv2:),*) rcut
    case('rmin')
        read(rtp(lv2:),*) rmin
    case('theta')
        read(rtp(lv2:),*) theta
    case('delta')
        read(rtp(lv2:),*) delta
    case('d_width')
        read(rtp(lv2:),*) d_width
    case('sigma_jitter')
        read(rtp(lv2:),*) sigma_jitter
    case('elements')
        read(rtp(lv2:),*) elements
    case('training')
        read(rtp(lv2:),*) ltrain
    case('testing')
        read(rtp(lv2:),*) ltest
    case default
        write(*,*) "Warning: '"//trim(cerr)//"' ignored."
    end select
    if (.not. allocated(elements)) allocate(elements(nspecies))
end do
close(60)

end subroutine read_input

subroutine u2l(ce)

character(len=*) :: ce
integer :: i

do i=1,len_trim(ce)
    if(ce(i:i)>'@'.and.ce(i:i)<'[') ce(i:i)=char(iachar(ce(i:i))+32)
end do
end subroutine u2l

SUBROUTINE read_structure(filename, at)
character(*),intent(in)                     :: filename
type(Structure),intent(inout),dimension(:)  :: at
 
integer                                     :: n_config
n_config = size(at)

open(2244,file=trim(adjustl(filename)))
read(2244,*)
data_c%nf = 0
do i = 1, n_config
!    print *, i
    read(2244,*)  na, nspecies
    data_c%nf = data_c%nf + 3*na
    call ini_structure(at(i), na, nspecies)
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
    call build_neighbor(at(i), elements)
enddo
data_c%ns = n_config*6
data_c%ne = n_config
data_c%nob = data_c%ne + data_c%nf + data_c%ns
close(2244)
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
END SUBROUTINE READ_ACSF

end module
