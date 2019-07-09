module io
use constants
implicit none
integer                                :: nsparse
integer                                :: nspecies, ninteraction
REAL(DP)                               :: theta, delta, d_width, sigma_jitter
REAL(DP)                               :: sigma_e, sigma_f, sigma_s
character(2),allocatable,dimension(:)  :: elements
REAL(DP)                               :: Rcut, Rmin
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


end module
