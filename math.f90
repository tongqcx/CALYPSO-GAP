module math
use constants
contains
function vectorlength(vc)
    real(8) :: vc(3),vectorlength
    vectorlength=sqrt(vc(1)**2+vc(2)**2+vc(3)**2)
end function

function recipvector(lat)
    real(8),intent(in) :: lat(:,:)
    real(8) :: recipvector(3,3)

    recipvector(:,1)=crossp(lat(:,2),lat(:,3))
    recipvector(:,2)=crossp(lat(:,3),lat(:,1))
    recipvector(:,3)=crossp(lat(:,1),lat(:,2))
    recipvector=recipvector/volume(lat)*pi*2._8

end function

function volume(lat)
    real(8),intent(in) :: lat(:,:)
    real(8) :: volume 

    volume=abs(sum(lat(:,1)*crossp(lat(:,2),lat(:,3))))

end function

function crossp(va,vb)
    real(8),intent(in) :: va(3),vb(3)
    real(8) :: crossp(3)

    crossp(1)=va(2)*vb(3)-va(3)*vb(2)
    crossp(2)=va(3)*vb(1)-va(1)*vb(3)
    crossp(3)=va(1)*vb(2)-va(2)*vb(1)
end function

!function MIV_num_columns( inchar ) result(num)
!IMPLICIT none
!character(*)                            :: inchar
!integer(i4b)                            :: num
!integer(i4b)                            :: length
!character(:),allocatable                :: ctp
!integer(i4b)                            :: i
!
!ctp = trim(adjustl(inchar))
!length = len(ctp)
!num = 1
!do i = 2,length-1
!    if( ctp(i:i) == ' ' .and. ctp(i+1:i+1) /= ' ' )  num = num+1
!end do
!return 
!END FUNCTION MIV_NUM_COLUMNS

end module math
