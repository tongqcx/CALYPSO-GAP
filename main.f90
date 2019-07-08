Program  test_nGAP
use struct
implicit none
integer            :: i,j
integer            :: na, nspecies
integer            :: nconfig
open(2211,file='config')
read(2211,*) nconfig
allocate(at(nconfig))
do i= 1, nconfig
    read(2211,*)  na, nspecies
    call ini_structure(at(i), na, nspecies)
    do j = 1,3
        read(2211,*) at(i)%lat(j,:)
    enddo
    read(2211,*) at(i)%stress(:)
    do j = 1, at(i)%natoms
        read(2211,*) at(i)%symbols(j), at(i)%pos(j,:), at(i)%force(j,:)
    enddo
    read(2211,*) at(i)%energy
    call build_neighbor(at(i))
enddo
end program

