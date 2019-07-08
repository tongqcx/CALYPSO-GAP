Program  test_nGAP
use struct
implicit none
integer            :: i,j
integer            :: na, nspecies
integer            :: nconfig
character(2)       :: element(2)
element(1) = 'B'
element(2) = 'C'
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
    call build_neighbor(at(i), element)
enddo
close(2211)
open(2211,file='11.dat')
    print *, at(i)%atom(1)%pos
    write(2211,*) at(i)%atom(1)%pos
    write(2211,*) at(i)%atom(1)%nneighbor
    do i = 1 , at(i)%atom(1)%nneighbor
        write(2211, *) at(i)%atom(1)%neighbor(1,i,:)
    enddo
close(2211)
open(2211,file='22.dat')
    write(2211,*) at(i)%atom(1)%pos
    write(2211,*) at(i)%atom(1)%nneighbor
    do i = 1 , at(i)%atom(1)%nneighbor
        write(2211, *) at(i)%atom(1)%neighbor(2,i,:)
    enddo
close(2211)
end program

