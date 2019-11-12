! GAP%dd is the only used variable in structure GAP
! DATA_C%nspecies is the only used variable in structure DATA_C

SUBROUTINE FGET_BOND(na, lat, elements, pos, rcut, min_bond)


INTEGER,intent(in)                             :: na
REAL(8),intent(in),dimension(3,3)              :: lat
INTEGER,intent(in),dimension(na)               :: elements
REAL(8),intent(in),dimension(na,3)             :: pos
REAL(8),intent(in)                             :: rcut
REAL(8),intent(out)                            :: min_bond
!local
REAL(8),dimension(3,3)                         :: recip_lat
REAL(8)                                        :: rmin
REAL(8)                                        :: dis
integer                                        :: max_neighbor
!REAL(8),allocatable,dimension(:,:,:)           :: neighbor
!INTEGER,allocatable,dimension(:)               :: neighbor_count
!REAL(8),allocatable,dimension(:)               :: weights
REAL(8),PARAMETER                              :: pi=3.141592654d0
REAL(8),dimension(3)                           :: xyz, dr
INTEGER,dimension(3)                           :: nabc
!INTEGER                                        :: nspecies, atom_number
!REAL(8)                                        :: atom_weight


!#############################################
! initial variables
rmin = 0.5d0
max_neighbor = 4000
min_bond = 10.0
!#############################################

natoms = size(pos,1)
!allocate(neighbor(natoms, max_neighbor, 6))
!allocate(neighbor_count(natoms))


recip_lat = recipvector(lat)
nabc(1)=ceiling(rcut*vectorlength(recip_lat(1,:))/pi/2)
nabc(2)=ceiling(rcut*vectorlength(recip_lat(2,:))/pi/2)
nabc(3)=ceiling(rcut*vectorlength(recip_lat(3,:))/pi/2)
!neighbor = 0.d0
!neighbor_count = 0
! this part could parallel with OpenMP
do i = 1, natoms
    do j = 1, natoms
        do n1 = -nabc(1), nabc(1)
            do n2 = -nabc(2), nabc(2)
                do n3 = -nabc(3), nabc(3)
                    if ((n1.eq.0).and.(n2.eq.0).and.(n3.eq.0).and.(i.eq.j)) cycle
                    xyz = pos(j,:) + dble(n1)*lat(1,:) + dble(n2)*lat(2,:) + dble(n3) * lat(3,:)
                    dr = pos(i,:) - xyz
                    dis = dsqrt(dr(1)**2 + dr(2)**2 + dr(3)**2)
                    if ( dis > rcut) cycle
                    if ( dis < rmin) then
                        print*, 'Warning: The distance of two atoms is very small',dis
                        !stop
                    endif

                    if (dis < min_bond)  min_bond = dis

                    !neighbor_count(i) = neighbor_count(i) + 1
                    !if (neighbor_count(i) > max_neighbor) then
                    !    print *, 'Atoms neighbor:', neighbor_count(i), 'large than max_neighbor',max_neighbor
                    !    print *, 'Please reset max_neighbor in nGAP/src/structure.f90'
                    !    stop
                    !endif
                    !neighbor(i,neighbor_count(i), 1:3) = xyz
                    !neighbor(i,neighbor_count(i), 4) = dis
                    !neighbor(i,neighbor_count(i), 5) = 0
                    !neighbor(i,neighbor_count(i), 6) = real(j)
                enddo
            enddo
        enddo
    enddo
enddo
!deallocate(neighbor)
!deallocate(neighbor_count)

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
    recipvector=recipvector/volume(lat)*pi*2.d0

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
END SUBROUTINE

