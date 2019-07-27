module struct
use constants
use math
implicit none
integer,parameter                      :: max_neighbor = 2000
type Atoms
    character(2)                       :: name
    integer                            :: atomic_number
    integer                            :: nneighbor
    integer,allocatable,dimension(:)   :: count
    real(8)                            :: pos(3)
    real(8)                            :: mass
    real(8),allocatable,dimension(:,:,:) :: neighbor
endtype Atoms

type Structure
    integer                                  :: natoms
    integer                                  :: nspecies
    type(Atoms),allocatable,dimension(:)     :: atom
    character(2),allocatable,dimension(:)    :: symbols
    integer,allocatable,dimension(:)         :: index
    real(DP),allocatable,dimension(:)        :: mlp_weights
    integer,allocatable,dimension(:,:)       :: pos_index
    real(DP),dimension(3,3)                  :: lat
    real(DP),dimension(3,3)                  :: recip_lat
    real(DP),dimension(:,:),allocatable      :: pos
    real(DP),dimension(:,:),allocatable      :: dpos
    real(DP)                                 :: sigma_e
    integer,dimension(:,:),allocatable       :: interaction_mat
    ! Properties
    real(DP)                                 :: volume
    real(DP),dimension(:,:),allocatable      :: force_ref, force_cal
    real(DP),dimension(6)                    :: stress_ref, stress_cal
    real(DP)                                 :: energy_ref, energy_cal
    real(DP),dimension(:),allocatable        :: atomic_energy
    real(DP),dimension(3,3)                  :: stress         
    ! for many-body descriptors ACSF
    REAL(DP),dimension(:,:),allocatable      :: xx, kk, ckm
    REAL(DP),dimension(:,:),allocatable      :: dedg
    REAL(DP),dimension(:,:,:,:),allocatable  :: dxdy, strs

endtype Structure
type(Structure),allocatable,dimension(:)   :: at
contains
SUBROUTINE INI_STRUCTURE(at)
type(Structure),intent(inout)  :: at

allocate(at%symbols(     at%natoms))
allocate(at%mlp_weights( at%natoms))
allocate(at%index(       at%natoms))
allocate(at%pos_index(at%nspecies,2))
allocate(at%atom(        at%natoms))
allocate(at%pos(         at%natoms,3))
allocate(at%dpos(        at%natoms,3))
allocate(at%force_ref(   at%natoms,3))
allocate(at%force_cal(   at%natoms,3))
allocate(at%atomic_energy(at%natoms))
END SUBROUTINE

!------------------------------------------------------

SUBROUTINE Build_neighbor(at, data_c)
type(Structure),intent(inout)  :: at
type(data_type),intent(in)     :: data_c
integer                        :: nabc(3)
real(DP)                       :: xyz(3), dr(3), dis
integer                        :: i, j, n1, n2, n3, count
character(2)                   :: temp
integer                        :: atom_index
!/////////////////////////////////////////////////////////////////////
!print*, 'data_c%nspecies',data_c%nspecies
do i = 1, at%natoms
    do j = 1, data_c%nspecies
        if (at%symbols(i) == data_c%elements(j)) at%index(i) = j
    enddo
enddo
at%pos_index = 0
atom_index = 1
at%pos_index(atom_index,1) = 1
temp = at%symbols(1)
do i = 2, at%natoms
    if (at%symbols(i) == temp) then
        at%pos_index(atom_index,2) = i
    else
        at%pos_index(atom_index,2) = i - 1
        atom_index = atom_index + 1
        at%pos_index(atom_index,1) = i
        temp = at%symbols(i)
    endif
enddo
!print*, 'at%pos_index', at%index
!stop
!////////////////////////////////////////////////////////////////////

!rcut = 9.d0
!rmin = 0.5d0

at%recip_lat = recipvector(at%lat)
nabc(1)=ceiling(data_c%rcut*vectorlength(at%recip_lat(1,:))/pi/2)
nabc(2)=ceiling(data_c%rcut*vectorlength(at%recip_lat(2,:))/pi/2)
nabc(3)=ceiling(data_c%rcut*vectorlength(at%recip_lat(3,:))/pi/2)
do i = 1, at%natoms

!------------------------------------------------
    allocate(at%atom(i)%neighbor(at%nspecies, max_neighbor, 5))
    allocate(at%atom(i)%count(at%nspecies))
    at%atom(i)%pos = at%pos(i,:)
    at%atom(i)%name = at%symbols(i)
    at%atom(i)%count = 0
    !count = 0
    do j = 1, at%natoms
        do n1 = -nabc(1), nabc(1)
            do n2 = -nabc(2), nabc(2)
                do n3 = -nabc(3), nabc(3)
                    if ((n1.eq.0).and.(n2.eq.0).and.(n3.eq.0).and.(i.eq.j)) cycle
                    xyz = at%pos(j,:) + dble(n1)*at%lat(1,:) + dble(n2)*at%lat(2,:) + dble(n3)*at%lat(3,:)
                    dr = at%pos(i,:) - xyz
                    dis = dsqrt(dr(1)**2 + dr(2)**2 + dr(3)**2)
                    if ( dis > data_c%rcut) cycle
                    if ( dis < data_c%rmin) then
                        print*, 'The distance of two atoms is very small'
                        stop
                    endif
                    at%atom(i)%count(at%index(j)) = at%atom(i)%count(at%index(j)) + 1
                    if (at%atom(i)%count(at%index(j)) > max_neighbor) then
                        print *, 'Reset max number of neighbor'
                        stop
                    endif
                    !at%atom(i)%nneighbor = count
                    at%atom(i)%neighbor(at%index(j),at%atom(i)%count(at%index(j)),1:3) = xyz
                    at%atom(i)%neighbor(at%index(j),at%atom(i)%count(at%index(j)),4) = dis
                    at%atom(i)%neighbor(at%index(j),at%atom(i)%count(at%index(j)),5) = real(j)
                enddo
            enddo
        enddo
    enddo
enddo

end SUBROUTINE

end module


