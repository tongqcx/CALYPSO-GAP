module struct
use constants
use math
implicit none
type Atoms
    character(2)                       :: name
    integer                            :: atomic_number
    real(8)                            :: pos
    real(8)                            :: mass
    real(8),allocatable,dimension(:,:,:) :: neighbor
endtype Atoms

type Structure
    integer                                :: natoms
    integer                                :: nspecies
    type(Atoms),allocatable,dimension(:)   :: atom
    character(2),allocatable,dimension(:)  :: symbols
    real(8),dimension(3,3)                 :: lat
    real(8),dimension(3,3)                 :: recip_lat
    real(8),dimension(:,:),allocatable     :: pos
    real(8),dimension(:,:),allocatable     :: dpos
    real(8),dimension(:,:),allocatable     :: force
    real(8),dimension(6)                   :: stress
    real(8)                                :: volume
    real(8)                                :: energy
endtype Structure
type(Structure),allocatable,dimension(:)   :: at
contains
SUBROUTINE INI_STRUCTURE(at, na, ns)
type(Structure),intent(inout)  :: at
integer,intent(in)             :: na, ns
at%natoms = na
at%nspecies = ns
allocate(at%symbols( at%natoms))
allocate(at%atom(    at%natoms))
allocate(at%pos(     at%natoms,3))
allocate(at%dpos(    at%natoms,3))
allocate(at%force(   at%natoms,3))
END SUBROUTINE

!------------------------------------------------------

SUBROUTINE Build_neighbor(at)
type(Structure),intent(inout)  :: at
real(DP)                       :: rcut
integer                        :: nabc(3)
rcut = 9.d0
at%recip_lat = recipvector(at%lat)
nabc(1)=ceiling(rcut*vectorlength(at%recip_lat(1,:))/pi/2)
nabc(2)=ceiling(rcut*vectorlength(at%recip_lat(2,:))/pi/2)
nabc(3)=ceiling(rcut*vectorlength(at%recip_lat(3,:))/pi/2)
print*, nabc
end SUBROUTINE



end module


