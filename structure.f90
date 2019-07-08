module structure
implicit none
type Atoms
    character(2)                       :: name
    integer                            :: atomic_number
    real(8)                            :: pos
    real(8)                            :: mass
    real(8),allocatable,dimension(:,:) :: neighbor
end type Atoms

type Structure
    integer                                :: Natoms
    character(2),allocatable,dimension(:)  :: symbols
    real(8),dimension(3,3)                 :: lat
    real(8),dimension(3,3)                 :: recip_lat
    real(8),dimension(:,:),allocatable     :: pos
    real(8),dimension(:,:),allocatable     :: dpos
    real(8),dimension(:,:),allocatable     :: force
    real(8),dimension(6)                   :: stress
    real(8)                                :: volume
end type Structure
end module


