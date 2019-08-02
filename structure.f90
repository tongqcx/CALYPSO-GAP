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

SUBROUTINE GET_RMSE(AT)
type(Structure),intent(inout),dimension(:)     :: at

!local
integer                                        :: i,j,k,ii
integer                                        :: nforce, n_config
n_config = size(at)
rmse_energy = 0.d0
rmse_force = 0.d0
rmse_stress = 0.d0
nforce = 0
do i = 1, n_config
    rmse_energy = rmse_energy + (at(i)%energy_cal/at(i)%natoms - at(i)%energy_ref/at(i)%natoms)**2
    do j = 1, at(i)%natoms
        do k = 1, 3
            nforce = nforce + 1
            rmse_force = rmse_force + (at(i)%force_cal(j,k) - at(i)%force_ref(j,k))**2
        enddo
    enddo
    at(i)%stress_ref = at(i)%stress_ref * (1.0/GPa2eVPang) * 10.0 / at(i)%volume
    do j = 1,6
        rmse_stress = rmse_stress + (at(i)%stress_cal(j)/10.0 - at(i)%stress_ref(j)/10.0)**2
    enddo
enddo
print *, 'RMSE ENERGY:', sqrt(rmse_energy/n_config)
print *, 'RMSE FORCE:', sqrt(rmse_force/nforce)
print *, 'RMSE STRESS Units GPa:', sqrt(rmse_stress/n_config/6.d0)
open(181,file="predited.datf")
do ii = 1, n_config
        write(181,*) "----------------------------------------------------"
        write(181,'(A9,X,I5,X,A30)') "Structure",ii
        write(181,*) "----------------------------------------------------"
        do i  =1,at(ii)%natoms
            do j = 1,3
            if (j==1) then
            write(181,'(I5,I5,X,A,A,X,3F15.6)') ii,i, "F","X",at(ii)%force_cal(i,j),at(ii)%force_ref(i,j),&
                                    abs(at(ii)%force_cal(i,j) - at(ii)%force_ref(i,j))
            elseif(j ==2 ) then
            write(181,'(I5,I5,X,A,A,X,3F15.6)') ii,i, "F","Y",at(ii)%force_cal(i,j),at(ii)%force_ref(i,j),&
                                    abs(at(ii)%force_cal(i,j) - at(ii)%force_ref(i,j))
            else
            write(181,'(I5,I5,X,A,A,X,3F15.6)') ii,i, "F","Z",at(ii)%force_cal(i,j),at(ii)%force_ref(i,j),&
                                    abs(at(ii)%force_cal(i,j) - at(ii)%force_ref(i,j))
            endif
            enddo
        enddo
        do j = 1,6
            if (j==1) then
            write(181,'(A2,X,3F20.5)') "XX",at(ii)%stress_cal(j),at(ii)%stress_ref(j),&
                                     abs(at(ii)%stress_cal(j) - at(ii)%stress_ref(j))
            elseif(j==2) then
            write(181,'(A2,X,3F20.5)') "XY",at(ii)%stress_cal(j),at(ii)%stress_ref(j),&
                                     abs(at(ii)%stress_cal(j) - at(ii)%stress_ref(j))
            elseif(j==3) then
            write(181,'(A2,X,3F20.5)') "XZ",at(ii)%stress_cal(j),at(ii)%stress_ref(j),&
                                     abs(at(ii)%stress_cal(j) - at(ii)%stress_ref(j))
            elseif(j==4) then
            write(181,'(A2,X,3F20.5)') "YY",at(ii)%stress_cal(j),at(ii)%stress_ref(j),&
                                     abs(at(ii)%stress_cal(j) - at(ii)%stress_ref(j))
            elseif(j==5) then
            write(181,'(A2,X,3F20.5)') "YZ",at(ii)%stress_cal(j),at(ii)%stress_ref(j),&
                                     abs(at(ii)%stress_cal(j) - at(ii)%stress_ref(j))
            else
            write(181,'(A2,X,3F20.5)') "ZZ",at(ii)%stress_cal(j),at(ii)%stress_ref(j),&
                                     abs(at(ii)%stress_cal(j) - at(ii)%stress_ref(j))
            endif
        enddo
        write(181,'(A6,X,I5,X,3F15.6)') "ENERGY",ii,at(ii)%energy_cal/at(ii)%natoms, at(ii)%energy_ref/at(ii)%natoms,&
               abs(at(ii)%energy_cal-at(ii)%energy_ref)/at(ii)%natoms
        !write(181,'(A10X3F15.8)') 'E-V:', at(i)%volume/at(i)%na, (at(i)%stress_cal(1) + at(i)%stress_cal(4) + at(i)%stress_cal(6))/30.0, &
        !at(i)%energy_cal/at(i)%na
enddo
close(181)
END SUBROUTINE


end module


