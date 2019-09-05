! GAP%dd is the only used variable in structure GAP
! DATA_C%nspecies is the only used variable in structure DATA_C

SUBROUTINE FCAR2WACSF(na, nf, lat, pos, xx, rcut)

!interface
!    SUBROUTINE  CAR2ACSF(neighbor, neighbor_count, xx)
!        !REAL(8), intent(in), dimension(:,:,:)    :: neighbor
!        !INTEGER, intent(in), dimension(:)        :: neighbor_count
!        !REAL(8), intent(out), dimension(:,:)     :: xx
!        REAL(8), dimension(:,:,:)    :: neighbor
!        INTEGER, dimension(:)        :: neighbor_count
!        REAL(8),  dimension(:,:)     :: xx
!    END SUBROUTINE
!end interface

!character(2), intent(in)               :: elements
INTEGER,intent(in)                     :: na, nf
REAL(8),intent(in),dimension(3,3)      :: lat
REAL(8),intent(in),dimension(na,3)      :: pos
REAL(8),intent(out),dimension(nf,na)     :: xx
REAL(8),intent(in)                     :: rcut
!local
REAL(8),dimension(3,3)                 :: recip_lat
REAL(8)                                :: rmin
REAL(8)                                :: dis
integer                                :: max_neighbor
REAL(8),allocatable,dimension(:,:,:)   :: neighbor
INTEGER,allocatable,dimension(:)       :: neighbor_count
REAL(8),PARAMETER                      :: pi=3.141592654d0
REAL(8),dimension(3)                   :: xyz, dr
INTEGER,dimension(3)                   :: nabc

!#############################################
! initial variables
rmin = 0.5d0
max_neighbor = 4000
!#############################################

natoms = size(pos,1)
allocate(neighbor(natoms, max_neighbor, 4))

recip_lat = recipvector(lat)
nabc(1)=ceiling(rcut*vectorlength(recip_lat(1,:))/pi/2)
nabc(2)=ceiling(rcut*vectorlength(recip_lat(2,:))/pi/2)
nabc(3)=ceiling(rcut*vectorlength(recip_lat(3,:))/pi/2)
neighbor = 0.d0
neighbor_count = 0
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
                    neighbor_count(i) = neighbor_count(i) + 1
                    if (count > max_neighbor) then
                        print *, 'Atoms neighbor:', count, 'large than max_neighbor',max_neighbor
                        print *, 'Please reset max_neighbor in nGAP/src/structure.f90'
                        stop
                    endif
                    neighbor(i,neighbor_count(i), 1:3) = xyz
                    neighbor(i,neighbor_count(i), 4) = dis
                enddo
            enddo
        enddo
    enddo
enddo
call car2acsf(neighbor, neighbor_count, xx)

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


SUBROUTINE CAR2ACSF(neighbor, neighbor_count, xx)

implicit real(8) (a-h,o-z)
REAL(8), intent(in), dimension(:,:,:)    :: neighbor
INTEGER, intent(in), dimension(:)        :: neighbor_count
REAL(8), intent(out), dimension(:,:)     :: xx

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
TYPE SF
INTEGER                                 :: ntype
REAL(8)                                :: alpha
REAL(8)                                :: cutoff
END TYPE SF

TYPE ACSF_type
INTEGER                                 :: nsf
REAL(8)                                :: global_cutoff
type(SF),dimension(:),allocatable       :: sf
END TYPE ACSF_type
TYPE(ACSF_type)                         :: ACSF
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


!local
REAL(8),dimension(3)                  :: xyz, xyz_j, xyz_k

natoms = size(neighbor,1)
xx = 0.d0
!allocate(xx(GAP%dd, natoms))
!allocate(at%dxdy(GAP%dd, at%natoms, at%natoms, 3))
!allocate(at%strs(3, 3, GAP%dd, at%natoms))
! @@@@ this three array must be initial
!at%dxdy = 0.d0
!at%strs = 0.d0
! @@@@
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
open(2244,file='neural.in')
read(2244,*)  acsf%global_cutoff
read(2244,*)  acsf%nsf
allocate(acsf%sf(acsf%nsf))
do i = 1, acsf%nsf
    read(2244,*) acsf%sf(i)%ntype, acsf%sf(i)%alpha, acsf%sf(i)%cutoff
enddo
close(2244)
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
rmin = 0.5d0
nnn = ACSF%nsf
do ii = 1, nnn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!G1 = SUM_j{exp(-alpha*rij**2)*fc(rij)}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (ACSF%sf(ii)%ntype.eq.1) then
        cutoff = ACSF%sf(ii)%cutoff
        alpha = ACSF%sf(ii)%alpha
        do i = 1, natoms
            ! ******************
            ! Be careful!!!
            ! i_type loop from 1 to data_c%nspecies not at%nspecies
            ! 2019.09.04 STUPID!!!
            ! ******************
            !do i_type = 1, data_c%nspecies
                do i_neighbor = 1, neighbor_count(i)
                    rij = neighbor(i, i_neighbor, 4)
                    if (rij.gt.cutoff) cycle
                    xyz = neighbor(i, i_neighbor, 1:3)
                    fcutij = 0.5d0 * (dcos(pi*rij/cutoff) + 1.d0)
                    xx(ii,i) = xx(ii,i) + dexp(-1.d0*alpha*rij**2)*fcutij
                enddo ! i_neighbor
        enddo ! i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! lamda = 1
! eta = 1
! G2 = SUM_jk{(1+lamda*costheta_ijk)^eta*
! exp(-alpha*(rij**2+rik**2+rjk**2))*fc(rij)*fc(rik)*fc(rjk)}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    elseif (ACSF%sf(ii)%ntype.eq.2) then
        cutoff = ACSF%sf(ii)%cutoff
        alpha = ACSF%sf(ii)%alpha
!        print*, 'cutoff',cutoff,'alpha',alpha
        do i = 1, natoms
!       lllll = 0
            do j_neighbor = 1, neighbor_count(i) 
                rij = neighbor(i, j_neighbor, 4)
                if (rij.gt.cutoff) cycle
                xyz_j = neighbor(i, j_neighbor, 1:3)
                fcutij=0.5d0*(dcos(pi*rij/cutoff)+1.d0)
                do k_neighbor = 1, neighbor_count(i)
                    ! ******************
                    ! Be careful
                    ! ******************
                    if (k_neighbor <= j_neighbor) cycle
                    rik = neighbor(i, k_neighbor,4)
                    if (rik.gt.cutoff) cycle
                    xyz_k = neighbor(i, k_neighbor,1:3)
                    fcutik=0.5d0*(dcos(pi*rik/cutoff)+1.d0)
                    rjk = (xyz_j(1) - xyz_k(1))**2 + (xyz_j(2) - xyz_k(2))**2 + (xyz_j(3) - xyz_k(3))**2
                    rjk = dsqrt(rjk)

                    if (rjk.gt.cutoff) cycle  ! CAUTAINS STUPID!!!
                    if (rjk < Rmin) then
                        print*, 'Rjk', rjk,' smaller than Rmin'
                        stop
                    endif
                    fcutjk=0.5d0*(dcos(pi*rjk/cutoff)+1.d0)
                    f=rjk**2 - rij**2 -rik**2
                    g=-2.d0*rij*rik
                    costheta=f/g
                    !!!!  2^1-eta (1+lamda coseta_ijk)^eta 
                    !!!!  eta = 1 lamda = +1.d0
                    costheta=1.d0 + costheta
                    expxyz=dexp(-alpha*(rij**2+rik**2+rjk**2))
                    xx(ii,i)=xx(ii,i)+costheta*expxyz*fcutij*fcutik*fcutjk
                enddo ! k_neighbor
            enddo ! j_neighbor
!       print*, 'lllll',lllll
        enddo ! i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!G3 = SUM_j{exp(-alpha*(rij-rshift)**2)*fc(rij)}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    elseif (ACSF%sf(ii)%ntype.eq.3) then
        cutoff = ACSF%sf(ii)%cutoff
        rshift = ACSF%sf(ii)%alpha
        alpha = 4.d0
        do i = 1, natoms
            do i_neighbor = 1, neighbor_count(i)
                rij = neighbor(i, i_neighbor,4)
                if (rij.gt.cutoff) cycle
                xyz = neighbor(i, i_neighbor,1:3)
                fcutij = 0.5d0 * (dcos(pi*rij/cutoff) + 1.d0)
                xx(ii,i)=xx(ii,i)+dexp(-1.d0*alpha*(rij-rshift)**2)*fcutij
            enddo ! i_neighbor
        enddo ! i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! lamda = -1.d0
! eta = 1
! G2 = SUM_jk{(1+lamda*costheta_ijk)^eta*
! exp(-alpha*(rij**2+rik**2+rjk**2))*fc(rij)*fc(rik)*fc(rjk)}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    elseif (ACSF%sf(ii)%ntype.eq.4) then
        cutoff = ACSF%sf(ii)%cutoff
        alpha = ACSF%sf(ii)%alpha
        do i = 1, natoms
            do j_neighbor = 1, neighbor_count(i)
                rij = neighbor(i, j_neighbor,4)
                if (rij.gt.cutoff) cycle
                xyz_j = neighbor(i, j_neighbor,1:3)
                fcutij=0.5d0*(dcos(pi*rij/cutoff)+1.d0)
                do k_neighbor = 1, neighbor_count(i)
                    if (k_neighbor <= j_neighbor) cycle
                    rik = neighbor(i,k_neighbor,4)
                    if (rik.gt.cutoff) cycle
                    xyz_k = neighbor(i, k_neighbor,1:3)
                    fcutik=0.5d0*(dcos(pi*rik/cutoff)+1.d0)
                    rjk = (xyz_j(1) - xyz_k(1))**2 + (xyz_j(2) - xyz_k(2))**2 + (xyz_j(3) - xyz_k(3))**2
                    rjk = dsqrt(rjk)
                    if (rjk.gt.cutoff) cycle  ! Be careful STUPID!!!
                    if (rjk < Rmin) then
                        print*, 'Rjk', rjk,' smaller than Rmin'
                        stop
                    endif
                    fcutjk=0.5d0*(dcos(pi*rjk/cutoff)+1.d0)

                    f=rjk**2 - rij**2 -rik**2
                    g=-2.d0*rij*rik
                    costheta=f/g
                    costheta=1.d0 - costheta  ! avoid negative values

                    expxyz=dexp(-alpha*(rij**2+rik**2+rjk**2))

                    xx(ii,i)=xx(ii,i)+costheta*expxyz*fcutij*fcutik*fcutjk

                enddo ! k_neighbor
            enddo ! j_neighbor
        enddo ! i
    else
        print *, 'Unknown function type',ii, ACSF%sf(ii)%ntype
    endif
enddo  ! types
!print*, 'CCC',at%xx(:,1)
END SUBROUTINE
END SUBROUTINE

