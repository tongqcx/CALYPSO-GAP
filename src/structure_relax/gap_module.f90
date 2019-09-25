!=====================================================
! 2019.09.18
! this module is a part of GAP source code, it could be used to 
! calculation the total energy, atomic force, and cell stress for given structure
! the neural.in and gap_parameters files should be contained in WorkDir
! Usages:
! call fgap_init()
! call fgap_calc()
!=====================================================
MODULE gap_module
integer, parameter                                         :: maxele = 107
double precision,parameter                                 :: delta = 1.0
double precision,parameter                                 :: GPa2eVPang =6.24219D-3
character(len=2), save                                     :: atsym(maxele)
data atsym/'H ','He','Li','Be','B ','C ','N ','O ','F ', &
           'Ne','Na','Mg','Al','Si','P ','S ','Cl','Ar', &
           'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co', &
           'Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr', &
           'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh', &
           'Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe', &
           'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu', &
           'Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf', &
           'Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl', &
           'Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th', &
           'Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es', &
           'Fm','Md','No','Lr','Rf','Ha','D ','X '/
!@ structure information
integer,         allocatable,dimension(:)                  :: gap_elements
integer,         allocatable,dimension(:)                  :: gap_natoms
double precision,            dimension(3,3)                :: gap_lat
double precision,allocatable,dimension(:,:)                :: gap_xyz
double precision,allocatable,dimension(:,:)                :: gap_dev
double precision                                           :: gap_energy, gap_variance, gap_volume
double precision,            dimension(6)                  :: gap_stress
character(2),    allocatable,dimension(:)                  :: gap_species
type Structure
    integer                                                :: na
    integer                                                :: ntype
    double precision,allocatable,dimension(:)              :: types
    double precision,            dimension(3,3)            :: lat
    double precision                                       :: volume
    double precision,allocatable,dimension(:,:)            :: xyz
    double precision,allocatable,dimension(:,:)            :: xx    !(nsf_max,natoms_max)
    double precision,allocatable,dimension(:,:)            :: kk    !(natoms_max,nsf_max)
    double precision,allocatable,dimension(:,:)            :: ckm   !(natoms_max,nsparseX_max)  ! k in the number of atoms
    double precision,allocatable,dimension(:)              :: e     !(natoms_max)
    double precision,allocatable,dimension(:)              :: covf  !(natoms_max)
    double precision                                       :: energy_ref, energy_cal
    double precision,dimension(6)                          :: stress_ref, stress_cal
    double precision,allocatable,dimension(:,:)            :: force_ref, force_cal
    double precision,allocatable,dimension(:,:,:,:)        :: dxdy  !(nsf_max,natoms_max,natoms_max,3)
    double precision,allocatable,dimension(:,:,:,:)        :: strs  !(3,3,nsf_max,natoms_max)
    double precision,allocatable,dimension(:,:)            :: dedg  !(natoms_max,nsf_max)
    double precision,allocatable,dimension(:,:)            :: dkdg  !(natoms_max,nsf_max,nsparseX_max)
end type
type(Structure)                                            :: at_

!@ Parameters for GAP
integer                                                    :: des_len
integer                                                    :: nsparseX
double precision, allocatable, dimension(:)                :: THETA ! des_len
double precision, allocatable, dimension(:,:)              :: MM    ! nsparseX, des_len
double precision, allocatable, dimension(:,:)              :: QMM   ! nsparseX, nsparseX
double precision, allocatable, dimension(:)                :: COEFF  ! nsparseX
double precision                                           :: Rcut
logical                                                    :: lgrad

private   at_test, des_len, nsparseX, THETA, MM, QMM, COEFF, Rcut, lgrad, GPa2eVPang, delta, maxele, atsym
private   GET_COV, FCAR2WACSF, write_array_2dim, set_elements

contains
SUBROUTINE  GET_COV(n,m,xx,yy,theta, cov)!{{{
    integer i,j,k
    integer,intent(in)       ::  n,m
    double precision, intent(in)       ::  xx(:,:)
    double precision, intent(in)       ::  yy(:,:)
    double precision, intent(in)       ::  theta(:)
    double precision, intent(out)      ::  cov(n,m)
    double precision                   ::  temp
    integer                  :: nf
    nf = size(xx,2)
    temp = 0.d0
    do i = 1,n
        do  j=1,m 
            temp = 0.d0
            do k = 1,nf
                temp = temp+((xx(i,k)-yy(j,k))/theta(k))**2
            enddo
            cov(i,j) = delta*exp(-0.5d0*temp)
        enddo
    enddo
END SUBROUTINE!}}}
FUNCTION    Det(matrix)!{{{
   double precision,  intent(in) :: Matrix(3,3)
   double precision              :: Det


   Det = Matrix(1,1)*(Matrix(2,2)*Matrix(3,3)-Matrix(2,3)*Matrix(3,2))&
      -Matrix(1,2)*(Matrix(2,1)*Matrix(3,3)-Matrix(2,3)*Matrix(3,1))&
      +Matrix(1,3)*(Matrix(2,1)*Matrix(3,2)-Matrix(3,1)*Matrix(2,2))

END FUNCTION!}}}
SUBROUTINE  FGAP_CALC(NA, SPECIES, LAT, POS, ENE, FORCE, STRESS, VARIANCE)!{{{

implicit none

integer, intent(in)                                                            :: NA
integer, intent(in),dimension(NA)                                              :: SPECIES
double precision, intent(in),dimension(3,3)                                    :: LAT
double precision, intent(in),dimension(NA,3)                                   :: POS


double precision, intent(out)                                                  :: ENE
double precision, intent(out)                                                  :: VARIANCE
double precision, intent(out),dimension(NA,3)                                  :: FORCE
double precision, intent(out),dimension(6)                                     :: STRESS

!-- local --
type(Structure)                                                                :: at
integer                                                                        :: i,j,k, k1, n
double precision                                                               :: rcut
logical                                                                        :: lgrad


rcut = 6.d0
lgrad = .true.
!--print*, Na,'Na'
at%na = NA
!if(.not.allocated(at%types))       allocate(at%types(at%na))
if(.not.allocated(at%e))           allocate(at%e(at%na))
if(.not.allocated(at%force_cal))   allocate(at%force_cal(at%na, 3))
if(.not.allocated(at%xyz))         allocate(at%xyz(at%na, 3))
if(.not.allocated(at%xx))          allocate(at%xx(des_len, at%na))
if(.not.allocated(at%kk))          allocate(at%kk(at%na, des_len))
if(.not.allocated(at%covf))        allocate(at%covf(at%na))
if(.not.allocated(at%ckm))         allocate(at%ckm(at%na, nsparseX))
if(.not.allocated(at%dxdy))        allocate(at%dxdy(des_len, at%na, at%na, 3))
if(.not.allocated(at%strs))        allocate(at%strs(3, 3, des_len, at%na))
if(.not.allocated(at%dedg))        allocate(at%dedg(at%na, des_len))

at%lat = LAT
do i = 1,at%na
    at%xyz(i,:) = POS(i,:)
enddo
call FCAR2WACSF(&
                 at%na, des_len,&
                 at%lat, SPECIES, &
                 at%xyz(1:at%na,:),&
                 at%xx(1:des_len, 1:at%na),&
                 at%dxdy(1:des_len,1:at%na,1:at%na,:),&
                 at%strs(:,:,1:des_len,1:at%na), &
                 rcut, lgrad)
do i = 1,at%na
    at%kk(i,1:des_len) = at%xx(1:des_len,i)
enddo
!call write_array_2dim(at%na, des_len, at%kk, 'kk.dat')

!--print*, 'get_cov'
at%ckm = 0.d0
call get_cov(at%na,&
             nsparseX,&
             at%kk(1:at%na,1:des_len),&
             mm(1:nsparseX,1:des_len),theta(1:des_len), &
             at%ckm(1:at%na,1:nsparseX))

!call write_array_2dim(at%na, nsparseX, at%ckm, 'ckm.dat')
!--print*, 'get energy'
at%e = matmul(at%ckm(1:at%na,1:nsparseX),coeff(1:nsparseX))
at%energy_cal = sum(at%e)

!call write_array_2dim(at%na, nsparseX, at%ckm(1:at%na,1:nsparseX), 'ckm.dat')
at%dedg = 0.d0
!--print*, 'get de/dd'
do i = 1 , at%na
    do j =  1,des_len
        do k = 1,nsparseX
            at%dedg(i,j) = at%dedg(i,j) - 1.d0*(at%kk(i,j)-mm(k,j))/theta(j)**2*at%ckm(i,k)*coeff(k)         
        enddo
    enddo
enddo
!call write_array_2dim(at%na, des_len, at%dedg(1:at%na, 1:des_len), 'dedg.dat')

!--print*, 'get force'
at%force_cal = 0.d0
do i = 1, at%na
    do n = 1, at%na
        do j = 1,3
            do k = 1, des_len
                at%force_cal(i,j) = at%force_cal(i,j) - at%dedg(n,k) * at%dxdy(k,n,i,j)
            enddo
        enddo
    enddo
enddo

!--print*, 'get stress'
at%stress_cal = 0.d0
k1 = 1
do i = 1, 3
    do j = i, 3
        do n = 1, at%na
            do k = 1, des_len
                at%stress_cal(k1) = at%stress_cal(k1) - at%dedg(n,k) * at%strs(i,j,k,n)
            enddo
        enddo
        k1 = k1 + 1
    enddo
enddo
at%volume = ABS(det(at%lat))
gap_volume = at%volume
at%stress_cal = at%stress_cal * (1.0/GPa2eVPang) / at%volume
!print*, at%stress_cal

!--print*, 'get variance'
VARIANCE = 0.d0
do i  = 1,at%na
    at%covf(i) = delta - dot_product(at%ckm(i,1:nsparseX),matmul(qmm,at%ckm(i,1:nsparseX)))
    VARIANCE = VARIANCE + at%covf(i)/at%na
enddo

ENE = at%energy_cal 
FORCE = at%force_cal(1:at%na,:)
STRESS = at%stress_cal
!print *, ENE
!print*, transpose(FORCE)
!print*, STRESS
!stop

END SUBROUTINE!}}}
SUBROUTINE  FGAP_INIT(nat)!{{{
implicit none
Integer, optional, intent(in)                 :: nat(:)

!-- local--
logical                                       :: alive
integer                                       :: i

inquire(file="gap_parameters",exist=alive)
if(.not.alive) then
    print*, "gap_parameters file does not exist!"
    stop
endif
open(111,file="gap_parameters") 
read(111,*) nsparseX,des_len
if (.not. allocated(MM)) allocate(MM(nsparseX, des_len))
if (.not. allocated(theta)) allocate(theta(des_len))
if (.not. allocated(QMM)) allocate(QMM(nsparseX, nsparseX))
if (.not. allocated(COEFF)) allocate(COEFF(nsparseX))
read(111,*) 
read(111,*) 
read(111,*) 
read(111,*) theta(1:des_len)
do i = 1,nsparseX
    read(111,*) MM(i,1:des_len)
enddo
!do i = 1,nsparseX
!    read(111,*) QMM(i,1:nsparseX)
!enddo
read(111,*) COEFF(1:nsparseX)
close(111)
if (present(nat)) call SET_ELEMENTS(nat)

END SUBROUTINE!}}}
SUBROUTINE  FCAR2WACSF(na, nf, lat, elements, pos, xx, dxdy, strs, rcut, lgrad)!{{{

INTEGER,intent(in)                             :: na, nf
REAL(8),intent(in),dimension(3,3)              :: lat
INTEGER,intent(in),dimension(na)               :: elements
REAL(8),intent(in),dimension(na,3)             :: pos
REAL(8),intent(out),dimension(nf,na)           :: xx
REAL(8),intent(out),dimension(nf, na, na, 3)   :: dxdy
REAL(8),intent(out),dimension(3, 3, nf, na)    :: strs
REAL(8),intent(in)                             :: rcut
LOGICAL,intent(in)                             :: lgrad
!local
REAL(8),dimension(3,3)                         :: recip_lat
REAL(8)                                        :: rmin, min_bond
REAL(8)                                        :: dis
integer                                        :: max_neighbor
REAL(8),allocatable,dimension(:,:,:)           :: neighbor
INTEGER,allocatable,dimension(:)               :: neighbor_count
REAL(8),allocatable,dimension(:)               :: weights
REAL(8),PARAMETER                              :: pi=3.141592654d0
REAL(8),dimension(3)                           :: xyz, dr
INTEGER,dimension(3)                           :: nabc
INTEGER                                        :: nspecies, atom_number
REAL(8)                                        :: atom_weight


!#############################################
! initial variables
rmin = 0.5d0
max_neighbor = 4000
!#############################################

natoms = size(pos,1)
allocate(neighbor(natoms, max_neighbor, 6))
allocate(neighbor_count(natoms))
allocate(weights(natoms))

open(2233, file='neural.in')
read(2233,*) nspecies
do i = 1, nspecies
    read(2233,*) atom_number, atom_weight
    do j = 1, natoms
        if (elements(j) == atom_number) weights(j) = atom_weight
    enddo
enddo
close(2233)


recip_lat = recipvector(lat)
nabc(1)=ceiling(rcut*vectorlength(recip_lat(1,:))/pi/2)
nabc(2)=ceiling(rcut*vectorlength(recip_lat(2,:))/pi/2)
nabc(3)=ceiling(rcut*vectorlength(recip_lat(3,:))/pi/2)
neighbor = 0.d0
neighbor_count = 0
min_bond = 10.d0
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
                    if ( dis < min_bond) then
                    !    print*, 'Warning: The distance of two atoms is very small',dis
                        !stop
                        min_bond = dis
                    endif
                    neighbor_count(i) = neighbor_count(i) + 1
                    if (neighbor_count(i) > max_neighbor) then
                        print *, 'Atoms neighbor:', neighbor_count(i), 'large than max_neighbor',max_neighbor
                        print *, 'Please reset max_neighbor in nGAP/src/structure.f90'
                        stop
                    endif
                    neighbor(i,neighbor_count(i), 1:3) = xyz
                    neighbor(i,neighbor_count(i), 4) = dis
                    neighbor(i,neighbor_count(i), 5) = weights(j)
                    neighbor(i,neighbor_count(i), 6) = real(j)
                enddo
            enddo
        enddo
    enddo
enddo
call car2acsf(natoms, max_neighbor, nf, pos, neighbor, neighbor_count, xx, dxdy, strs, lgrad)
deallocate(weights)
deallocate(neighbor)
deallocate(neighbor_count)

if (min_bond < rmin) write(*,'(A35,X,F10.3)'), 'CAUTIONS: the shortest bond length:', min_bond

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

SUBROUTINE CAR2ACSF(NA, max_neighbor, nf, pos, neighbor, neighbor_count, xx, dxdy, strs, lgrad)

implicit real(8) (a-h,o-z)
INTEGER, intent(in)                                      :: NA, max_neighbor, NF
REAL(8), intent(in), dimension(NA,3)                     :: pos
REAL(8), intent(in), dimension(NA,max_neighbor,6)        :: neighbor
INTEGER, intent(in), dimension(NA)                       :: neighbor_count
REAL(8), intent(out), dimension(NF, NA)                  :: xx
REAL(8), intent(out), dimension(NF,NA,NA,3)              :: dxdy
REAL(8), intent(out), dimension(3,3,NF,NA)               :: strs
LOGICAL, intent(in)                                      :: lgrad

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
TYPE SF
INTEGER                                      :: ntype
REAL(8)                                      :: alpha
REAL(8)                                      :: cutoff
END TYPE SF

TYPE ACSF_type
INTEGER                                      :: nsf
REAL(8)                                      :: global_cutoff
type(SF),dimension(:),allocatable            :: sf
END TYPE ACSF_type
TYPE(ACSF_type)                              :: ACSF
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


!local
REAL(8),PARAMETER                            :: pi=3.141592654d0
REAL(8),dimension(3)                         :: xyz, xyz_j, xyz_k
logical                                      :: alive
INTEGER                                      :: nspecies
REAL(8)                                      :: weights, weights_j, weights_k
REAL(8)                                      :: rij, fcutij, rik, fcutik, rjk, fcutjk

natoms = size(neighbor,1)
!print*, 'NA, NF' , NA, NF
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
inquire(file="neural.in",exist=alive)
if(.not.alive) then
    print*, "neural.in file does not exist!"
    stop
endif
open(2244,file='neural.in')
read(2244,*)  nspecies
do i = 1, nspecies
    read(2244, *)
enddo
!read(2244,*)  acsf%global_cutoff
read(2244,*)  acsf%nsf
allocate(acsf%sf(acsf%nsf))
do i = 1, acsf%nsf
    read(2244,*) acsf%sf(i)%ntype, acsf%sf(i)%alpha, acsf%sf(i)%cutoff
enddo
close(2244)
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
rmin = 0.5d0
nnn = ACSF%nsf

xx = 0.d0
dxdy = 0.d0
strs = 0.d0

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
                    weights = neighbor(i, i_neighbor, 5)
                    n = int(neighbor(i, i_neighbor, 6))
                    fcutij = 0.5d0 * (dcos(pi*rij/cutoff) + 1.d0)
                    xx(ii,i) = xx(ii,i) + dexp(-1.d0*alpha*rij**2)*fcutij
                    xx(ii + nnn, i) = xx(ii + nnn, i) + dexp(-1.d0*alpha*rij**2)*fcutij * weights !!!!!!! 

                    if (lgrad) then
                        deltaxj = -1.d0*(pos(i, 1) - xyz(1))
                        deltayj = -1.d0*(pos(i, 2) - xyz(2))
                        deltazj = -1.d0*(pos(i, 3) - xyz(3))
                        drijdxi = -1.d0*deltaxj/rij
                        drijdyi = -1.d0*deltayj/rij
                        drijdzi = -1.d0*deltazj/rij
                        drijdxj = -1.d0*drijdxi
                        drijdyj = -1.d0*drijdyi
                        drijdzj = -1.d0*drijdzi
                        temp1=0.5d0*(-dsin(pi*rij/cutoff))*(pi/cutoff)
                        dfcutijdxi=temp1*drijdxi
                        dfcutijdyi=temp1*drijdyi
                        dfcutijdzi=temp1*drijdzi
                        dfcutijdxj=-1.d0*dfcutijdxi
                        dfcutijdyj=-1.d0*dfcutijdyi
                        dfcutijdzj=-1.d0*dfcutijdzi
                        !dxx/dx
                        temp1=-2.d0*alpha*rij*dexp(-1.d0*alpha*rij**2)*fcutij
                        temp2= dexp(-1.d0*alpha*rij**2)

                        dxdy(ii,i,i,1)=dxdy(ii,i,i,1)+(drijdxi*temp1 + temp2*dfcutijdxi)
                        dxdy(ii+nnn,i,i,1)=dxdy(ii+nnn,i,i,1) + (drijdxi*temp1+ temp2*dfcutijdxi) * weights
                
                        temp3=drijdxj*temp1 + temp2*dfcutijdxj
                        dxdy(ii,i,n,1)=dxdy(ii,i,n,1)+temp3
                
                        temp4=temp3*weights
                        dxdy(ii + nnn,i,n,1)=dxdy(ii + nnn,i,n,1)+temp4
                
                        strs(1,1,ii,i)=strs(1,1,ii,i)+deltaxj*temp3
                        strs(2,1,ii,i)=strs(2,1,ii,i)+deltayj*temp3
                        strs(3,1,ii,i)=strs(3,1,ii,i)+deltazj*temp3
                
                        strs(1,1,ii+nnn,i)=strs(1,1,ii+nnn,i)+deltaxj*temp4
                        strs(2,1,ii+nnn,i)=strs(2,1,ii+nnn,i)+deltayj*temp4
                        strs(3,1,ii+nnn,i)=strs(3,1,ii+nnn,i)+deltazj*temp4
                        !dxx/dy
                        dxdy(ii,i,i,2)=dxdy(ii,i,i,2)+(drijdyi*temp1+temp2*dfcutijdyi)
                        dxdy(ii+nnn,i,i,2)=dxdy(ii+nnn,i,i,2)+(drijdyi*temp1+temp2*dfcutijdyi)*weights
                        temp3= drijdyj*temp1 + temp2*dfcutijdyj
                        dxdy(ii,i,n,2)=dxdy(ii,i,n,2)+temp3
                        temp4=temp3 * weights
                        dxdy(ii + nnn,i,n,2)=dxdy(ii + nnn ,i,n,2)+temp4
                
                        strs(1,2,ii,i)=strs(1,2,ii,i)+deltaxj*temp3
                        strs(2,2,ii,i)=strs(2,2,ii,i)+deltayj*temp3
                        strs(3,2,ii,i)=strs(3,2,ii,i)+deltazj*temp3
                
                        strs(1,2,ii + nnn,i)=strs(1,2,ii + nnn,i)+deltaxj*temp4
                        strs(2,2,ii + nnn,i)=strs(2,2,ii + nnn,i)+deltayj*temp4
                        strs(3,2,ii + nnn,i)=strs(3,2,ii + nnn,i)+deltazj*temp4
                        !dxx/dz
                        dxdy(ii,i,i,3)=dxdy(ii,i,i,3)+&
                               (drijdzi*temp1&
                              + temp2*dfcutijdzi)
                        dxdy(ii + nnn,i,i,3)=dxdy(ii + nnn,i,i,3)+&
                               (drijdzi*temp1&
                              + temp2*dfcutijdzi)*weights
                        temp3=drijdzj*temp1 + temp2*dfcutijdzj
                        dxdy(ii,i,n,3)=dxdy(ii,i,n,3)+temp3
                        temp4=temp3*weights
                        dxdy(ii + nnn,i,n,3)=dxdy(ii + nnn,i,n,3)+temp4
                
                        strs(1,3,ii,i)=strs(1,3,ii,i)+deltaxj*temp3
                        strs(2,3,ii,i)=strs(2,3,ii,i)+deltayj*temp3
                        strs(3,3,ii,i)=strs(3,3,ii,i)+deltazj*temp3
                
                        strs(1,3,ii + nnn,i)=strs(1,3,ii + nnn,i)+deltaxj*temp4
                        strs(2,3,ii + nnn,i)=strs(2,3,ii + nnn,i)+deltayj*temp4
                        strs(3,3,ii + nnn,i)=strs(3,3,ii + nnn,i)+deltazj*temp4
                    endif ! lgrad
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
                weights_j = neighbor(i, j_neighbor, 5)
                n = int(neighbor(i, j_neighbor, 6))
                if (lgrad) then
                    deltaxj = -1.d0*(pos(i, 1) - xyz_j(1))
                    deltayj = -1.d0*(pos(i, 2) - xyz_j(2))
                    deltazj = -1.d0*(pos(i, 3) - xyz_j(3))
                    drijdxi = -1.d0*deltaxj/rij
                    drijdyi = -1.d0*deltayj/rij
                    drijdzi = -1.d0*deltazj/rij
                    drijdxj = -1.d0*drijdxi
                    drijdyj = -1.d0*drijdyi
                    drijdzj = -1.d0*drijdzi
                    drijdxk = 0.d0
                    drijdyk = 0.d0
                    drijdzk = 0.d0
                    temp1=0.5d0*(-dsin(pi*rij/cutoff))*(pi/cutoff)
                    dfcutijdxi=temp1*drijdxi
                    dfcutijdyi=temp1*drijdyi
                    dfcutijdzi=temp1*drijdzi
                    dfcutijdxj=-1.d0*dfcutijdxi
                    dfcutijdyj=-1.d0*dfcutijdyi
                    dfcutijdzj=-1.d0*dfcutijdzi
                    dfcutijdxk=0.0d0
                    dfcutijdyk=0.0d0
                    dfcutijdzk=0.0d0
                endif
                do k_neighbor = 1, neighbor_count(i)
                    ! ******************
                    ! Be careful
                    ! ******************
                    if (k_neighbor <= j_neighbor) cycle
                    rik = neighbor(i, k_neighbor,4)
                    if (rik.gt.cutoff) cycle
                    xyz_k = neighbor(i, k_neighbor,1:3)
                    weights_k = neighbor(i, k_neighbor,5)
                    m = int(neighbor(i, k_neighbor,6))
                    fcutik=0.5d0*(dcos(pi*rik/cutoff)+1.d0)
                    if (lgrad) then
                        deltaxk = -1.d0*(pos(i, 1) - xyz_k(1))
                        deltayk = -1.d0*(pos(i, 2) - xyz_k(2))
                        deltazk = -1.d0*(pos(i, 3) - xyz_k(3))
                        drikdxi = -deltaxk/rik
                        drikdyi = -deltayk/rik
                        drikdzi = -deltazk/rik
                        drikdxk = -1.d0*drikdxi
                        drikdyk = -1.d0*drikdyi
                        drikdzk = -1.d0*drikdzi
                        drikdxj = 0.d0
                        drikdyj = 0.d0
                        drikdzj = 0.d0
                        temp1=0.5d0*(-dsin(pi*rik/cutoff))*(pi/cutoff)
                        dfcutikdxi=temp1*drikdxi
                        dfcutikdyi=temp1*drikdyi
                        dfcutikdzi=temp1*drikdzi
                        dfcutikdxj=0.0d0
                        dfcutikdyj=0.0d0
                        dfcutikdzj=0.0d0
                        dfcutikdxk=-1.d0*dfcutikdxi
                        dfcutikdyk=-1.d0*dfcutikdyi
                        dfcutikdzk=-1.d0*dfcutikdzi
                    endif
                    rjk = (xyz_j(1) - xyz_k(1))**2 + (xyz_j(2) - xyz_k(2))**2 + (xyz_j(3) - xyz_k(3))**2
                    rjk = dsqrt(rjk)

                    if (rjk.gt.cutoff) cycle  ! CAUTAINS STUPID!!!
                    if (rjk < Rmin) then
                    !    print*, 'Rjk', rjk,' smaller than Rmin'
                    !    stop
                    endif
                    fcutjk=0.5d0*(dcos(pi*rjk/cutoff)+1.d0)
                    if (lgrad) then
                        drjkdxj = (xyz_j(1) - xyz_k(1))/rjk
                        drjkdyj = (xyz_j(2) - xyz_k(2))/rjk
                        drjkdzj = (xyz_j(3) - xyz_k(3))/rjk
                        drjkdxk = -1.d0*drjkdxj
                        drjkdyk = -1.d0*drjkdyj
                        drjkdzk = -1.d0*drjkdzj
                        drjkdxi = 0.d0
                        drjkdyi = 0.d0
                        drjkdzi = 0.d0
                        temp1=0.5d0*(-dsin(pi*rjk/cutoff))*(pi/cutoff)
                        dfcutjkdxj=temp1*drjkdxj
                        dfcutjkdyj=temp1*drjkdyj
                        dfcutjkdzj=temp1*drjkdzj
                        dfcutjkdxk=-1.d0*dfcutjkdxj
                        dfcutjkdyk=-1.d0*dfcutjkdyj
                        dfcutjkdzk=-1.d0*dfcutjkdzj
                        dfcutjkdxi=0.0d0
                        dfcutjkdyi=0.0d0
                        dfcutjkdzi=0.0d0
                    endif
                    f=rjk**2 - rij**2 -rik**2
                    g=-2.d0*rij*rik
                    costheta=f/g
                    !!!!  2^1-eta (1+lamda coseta_ijk)^eta 
                    !!!!  eta = 1 lamda = +1.d0
                    costheta=1.d0 + costheta
                    if (lgrad) then
                        dfdxi=-2.d0*rij*drijdxi - 2.d0*rik*drikdxi
                        dfdyi=-2.d0*rij*drijdyi - 2.d0*rik*drikdyi
                        dfdzi=-2.d0*rij*drijdzi - 2.d0*rik*drikdzi

                        dfdxj=2.d0*rjk*drjkdxj - 2.d0*rij*drijdxj
                        dfdyj=2.d0*rjk*drjkdyj - 2.d0*rij*drijdyj
                        dfdzj=2.d0*rjk*drjkdzj - 2.d0*rij*drijdzj

                        dfdxk=2.d0*rjk*drjkdxk - 2.d0*rik*drikdxk
                        dfdyk=2.d0*rjk*drjkdyk - 2.d0*rik*drikdyk
                        dfdzk=2.d0*rjk*drjkdzk - 2.d0*rik*drikdzk

                        dgdxi=-2.d0*(drijdxi*rik + rij*drikdxi)
                        dgdyi=-2.d0*(drijdyi*rik + rij*drikdyi)
                        dgdzi=-2.d0*(drijdzi*rik + rij*drikdzi)

                        dgdxj=-2.d0*drijdxj*rik
                        dgdyj=-2.d0*drijdyj*rik
                        dgdzj=-2.d0*drijdzj*rik

                        dgdxk=-2.d0*rij*drikdxk
                        dgdyk=-2.d0*rij*drikdyk
                        dgdzk=-2.d0*rij*drikdzk

                        temp1=1.d0/g**2
                        dcosthetadxi=(dfdxi*g - f*dgdxi)*temp1
                        dcosthetadyi=(dfdyi*g - f*dgdyi)*temp1
                        dcosthetadzi=(dfdzi*g - f*dgdzi)*temp1
                        dcosthetadxj=(dfdxj*g - f*dgdxj)*temp1
                        dcosthetadyj=(dfdyj*g - f*dgdyj)*temp1
                        dcosthetadzj=(dfdzj*g - f*dgdzj)*temp1
                        dcosthetadxk=(dfdxk*g - f*dgdxk)*temp1
                        dcosthetadyk=(dfdyk*g - f*dgdyk)*temp1
                        dcosthetadzk=(dfdzk*g - f*dgdzk)*temp1
                    endif
                    expxyz=dexp(-alpha*(rij**2+rik**2+rjk**2))
                    if (lgrad) then
                        temp1=-alpha*2.0d0*expxyz
                        dexpxyzdxi=(rij*drijdxi+rik*drikdxi+rjk*drjkdxi)*temp1
                        dexpxyzdyi=(rij*drijdyi+rik*drikdyi+rjk*drjkdyi)*temp1
                        dexpxyzdzi=(rij*drijdzi+rik*drikdzi+rjk*drjkdzi)*temp1
                        dexpxyzdxj=(rij*drijdxj+rik*drikdxj+rjk*drjkdxj)*temp1
                        dexpxyzdyj=(rij*drijdyj+rik*drikdyj+rjk*drjkdyj)*temp1
                        dexpxyzdzj=(rij*drijdzj+rik*drikdzj+rjk*drjkdzj)*temp1
                        dexpxyzdxk=(rij*drijdxk+rik*drikdxk+rjk*drjkdxk)*temp1
                        dexpxyzdyk=(rij*drijdyk+rik*drikdyk+rjk*drjkdyk)*temp1
                        dexpxyzdzk=(rij*drijdzk+rik*drikdzk+rjk*drjkdzk)*temp1
                    endif
                    xx(ii,i)=xx(ii,i)+costheta*expxyz*fcutij*fcutik*fcutjk
                    xx(ii + nnn,i)=xx(ii + nnn,i)+&
                    costheta*expxyz*fcutij*fcutik*fcutjk*weights_j*weights_k

                    if (lgrad) then
                        temp1=(dcosthetadxi*expxyz*fcutij*fcutik*fcutjk&
                              +costheta*dexpxyzdxi*fcutij*fcutik*fcutjk&
                              +costheta*expxyz*dfcutijdxi*fcutik*fcutjk&
                              +costheta*expxyz*fcutij*dfcutikdxi*fcutjk&
                              +costheta*expxyz*fcutij*fcutik*dfcutjkdxi)
                        temp2=(dcosthetadxj*expxyz*fcutij*fcutik*fcutjk&
                              +costheta*dexpxyzdxj*fcutij*fcutik*fcutjk&
                              +costheta*expxyz*dfcutijdxj*fcutik*fcutjk&
                              +costheta*expxyz*fcutij*dfcutikdxj*fcutjk&
                              +costheta*expxyz*fcutij*fcutik*dfcutjkdxj)
                        temp3=(dcosthetadxk*expxyz*fcutij*fcutik*fcutjk&
                              +costheta*dexpxyzdxk*fcutij*fcutik*fcutjk&
                              +costheta*expxyz*dfcutijdxk*fcutik*fcutjk&
                              +costheta*expxyz*fcutij*dfcutikdxk*fcutjk&
                              +costheta*expxyz*fcutij*fcutik*dfcutjkdxk)
                        temp4 = temp1 * weights_j * weights_k
                        temp5 = temp2 * weights_j * weights_k
                        temp6 = temp3 * weights_j * weights_k
                        dxdy(ii,i,i,1)=dxdy(ii,i,i,1)+temp1
                        dxdy(ii,i,n,1)=dxdy(ii,i,n,1)+temp2
                        dxdy(ii,i,m,1)=dxdy(ii,i,m,1)+temp3
                        dxdy(ii + nnn,i,i,1)=dxdy(ii + nnn,i,i,1)+temp4
                        dxdy(ii + nnn,i,n,1)=dxdy(ii + nnn,i,n,1)+temp5
                        dxdy(ii + nnn,i,m,1)=dxdy(ii + nnn,i,m,1)+temp6

                        strs(1,1,ii,i)=strs(1,1,ii,i)+deltaxj*temp2+deltaxk*temp3
                        strs(2,1,ii,i)=strs(2,1,ii,i)+deltayj*temp2+deltayk*temp3
                        strs(3,1,ii,i)=strs(3,1,ii,i)+deltazj*temp2+deltazk*temp3
                        strs(1,1,ii + nnn,i)=strs(1,1,ii + nnn,i)+deltaxj*temp5+deltaxk*temp6
                        strs(2,1,ii + nnn,i)=strs(2,1,ii + nnn,i)+deltayj*temp5+deltayk*temp6
                        strs(3,1,ii + nnn,i)=strs(3,1,ii + nnn,i)+deltazj*temp5+deltazk*temp6
                        ! dxxii/dy_i
                        temp1=(dcosthetadyi*expxyz*fcutij*fcutik*fcutjk&
                              +costheta*dexpxyzdyi*fcutij*fcutik*fcutjk&
                              +costheta*expxyz*dfcutijdyi*fcutik*fcutjk&
                              +costheta*expxyz*fcutij*dfcutikdyi*fcutjk&
                              +costheta*expxyz*fcutij*fcutik*dfcutjkdyi)
                        temp2=(dcosthetadyj*expxyz*fcutij*fcutik*fcutjk&
                              +costheta*dexpxyzdyj*fcutij*fcutik*fcutjk&
                              +costheta*expxyz*dfcutijdyj*fcutik*fcutjk&
                              +costheta*expxyz*fcutij*dfcutikdyj*fcutjk&
                              +costheta*expxyz*fcutij*fcutik*dfcutjkdyj)
                        temp3=(dcosthetadyk*expxyz*fcutij*fcutik*fcutjk&
                              +costheta*dexpxyzdyk*fcutij*fcutik*fcutjk&
                              +costheta*expxyz*dfcutijdyk*fcutik*fcutjk&
                              +costheta*expxyz*fcutij*dfcutikdyk*fcutjk&
                              +costheta*expxyz*fcutij*fcutik*dfcutjkdyk)
                        temp4 = temp1 * weights_j * weights_k
                        temp5 = temp2 * weights_j * weights_k
                        temp6 = temp3 * weights_j * weights_k
                        dxdy(ii,i,i,2)=dxdy(ii,i,i,2)+temp1
                        dxdy(ii,i,n,2)=dxdy(ii,i,n,2)+temp2
                        dxdy(ii,i,m,2)=dxdy(ii,i,m,2)+temp3

                        dxdy(ii + nnn,i,i,2)=dxdy(ii + nnn,i,i,2)+temp4
                        dxdy(ii + nnn,i,n,2)=dxdy(ii + nnn,i,n,2)+temp5
                        dxdy(ii + nnn,i,m,2)=dxdy(ii + nnn,i,m,2)+temp6
                        strs(1,2,ii,i)=strs(1,2,ii,i)+deltaxj*temp2+deltaxk*temp3
                        strs(2,2,ii,i)=strs(2,2,ii,i)+deltayj*temp2+deltayk*temp3
                        strs(3,2,ii,i)=strs(3,2,ii,i)+deltazj*temp2+deltazk*temp3

                        strs(1,2,ii + nnn,i)=strs(1,2,ii + nnn,i)+deltaxj*temp5+deltaxk*temp6
                        strs(2,2,ii + nnn,i)=strs(2,2,ii + nnn,i)+deltayj*temp5+deltayk*temp6
                        strs(3,2,ii + nnn,i)=strs(3,2,ii + nnn,i)+deltazj*temp5+deltazk*temp6
                        ! dxxii/dz_i
                        temp1=(dcosthetadzi*expxyz*fcutij*fcutik*fcutjk&
                              +costheta*dexpxyzdzi*fcutij*fcutik*fcutjk&
                              +costheta*expxyz*dfcutijdzi*fcutik*fcutjk&
                              +costheta*expxyz*fcutij*dfcutikdzi*fcutjk&
                              +costheta*expxyz*fcutij*fcutik*dfcutjkdzi)
                        temp2=(dcosthetadzj*expxyz*fcutij*fcutik*fcutjk&
                              +costheta*dexpxyzdzj*fcutij*fcutik*fcutjk&
                              +costheta*expxyz*dfcutijdzj*fcutik*fcutjk&
                              +costheta*expxyz*fcutij*dfcutikdzj*fcutjk&
                              +costheta*expxyz*fcutij*fcutik*dfcutjkdzj)
                        temp3=(dcosthetadzk*expxyz*fcutij*fcutik*fcutjk&
                              +costheta*dexpxyzdzk*fcutij*fcutik*fcutjk&
                              +costheta*expxyz*dfcutijdzk*fcutik*fcutjk&
                              +costheta*expxyz*fcutij*dfcutikdzk*fcutjk&
                              +costheta*expxyz*fcutij*fcutik*dfcutjkdzk)
                        temp4 = temp1 * weights_j * weights_k
                        temp5 = temp2 * weights_j * weights_k
                        temp6 = temp3 * weights_j * weights_k
                        dxdy(ii,i,i,3)=dxdy(ii,i,i,3)+temp1
                        dxdy(ii,i,n,3)=dxdy(ii,i,n,3)+temp2
                        dxdy(ii,i,m,3)=dxdy(ii,i,m,3)+temp3

                        dxdy(ii + nnn,i,i,3)=dxdy(ii + nnn,i,i,3)+temp4
                        dxdy(ii + nnn,i,n,3)=dxdy(ii + nnn,i,n,3)+temp5
                        dxdy(ii + nnn,i,m,3)=dxdy(ii + nnn,i,m,3)+temp6
                        strs(1,3,ii,i)=strs(1,3,ii,i)+deltaxj*temp2+deltaxk*temp3
                        strs(2,3,ii,i)=strs(2,3,ii,i)+deltayj*temp2+deltayk*temp3
                        strs(3,3,ii,i)=strs(3,3,ii,i)+deltazj*temp2+deltazk*temp3

                        strs(1,3,ii + nnn,i)=strs(1,3,ii + nnn,i)+deltaxj*temp5+deltaxk*temp6
                        strs(2,3,ii + nnn,i)=strs(2,3,ii + nnn,i)+deltayj*temp5+deltayk*temp6
                        strs(3,3,ii + nnn,i)=strs(3,3,ii + nnn,i)+deltazj*temp5+deltazk*temp6
                    endif
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
                weights = neighbor(i, i_neighbor, 5)
                n = int(neighbor(i, i_neighbor, 6))
                fcutij = 0.5d0 * (dcos(pi*rij/cutoff) + 1.d0)
                if (lgrad) then
                    deltaxj = -1.d0*(pos(i, 1) - xyz(1))
                    deltayj = -1.d0*(pos(i, 2) - xyz(2))
                    deltazj = -1.d0*(pos(i, 3) - xyz(3))
                    drijdxi = -1.d0*deltaxj/rij
                    drijdyi = -1.d0*deltayj/rij
                    drijdzi = -1.d0*deltazj/rij
                    drijdxj = -1.d0*drijdxi
                    drijdyj = -1.d0*drijdyi
                    drijdzj = -1.d0*drijdzi
                    temp1=0.5d0*(-dsin(pi*rij/cutoff))*(pi/cutoff)
                    dfcutijdxi=temp1*drijdxi
                    dfcutijdyi=temp1*drijdyi
                    dfcutijdzi=temp1*drijdzi
                    dfcutijdxj=-1.d0*dfcutijdxi
                    dfcutijdyj=-1.d0*dfcutijdyi
                    dfcutijdzj=-1.d0*dfcutijdzi
                endif
                xx(ii,i)=xx(ii,i)+dexp(-1.d0*alpha*(rij-rshift)**2)*fcutij
                xx(ii + nnn,i)=xx(ii + nnn,i)+dexp(-1.d0*alpha*(rij-rshift)**2)*fcutij*weights
                if (lgrad) then
                    temp1=-2.d0*alpha*(rij-rshift)
                    temp2=dexp(-1.d0*alpha*(rij-rshift)**2)
                    ! dxx/dx
                    dxdy(ii,i,i,1)=dxdy(ii,i,i,1)+&
                           (temp1*drijdxi*temp2*fcutij&
                          + temp2*dfcutijdxi)

                    dxdy(ii + nnn,i,i,1)=dxdy(ii + nnn,i,i,1)+&
                           (temp1*drijdxi*temp2*fcutij&
                          + temp2*dfcutijdxi)*weights
                    temp3=temp1*drijdxj*temp2*fcutij + temp2*dfcutijdxj
                    dxdy(ii,i,n,1)=dxdy(ii,i,n,1)+temp3
                    temp4 = temp3 * weights
                    dxdy(ii + nnn,i,n,1)=dxdy(ii + nnn,i,n,1)+temp4
                    strs(1,1,ii,i)=strs(1,1,ii,i)+deltaxj*temp3
                    strs(2,1,ii,i)=strs(2,1,ii,i)+deltayj*temp3
                    strs(3,1,ii,i)=strs(3,1,ii,i)+deltazj*temp3
                    strs(1,1,ii + nnn,i)=strs(1,1,ii + nnn,i)+deltaxj*temp4
                    strs(2,1,ii + nnn,i)=strs(2,1,ii + nnn,i)+deltayj*temp4
                    strs(3,1,ii + nnn,i)=strs(3,1,ii + nnn,i)+deltazj*temp4
                    ! dxx/dy
                    dxdy(ii,i,i,2)=dxdy(ii,i,i,2)+&
                           (temp1*drijdyi*temp2*fcutij&
                          + temp2*dfcutijdyi)
                    dxdy(ii + nnn,i,i,2)=dxdy(ii + nnn,i,i,2)+&
                           (temp1*drijdyi*temp2*fcutij&
                          + temp2*dfcutijdyi)*weights
                    temp3= temp1*drijdyj*temp2*fcutij + temp2*dfcutijdyj
                    dxdy(ii,i,n,2)=dxdy(ii,i,n,2)+temp3
                    temp4 = temp3 * weights
                    dxdy(ii + nnn,i,n,2)=dxdy(ii + nnn,i,n,2)+temp4
                    strs(1,2,ii,i)=strs(1,2,ii,i)+deltaxj*temp3
                    strs(2,2,ii,i)=strs(2,2,ii,i)+deltayj*temp3
                    strs(3,2,ii,i)=strs(3,2,ii,i)+deltazj*temp3

                    strs(1,2,ii + nnn,i)=strs(1,2,ii + nnn,i)+deltaxj*temp4
                    strs(2,2,ii + nnn,i)=strs(2,2,ii + nnn,i)+deltayj*temp4
                    strs(3,2,ii + nnn,i)=strs(3,2,ii + nnn,i)+deltazj*temp4
                    ! dxx/dz
                    dxdy(ii,i,i,3)=dxdy(ii,i,i,3)+&
                           (temp1*drijdzi*temp2*fcutij&
                          + temp2*dfcutijdzi)
                    dxdy(ii + nnn,i,i,3)=dxdy(ii + nnn,i,i,3)+&
                           (temp1*drijdzi*temp2*fcutij&
                          + temp2*dfcutijdzi)*weights
                    temp3=temp1*drijdzj*temp2*fcutij + temp2*dfcutijdzj
                    dxdy(ii,i,n,3)=dxdy(ii,i,n,3)+temp3
                    temp4 = temp3 * weights
                    dxdy(ii + nnn,i,n,3)=dxdy(ii + nnn,i,n,3)+temp4
                    strs(1,3,ii,i)=strs(1,3,ii,i)+deltaxj*temp3
                    strs(2,3,ii,i)=strs(2,3,ii,i)+deltayj*temp3
                    strs(3,3,ii,i)=strs(3,3,ii,i)+deltazj*temp3

                    strs(1,3,ii + nnn,i)=strs(1,3,ii + nnn,i)+deltaxj*temp4
                    strs(2,3,ii + nnn,i)=strs(2,3,ii + nnn,i)+deltayj*temp4
                    strs(3,3,ii + nnn,i)=strs(3,3,ii + nnn,i)+deltazj*temp4
                endif
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
                weights_j = neighbor(i, j_neighbor, 5)
                n = int(neighbor(i, j_neighbor, 6))
                fcutij=0.5d0*(dcos(pi*rij/cutoff)+1.d0)
                if (lgrad) then
                    deltaxj = -1.d0*(pos(i, 1) - xyz_j(1))
                    deltayj = -1.d0*(pos(i, 2) - xyz_j(2))
                    deltazj = -1.d0*(pos(i, 3) - xyz_j(3))
                    drijdxi = -1.d0*deltaxj/rij
                    drijdyi = -1.d0*deltayj/rij
                    drijdzi = -1.d0*deltazj/rij
                    drijdxj = -1.d0*drijdxi
                    drijdyj = -1.d0*drijdyi
                    drijdzj = -1.d0*drijdzi
                    drijdxk = 0.d0
                    drijdyk = 0.d0
                    drijdzk = 0.d0
                    temp1=0.5d0*(-dsin(pi*rij/cutoff))*(pi/cutoff)
                    dfcutijdxi=temp1*drijdxi
                    dfcutijdyi=temp1*drijdyi
                    dfcutijdzi=temp1*drijdzi
                    dfcutijdxj=-1.d0*dfcutijdxi
                    dfcutijdyj=-1.d0*dfcutijdyi
                    dfcutijdzj=-1.d0*dfcutijdzi
                    dfcutijdxk=0.0d0
                    dfcutijdyk=0.0d0
                    dfcutijdzk=0.0d0
                endif
                do k_neighbor = 1, neighbor_count(i)
                    if (k_neighbor <= j_neighbor) cycle
                    rik = neighbor(i,k_neighbor,4)
                    if (rik.gt.cutoff) cycle
                    xyz_k = neighbor(i, k_neighbor,1:3)
                    weights_k = neighbor(i, k_neighbor,5)
                    m = int(neighbor(i, k_neighbor,6))
                    fcutik=0.5d0*(dcos(pi*rik/cutoff)+1.d0)

                    if (lgrad) then
                        deltaxk = -1.d0*(pos(i, 1) - xyz_k(1))
                        deltayk = -1.d0*(pos(i, 2) - xyz_k(2))
                        deltazk = -1.d0*(pos(i, 3) - xyz_k(3))
                        drikdxi = -deltaxk/rik
                        drikdyi = -deltayk/rik
                        drikdzi = -deltazk/rik
                        drikdxk = -1.d0*drikdxi
                        drikdyk = -1.d0*drikdyi
                        drikdzk = -1.d0*drikdzi
                        drikdxj = 0.d0
                        drikdyj = 0.d0
                        drikdzj = 0.d0
                        temp1=0.5d0*(-dsin(pi*rik/cutoff))*(pi/cutoff)
                        dfcutikdxi=temp1*drikdxi
                        dfcutikdyi=temp1*drikdyi
                        dfcutikdzi=temp1*drikdzi
                        dfcutikdxj=0.0d0
                        dfcutikdyj=0.0d0
                        dfcutikdzj=0.0d0
                        dfcutikdxk=-1.d0*dfcutikdxi
                        dfcutikdyk=-1.d0*dfcutikdyi
                        dfcutikdzk=-1.d0*dfcutikdzi
                    endif
                    rjk = (xyz_j(1) - xyz_k(1))**2 + (xyz_j(2) - xyz_k(2))**2 + (xyz_j(3) - xyz_k(3))**2
                    rjk = dsqrt(rjk)

                    if (rjk.gt.cutoff) cycle  ! Be careful STUPID!!!
                    if (rjk < Rmin) then
                    !    print*, 'Rjk', rjk,' smaller than Rmin'
                    !    stop
                    endif
                    fcutjk=0.5d0*(dcos(pi*rjk/cutoff)+1.d0)
                    if (lgrad) then
                        drjkdxj = (xyz_j(1) - xyz_k(1))/rjk
                        drjkdyj = (xyz_j(2) - xyz_k(2))/rjk
                        drjkdzj = (xyz_j(3) - xyz_k(3))/rjk
                        drjkdxk = -1.d0*drjkdxj
                        drjkdyk = -1.d0*drjkdyj
                        drjkdzk = -1.d0*drjkdzj
                        drjkdxi = 0.d0
                        drjkdyi = 0.d0
                        drjkdzi = 0.d0
                        temp1=0.5d0*(-dsin(pi*rjk/cutoff))*(pi/cutoff)
                        dfcutjkdxj=temp1*drjkdxj
                        dfcutjkdyj=temp1*drjkdyj
                        dfcutjkdzj=temp1*drjkdzj
                        dfcutjkdxk=-1.d0*dfcutjkdxj
                        dfcutjkdyk=-1.d0*dfcutjkdyj
                        dfcutjkdzk=-1.d0*dfcutjkdzj
                        dfcutjkdxi=0.0d0
                        dfcutjkdyi=0.0d0
                        dfcutjkdzi=0.0d0
                    endif

                    f=rjk**2 - rij**2 -rik**2
                    g=-2.d0*rij*rik
                    costheta=f/g
                    costheta=1.d0 - costheta  ! avoid negative values
                    if (lgrad) then
                        dfdxi=-2.d0*rij*drijdxi - 2.d0*rik*drikdxi
                        dfdyi=-2.d0*rij*drijdyi - 2.d0*rik*drikdyi
                        dfdzi=-2.d0*rij*drijdzi - 2.d0*rik*drikdzi

                        dfdxj=2.d0*rjk*drjkdxj - 2.d0*rij*drijdxj
                        dfdyj=2.d0*rjk*drjkdyj - 2.d0*rij*drijdyj
                        dfdzj=2.d0*rjk*drjkdzj - 2.d0*rij*drijdzj

                        dfdxk=2.d0*rjk*drjkdxk - 2.d0*rik*drikdxk
                        dfdyk=2.d0*rjk*drjkdyk - 2.d0*rik*drikdyk
                        dfdzk=2.d0*rjk*drjkdzk - 2.d0*rik*drikdzk

                        dgdxi=-2.d0*(drijdxi*rik + rij*drikdxi)
                        dgdyi=-2.d0*(drijdyi*rik + rij*drikdyi)
                        dgdzi=-2.d0*(drijdzi*rik + rij*drikdzi)

                        dgdxj=-2.d0*drijdxj*rik
                        dgdyj=-2.d0*drijdyj*rik
                        dgdzj=-2.d0*drijdzj*rik

                        dgdxk=-2.d0*rij*drikdxk
                        dgdyk=-2.d0*rij*drikdyk
                        dgdzk=-2.d0*rij*drikdzk

                        temp1=1.d0/g**2
                        !!!! Be careful costheta = 1.d0 - costheta 2019.07.25
                        dcosthetadxi=-1.d0 * (dfdxi*g - f*dgdxi)*temp1  
                        dcosthetadyi=-1.d0 * (dfdyi*g - f*dgdyi)*temp1 
                        dcosthetadzi=-1.d0 * (dfdzi*g - f*dgdzi)*temp1 
                        dcosthetadxj=-1.d0 * (dfdxj*g - f*dgdxj)*temp1 
                        dcosthetadyj=-1.d0 * (dfdyj*g - f*dgdyj)*temp1 
                        dcosthetadzj=-1.d0 * (dfdzj*g - f*dgdzj)*temp1 
                        dcosthetadxk=-1.d0 * (dfdxk*g - f*dgdxk)*temp1 
                        dcosthetadyk=-1.d0 * (dfdyk*g - f*dgdyk)*temp1 
                        dcosthetadzk=-1.d0 * (dfdzk*g - f*dgdzk)*temp1 
                    endif

                    expxyz=dexp(-alpha*(rij**2+rik**2+rjk**2))

                    xx(ii,i)=xx(ii,i)+costheta*expxyz*fcutij*fcutik*fcutjk
                    xx(ii + nnn,i)=xx(ii + nnn,i)+&
                    costheta*expxyz*fcutij*fcutik*fcutjk*weights_j*weights_k
                    if (lgrad) then
                        temp1=-alpha*2.0d0*expxyz
                        dexpxyzdxi=(rij*drijdxi+rik*drikdxi+rjk*drjkdxi)*temp1
                        dexpxyzdyi=(rij*drijdyi+rik*drikdyi+rjk*drjkdyi)*temp1
                        dexpxyzdzi=(rij*drijdzi+rik*drikdzi+rjk*drjkdzi)*temp1
                        dexpxyzdxj=(rij*drijdxj+rik*drikdxj+rjk*drjkdxj)*temp1
                        dexpxyzdyj=(rij*drijdyj+rik*drikdyj+rjk*drjkdyj)*temp1
                        dexpxyzdzj=(rij*drijdzj+rik*drikdzj+rjk*drjkdzj)*temp1
                        dexpxyzdxk=(rij*drijdxk+rik*drikdxk+rjk*drjkdxk)*temp1
                        dexpxyzdyk=(rij*drijdyk+rik*drikdyk+rjk*drjkdyk)*temp1
                        dexpxyzdzk=(rij*drijdzk+rik*drikdzk+rjk*drjkdzk)*temp1
                        temp1=(dcosthetadxi*expxyz*fcutij*fcutik*fcutjk&
                              +costheta*dexpxyzdxi*fcutij*fcutik*fcutjk&
                              +costheta*expxyz*dfcutijdxi*fcutik*fcutjk&
                              +costheta*expxyz*fcutij*dfcutikdxi*fcutjk&
                              +costheta*expxyz*fcutij*fcutik*dfcutjkdxi)
                        temp2=(dcosthetadxj*expxyz*fcutij*fcutik*fcutjk&
                              +costheta*dexpxyzdxj*fcutij*fcutik*fcutjk&
                              +costheta*expxyz*dfcutijdxj*fcutik*fcutjk&
                              +costheta*expxyz*fcutij*dfcutikdxj*fcutjk&
                              +costheta*expxyz*fcutij*fcutik*dfcutjkdxj)
                        temp3=(dcosthetadxk*expxyz*fcutij*fcutik*fcutjk&
                              +costheta*dexpxyzdxk*fcutij*fcutik*fcutjk&
                              +costheta*expxyz*dfcutijdxk*fcutik*fcutjk&
                              +costheta*expxyz*fcutij*dfcutikdxk*fcutjk&
                              +costheta*expxyz*fcutij*fcutik*dfcutjkdxk)
                        temp4 = temp1 * weights_j * weights_k
                        temp5 = temp2 * weights_j * weights_k
                        temp6 = temp3 * weights_j * weights_k
                        dxdy(ii,i,i,1)=dxdy(ii,i,i,1)+temp1
                        dxdy(ii,i,n,1)=dxdy(ii,i,n,1)+temp2
                        dxdy(ii,i,m,1)=dxdy(ii,i,m,1)+temp3
                        dxdy(ii + nnn,i,i,1)=dxdy(ii + nnn,i,i,1)+temp4
                        dxdy(ii + nnn,i,n,1)=dxdy(ii + nnn,i,n,1)+temp5
                        dxdy(ii + nnn,i,m,1)=dxdy(ii + nnn,i,m,1)+temp6

                        strs(1,1,ii,i)=strs(1,1,ii,i)+deltaxj*temp2+deltaxk*temp3
                        strs(2,1,ii,i)=strs(2,1,ii,i)+deltayj*temp2+deltayk*temp3
                        strs(3,1,ii,i)=strs(3,1,ii,i)+deltazj*temp2+deltazk*temp3
                        strs(1,1,ii + nnn,i)=strs(1,1,ii + nnn,i)+deltaxj*temp5+deltaxk*temp6
                        strs(2,1,ii + nnn,i)=strs(2,1,ii + nnn,i)+deltayj*temp5+deltayk*temp6
                        strs(3,1,ii + nnn,i)=strs(3,1,ii + nnn,i)+deltazj*temp5+deltazk*temp6
                        ! dxxii/dy_i
                        temp1=(dcosthetadyi*expxyz*fcutij*fcutik*fcutjk&
                              +costheta*dexpxyzdyi*fcutij*fcutik*fcutjk&
                              +costheta*expxyz*dfcutijdyi*fcutik*fcutjk&
                              +costheta*expxyz*fcutij*dfcutikdyi*fcutjk&
                              +costheta*expxyz*fcutij*fcutik*dfcutjkdyi)
                        temp2=(dcosthetadyj*expxyz*fcutij*fcutik*fcutjk&
                              +costheta*dexpxyzdyj*fcutij*fcutik*fcutjk&
                              +costheta*expxyz*dfcutijdyj*fcutik*fcutjk&
                              +costheta*expxyz*fcutij*dfcutikdyj*fcutjk&
                              +costheta*expxyz*fcutij*fcutik*dfcutjkdyj)
                        temp3=(dcosthetadyk*expxyz*fcutij*fcutik*fcutjk&
                              +costheta*dexpxyzdyk*fcutij*fcutik*fcutjk&
                              +costheta*expxyz*dfcutijdyk*fcutik*fcutjk&
                              +costheta*expxyz*fcutij*dfcutikdyk*fcutjk&
                              +costheta*expxyz*fcutij*fcutik*dfcutjkdyk)
                        temp4 = temp1 * weights_j * weights_k
                        temp5 = temp2 * weights_j * weights_k
                        temp6 = temp3 * weights_j * weights_k
                        dxdy(ii,i,i,2)=dxdy(ii,i,i,2)+temp1
                        dxdy(ii,i,n,2)=dxdy(ii,i,n,2)+temp2
                        dxdy(ii,i,m,2)=dxdy(ii,i,m,2)+temp3

                        dxdy(ii + nnn,i,i,2)=dxdy(ii + nnn,i,i,2)+temp4
                        dxdy(ii + nnn,i,n,2)=dxdy(ii + nnn,i,n,2)+temp5
                        dxdy(ii + nnn,i,m,2)=dxdy(ii + nnn,i,m,2)+temp6
                        strs(1,2,ii,i)=strs(1,2,ii,i)+deltaxj*temp2+deltaxk*temp3
                        strs(2,2,ii,i)=strs(2,2,ii,i)+deltayj*temp2+deltayk*temp3
                        strs(3,2,ii,i)=strs(3,2,ii,i)+deltazj*temp2+deltazk*temp3

                        strs(1,2,ii + nnn,i)=strs(1,2,ii + nnn,i)+deltaxj*temp5+deltaxk*temp6
                        strs(2,2,ii + nnn,i)=strs(2,2,ii + nnn,i)+deltayj*temp5+deltayk*temp6
                        strs(3,2,ii + nnn,i)=strs(3,2,ii + nnn,i)+deltazj*temp5+deltazk*temp6
                        ! dxxii/dz_i
                        temp1=(dcosthetadzi*expxyz*fcutij*fcutik*fcutjk&
                              +costheta*dexpxyzdzi*fcutij*fcutik*fcutjk&
                              +costheta*expxyz*dfcutijdzi*fcutik*fcutjk&
                              +costheta*expxyz*fcutij*dfcutikdzi*fcutjk&
                              +costheta*expxyz*fcutij*fcutik*dfcutjkdzi)
                        temp2=(dcosthetadzj*expxyz*fcutij*fcutik*fcutjk&
                              +costheta*dexpxyzdzj*fcutij*fcutik*fcutjk&
                              +costheta*expxyz*dfcutijdzj*fcutik*fcutjk&
                              +costheta*expxyz*fcutij*dfcutikdzj*fcutjk&
                              +costheta*expxyz*fcutij*fcutik*dfcutjkdzj)
                        temp3=(dcosthetadzk*expxyz*fcutij*fcutik*fcutjk&
                              +costheta*dexpxyzdzk*fcutij*fcutik*fcutjk&
                              +costheta*expxyz*dfcutijdzk*fcutik*fcutjk&
                              +costheta*expxyz*fcutij*dfcutikdzk*fcutjk&
                              +costheta*expxyz*fcutij*fcutik*dfcutjkdzk)
                        temp4 = temp1 * weights_j * weights_k
                        temp5 = temp2 * weights_j * weights_k
                        temp6 = temp3 * weights_j * weights_k
                        dxdy(ii,i,i,3)=dxdy(ii,i,i,3)+temp1
                        dxdy(ii,i,n,3)=dxdy(ii,i,n,3)+temp2
                        dxdy(ii,i,m,3)=dxdy(ii,i,m,3)+temp3

                        dxdy(ii + nnn,i,i,3)=dxdy(ii + nnn,i,i,3)+temp4
                        dxdy(ii + nnn,i,n,3)=dxdy(ii + nnn,i,n,3)+temp5
                        dxdy(ii + nnn,i,m,3)=dxdy(ii + nnn,i,m,3)+temp6
                        strs(1,3,ii,i)=strs(1,3,ii,i)+deltaxj*temp2+deltaxk*temp3
                        strs(2,3,ii,i)=strs(2,3,ii,i)+deltayj*temp2+deltayk*temp3
                        strs(3,3,ii,i)=strs(3,3,ii,i)+deltazj*temp2+deltazk*temp3

                        strs(1,3,ii + nnn,i)=strs(1,3,ii + nnn,i)+deltaxj*temp5+deltaxk*temp6
                        strs(2,3,ii + nnn,i)=strs(2,3,ii + nnn,i)+deltayj*temp5+deltayk*temp6
                        strs(3,3,ii + nnn,i)=strs(3,3,ii + nnn,i)+deltazj*temp5+deltazk*temp6
                    endif ! lgrad
                enddo ! k_neighbor
            enddo ! j_neighbor
        enddo ! i
    else
        print *, 'Unknown function type',ii, ACSF%sf(ii)%ntype
    endif
enddo  ! types
END SUBROUTINE
END SUBROUTINE!}}}
SUBROUTINE  write_array_2dim(n,m, a,name)!{{{
REAL(8),intent(in),dimension(n,m)       :: a
character(*),intent(in)                :: name
integer                                  :: i,j
open(2244,file=trim(adjustl(name)))
do i = 1, n
    do j = 1, m
        write(2244,'(F20.10,$)') a(i,j)
    enddo
    write(2244,*)
enddo
close(2244)
END SUBROUTINE!}}}
SUBROUTINE  SET_ELEMENTS(nat)!{{{
Integer, intent(in)    :: nat(:)
Integer                :: temp_mat(10)
Integer                :: temp_name(10)
integer                :: na, xx, k, index, i
na = size(nat)
temp_mat = 0
xx = nat(1)
k = 1
index = 1
temp_name(1) = nat(1)
do i =  2,na
    if (nat(i) == xx) then
        k = k + 1
    else 
       xx = nat(i)
       k = 1
       index = index + 1
       temp_name(index) = nat(i)
    endif
    temp_mat(index) = k
enddo
allocate(gap_natoms(index))
allocate(gap_species(index))
do i = 1, index
    gap_natoms(i) = temp_mat(i) 
    gap_species(i) = atsym(temp_name(i))
    !call get_element_name(temp_name(i), gap_species(i))
enddo
END SUBROUTINE!}}}

END MODULE

