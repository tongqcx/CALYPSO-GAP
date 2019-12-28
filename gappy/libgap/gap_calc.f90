SUBROUTINE  FGAP_CALC(NA, SPECIES, LAT, POS, &
                     ENE, FORCE, STRESS, VARIANCE, &
                     nsparseX,des_len, &
                     theta, MM, qmm, coeff, &
                     Rcut, lgrad)

implicit none
double precision,parameter        :: delta = 1.0
double precision,parameter        :: GPa2eVPang =6.24219D-3
type Structure
    integer                                                                    :: na
    integer                                                                    :: ntype
    double precision                                                           :: lat(3,3)
    double precision                                                           :: volume
    double precision, allocatable, dimension(:,:)                              :: xyz
    double precision, allocatable, dimension(:,:)                              :: xx
    double precision, allocatable, dimension(:,:)                              :: kk
    double precision, allocatable, dimension(:,:)                              :: ckm  ! k in the number of atoms
    double precision, allocatable, dimension(:)                                :: e
    double precision                                                           :: energy_ref, energy_cal
    double precision, dimension(6)                                             :: stress_ref, stress_cal
    double precision, allocatable, dimension(:,:)                              :: force_ref, force_cal
    double precision, allocatable, dimension(:,:,:,:)                          :: dxdy
    double precision, allocatable, dimension(:,:,:,:)                          :: strs
    double precision, allocatable, dimension(:,:)                              :: dedg
    double precision, allocatable, dimension(:,:,:)                            :: dkdg
end type

integer, intent(in)                                                            :: NA
integer, intent(in),dimension(NA)                                              :: SPECIES
double precision, intent(in),dimension(3,3)                                    :: LAT
double precision, intent(in),dimension(NA,3)                                   :: POS

integer, intent(in)                                                            :: des_len
integer, intent(in)                                                            :: nsparseX
double precision, intent(in),dimension(des_len)                                :: THETA
double precision, intent(in),dimension(nsparseX,des_len)                       :: MM
double precision, intent(in),dimension(nsparseX,nsparseX)                      :: QMM
double precision, intent(in),dimension(nsparseX)                               :: COEFF
double precision, intent(in)                                                   :: Rcut
logical         , intent(in)                                                   :: lgrad

double precision, intent(out)                                                  :: ENE
double precision, intent(out)                                                  :: VARIANCE
double precision, intent(out),dimension(NA,3)                                  :: FORCE
double precision, intent(out),dimension(6)                                     :: STRESS

!-- local --
type(Structure)                                                                :: at
integer                                                                        :: i,j,k, k1, n
!######################################################################
REAL(8),dimension(3,3)                         :: recip_lat
REAL(8)                                        :: rmin
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
INTEGER                                        :: natoms, n1, n2, n3
!real(8)                               :: vectorlength, recipvector, volume, crossp
!######################################################################
rmin = 0.5d0
max_neighbor = 1000

natoms = size(pos,1)
allocate(neighbor(natoms, max_neighbor, 6))
allocate(neighbor_count(natoms))
allocate(weights(natoms))

open(2233, file='gap_parameters')
read(2233,*) nspecies
do i = 1, nspecies
    read(2233,*) atom_number, atom_weight
    do j = 1, natoms
        if (species(j) == atom_number) weights(j) = atom_weight
    enddo
enddo
close(2233)

recip_lat = recipvector(lat)
nabc(1)=ceiling(rcut*vectorlength(recip_lat(1,:))/pi/2)
nabc(2)=ceiling(rcut*vectorlength(recip_lat(2,:))/pi/2)
nabc(3)=ceiling(rcut*vectorlength(recip_lat(3,:))/pi/2)
neighbor = 0.d0
neighbor_count = 0
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

allocate(at%xx(des_len, natoms))
allocate(at%dxdy(des_len, natoms, natoms, 3))
allocate(at%strs(3, 3, des_len, natoms))
call car2acsf(natoms, max_neighbor, des_len, pos, neighbor, neighbor_count, &
                    at%xx, &
                    at%dxdy, &
                    at%strs, lgrad)
deallocate(weights)
deallocate(neighbor)
deallocate(neighbor_count)
!###################################################

allocate(at%kk(natoms, des_len))  !&&&&&&&
at%na = natoms
do i = 1,at%na
    at%kk(i,:) = at%xx(:,i)
enddo
!call write_array_2dim(at%na, des_len, at%kk, 'kk.dat')

allocate(at%ckm(natoms, nsparseX))  !&&&&&&&
!--print*, 'get_cov'
call get_cov(at%na,&
             nsparseX,&
             at%kk,&
             mm(1:nsparseX,1:des_len),theta(1:des_len), &
             at%ckm)

allocate(at%e(natoms))   !&&&&&&&&&&&&&&
! at%e must be initial!! stupid 2019.10.25
!--print*, 'get energy'
at%e = 0.d0
at%e = matmul(at%ckm,coeff(1:nsparseX))
at%energy_cal = sum(at%e)

!call write_array_2dim(at%na, nsparseX, at%ckm(1:at%na,1:nsparseX), 'ckm.dat')
allocate(at%dedg(natoms, des_len))    !&&&&&&&&&&&&&
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
!&&&&&&&&&&&&&&&&&&&&&&
deallocate(at%ckm, at%kk, at%e)
!&&&&&&&&&&&&&&&&&&&&&&

allocate(at%force_cal(natoms, 3))  !&&&&&&&&&&&&&&&&&
at%force_cal = 0.d0
!--print*, 'get force'
!!$OMP parallel private(i, n, j, k)
!!$OMP DO
do i = 1, at%na
    do n = 1, at%na
        do j = 1,3
            do k = 1, des_len
                at%force_cal(i,j) = at%force_cal(i,j) - at%dedg(n,k) * at%dxdy(k,n,i,j)
            enddo
        enddo
    enddo
enddo
!!$OMP END PARALLEL 

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
at%volume = ABS(det(lat))
!at%stress_cal = at%stress_cal * (1.0/GPa2eVPang) * 10.0 / at%volume
at%stress_cal = at%stress_cal * (1.0/GPa2eVPang) / at%volume

!--print*, 'get variance'
VARIANCE = 0.d0
!do i  = 1,at%na
!    at%covf(i) = delta - dot_product(at%ckm(i,1:nsparseX),matmul(qmm,at%ckm(i,1:nsparseX)))
!    VARIANCE = VARIANCE + at%covf(i)/at%na
!enddo
!&&&&&&&&&&&&&&&&&&&&&&
deallocate(at%dedg)
!&&&&&&&&&&&&&&&&&&&&&&

ENE = at%energy_cal
FORCE = at%force_cal(1:at%na,:)
!&&&&&&&&&&&&&&&&&&&&&&
deallocate(at%force_cal)
!&&&&&&&&&&&&&&&&&&&&&&

STRESS(1) = at%stress_cal(1)
STRESS(2) = at%stress_cal(4)
STRESS(3) = at%stress_cal(6)
STRESS(4) = at%stress_cal(2)
STRESS(5) = at%stress_cal(5)
STRESS(6) = at%stress_cal(3)

!print *, 'ene',ENE
!do i = 1,natoms
!    print*, FORCE(i,:)
!enddo
!print*, 'stress',STRESS
!print*, '//////////GAP/////////////'
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

SUBROUTINE GET_COV(n,m,xx,yy,theta, cov)!{{{
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

FUNCTION Det(matrix)!{{{
   double precision,  intent(in) :: Matrix(3,3)
   double precision              :: Det


   Det = Matrix(1,1)*(Matrix(2,2)*Matrix(3,3)-Matrix(2,3)*Matrix(3,2))&
      -Matrix(1,2)*(Matrix(2,1)*Matrix(3,3)-Matrix(2,3)*Matrix(3,1))&
      +Matrix(1,3)*(Matrix(2,1)*Matrix(3,2)-Matrix(3,1)*Matrix(2,2))

END FUNCTION!}}}

END SUBROUTINE FGAP_CALC

SUBROUTINE  FGAP_READ(nsparseX, des_len, theta, MM, invcmm, coeff)

implicit none
integer, parameter                  :: nsf_max = 100
integer, parameter                  :: nsparseX_max = 4000

integer, intent(out)                :: nsparseX
integer, intent(out)                :: des_len
double precision, intent(out)       :: theta(nsf_max)
double precision, intent(out)       :: MM(nsparseX_max, nsf_max)
double precision, intent(out)       :: invcmm(nsparseX_max, nsparseX_max)
double precision, intent(out)       :: coeff(nsparseX_max)

!-- local--
logical                             :: alive
integer                             :: nspecies, nsf
!double precision                    :: mean_e,cov_e, &
!                                       mean_f,cov_f, &
!                                       mean_s,cov_s
integer                             :: i
    inquire(file="gap_parameters",exist=alive)
    if(.not.alive) then
        print*, "gap_parameters file does not exist!"
        stop
    endif
    open(111,file="gap_parameters") 

    read(111,*) nspecies
    do i = 1, nspecies
        read(111,*)
    enddo

    read(111,*) nsf
    do i = 1, nsf
        read(111,*)
    enddo

    read(111,*) nsparseX,des_len
    if (nsparseX > nsparseX_max) then
        print *, 'The siez of sparse set large than nsparseX_max=',nsparseX_max, &
                   "and you should modify it in gap_init.f90"
        stop
    endif
    if (des_len > nsf_max) then
        print *, 'The length of descriptors large than nsf_max=',nsf_max, &
                   "and you should modify it in gap_init.f90"
        stop
    endif
    read(111,*) 
    read(111,*) 
    read(111,*) 
    read(111,*) theta(1:des_len)
    do i = 1,nsparseX
        read(111,*) MM(i,1:des_len)
    enddo
    !do i = 1,nsparseX
    !    read(111,*) INVCMM(i,1:nsparseX)
    !enddo
    INVCMM = 0.d0
    read(111,*) COEFF(1:nsparseX)
    close(111)
END SUBROUTINE

