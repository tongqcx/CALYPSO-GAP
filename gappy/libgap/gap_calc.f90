SUBROUTINE  FGAP_CALC(NA, SPECIES, LAT, POS, &
                     ENE, FORCE, STRESS, VARIANCE, &
                     nsparseX,des_len, &
                     theta, MM, qmm, coeff)
USE GAP_INIT

implicit none

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

double precision, intent(out)                                                  :: ENE
double precision, intent(out)                                                  :: VARIANCE
double precision, intent(out),dimension(NA,3)                                  :: FORCE
double precision, intent(out),dimension(6)                                     :: STRESS

!-- local --
type(Structure)                                                                :: at
integer                                                                        :: i,j,k, k1, n
double precision                                                               :: Rcut
logical                                                                        :: lgrad


!--print*, Na,'Na'
if (Na > natoms_max) then
    print *, 'The number of atom in structure large than natoms_max=',natoms_max, &
               "and you should modify it in gap_init.f90"
    stop
endif
Rcut = 6.0
lgrad = .true.
at%na = NA
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
print*, 'BEGIN gal_calc'
do i = 1,at%na
    at%kk(i,1:des_len) = at%xx(1:des_len,i)
enddo

!--print*, 'get_cov'
call get_cov(at%na,&
             nsparseX,&
             at%kk(1:at%na,1:des_len),&
             mm(1:nsparseX,1:des_len),theta(1:des_len), &
             at%ckm(1:at%na,1:nsparseX))

!--print*, 'get energy'
at%e = matmul(at%ckm(1:at%na,1:nsparseX),coeff(1:nsparseX))
at%energy_cal = sum(at%e)

at%dedg = 0.d0
!--print*, 'get de/dd'
do i = 1 , at%na
    do j =  1,des_len
        do k = 1,nsparseX
            at%dedg(i,j) = at%dedg(i,j) - 1.d0*(at%kk(i,j)-mm(k,j))/theta(j)**2*at%ckm(i,k)*coeff(k)         
        enddo
    enddo
enddo

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
!at%stress_cal = at%stress_cal * (1.0/GPa2eVPang) * 10.0 / at%volume

!--print*, 'get variance'
!VARIANCE = 0.d0
!do i  = 1,at%na
!    at%covf(i) = delta - dot_product(at%ckm(i,1:nsparseX),matmul(qmm,at%ckm(i,1:nsparseX)))
!    VARIANCE = VARIANCE + at%covf(i)/at%na
!enddo

ENE = at%energy_cal
FORCE = at%force_cal(1:at%na,:)
STRESS = at%stress_cal

END SUBROUTINE

SUBROUTINE  FGAP_READ(nsparseX, des_len, theta, MM, invcmm, coeff)

USE GAP_INIT, only: natoms_max, nsf_max, nsparseX_max
implicit none

integer, intent(out)                :: nsparseX
integer, intent(out)                :: des_len
double precision, intent(out)       :: theta(nsf_max)
double precision, intent(out)       :: MM(nsparseX_max, nsf_max)
double precision, intent(out)       :: invcmm(nsparseX_max, nsparseX_max)
double precision, intent(out)       :: coeff(nsparseX_max)

!-- local--
logical                             :: alive
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
    do i = 1,nsparseX
        read(111,*) INVCMM(i,1:nsparseX)
    enddo
    read(111,*) COEFF(1:nsparseX)
    close(111)
END SUBROUTINE

