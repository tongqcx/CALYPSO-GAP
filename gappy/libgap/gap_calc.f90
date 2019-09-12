SUBROUTINE  FGAP_CALC(NA, SPECIES, LAT, POS, &
                     ENE, FORCE, STRESS, VARIANCE, &
                     sparse_cov_mat_len,des_len, &
                     theta, MM, qmm, coeff)
USE GAP_INIT

implicit none

integer, intent(in)                                                            :: NA
integer, intent(in),dimension(NA)                                              :: SPECIES
double precision, intent(in),dimension(3,3)                                    :: LAT
double precision, intent(in),dimension(NA,3)                                   :: POS

integer, intent(in)                                                            :: des_len
integer, intent(in)                                                            :: sparse_cov_mat_len
double precision, intent(in),dimension(des_len)                                :: THETA
double precision, intent(in),dimension(sparse_cov_mat_len,des_len)             :: MM
double precision, intent(in),dimension(sparse_cov_mat_len,sparse_cov_mat_len)  :: QMM
double precision, intent(in),dimension(sparse_cov_mat_len)                     :: COEFF

double precision, intent(out)                                                  :: ENE
double precision, intent(out)                                                  :: VARIANCE
double precision, intent(out),dimension(NA,3)                                  :: FORCE
double precision, intent(out),dimension(6)                                     :: STRESS

!-- local --
type(Structure)                                                                :: at_test
integer                                                                        :: i,j,k
logical,parameter                                                              :: lcrystal = .true.

!--print*, Na,'Na'
if (Na > tna) then
    print *, 'The number of atom in structure large than tna=',tna, &
               "and you should modify it in gap_init.f90"
    stop
endif
at_test%na = NA
at_test%lat = LAT
do i = 1,at_test%na
    at_test%xyz(i,:) = POS(i,:)
    call get_ele_weights(SPECIES(i),at_test%types(i))
enddo
call functions(1,at_test%na,des_len,&
                 at_test%types(1:at_test%na), &
                 at_test%lat,at_test%xyz(1:at_test%na,:),&
                 at_test%xx(1:des_len,1,1:at_test%na),&
                 .true.,lcrystal,&
                 at_test%dxdy(1:des_len,1:at_test%na,1:at_test%na,:),&
                 at_test%strs(:,:,1:des_len,1:at_test%na))
do i = 1,at_test%na
    at_test%kk(i,1:des_len) = at_test%xx(1:des_len,1,i)
enddo

!--print*, 'get_cov'
call get_cov(at_test%na,&
             sparse_cov_mat_len,&
             at_test%kk(1:at_test%na,1:des_len),&
             mm(1:sparse_cov_mat_len,1:des_len),theta(1:des_len), &
             at_test%ckm(1:at_test%na,1:sparse_cov_mat_len))

!--print*, 'get energy'
at_test%e = matmul(at_test%ckm(1:at_test%na,1:sparse_cov_mat_len),coeff(1:sparse_cov_mat_len))

at_test%energy_cal = sum(at_test%e)
at_test%force_cal = 0.d0
at_test%stress_cal = 0.d0
at_test%derf = 0.d0
at_test%ders = 0.d0

!--print*, 'get de/dd'
do i = 1 , at_test%na
    do j =  1,des_len
        do k = 1,sparse_cov_mat_len
            at_test%derf(i,j) = at_test%derf(i,j) - 1.d0*(at_test%kk(i,j)-mm(k,j))/theta(j)**2*at_test%ckm(i,k)*coeff(k)         
            at_test%ders(i,j) =  at_test%derf(i,j)       
        enddo
    enddo
enddo

!--print*, 'get force'
call nnpotf(at_test%derf(1:at_test%na,1:des_len),&
            at_test%force_cal(1:at_test%na,:),&
            at_test%dxdy(1:des_len,1:at_test%na,1:at_test%na,:),&
            at_test%na,des_len,.true.)

!--print*, 'get stress'
call nnpots(at_test%ders(1:at_test%na,1:des_len),&
            at_test%strs(:,:,1:des_len,1:at_test%na),&
            at_test%stress_cal,&
            at_test%na,des_len,.true.)
at_test%volume = ABS(det(at_test%lat))
at_test%stress_cal = at_test%stress_cal * (1.0/GPa2eVPang) * 10.0 / at_test%volume

!--print*, 'get variance'
VARIANCE = 0.d0
do i  = 1,at_test%na
    at_test%covf(i) = delta - dot_product(at_test%ckm(i,1:sparse_cov_mat_len),matmul(qmm,at_test%ckm(i,1:sparse_cov_mat_len)))
    VARIANCE = VARIANCE + at_test%covf(i)/at_test%na
enddo

ENE = at_test%energy_cal
FORCE = at_test%force_cal(1:at_test%na,:)
STRESS = at_test%stress_cal

END SUBROUTINE

SUBROUTINE  FGAP_READ(sparse_cov_mat_len, des_len, theta, MM, invcmm, coeff)

USE GAP_INIT, only: tna,tnf,tcc
implicit none

integer, intent(out)                :: sparse_cov_mat_len
integer, intent(out)                :: des_len
double precision, intent(out)       :: theta(tnf)
double precision, intent(out)       :: MM(tcc, tnf)
double precision, intent(out)       :: invcmm(tcc, tcc)
double precision, intent(out)       :: coeff(tcc)

!-- local--
logical                             :: alive
double precision                    :: mean_e,cov_e, &
                                       mean_f,cov_f, &
                                       mean_s,cov_s
integer                             :: i
    inquire(file="gap_parameters",exist=alive)
    if(.not.alive) then
        print*, "gap_parameters file does not exist!"
        stop
    endif
    open(111,file="gap_parameters") 
    read(111,*) sparse_cov_mat_len,des_len
    if (sparse_cov_mat_len > tcc) then
        print *, 'The siez of sparse set large than tcc=',tcc, &
                   "and you should modify it in gap_init.f90"
        stop
    endif
    if (des_len > tnf) then
        print *, 'The length of descriptors large than tnf=',tnf, &
                   "and you should modify it in gap_init.f90"
        stop
    endif
    read(111,*) mean_e,cov_e
    read(111,*) mean_f,cov_f
    read(111,*) mean_s,cov_s
    read(111,*) theta(1:des_len)
    do i = 1,sparse_cov_mat_len
        read(111,*) MM(i,1:des_len)
    enddo
    do i = 1,sparse_cov_mat_len
        read(111,*) INVCMM(i,1:sparse_cov_mat_len)
    enddo
    read(111,*) COEFF(1:sparse_cov_mat_len)
    close(111)
END SUBROUTINE

