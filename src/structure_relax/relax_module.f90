! this module perform structure optimization task
! 2019.09.19
module relax_module
implicit none


private   struct2relaxv, relaxv2struct, lat2matrix, upper, lower, lat_inv
contains
SUBROUTINE  relax_main(NA, SPECIES, LAT, POS, EXTSTRESS)!{{{
implicit none
double precision, intent(in)                  :: NA
double precision, intent(in),dimension(NA)    :: SPECIES
double precision, intent(in),dimension(3,3)   :: LAT
double precision, intent(in),dimension(Na,3)  :: POS
double precision, intent(in),dimension(6)     :: EXTSTRESS

!local
double precision                              :: ENE, VARIANCE
double precision, allocatable,dimension(:,:)  :: FORCE
double precision,             dimension(6)    :: STRESS

!variables for lbfgs
integer                                       :: n, m, iprint
double precision, parameter                   :: factor = 1.0d7
double precision, parameter                   :: pgtol = 1.0d-5
character(len=60)                             :: task, csave
logical                                       :: lsave(4)
integer                                       :: isave(44)
double precision                              :: f
double precision                              :: dsave(29)
integer,  allocatable,dimension(:)            :: nbd, iwa
double precision, allocatable,dimension(:)    :: x, l, u, g, wa


n = 3*NA + 6
m = 5
iprint = 1
if (.not. allocate(FORCE))  allocate(FORCE(NA, 3))
allocate ( nbd(n), x(n), l(n), u(n), g(n) )
allocate ( iwa(3*n) )
allocate ( wa(2*m*n + 5*n + 11*m*m + 8*m) )

task = 'START'
call  FGAP_CALC(NA, SPECIES, LAT, POS, ENE, FOECE, STRESS, VARIANCE)

! begin lbfgs loop
do while(task(1:2).eq.'FG'.or.task.eq.'NEW_X'.or.task.eq.'START') 
 
    call struct2relaxv(NA, LAT, POS, ENE, FORCE, STRESS, EXTSTRESS, x, f, g)
    call setulb ( n, m, x, l, u, nbd, f, g, factr, pgtol, &
                       wa, iwa, task, iprint,&
                       csave, lsave, isave, dsave )
    if (task(1:2) .eq. 'FG') then
        call relaxv2struct(x, LAT, POS)
        call  FGAP_CALC(NA, SPECIES, LAT, POS, ENE, FOECE, STRESS, VARIANCE)
    endif
enddo
! end of lbfgs loop

END SUBROUTINE!}}}

SUBROUTINE  struct2relaxv(NA, LAT, POS, ENE, FORCE, STRESS, EXTSTRESS, n, x, f, g)!{{{
implicit none

INTEGER         , intent(in)                           :: NA
double precision, intent(in),dimension(3,3)            :: LAT
double precision, intent(in),dimension(NA,3)           :: POS, FORCE
double precision, intent(in)                           :: ENE
double precision, intent(in),dimension(6)              :: STRESS, EXTSTRESS
INTEGER         , intent(in)                           :: n
double precision, intent(inout),dimension(n)           :: x, g
double precision, intent(inout) i                      :: f
! local
integer                                                :: i,j,k
double precision, dimension(6)                         :: cellp, strderv, cellderv
double precision, allocatable, dimension(:,:)          :: POS_FRAC, FORCE_FRAC

if (.not. allocate(FORCE_FRAC)) allocate(FORCE_FRAC(NA, 3))
if (.not. allocate(POS_FRAC)) allocate(POS_FRAC(NA, 3))

CALL LAT2MATRIX(cellp, LAT, 2)
CALL CART2FRAC(NA, LAT, POS, POS_FRAC)
FORCE_FRAC = matmul(FORCE, transpose(LAT))
strderv(1) = STRESS(1)
strderv(2) = STRESS(4)
strderv(3) = STRESS(6)
strderv(4) = STRESS(5)
strderv(5) = STRESS(3)
strderv(6) = STRESS(2)
strderv = strderv - EXTSTRESS
! calculating the dev of strain to cell parameters
CALL DEVCELL(strderv, cellderv)
k = 0
do i = 1, 6
    k = k + 1
    xc(k) = cellp(i)
    gc(k) = cellderv(i)
enddo
do i = 1, NA
    do j = 1, 3
        k = k + 1
        xc(k) = POS_FRAC(i,j)
        gc(k) = FORCE_FRAC(i,j)
    enddo
enddo
f = ENE
END SUBROUTINE!}}}

SUBROUTINE  relaxv2struct(n, xc, NA, LAT, POS)!{{{
implicit none

INTEGER         , intent(in)                           :: n
double precision, intent(in),dimension(n)              :: xc
INTEGER         , intent(in)                           :: NA
double precision, intent(inout),dimension(3,3)         :: LAT
double precision, intent(inout),dimension(NA,3)        :: POS

! local
integer                                                :: i,j
double precision, dimension(6)                         :: cellp
cellp(1:6) = xc(1:6)
CALL LAT2MATRIX(cellp, LAT, 1)
do i = 1, NA
    do j = 1, 3
        POS_FRAC(i,j) = xc(6 + 3*(i - 1) + j)
    enddo
enddo
CALL FRAC2CART(NA, LAT, POS_FRAC, POS)
END SUBROUTINE!}}}
    
SUBROUTINE  CART2FRAC(NA, LAT, POS, POS_FCAR)!{{{
implicit none

INTEGER         , intent(in)                           :: NA
double precision, intent(in),dimension(3,3)            :: LAT
double precision, intent(in),dimension(NA,3)           :: POS
double precision, intent(inout),dimension(NA,3)        :: POS_FRAC

! local
integer                                                :: i,j
double precision, dimension(3,3)                       :: INV_LAT
call lat_inv(LAT, INV_LAT)
POS_FRAC = matmul(POS, INV_LAT)
do i = 1, NA
    do j = 1,3
        if (POS_FRAC(i,j) < 0.d0) POS_FRAC(i,j) = POS_FRAC(i,j) + 1.d0
        if (POS_FRAC(i,j) > 1.d0) POS_FRAC(i,j) = POS_FRAC(i,j) - 1.d0
    enddo
enddo
END SUBROUTINE CART2FRAC!}}}

SUBROUTINE  FRAC2CART(NA, LAT, POS_FRAC, POS)!{{{
implicit none

INTEGER         , intent(in)                           :: NA
double precision, intent(in),dimension(3,3)            :: LAT
double precision, intent(in),dimension(NA,3)           :: POS_FRAC
double precision, intent(inout),dimension(NA,3)        :: POS
do i = 1, NA
    do j = 1,3
        if (POS_FRAC(i,j) < 0.d0) POS_FRAC(i,j) = POS_FRAC(i,j) + 1.d0
        if (POS_FRAC(i,j) > 1.d0) POS_FRAC(i,j) = POS_FRAC(i,j) - 1.d0
    enddo
enddo
POS = matmul(POS_FRAC, LAT)
END SUBROUTINE FRAC2CART!}}}
        

subroutine lat2matrix(lat,matrix,iflag)!{{{
implicit none

! if iflag==1, abc2matrix; iflag==2. matrix2abc
integer,          intent(in)        :: iflag 
double precision, intent(inout)     :: lat(6),matrix(3,3)

!local parameters
double precision                    :: ra,rb,rc,&
                                       cosinea, cosineb,cosinec,&
                                       anglea,angleb,anglec

if (iflag==1) then
   matrix=0.0
   matrix(1,1) = lat(1)
   matrix(2,1) = lat(2)*cos(lat(6))
   matrix(2,2) = lat(2)*sin(lat(6))
   matrix(3,1) = lat(3)*cos(lat(5))
   matrix(3,2) = lat(3)*cos(lat(4))*sin(lat(6))-((lat(3)*cos(lat(5))&
   -lat(3)*cos(lat(4))*cos(lat(6)))/tan(lat(6)))
   matrix(3,3) = sqrt(lat(3)**2 -matrix(3,1)**2 - matrix(3,2)**2)
else
   lat=0.0
   ra=sqrt(matrix(1,1)**2+matrix(1,2)**2+matrix(1,3)**2)
   rb=sqrt(matrix(2,1)**2+matrix(2,2)**2+matrix(2,3)**2)
   rc=sqrt(matrix(3,1)**2+matrix(3,2)**2+matrix(3,3)**2)
   cosinea=(matrix(2,1)*matrix(3,1)+matrix(2,2)*matrix(3,2)+matrix(2,3)*matrix(3,3))/rb/rc
   cosineb=(matrix(1,1)*matrix(3,1)+matrix(1,2)*matrix(3,2)+matrix(1,3)*matrix(3,3))/ra/rc
   cosinec=(matrix(1,1)*matrix(2,1)+matrix(1,2)*matrix(2,2)+matrix(1,3)*matrix(2,3))/ra/rb
   anglea=acos(cosinea)
   angleb=acos(cosineb)
   anglec=acos(cosinec)
   lat(1)=ra
   lat(2)=rb
   lat(3)=rc
   lat(4)=anglea
   lat(5)=angleb
   lat(6)=anglec
endif
end subroutine lat2matrix!}}}

subroutine upper(matrix1,matrix2)!{{{
implicit none
real(8) :: matrix1(:,:)
real(8) :: matrix2(:,:)
integer :: i,j,m,n
real :: a 
m=size(matrix1,1)
n=size(matrix1,2)
do i=1,n
   do j=i+1,n
	  a=matrix1(j,i)/matrix1(i,i)
	  matrix1(j,:)=matrix1(j,:)-a*matrix1(i,:)
	  matrix2(j,:)=matrix2(j,:)-a*matrix2(i,:)
   end do
end do 
end subroutine upper!}}}

subroutine lower(matrix1,matrix2)!{{{
implicit none

real(8) :: matrix1(:,:)
real(8) :: matrix2(:,:)
integer :: i,j,m,n
real :: a
m=size(matrix1,1)
n=size(matrix1,2)
do i=n,2,-1
   do j=i-1,1,-1
	  a=matrix1(j,i)/matrix1(i,i)
	  matrix1(j,:)=matrix1(j,:)-a*matrix1(i,:)
	  matrix2(j,:)=matrix2(j,:)-a*matrix2(i,:)
   end do
end do 
end subroutine!}}}

subroutine lat_inv(matrix3,matrix2)!{{{
implicit none

real(8) :: matrix1(3,3)
real(8) :: matrix2(3,3),matrix3(3,3)
integer :: i,j
matrix1=matrix3
do i=1,3
   do j=1,3
	  if(i==j)then
		 matrix2(i,j)=1
	  else
		 matrix2(i,j)=0
	  end if
   end do 
end do
call upper(matrix1,matrix2)
call lower(matrix1,matrix2)
do i=1,3
   matrix2(i,:)=matrix2(i,:)/matrix1(i,i)
end do
end subroutine!}}}

END MODULE
