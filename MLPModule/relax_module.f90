! this module perform structure optimization task
! 2019.09.19
module relax_module
use read_module, ONLY: atsym
use gap_module
implicit none


private   struct2relaxv, relaxv2struct, lat2matrix, inv_33
contains

!CONJ
SUBROUTINE  relax_main_conj(NA, SPECIES, LAT, POS, EXTSTRESS, maxcycle, ftol, gtol)
implicit none
INTEGER,          intent(in)                     :: NA
INTEGER,          intent(in),dimension(NA)       :: SPECIES
double precision, intent(inout),dimension(3,3)   :: LAT
double precision, intent(inout),dimension(Na,3)  :: POS
double precision, intent(in),dimension(6)        :: EXTSTRESS
INTEGER,          intent(in)                     :: maxcycle
double precision, intent(in)                     :: ftol, gtol

!local
double precision                              :: ENE, VARIANCE
double precision, allocatable,dimension(:,:)  :: FORCE
double precision,             dimension(6)    :: STRESS

!variables for lbfgs
integer                                       :: n, nmin
double precision, allocatable,dimension(:)    :: x, g, xlast
double precision                              :: f
double precision, allocatable,dimension(:)    :: pvect, gg, glast
double precision                              :: pnorm, gsca, ggg, dggg, gam, pnlast, alp
logical                                       :: okf
integer                                       :: imode, jcyc, iflag

double precision                              :: f_bak, d, volume, df, norm_g
integer                                       :: i
double precision, allocatable,dimension(:)    :: x_save, dmax
logical                                       :: lfirst, lerror
integer                                       :: tt1, tt2
integer                                       :: lm_err
!double precision                              :: ftol, gtol



n = 3*NA + 6
if (.not. allocated(FORCE))  allocate(FORCE(NA, 3))
allocate (x(n), g(n), pvect(n), gg(n), glast(n), xlast(n))
print*, '***************************************************'
print*, '***            GEOMETRY OPTIMIZATION            ***' 
print*, '***************************************************'
write(*, '(A17, I8)'), 'Number of atoms =',NA
write(*, '(A21, I8)'), 'Number of variables =',n
write(*, '(A27, F8.3)'), 'Pressure of configuration =', SUM(EXTSTRESS(1:3))/3.d0
write(*, '(A10, F15.10)'), 'ftol =', ftol
write(*, '(A10, F15.10)'), 'gtol =', gtol

gg = 0.d0
nmin = 1
gsca = 0.001d0
volume = abs(det(lat))
pnorm = 1.d0/volume**(1.d0/3.d0)
pnlast = pnorm
jcyc = 0
alp = 1.d0
imode = 1
lm_err = 0
!ftol = 0.0005d0
!gtol = 0.005d0

lfirst = .true.
lerror = .false.
call  FGAP_INIT()
call  FGAP_CALC(NA, SPECIES, LAT, POS, ENE, FORCE, STRESS, VARIANCE)
call  struct2relaxv(NA, LAT, POS, ENE, FORCE, STRESS, EXTSTRESS, n, x, f, g)

f_bak = f
! begin lbfgs loop
CALL  SYSTEM_CLOCK(tt1)
do while(.true.) 
    okf = .true.
    if (lfirst) then
        do i = 1, n
            xlast(i) = x(i)
            glast(i) = -1.d0 * gsca * g(i)
            pvect(i) = glast(i)
            lfirst = .false.
        enddo
    endif
    ggg = 0.d0
    dggg = 0.d0
    do i = 1, n
        ggg = ggg + glast(i)*glast(i)
        dggg = dggg + (glast(i) + gsca*g(i)) * g(i) * gsca
    enddo
    gam = dggg/ggg
    !print*, 'gam',gam
    do i = 1, n
        xlast(i) = x(i)
        glast(i) = -1.d0 * gsca * g(i)
        pvect(i) = glast(i) + gam * pvect(i)
    enddo
    pnorm = gnorm(n, pvect) * n
    if (pnorm > 1.5d0 * pnlast) then
        pvect = pvect * 1.5d0* pnlast/pnorm
        pnorm = 1.5d0 * pnlast
    endif
    pnlast = pnorm
    !print *, 'pnorm', pnorm
    !print*, 'x'
    !print*, x
    !print*, 'pvect', pvect

    call olinmin(x, alp, pvect, n, nmin, f, okf, gg, imode, NA, SPECIES, EXTSTRESS)
    if (.not. okf) then
        print *, 'can not locate minimum, use default step'
        alp = 1.d0
        lm_err = lm_err + 1
        if (lm_err > 3) then
            print*, '**Optimization failure**'
            lerror = .true.
            !call write_fake_vasp(NA, SPECIES)
            !exit
        endif
    endif
        
    call funct(iflag, n, x, f, g, NA, SPECIES, EXTSTRESS)
    jcyc = jcyc + 1
    df = f - f_bak
    norm_g = gnorm(n, g)
    CALL  SYSTEM_CLOCK(tt2)
    write(*,'(''  Cycle: '',i6,''  Energy:'',f17.6,''  DE:'', f17.6,''  Gnorm:'',f14.6, ''  CPU:'',f8.3)') jcyc, f, df, norm_g, (tt2-tt1)/10000.0
    f_bak = f

    call relaxv2struct(n, x, NA, LAT, POS)
    call FGAP_CALC(NA, SPECIES, LAT, POS, ENE, FORCE, STRESS, VARIANCE)
    call write_vasp(NA, SPECIES, LAT, POS, ENE, FORCE, STRESS, EXTSTRESS)

    if (abs(df) < ftol) then
        write(*, *) '** Energy convergence**',ftol
        exit
    endif
    if (norm_g < gtol) then
        print*, '**Gtol satisfied**',gtol
        exit
    endif
    if (jcyc == maxcycle) then
        print*, '**Reach maxcycle**', maxcycle
        exit
    endif
 
enddo
print*, '' 
print*, 'END OF STRUCTURE OPTIMIZATION'
print*, '' 

END SUBROUTINE

SUBROUTINE funct(iflag, n, xc, fc, gc, NA, SPECIES, EXTSTRESS)
implicit none

integer, intent(inout)                       :: iflag
integer, intent(in)                          :: n
integer, intent(in)                          :: NA
integer, intent(in),           dimension(NA) :: SPECIES
double precision,intent(in),   dimension(6)  :: EXTSTRESS
double precision,intent(inout),dimension(n)  :: xc, gc
double precision,intent(inout)               :: fc

! local 
double precision, allocatable, dimension(:,:) :: POS, FORCE
double precision,              dimension(3,3) :: LAT
double precision,              dimension(6)   :: STRESS
double precision                              :: ENE, VARIANCE

if (.not. allocated(POS)) allocate(POS(NA, 3))
if (.not. allocated(FORCE)) allocate(FORCE(NA, 3))
call  relaxv2struct(n, xc, NA, LAT, POS)
call  FGAP_CALC(NA, SPECIES, LAT, POS, ENE, FORCE, STRESS, VARIANCE)
call  struct2relaxv(NA, LAT, POS, ENE, FORCE, STRESS, EXTSTRESS, n, xc, fc, gc)
END SUBROUTINE

SUBROUTINE  struct2relaxv(NA, LAT, POS, ENE, FORCE, STRESS, EXTSTRESS, n, xc, fc, gc)
implicit none

INTEGER         , intent(in)                           :: NA
double precision,parameter                             :: cfactor = 6.241460893d-3
double precision, intent(inout),dimension(3,3)            :: LAT
double precision, intent(in),dimension(NA,3)           :: POS, FORCE
double precision, intent(in)                           :: ENE
double precision, intent(in),dimension(6)              :: STRESS, EXTSTRESS
INTEGER         , intent(in)                           :: n
double precision, intent(inout),dimension(n)           :: xc, gc
double precision, intent(inout)                        :: fc
! local
integer                                                :: i,j,k
double precision, dimension(6)                         :: cellp, strderv, cellderv
double precision, allocatable, dimension(:,:)          :: POS_FRAC, FORCE_FRAC

if (.not. allocated(FORCE_FRAC)) allocate(FORCE_FRAC(NA, 3))
if (.not. allocated(POS_FRAC)) allocate(POS_FRAC(NA, 3))

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
CALL CELLDRV(transpose(LAT), cellp, strderv, cellderv)
!print*, 'strderv'
!print*, strderv
!print*, 'cellderv'
!print*, cellderv
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
fc = ENE + SUM(EXTSTRESS(1:3))/3.d0 * ABS(DET(LAT)) * cfactor

!gc = gc * -1.d0
!print *, 'fc',fc
!stop
END SUBROUTINE

SUBROUTINE  relaxv2struct(n, xc, NA, LAT, POS)
implicit none

INTEGER         , intent(in)                           :: n
double precision, intent(in),dimension(n)              :: xc
INTEGER         , intent(in)                           :: NA
double precision, intent(out),dimension(3,3)         :: LAT
double precision, intent(out),dimension(NA,3)        :: POS

! local
integer                                                :: i,j
double precision, dimension(6)                         :: cellp
double precision, allocatable, dimension(:,:)          :: POS_FRAC

if (.not. allocated(POS_FRAC)) allocate(POS_FRAC(NA, 3))
cellp(1:6) = xc(1:6)
CALL LAT2MATRIX(cellp, LAT, 1)
do i = 1, NA
    do j = 1, 3
        POS_FRAC(i,j) = xc(6 + 3*(i - 1) + j)
    enddo
enddo
CALL FRAC2CART(NA, LAT, POS_FRAC, POS)
END SUBROUTINE
    
SUBROUTINE  CART2FRAC(NA, LAT, POS, POS_FRAC)
implicit none

INTEGER         , intent(in)                           :: NA
double precision, intent(in),dimension(3,3)            :: LAT
double precision, intent(in),dimension(NA,3)           :: POS
double precision, intent(inout),dimension(NA,3)        :: POS_FRAC

! local
integer                                                :: i,j
double precision, dimension(3,3)                       :: INV_LAT
!call lat_inv(LAT, INV_LAT)
INV_LAT = inv_33(LAT)
POS_FRAC = matmul(POS, INV_LAT)
do i = 1, NA
    do j = 1,3
        if (POS_FRAC(i,j) < 0.d0) POS_FRAC(i,j) = POS_FRAC(i,j) + 1.d0
        if (POS_FRAC(i,j) > 1.d0) POS_FRAC(i,j) = POS_FRAC(i,j) - 1.d0
    enddo
enddo
END SUBROUTINE CART2FRAC

SUBROUTINE  FRAC2CART(NA, LAT, POS_FRAC, POS)
implicit none

INTEGER         , intent(in)                           :: NA
double precision, intent(in),dimension(3,3)            :: LAT
double precision, intent(inout),dimension(NA,3)        :: POS_FRAC
double precision, intent(inout),dimension(NA,3)        :: POS
! loca
integer                                                :: i,j
do i = 1, NA
    do j = 1,3
        if (POS_FRAC(i,j) < 0.d0) POS_FRAC(i,j) = POS_FRAC(i,j) + 1.d0
        if (POS_FRAC(i,j) > 1.d0) POS_FRAC(i,j) = POS_FRAC(i,j) - 1.d0
    enddo
enddo
POS = matmul(POS_FRAC, LAT)
END SUBROUTINE FRAC2CART
        

subroutine lat2matrix(lat,matrix,iflag)!{{{
implicit none

! if iflag==1, abc2matrix; iflag==2. matrix2abc
integer,          intent(in)        :: iflag 
double precision, intent(inout)     :: lat(6),matrix(3,3)

!local parameters
double precision                    :: ra,rb,rc,&
                                       cosinea, cosineb,cosinec,&
                                       anglea,angleb,anglec
double precision, parameter         :: radtodeg = 57.29577951d0
double precision, parameter         :: degtorad = 1.0/radtodeg

if (iflag==1) then
   lat(4:6) = lat(4:6) * degtorad
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
   lat(4)=anglea * radtodeg
   lat(5)=angleb * radtodeg
   lat(6)=anglec * radtodeg
endif
end subroutine lat2matrix!}}}

function gnorm(n,x)!{{{
implicit none

integer, intent(in)                        :: n
double precision, intent(in), dimension(n) :: x
double precision                           :: gnorm

! local 
integer                                    :: i
gnorm = 0.d0
do i = 1, n
    gnorm = gnorm + x(i)**2
enddo
gnorm = dsqrt(gnorm)/n
end function gnorm!}}}

SUBROUTINE WRITE_VASP(NA, SPECIES, LAT, POS, ENE, FORCE, STRESS, EXTSTRESS)
IMPLICIT NONE

integer, intent(in)                         :: NA
integer, intent(in), dimension(NA)          :: SPECIES
double precision,intent(in),dimension(3,3)  :: LAT
double precision,intent(in),dimension(NA,3) :: POS, FORCE
double precision,intent(in),dimension(6)    :: STRESS, EXTSTRESS
double precision,intent(in)                 :: ENE
! loca
integer                                     :: NS
integer, allocatable, dimension(:)          :: nele, ele_number
character,allocatable, dimension(:)         :: ele
double precision,allocatable,dimension(:,:) :: POS_FRAC
integer                                     :: i, j, k, i_temp, count
double precision                            :: press, enthalpy
double precision,parameter                  :: cfactor = 6.241460893d-3

i_temp = SPECIES(1)
NS = 1
do i = 2, NA
    if (SPECIES(i) .ne. i_temp) then
        NS = NS + 1
        i_temp = SPECIES(i)
    endif
enddo

if (.not. allocated(ele)) allocate(ele(ns))
if (.not. allocated(nele)) allocate(nele(ns))
if (.not. allocated(ele_number)) allocate(ele_number(ns))
if (.not. allocated(pos_frac)) allocate(pos_frac(na, 3))

ele_number(1) = SPECIES(1)
ele(1) = atsym(ele_number(1))
i_temp = SPECIES(1)
k = 1
count = 1
do i = 2, NA
    if (SPECIES(i) .ne. i_temp) then
        i_temp = SPECIES(i)
        k = k + 1
        ele_number(k) = i_temp
        ele(k) = atsym(ele_number(k))
    endif
enddo

nele = 0

do i = 1, NS
    do j = 1, NA
        if (SPECIES(j) == ele_number(i)) nele(i) = nele(i) + 1
    enddo
enddo

call CART2FRAC(NA, LAT, POS, POS_FRAC)
!@ write CONTCAR 
open(1123, file = 'CONTCAR')
write(1123,*) 'written by GAP'
write(1123,*) '1.0'
DO i = 1, 3
    write(1123,'(3F20.10)') LAT(i,:)
ENDDO
write(1123,'(3A5)') ele
write(1123,'(3I5)') nele
write(1123,'(A6)') 'Direct'
DO i = 1, NA
    write(1123,'(3F20.10)') POS_FRAC(i,:)
ENDDO
close(1123)

press = SUM(EXTSTRESS(1:3))/3.0
enthalpy = ene + press * ABS(DET(LAT)) * cfactor

open(212, file = 'OUTCAR')
do i = 1, ns
    write(212,*) 'VRHFIN ='//trim(ele(i))
enddo
write(212,*) 'ions per type', ' = ', nele
write(212,*) 'direct lattice vectors'
write(212,*) 'volume of cell :'
write(212,*) 'Direction    XX          YY          ZZ          XY          YZ          ZX'
write(212,'(A6, 6F15.6)') 'in kB', stress(1)*10, stress(4)*10, stress(6)*10, stress(2)*10, stress(5)*10, stress(3)*10
write(212,'(A20,F15.6, A5, A20, F15.6, A5)') 'external pressure =', (stress(1) + stress(4) + stress(6))*10.0/3.0 - press*10.0, &
'kB', 'Pullay stress =', press*10.0 ,'kB'
write(212,*) 'volume of cell :', abs(det(lat))
write(212,*) 'direct lattice vectors'
do i = 1,3
   write(212,'(6F15.6)') lat(i,:)
enddo
write(212,*) 'POSITION                                       TOTAL-FORCE(eV/Angst)'
write(212,*) '-------------------------------------------------------------------'
do i = 1, na
   write(212,'(6F15.6)') pos(i,:), force(i,:)
enddo
write(212,*) '-------------------------------------------------------------------'
write(212,'(A31, 2F15.6)') 'energy  without entropy=', ene, ene/na
write(212,'(A30, 2F15.6)') 'enthalpy is  TOTEN    =', enthalpy, enthalpy/na
close(212)

deallocate(ele, nele, ele_number, pos_frac)

END SUBROUTINE

SUBROUTINE WRITE_FAKE_VASP(NA, SPECIES)
IMPLICIT NONE

integer, intent(in)                         :: NA
integer, intent(in), dimension(NA)          :: SPECIES
! loca
integer                                     :: NS
integer, allocatable, dimension(:)          :: nele, ele_number
character,allocatable, dimension(:)         :: ele
integer                                     :: i, j, k, i_temp, count
double precision                            :: fake_cp, fake_fp, fake_e

fake_cp = 1.d0
fake_fp = 0.d0
fake_e = 9999.d0
i_temp = SPECIES(1)
NS = 1
do i = 2, NA
    if (SPECIES(i) .ne. i_temp) then
        NS = NS + 1
        i_temp = SPECIES(i)
    endif
enddo
if (.not. allocated(ele)) allocate(ele(ns))
if (.not. allocated(nele)) allocate(nele(ns))
if (.not. allocated(ele_number)) allocate(ele_number(ns))
ele_number(1) = SPECIES(1)
ele(1) = atsym(ele_number(1))
i_temp = SPECIES(1)
k = 1
count = 1
do i = 2, NA
    if (SPECIES(i) .ne. i_temp) then
        i_temp = SPECIES(i)
        k = k + 1
        ele_number(k) = i_temp
        ele(k) = atsym(ele_number(k))
    endif
enddo

nele = 0

do i = 1, NS
    do j = 1, NA
        if (SPECIES(j) == ele_number(i)) nele(i) = nele(i) + 1
    enddo
enddo
open(1123, file = 'CONTCAR')
write(1123,*) 'written by GAP'
write(1123,*) '1.0'
DO i = 1, 3
    write(1123,'(3F20.10)') fake_cp, fake_cp, fake_cp
ENDDO
write(1123,'(3A5)') ele
write(1123,'(3I5)') nele
write(1123,'(A9)') 'Cartesian'
DO i = 1, NA
    write(1123,'(3F20.10)') fake_fp, fake_fp, fake_fp
ENDDO
close(1123)


open(212, file = 'OUTCAR')
do i = 1, ns
    write(212,*) 'VRHFIN ='//trim(ele(i))
enddo
write(212,*) 'ions per type', ' = ', nele
write(212,*) 'direct lattice vectors'
write(212,*) 'volume of cell :'
write(212,*) 'Direction    XX          YY          ZZ          XY          YZ          ZX'
write(212,'(A6, 6F15.6)') 'in kB',  fake_fp, fake_fp, fake_fp, fake_fp, fake_fp, fake_fp
write(212,'(A20,F15.6, A5, A20, F15.6, A5)') 'external pressure =', fake_fp, &
'kB', 'Pullay stress =', fake_fp ,'kB'
write(212,*) 'volume of cell :', fake_fp
write(212,*) 'direct lattice vectors'
do i = 1,3
   write(212,'(6F15.6)') fake_cp, fake_cp, fake_cp
enddo
write(212,*) 'POSITION                                       TOTAL-FORCE(eV/Angst)'
write(212,*) '-------------------------------------------------------------------'
do i = 1, na
   write(212,'(6F15.6)') fake_cp, fake_cp, fake_cp, fake_cp * 10.0, fake_cp * 10.0, fake_cp * 10.0
enddo
write(212,*) '-------------------------------------------------------------------'
write(212,'(A31, F15.6)') 'energy  without entropy=', fake_e
write(212,'(A30, F15.6)') 'enthalpy is  TOTEN    =', fake_e
close(212)

deallocate(ele, nele, ele_number)
END SUBROUTINE

function inv_33(M)
double precision, intent(in) :: M(3,3)
double precision             :: Inv(3,3),inv_33(3,3)
double precision             :: d
d = Det(M)
inv(1,1)=M(2,2)*M(3,3) - M(2,3)*M(3,2)
inv(2,1)=M(2,3)*M(3,1) - M(2,1)*M(3,3)
inv(3,1)=M(2,1)*M(3,2) - M(2,2)*M(3,1)
inv(1,2)=M(1,3)*M(3,2) - M(1,2)*M(3,3)
inv(2,2)=M(1,1)*M(3,3) - M(1,3)*M(3,1)
inv(3,2)=M(1,2)*M(3,1) - M(1,1)*M(3,2)
inv(1,3)=M(1,2)*M(2,3) - M(1,3)*M(2,2)
inv(2,3)=M(1,3)*M(2,1) - M(1,1)*M(2,3)
inv(3,3)=M(1,1)*M(2,2) - M(1,2)*M(2,1)
inv_33=inv/d
end function

FUNCTION rmse_force(NA, FORCE)
implicit none

integer,         intent(in)                  :: NA
double precision,intent(in), dimension(NA,3) :: FORCE
double precision                             :: rmse_force
! local
integer                                      :: i,j
rmse_force = 0.d0
do i = 1, na
    do j = 1,3
        rmse_force = rmse_force + FORCE(i,j)**2/(3*NA)
    enddo
enddo
rmse_force = dsqrt(rmse_force)
END FUNCTION

END MODULE
