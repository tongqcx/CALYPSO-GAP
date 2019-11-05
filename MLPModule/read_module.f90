module read_module
implicit none
INTEGER                                         :: NA
INTEGER, ALLOCATABLE, DIMENSION(:)              :: SPECIES
double precision, DIMENSION(3,3)                :: LAT
double precision, DIMENSION(6)                  :: STRESS
double precision, ALLOCATABLE, DIMENSION(:,:)   :: POS
double precision                                :: pressure
double precision                                :: gtol,ftol
integer                                         :: maxcycle
integer, parameter                                        :: maxele = 107
character(len=2), save                                    :: atsym(maxele)
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

private   ::  maxele
contains
SUBROUTINE  READ_POSCAR()
!INTEGER, intent(inout)                                    :: NA
!INTEGER, intent(inout), DIMENSION(:)                      :: SPECIES
!double precision, intent(inout), DIMENSION(3,3)           :: LAT
!double precision, intent(inout), DIMENSION(:,:)           :: POS

! local
integer                                                   :: i
character(200)                                            :: c_temp
double precision                                          :: lat_scale

open(2233, file = 'POSCAR')
read(2233,*) lat_scale
do i = 1, 3
    read(2233,*) LAT(i,:)
enddo
read(2233,*)  c_temp
read(2233,*)  
do i = 1, NA
    read(2233,*) POS(i,:)
enddo
END SUBROUTINE READ_POSCAR

SUBROUTINE READ_GAP()

! local
integer                                                   :: i,j
character(2)                                              :: ctemp
integer                                                   :: ns
double precision                                          :: pressure

pressure = 100.0
open(2233, file = 'config')
READ(2233, *) 
READ(2233, *) NA, NS
if (.not. allocated(POS))     allocate(POS(NA,3))
if (.not. allocated(SPECIES)) allocate(SPECIES(NA))
do i = 1, 3
    read(2233,*) LAT(i,:)
enddo
READ(2233,*)
do i = 1, NA
    read(2233, *) ctemp, POS(i,:)
    do j = 1, size(atsym)
        if (ctemp == atsym(j)) SPECIES(i) = j
    enddo
enddo
READ(2233,*)
STRESS = 0.d0
do i = 1, 3
    STRESS(i) = pressure
enddo
END SUBROUTINE

SUBROUTINE READ_GULP()
integer                                                   :: i,j
character(2)                                              :: ctemp
integer                                                   :: ns
character(100)                                            :: line
character(4)                                              :: core
integer                                                   :: ios

open(2233, file = 'gulpinput')
READ(2233,*)
READ(2233,*)
do i = 1, 3
    read(2233,*) LAT(i,:)
enddo
READ(2233,*)
NA = 0
DO while (.true.)
    read(2233,*) line
    if (index(line, 'space') .ne. 0) exit
    NA = NA + 1
ENDDO
if (.not. allocated(POS))     allocate(POS(NA,3))
if (.not. allocated(SPECIES)) allocate(SPECIES(NA))
REWIND(2233)
do i = 1, 6
    read(2233,*)
enddo

do i = 1, NA
    read(2233, *) ctemp, core, POS(i,:)
    do j = 1, size(atsym)
        if (ctemp == atsym(j)) SPECIES(i) = j
    enddo
enddo

! default value
pressure = 0.d0
maxcycle = 50
gtol = 0.005d0
ftol = 0.0005d0

DO while (.true.)
    read(2233,'(A100)', iostat=ios) line
    if (ios .ne. 0) exit
    if (index(line, 'pressure') .ne. 0) read(line,*)  ctemp, pressure
    if (index(line, 'maxcycle') .ne. 0) read(line,*)  ctemp, maxcycle
    if (index(line, 'gtol') .ne. 0) read(line,*)  ctemp, gtol
    if (index(line, 'ftol') .ne. 0) read(line,*)  ctemp, ftol
ENDDO
!print*, 'pressure', pressure
!print*, 'maxcycle', maxcycle
!print*, transpose(lat)
!print*, transpose(pos)
!print*, species
STRESS = 0.d0
do i = 1, 3
    STRESS(i) = pressure
enddo
!print*, STRESS
END SUBROUTINE

end module
