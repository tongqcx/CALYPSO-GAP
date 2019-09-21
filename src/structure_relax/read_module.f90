module read_module
implicit none

contains
SUBROUTINE  READ_POSCAR(NA, SPECIES, LAT, POS)
INTEGER, intent(inout)                                    :: NA
INTEGER, intent(inout), DIMENSION(:)                      :: SPECIES
double precision, intent(inout), DIMENSION(3,3)           :: LAT
double precision, intent(inout), DIMENSION(:,:)           :: POS

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

SUBROUTINE READ_GAP(NA, SPECIES, LAT, POS)
INTEGER, intent(inout)                                    :: NA
INTEGER, intent(inout), DIMENSION(:)                      :: SPECIES
double precision, intent(inout), DIMENSION(3,3)           :: LAT
double precision, intent(inout), DIMENSION(:,:)           :: POS

! local
integer                                                   :: i,j
integer, parameter                                        :: maxele = 107
character(2)                                              :: ctemp
character(len=2), save                                    :: atsym(maxele)
integer                                                   :: ns
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
    do j = 1, size(asym)
        if (ctemp == asym(j)) SPECIES(i) = j
    enddo
enddo
READ(2233,*)
END SUBROUTINE
end module
