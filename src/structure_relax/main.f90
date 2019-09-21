Program  test_relax
use read_module
use gap_module
use relax_module


INTEGER                                         :: NA
INTEGER, ALLOCATABLE, DIMENSION(NA)             :: SPECIES
double precision, ALLOCATABLE, DIMENSION(3,3)   :: LAT
double precision, ALLOCATABLE, DIMENSION(:,:)   :: POS
double precision, dimension(6)                  :: EXTSTRESS
data extstress/100.0, 100.0, 100.0, 0.0, 0.0, 0.0/
call read_gap(NA, SPECIES, LAT, POS)
call relax_main(NA, SPECIES, LAT, POS, EXTSTRESS)
end subroutine
