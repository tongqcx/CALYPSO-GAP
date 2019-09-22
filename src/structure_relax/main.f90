Program  test_relax
use read_module
use gap_module
use relax_module


double precision, dimension(6)                  :: EXTSTRESS
data extstress/100.0, 100.0, 100.0, 0.0, 0.0, 0.0/
call read_gap()
!print *, SPECIES
!print *, transpose(lat)
!print *, transpose(pos)
call relax_main(NA, SPECIES, LAT, POS, EXTSTRESS)
end program
