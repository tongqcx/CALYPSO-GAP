Program  test_relax
use read_module
use gap_module
use relax_module


call read_gulp()
!call read_gap()
!print *, SPECIES
!print *, transpose(lat)
!print *, transpose(pos)
call relax_main_conj(NA, SPECIES, LAT, POS, STRESS, maxcycle)
end program
