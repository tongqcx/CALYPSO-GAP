! this module perform structure optimization task
! 2019.09.19
module relax_module
implicit none



private   struct2relaxv, relaxv2struct
contains
SUBROUTINE  relax_main(NA, SPECIES, LAT, POS)!{{{
implicit none
double precision, intent(in)                  :: NA
double precision, intent(in),dimension(NA)    :: SPECIES
double precision, intent(in),dimension(3,3)   :: LAT
double precision, intent(in),dimension(Na,3)  :: POS

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
m = 
iprint = 1
if (.not. allocate(FORCE))  allocate(FORCE(NA, 3))
allocate ( nbd(n), x(n), l(n), u(n), g(n) )
allocate ( iwa(3*n) )
allocate ( wa(2*m*n + 5*n + 11*m*m + 8*m) )

task = 'START'
call  FGAP_CALC(NA, SPECIES, LAT, POS, ENE, FOECE, STRESS, VARIANCE)

! begin lbfgs loop
do while(task(1:2).eq.'FG'.or.task.eq.'NEW_X'.or.task.eq.'START') 
 
    call struct2relaxv(NA, LAT, POS, ENE, FORCE, STRESS, x, f, g)
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

SUBROUTINE  struct2relaxv(NA, LAT, POS, ENE, FORCE, STRESS, n, x, f, g)!{{{
implicit none

INTEGER         , intent(in)                           :: NA
double precision, intent(in),dimension(3,3)            :: LAT
double precision, intent(in),dimension(NA,3)           :: POS, FORCE
double precision, intent(in)                           :: ENE
double precision, intent(in),dimension(6)              :: STRESS
INTEGER         , intent(in)                           :: n
double precision, intent(inout),dimension(n)           :: x, g
double precision, intent(inout) i                      :: f
END SUBROUTINE!}}}

SUBROUTINE  relaxv2struct(n, x, NA, LAT, POS)!{{{
implicit none

INTEGER         , intent(in)                           :: n
double precision, intent(in),dimension(n)              :: x
INTEGER         , intent(in)                           :: NA
double precision, intent(inout),dimension(3,3)         :: LAT
double precision, intent(inout),dimension(NA,3)        :: POS
END SUBROUTINE!}}}
    
END MODULE
