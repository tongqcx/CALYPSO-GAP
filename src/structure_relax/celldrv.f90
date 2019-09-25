  subroutine celldrv(lat, cellp, strderv,cellderv)
!
!  Converts the derivatives in terms of the strains into 
!  cell parameter derivatives.
!
!  11/07 Unused variables cleaned up
!   4/08 lgrad2 removed from arguments as this is unused
!
!  Conditions of use:
!
!  GULP is available free of charge to academic institutions
!  and non-commerical establishments only. Copies should be
!  obtained from the author only and should not be distributed
!  in any form by the user to a third party without the express
!  permission of the author. This notice applies to all parts
!  of the program, except any library routines which are
!  distributed with the code for completeness. All rights for
!  such routines remain with the original distributor.
!
!  No claim is made that this program is free from errors and
!  no liability will be accepted for any loss or damage that
!  may result. The user is responsible for checking the validity
!  of their results.
!
!  Copyright Curtin University 2008
!
!  Julian Gale, NRI, Curtin University, April 2008
!
!  use current
!  use derivatives, only : cderv
  implicit none
  integer, parameter :: dp  = kind(1.0d0)
  integer, parameter :: i4  = selected_int_kind(5)
  
!
!  Passed variables
  real(dp),    intent(in)  :: lat(3,3)
  real(dp),    intent(in)  :: cellp(6)
  real(dp),    intent(out) :: cellderv(6)
  real(dp),    intent(in)  :: strderv(6)
!
!  Local variables
!
  integer(i4)              :: i
  integer(i4)              :: j
  real(dp)                 :: cderv(6,6)                       ! Second derivatives - cell-parameter - cell-parameter
!
!  Initialise cell derivative vector
!
  do i = 1,6
    cellderv(i) = 0.0_dp
  enddo
  call setcellderv3D(lat, cellp, cderv, .false.,.false.,.true.)
!
!  Calculate first derivatives of energy with respect to cell parameters
!
!print*, 'celldrv: strderv', strderv
!stop
  do i = 1,6
      do j = 1,6
        cellderv(j) = cellderv(j) + strderv(i)*cderv(j,i)
      enddo
  enddo
!
  return
  end
