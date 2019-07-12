module constants
INTEGER,  PARAMETER      :: DP = KIND(1.0D0)
INTEGER,  PARAMETER      :: QP = KIND(1.0D0)
REAL(DP), PARAMETER      :: pi=3.141592653589793238462643383279502884197_DP
REAL(DP), PARAMETER      :: ene_cons = -6.0d0
!  *************  Parameters of GPR
integer                                :: nsparse
integer                                :: nspecies, ninteraction
REAL(DP)                               :: theta, delta, d_width, sigma_jitter
REAL(DP)                               :: sigma_e, sigma_f, sigma_s
character(2),allocatable,dimension(:)  :: elements
REAL(DP)                               :: Rcut, Rmin
REAL(DP)                               :: RMSE_ENERGY, RMSE_FORCE, RMSE_STRESS
INTEGER                                :: nforce
end module constants
