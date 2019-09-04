module GPR_BASE
use constants
use io
use math
use linearalgebra
use struct

!interface covariance
!    module procedure covariance_2B, covariance_MB
!end interface covariance
!interface INI_GAP
!    module procedure INI_GAP_2B, INI_GAP_MB
!end interface INI_GAP

real(dp), dimension(:,:), allocatable :: c_subYY_sqrtInverseLambda
real(dp), dimension(:,:), allocatable :: factor_c_subYsubY
real(dp), dimension(:,:), allocatable :: a
real(dp), dimension(:),   allocatable :: globalY
type(LA_Matrix)                       :: LA_c_subYsubY, LA_q_subYsubY
integer                               :: n_globalSparseX , n_globalY
integer                               :: error



contains
SUBROUTINE GPR(cmm, cmo, lamdaobe, alpha)
implicit none
real(dp),intent(in),dimension(:,:)      :: cmm
real(dp),intent(in),dimension(:,:)      :: cmo
real(dp),intent(in),dimension(:)        :: lamdaobe
real(dp),intent(out),dimension(:)       :: alpha
!--- local ---
integer                                 :: i,j
n_globalSparseX = size(cmo, 1)
n_globalY = size(cmo, 2)
!==========================================================================
allocate(factor_c_subYsubY(n_globalSparseX,n_globalSparseX))
allocate(a(n_globalY+n_globalSparseX,n_globalSparseX))
allocate(globalY(n_globalY+n_globalSparseX))

globalY = 0.d0
call LA_Matrix_initialise(LA_c_subYsubY,cmm)
call LA_Matrix_Factorise(LA_c_subYsubY,factor_c_subYsubY,error=error)
call LA_Matrix_finalise(LA_c_subYsubY)


do i = 1, n_globalSparseX-1
   do j = i+1, n_globalSparseX
      factor_c_subYsubY(j,i) = 0.0_qp
   end do
end do

a(1:n_globalY,:) = transpose(cmo(:,:))
a(n_globalY+1:,:) = factor_c_subYsubY
do i = 1, n_globalY
    globalY(i) = lamdaobe(i)
enddo
call LA_matrix_initialise(LA_q_subYsubY,a)
call LA_Matrix_QR_Solve_Vector(LA_q_subYsubY,globalY,alpha)
call LA_matrix_finalise(LA_q_subYsubY)
deallocate(factor_c_subYsubY)
deallocate(a)
deallocate(globalY)
END subroutine

FUNCTION  covariance_2B(delta, x, y, theta)
implicit none
REAL(DP),intent(in)  :: delta
real(DP),intent(in)  :: x
real(DP),intent(in)  :: y
REAL(DP),intent(in)  :: theta
real(DP)             :: covariance_2B

!integer  i 
REAL(DP)             :: fc_i, fc_j
fc_i = fcut_ij(x)
fc_j = fcut_ij(y)
covariance_2B = 0.d0
covariance_2B = covariance_2B + ((x-y)/theta)**2
covariance_2B = delta**2*exp(-0.5d0*covariance_2B) * fc_i * fc_j
END FUNCTION covariance_2B

FUNCTION  covariance_MB(delta, x,y, theta)
implicit none
real(DP),intent(in)  :: delta
real(DP),intent(in)  :: x(:)
real(DP),intent(in)  :: y(:)
real(DP),intent(in)  :: theta(:)
real(DP)             :: covariance_MB

integer  i 
!REAL(DP)            :: fc_i, fc_j
covariance_MB = 0.d0
do i = 1, size(x)
    covariance_MB = covariance_MB + ((x(i)-y(i))/theta(i))**2
enddo
covariance_MB = delta**2*exp(-0.5d0*covariance_MB) 
END FUNCTION covariance_MB

FUNCTION  dcovdx_2B(delta, x, y, theta)
implicit none
REAL(DP),intent(in)  :: delta
real(DP),intent(in)  :: x
real(DP),intent(in)  :: y
REAL(DP),intent(in)  :: theta
real(DP)             :: dcovdx_2B
!integer  i 
REAL(DP)            :: fc_i, fc_j, dfc_i, exp_part

dcovdx_2B = 0.d0
exp_part = 0.d0
fc_i = fcut_ij(x)
fc_j = fcut_ij(y)
dfc_i = dfcut_ij(x)
exp_part = exp_part + ((x-y)/theta)**2
exp_part = delta**2 * exp(-0.5d0 * exp_part)
dcovdx_2B = exp_part * -1.d0 * (x-y)/theta**2
dcovdx_2B = dcovdx_2B * fc_i + exp_part * dfc_i
dcovdx_2B = dcovdx_2B * fc_j
END FUNCTION dcovdx_2B

subroutine matmuldiag(x,y)
real(8),intent(in)   :: x(:)
real(8),intent(inout) :: y(:,:)
do i = 1,size(x)
    y(i,:) = x(i)*y(i,:)
enddo
end subroutine matmuldiag

subroutine matmuldiag_T(y,x)
real(8),intent(in)   :: x(:)
real(8),intent(inout) :: y(:,:)
do i = 1,size(x)
    y(:,i) = x(i)*y(:,i)
enddo
end subroutine matmuldiag_T

SUBROUTINE GAP_SET_THETA(x, theta)
implicit none

REAL(DP),intent(in),dimension(:,:)   :: x
REAL(DP),intent(inout),dimension(:)  :: theta
! local
integer                              :: i,j,n,m
REAL(DP),allocatable,dimension(:)    :: t

theta = 0.d0
n = size(x,1)
m = size(x,2)
allocate(t(n))
do i = 1, m
    do j = 1, n
        t(j) = x(j,i)
    enddo
print*, t
    theta(i) = my_cov(t)
print *, theta(i)
stop
enddo
deallocate(t)
END SUBROUTINE

FUNCTION  MY_COV(t)
implicit none
real(DP)                          :: my_cov
real(DP),intent(in),dimension(:)  :: t
real(DP),allocatable,dimension(:) :: t2
integer                           :: i,n
n = size(t)
allocate(t2(n))
do i = 1,n
    t2(i) = t(i)*t(i)
enddo
my_cov = sum(t2)/n - (sum(t)/n)**2
my_cov = my_cov + 1.0d0
deallocate(t2)
END FUNCTION

SUBROUTINE GET_COV(delta, xx, yy, theta, cov)
integer i,j,k
REAL(DP),intent(in)                 :: delta
real(DP),intent(in),dimension(:,:)  :: xx, yy
REAL(DP),intent(in),dimension(:)    :: theta
real(DP),intent(out),dimension(:,:) :: cov
!local

real(DP),allocatable,dimension(:,:) :: temp
integer                             :: nf
n = size(xx,1)
m = size(yy,1)
allocate(temp(n,m))
nf = size(xx,2)
temp = 0.d0

do i = 1,n
    do  j=1,m 
        do k = 1,nf
            temp(i,j) = temp(i,j)+((xx(i,k)-yy(j,k))/theta(k))**2
        enddo
        cov(i,j) = delta * exp(-0.5d0*temp(i,j))
    enddo
enddo
deallocate(temp)
END SUBROUTINE

  subroutine cur_decomposition(this, index_out, rank, n_iter)
    ! based on 10.1073/pnas.0803205106

    real(dp), intent(in), dimension(:,:) :: this
    integer, dimension(:), intent(out) :: index_out
    integer, intent(in), optional :: rank, n_iter

    integer :: n
    integer :: expected_columns
    integer :: my_n_iter, my_rank
    type(LA_Matrix) :: LA_this
    real(dp), allocatable, dimension(:) :: p, s, p_minus_ran_uniform
    real(dp), allocatable, dimension(:,:) :: v
    integer :: j, l
    integer, allocatable, dimension(:), target :: p_index
    integer, pointer, dimension(:) :: tmp_index_out => null()
    real(dp), allocatable, dimension(:,:) :: C, Cp
    real(dp) :: err, min_err
    integer :: error
    real(dp) :: temp_random

    expected_columns = size(index_out)

    if( expected_columns <= 0 ) then
!       call print_warning("cur_decomposition: called with expected_columns "//expected_columns//", can't be zero or less")
       return
    endif

    call LA_Matrix_initialise(LA_this,this)

    my_n_iter = optional_default(1, n_iter)

    if (present(rank)) then
       call LA_Matrix_SVD_Allocate(LA_this,v=v,error=error)
!       HANDLE_ERROR(error)
       call LA_Matrix_SVD(LA_this,v=v,error=error)
!       HANDLE_ERROR(error)
       my_rank = rank
    else
       call LA_Matrix_SVD_Allocate(LA_this,s=s,v=v,error=error)
!       HANDLE_ERROR(error)
       call LA_Matrix_SVD(LA_this,s=s,v=v,error=error)
!       HANDLE_ERROR(error)
       my_rank = count(s > TOL_SVD) / 2
    endif

    n = size(v,1)
    allocate(p(n), p_minus_ran_uniform(n), p_index(n))
    allocate( C(size(this,1),expected_columns), Cp(expected_columns,size(this,1)) )

    p = sum(v(:,1:my_rank)**2, dim=2)
    p = p * expected_columns
    p = p / my_rank
    p = min(p,1.0_dp)

    if(my_n_iter <= 0) then ! do not do probabilistic selection of columns
       p_index = (/(j, j=1,n )/)
       p_minus_ran_uniform = -p
       call heap_sort(p_minus_ran_uniform,i_data=p_index)
       index_out = p_index(1:expected_columns)
    else
       min_err = huge(1.0_dp)
       do l = 1, my_n_iter

          ! randomly select columns according to the probabilities
          do j = 1, n
             temp_random = ran_uniform()
!             print*, 'random', temp_random
             p_minus_ran_uniform(j) = temp_random - p(j)
             p_index(j) = j ! initialise index array
          end do

          call heap_sort(p_minus_ran_uniform,i_data=p_index)
          tmp_index_out => p_index(1:expected_columns)

          C = this(:,tmp_index_out)
          ! pinv: Moore-Penrose pseudo-inverse
          call pseudo_inverse(C,Cp)
          err = sum( (this - ( C .mult. Cp .mult. this))**2 )

!          call print("cur_decomposition: iteration: "//l//", error: "//err)
          if(err < min_err) then        ! this happens at least once
             index_out = tmp_index_out
             min_err = err
          endif

       end do
    endif

    call LA_Matrix_finalise(LA_this)

    tmp_index_out => null()
    if(allocated(s)) deallocate(s)
    if(allocated(v)) deallocate(v)
    if(allocated(p)) deallocate(p)
    if(allocated(p_minus_ran_uniform)) deallocate(p_minus_ran_uniform)
    if(allocated(p_index)) deallocate(p_index)
    if(allocated(C)) deallocate(C)
    if(allocated(Cp)) deallocate(Cp)

  end subroutine cur_decomposition

SUBROUTINE SET_GPR_OB(AT, DATA_C)
type(Structure),intent(inout),dimension(:)  :: at
type(DATA_TYPE),intent(inout)               :: data_c

n_config = size(at)
kf = 0
do i = 1, n_config
    kf = kf + 1

    if (data_c%ltrain_2b .and. .not. data_c%ltrain_mb) then
        data_c%obe(i) = at(i)%energy_ref - at(i)%natoms * ene_cons
    endif

    if (data_c%ltrain_mb .and. .not. data_c%ltrain_2b) then
        data_c%ob(kf) = at(i)%energy_ref
        do j = 1, at(i)%natoms
            do k1 = 1, 3
                kf = kf + 1
                data_c%ob(kf) = at(i)%force_ref(j,k1)
            enddo
        enddo
        do j = 1, 6
            kf = kf + 1
            data_c%ob(kf) = at(i)%stress_ref(j)
        enddo
    endif
enddo
print*, 'Number of structures/energy/stress:', data_c%ne, data_c%ne, data_c%ns
print*, 'Number of atoms/forces:', data_c%natoms, data_c%nf
print*, 'Number of Observable variables:',data_c%nob
ENDSUBROUTINE SET_GPR_OB

SUBROUTINE SET_GPR_OB_MB(AT, DATA_C)
type(Structure),intent(in),dimension(:)     :: at
type(DATA_TYPE),intent(inout)               :: data_c

n_config = size(at)
if (.not. allocated(data_c%ob)) allocate(data_c%ob(data_c%nob))
kf = 0
do i = 1, n_config
    kf = kf + 1
    data_c%ob(kf) = at(i)%energy_ref - at(i)%energy_cal_2b
    do j = 1, at(i)%natoms
        do k1 = 1, 3
            kf = kf + 1
            data_c%ob(kf) = at(i)%force_ref(j,k1) - at(i)%force_cal_2b(j,k1)
        enddo
    enddo
    do j = 1, 6
        kf = kf + 1
        data_c%ob(kf) = at(i)%stress_ref(j) - at(i)%stress_cal_2b(j)
    enddo
enddo
ENDSUBROUTINE SET_GPR_OB_MB

SUBROUTINE SET_GPR_OB_2B(AT, DATA_C)
type(Structure),intent(in),dimension(:)     :: at
type(DATA_TYPE),intent(inout)               :: data_c

n_config = size(at)
data_c%ne = n_config
if (.not. allocated(data_c%obe)) allocate(data_c%obe(data_c%ne))
kf = 0
do i = 1, n_config
    data_c%obe(i) = at(i)%energy_ref - at(i)%energy_cal_mb - at(i)%natoms * ene_cons
enddo
ENDSUBROUTINE SET_GPR_OB_2B

END module
