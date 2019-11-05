module math
use gap_constants

interface  write_array
    module procedure write_array_2dim, write_array_1dim
end interface write_array

contains

function vectorlength(vc)
    real(8) :: vc(3),vectorlength
    vectorlength=sqrt(vc(1)**2+vc(2)**2+vc(3)**2)
end function

function recipvector(lat)
    real(8),intent(in) :: lat(:,:)
    real(8) :: recipvector(3,3)

    recipvector(:,1)=crossp(lat(:,2),lat(:,3))
    recipvector(:,2)=crossp(lat(:,3),lat(:,1))
    recipvector(:,3)=crossp(lat(:,1),lat(:,2))
    recipvector=recipvector/volume(lat)*pi*2._8

end function

function volume(lat)
    real(8),intent(in) :: lat(:,:)
    real(8) :: volume 

    volume=abs(sum(lat(:,1)*crossp(lat(:,2),lat(:,3))))

end function

function crossp(va,vb)
    real(8),intent(in) :: va(3),vb(3)
    real(8) :: crossp(3)

    crossp(1)=va(2)*vb(3)-va(3)*vb(2)
    crossp(2)=va(3)*vb(1)-va(1)*vb(3)
    crossp(3)=va(1)*vb(2)-va(2)*vb(1)
end function

!function MIV_num_columns( inchar ) result(num)
!IMPLICIT none
!character(*)                            :: inchar
!integer(i4b)                            :: num
!integer(i4b)                            :: length
!character(:),allocatable                :: ctp
!integer(i4b)                            :: i
!
!ctp = trim(adjustl(inchar))
!length = len(ctp)
!num = 1
!do i = 2,length-1
!    if( ctp(i:i) == ' ' .and. ctp(i+1:i+1) /= ' ' )  num = num+1
!end do
!return 
!END FUNCTION MIV_NUM_COLUMNS

FUNCTION fcut_ij(r)
real(DP)   :: fcut_ij, r
if (r <= DATA_C%rcut - DATA_C%d_width ) then
    fcut_ij = 1.d0
elseif ( r > DATA_C%rcut) then
    fcut_ij = 0.d0
else
    fcut_ij = ((cos(pi*(r - DATA_C%rcut + DATA_C%d_width)/DATA_C%d_width) + 1))/2
endif
return
END FUNCTION fcut_ij

FUNCTION dfcut_ij(r)
real(DP)   :: dfcut_ij, r
if (r <= DATA_C%rcut - DATA_C%d_width ) then
    dfcut_ij = 0.d0
elseif ( r > DATA_C%rcut) then
    dfcut_ij = 0.d0
else
    dfcut_ij = -0.5d0 * pi * sin(pi*(r - DATA_C%rcut + DATA_C%d_width)/DATA_C%d_width) / DATA_C%d_width 
endif
return
END FUNCTION dfcut_ij

SUBROUTINE  write_array_2dim(a,name)
REAL(DP),intent(in),dimension(:,:)       :: a
character(*),intent(in)                :: name
integer                                  :: i,j
open(2244,file=trim(adjustl(name)))
do i = 1, size(a,1)
    do j = 1,size(a,2)
        write(2244,'(F20.10,$)') a(i,j)
    enddo
    write(2244,*)
enddo
close(2244)
END SUBROUTINE

SUBROUTINE  write_array_1dim(a,name)
REAL(DP),intent(in),dimension(:)       :: a
character(*),intent(in)                :: name
integer                                  :: i,j
open(2244,file=trim(adjustl(name)))
do i = 1, size(a)
        write(2244,'(F20.10)') a(i)
enddo
close(2244)
END SUBROUTINE

FUNCTION length(x,theta)
    real(8)    :: length
    real(8)    :: x(:)
    real(8)    :: theta(:)
    integer    :: i
    length = 0.d0
    do i = 1,size(x)
        length = length + (x(i)/theta(i))**2
    enddo
    length = sqrt(length)
END FUNCTION

FUNCTION MY_DET(X)
    REAL(DP)                 :: MY_DET
    integer,allocatable      ::  ipiv(:)
    REAL(DP)                 :: X(:,:)
    REAL(DP)                 :: RES
    INTEGER                  :: I,J,N
    N = SIZE(X,1)
    allocate(ipiv(n))
    call dgetrf(n,n,x,n,ipiv,info)
    res = 1.d0
    do i=1,n
        if (ipiv(i).ne.i)  then
            res = -1.d0*res*x(i,i)
        else
            res = 1.d0*res*x(i,i)
        endif
    enddo
    deallocate(ipiv)
    MY_DET = RES
END FUNCTION

SUBROUTINE MY_INVERSE(n,x,ix)
    integer                 :: n,lda,info,lwork
    real(8)                 :: x(n,n),ix(n,n)
    integer,allocatable     :: ipiv(:)
    real(8),allocatable     :: work(:)
    real(8),allocatable     :: e(:,:)
    real(8),allocatable     :: mat_x_ix(:,:)
    logical                 :: icheck
    allocate(ipiv(n))
    allocate(e(n,n))
    do  i =1,n
        do j =1,n
            if (i==j) then
                e(i,j) = 1.d0
            else
                e(i,j) = 0.d0
            endif
        enddo
    enddo
    ix = x
    lwork = n
    lda = n
    allocate(work(lwork))
    print*,"INV BEGIN",n
    !call get_time(initime)
    call dgetrf(n,n,ix,lda,ipiv,info)
    if ( info/=0 )  then
       print*, "Error in dgetrf!"
       write(14,*) "Error in dgetrf!"
       stop
    endif
    !call dgetri(n,ix,lda,ipiv,work,lwork,info)
    call dgetrs('N',n,n,ix,lda,ipiv,e,n,info)
    if ( info/=0 )  then
       print*, "Error in dgetrf!"
       write(14,*) "Error in dgetrf!"
       stop
    endif
    !call get_time(fintime)
    print*, "INV FINISHED" 
    
    deallocate(ipiv)
    deallocate(work)
    ix = e
    deallocate(e)
    allocate(mat_x_ix(n,n))
    call matmul_real(x, ix, 'N', 'N', mat_x_ix)
    call check_matrix(mat_x_ix)
    !if (abs(matrix_mean(mat_x_ix))-1.d0) > 1.0E-1)  then
    !    print*, "inv ERROR"
    !    write(14,*) "inv ERROR"
    !    stop
    !endif
    deallocate(mat_x_ix)
END SUBROUTINE

SUBROUTINE check_matrix(x)
    real(8),intent(in),dimension(:,:)      :: x
    !logical                                :: check_matrix
    !local
    integer                                :: i,j, n,m, icount
    real(8)                                :: error, iprecision, max_error
    icount = 0
    iprecision = 1.d-4
    max_error = 0.d0
    n = size(x, 1)
    m = size(x, 2)
    do i = 1, n
        do j = 1, m
            if (i == j) then
                error = abs(x(i,j) - 1.d0)
                max_error = max(max_error, error)
                if (error > iprecision)  icount = icount + 1
            else
                error = abs(x(i,j) - 0.d0)
                max_error = max(max_error, error)
                if (error > iprecision)  icount = icount + 1
            endif
        enddo
    enddo
    if (icount > 0) then
        write(*,'(A48,F10.5,I10)'), 'THERE IS AN UNEXPECTED ERROR IN MATRIX INVERSION', max_error, icount
    !    check_matrix = .FALSE.
    endif
END SUBROUTINE
!######################################################
!     Author: QiangXu  ATLAS2.0
!     Date:   2018.04.14  10:18
!######################################################
SUBROUTINE matmul_real(matA,matB,opA,opB,matC)
   IMPLICIT NONE
   REAL(8),INTENT(IN)  :: matA(:,:),matB(:,:)
   REAL(8),INTENT(OUT) :: matC(:,:)
   CHARACTER(1),INTENT(IN) :: opA,opB
   !
   INTEGER(4) :: LDA,LDB,LDC,M,N,K
   REAL(8)  :: alpha=1.d0  &
           &    ,  beta=0.d0
   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !C
   matC(:,:)=0.d0
   LDC=SIZE(matC,1)
   M=LDC
   N=SIZE(matC,2)
   !A
   IF(opA=='N'.OR.opA=='n')THEN
      K=SIZE(matA,2)
      LDA=MAX(1,M)
   ELSE
      K=SIZE(matA,1)
      LDA=MAX(1,K)
   ENDIF
   !B
   IF(opB=='N'.OR.opB=='n')THEN
      LDB=MAX(1,K)
   ELSE 
      LDB=MAX(1,N)
   ENDIF
   !
   CALL DGEMM(opA,opB,M,N,K,alpha,matA,LDA,matB,LDB,beta,matC,LDC)
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ENDSUBROUTINE matmul_real

end module math
