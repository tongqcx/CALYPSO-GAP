! H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! H0 X
! H0 X   libAtoms+QUIP: atomistic simulation library
! H0 X
! H0 X   Portions of this code were written by
! H0 X     Albert Bartok-Partay, Silvia Cereda, Gabor Csanyi, James Kermode,
! H0 X     Ivan Solt, Wojciech Szlachta, Csilla Varnai, Steven Winfield.
! H0 X
! H0 X   Copyright 2006-2010.
! H0 X
! H0 X   These portions of the source code are released under the GNU General
! H0 X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
! H0 X
! H0 X   If you would like to license the source code under different terms,
! H0 X   please contact Gabor Csanyi, gabor@csanyi.net
! H0 X
! H0 X   Portions of this code were written by Noam Bernstein as part of
! H0 X   his employment for the U.S. Government, and are not subject
! H0 X   to copyright in the USA.
! H0 X
! H0 X
! H0 X   When using this software, please cite the following reference:
! H0 X
! H0 X   http://www.libatoms.org
! H0 X
! H0 X  Additional contributions by
! H0 X    Alessio Comisso, Chiara Gattinoni, and Gianpietro Moras, Qunchao Tong
! H0 X
! H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
module linearalgebra
use constants
  integer, parameter :: NOT_FACTORISED = 0
  integer, parameter :: CHOLESKY       = 1
  integer, parameter :: QR             = 2
  real(dp), public, parameter :: TOL_SVD = 1e-13_dp
  integer, parameter, private      :: SYSTEM_STRING_LENGTH = 1024 !max line length read
  integer, parameter, private      :: SYSTEM_STRING_LENGTH_LONG = 102400 !max line length read
  logical, parameter, private      :: use_intrinsic_blas = .false. 
  interface check_size
     module procedure check_size_int_dim1, check_size_int_dim1_s, check_size_int_dim2  
     module procedure check_size_real_dim1, check_size_real_dim1_s, check_size_real_dim2, check_size_real_dim3
     module procedure check_size_complex_dim1, check_size_complex_dim1_s, check_size_complex_dim2
     module procedure check_size_log_dim1, check_size_log_dim1_s, check_size_log_dim2  
  end interface check_size

  interface heap_sort
    module procedure heap_sort_i, heap_sort_r, heap_sort_r_2dim, heap_sort_i_2dim
  end interface heap_sort

  interface operator(.feq.)
!     module procedure real_feq,complex_feq,matrix_feq,vector_feq
     module procedure real_feq,complex_feq
  end interface

  interface operator(.multd.)
     module procedure matrix_product_vect_asdiagonal_dd, matrix_product_vect_asdiagonal_zz
  end interface

  interface operator(.mult.)
     module procedure  matrix_product_vect, matrix_product_ddd, matrix_product_int_vect, matrix_product_zzz
!     module procedure matrix_product_int_mat
  end interface
  interface optional_default
    module procedure optional_default_l, optional_default_i, optional_default_r
    module procedure optional_default_c, optional_default_ca, optional_default_z
    module procedure optional_default_ia, optional_default_ra
  end interface 

  interface matrix_product_vect_asdiagonal_sub
    module procedure matrix_product_vect_asdiagonal_sub_ddd
    module procedure matrix_product_vect_asdiagonal_sub_zzd
    module procedure matrix_product_vect_asdiagonal_sub_zdz
    module procedure matrix_product_vect_asdiagonal_sub_zzz
    module procedure vect_asdiagonal_product_matrix_sub_ddd
    module procedure vect_asdiagonal_product_matrix_sub_zzd
    module procedure vect_asdiagonal_product_matrix_sub_zdz
    module procedure vect_asdiagonal_product_matrix_sub_zzz
  end interface matrix_product_vect_asdiagonal_sub

  interface matrix_product_sub
    module procedure matrix_product_sub_ddd, matrix_product_sub_zzz, matrix_vector_product_sub_ddd
  end interface

  type LA_Matrix
     real(qp), dimension(:,:), allocatable :: matrix, factor
     real(qp), dimension(:), allocatable :: s, tau
     integer :: n, m
     logical :: initialised = .false.
     logical :: equilibrated = .false.
     integer :: factorised = NOT_FACTORISED
  endtype LA_Matrix
  contains
  subroutine LA_Matrix_Initialise(this,matrix)

     type(LA_Matrix), intent(inout) :: this
     real(qp), dimension(:,:), intent(in) :: matrix
     integer i,j

     !if(this%initialised) call finalise(this)

     this%n = size(matrix,1)
     this%m = size(matrix,2)
     
     allocate(this%matrix(this%n,this%m), this%factor(this%n,this%m), this%s(this%n), &
     this%tau(this%m) )

     this%matrix = matrix
     this%initialised = .true.
     !print*, 'NNN>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
     !do i =1 ,this%n
     !    do j =1 ,this%m
     !    write(*,'(F20.10,$)'), this%matrix(i,j)
     !    enddo
     !    write(*,*)
     !enddo
     !print*, 'NNN>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'

  endsubroutine LA_Matrix_Initialise

  subroutine LA_Matrix_Finalise(this)

     type(LA_Matrix), intent(inout) :: this

     if(.not. this%initialised) return

     this%n = 0
     this%m = 0
     if(allocated(this%matrix) ) deallocate(this%matrix)
     if(allocated(this%factor) ) deallocate(this%factor)
     if(allocated(this%s) ) deallocate(this%s)
     if(allocated(this%tau) ) deallocate(this%tau)
     this%initialised = .false.
     this%equilibrated = .false.
     this%factorised = NOT_FACTORISED

  endsubroutine LA_Matrix_Finalise

  subroutine LA_Matrix_Factorise(this,factor,error)

     type(LA_Matrix), intent(inout) :: this
     real(qp), dimension(:,:), intent(out), optional :: factor
     integer, optional, intent(out) :: error

     integer :: i, j, info
     real(dp) :: scond, amax

     !INIT_ERROR(error)

     if(.not. this%initialised) then
        print* ,'LA_Matrix_Factorise: object not initialised'
     endif

     if( this%n /= this%m ) then
        print*, 'LA_Matrix_Factorise: matrix not square'
     endif

     this%s = 1.0_qp

     do i = 1, this%n
        this%s(i) = 1.0_qp / sqrt(this%matrix(i,i))
     enddo
     scond = maxval(this%s) / minval(this%s)
     amax = maxval(this%matrix)
     !print*, 'scond',scond, amax

     this%equilibrated = ( scond < 0.1_qp )

     if( this%equilibrated ) then
        do i = 1, this%n
           this%factor(:,i) = this%matrix(:,i)*this%s(:)*this%s(i)
        enddo
     else
        this%factor = this%matrix
     endif

!#if defined(HAVE_QP) || defined(ALBERT_LAPACK)
!print*, 'HAVE QP'
!     call qpotrf(this%factor,info)
!#else
!print*, 'not have HAVE QP'
     call dpotrf('L', this%n, this%factor, this%n, info)
     print*, 'info of dpotrf', info
     do i = 2, this%n
        do j = 1, i
           this%factor(j,i) = this%factor(i,j)
        enddo
     enddo
!#endif

     if( info /= 0 ) then
        print*, 'LA_Matrix_Factorise: cannot factorise, error: '
     endif

     if( present(factor) ) then
        if( this%equilibrated ) then
           factor = 0.0_qp
           do i = 1, this%n
              do j = 1, this%n
                 factor(j,i) = this%factor(j,i) / this%s(i)
              enddo
           enddo
        else
           factor = this%factor
        endif
     endif

     this%factorised = CHOLESKY
        
  endsubroutine LA_Matrix_Factorise


  subroutine LA_Matrix_QR_Factorise(this,q,r,error)
     type(LA_Matrix), intent(inout) :: this         
     real(qp), dimension(:,:), intent(out), optional :: q, r
     integer, intent(out), optional :: error

     integer :: info
!#ifndef HAVE_QP
     real(dp), dimension(:), allocatable :: work
     integer :: lwork
!#endif

!     INIT_ERROR(error)

     this%factor = this%matrix

!#if defined(HAVE_QP) || defined(ALBERT_LAPACK)
!     call qgeqrf(this%factor,this%tau,info)
!#else

     allocate(work(1))
     lwork = -1
     call dgeqrf(this%n, this%m, this%factor, this%n, this%tau, work, lwork, info)
     lwork = nint(work(1))
     deallocate(work)

     allocate(work(lwork))
     call dgeqrf(this%n, this%m, this%factor, this%n, this%tau, work, lwork, info)
     deallocate(work)
!#endif

     if( info /= 0 ) then
        !RAISE_ERROR('LA_Matrix_QR_Factorise: '//(-info)//'-th parameter had an illegal value.',error)
     endif

     this%factorised = QR

     if( present(q) .or. present(r) ) call LA_Matrix_GetQR(this,q,r,error)

  endsubroutine LA_Matrix_QR_Factorise

  subroutine LA_Matrix_GetQR(this,q,r,error)
     type(LA_Matrix), intent(inout) :: this         
     real(qp), dimension(:,:), intent(out), optional :: q, r
     integer, intent(out), optional :: error

     integer :: j, info
!#ifndef HAVE_QP
     real(dp), dimension(:), allocatable :: work
     integer :: lwork
!#endif

!     INIT_ERROR(error)

     if( this%factorised /= QR ) then
        print*, 'LA_Matrix_GetQR: not QR-factorised, call LA_Matrix_QR_Factorise first.'
     endif

     if(present(q)) then
        call check_size('q', q, (/this%n,this%m/),'LA_Matrix_GetQR',error)
        q = this%factor

!#if defined(HAVE_QP) || defined(ALBERT_LAPACK)
!        call qorgqr(this%factor,this%tau,q,info)
!#else
        allocate(work(1))
        lwork = -1
        call dorgqr(this%n, this%m, this%m, q, this%n, this%tau, work, lwork, info)
        lwork = nint(work(1))
        deallocate(work)

        allocate(work(lwork))
        call dorgqr(this%n, this%m, this%m, q, this%n, this%tau, work, lwork, info)
        deallocate(work)
!#endif
     endif

     if(present(r)) then
        call check_size('r', r, (/this%m,this%m/),'LA_Matrix_GetQR',error)
        r = this%factor(1:this%m,1:this%m)
        do j = 1, this%m - 1
           r(j+1:this%m,j) = 0.0_dp
        enddo
     endif

     if( info /= 0 ) then
        !print *, 'LA_Matrix_QR_Factorise: '//(str(info))//'-th parameter had an illegal value.'
        print *, 'LA_Matrix_QR_Factorise: -th parameter had an illegal value.'
     endif

  endsubroutine LA_Matrix_GetQR

  subroutine LA_Matrix_QR_Solve_Matrix(this,matrix,result,error)
     type(LA_Matrix), intent(inout) :: this
     real(qp), dimension(:,:), intent(in) :: matrix
     real(qp), dimension(:,:), intent(out) :: result
     integer, intent(out), optional :: error

     real(qp), dimension(:,:), allocatable :: my_result
     integer :: info, i, j, n, o
!#ifndef HAVE_QP
     real(dp), dimension(:), allocatable :: work
     integer :: lwork
!#endif

!     INIT_ERROR(error)

     if(this%factorised == NOT_FACTORISED) then
        call LA_Matrix_QR_Factorise(this,error=error)
     elseif(this%factorised /= QR) then
        print*,  'LA_Matrix_QR_Solve_Matrix: matrix not QR-factorised'
     endif

     n = size(matrix,1)
     o = size(matrix,2)
     call check_size('result', result, (/this%m,o/),'LA_Matrix_QR_Solve_Matrix',error)

     if( n /= this%n ) then
        !RAISE_ERROR('LA_Matrix_QR_Solve_Matrix: dimensions of Q and matrix do not match.',error)
     endif

     allocate(my_result(n,o))
     my_result = matrix
!#if defined(HAVE_QP) || defined(ALBERT_LAPACK)
!     call qormqr(this%factor, this%tau, my_result, info)
!#else
     lwork = -1
     allocate(work(1))
     call dormqr('L', 'T', this%n, o, this%m, this%factor, this%n, this%tau, my_result, this%n, work, lwork, info)
     lwork = nint(work(1))
     deallocate(work)

     allocate(work(lwork))
     call dormqr('L', 'T', this%n, o, this%m, this%factor, this%n, this%tau, my_result, this%n, work, lwork, info)
     deallocate(work)
!#endif

     if( info /= 0 ) then
        !RAISE_ERROR('LA_Matrix_QR_QR_Solve_Matrix: '//(-info)//'-th parameter had an illegal value.',error)
     endif

     do i = 1, o
        do j = this%m, 2, -1
           my_result(j,i) = my_result(j,i)/this%factor(j,j)
           my_result(1:j-1,i) = my_result(1:j-1,i) - my_result(j,i)*this%factor(1:j-1,j)
        enddo
        my_result(1,i) = my_result(1,i) / this%factor(1,1)
     enddo

     result = my_result(1:this%m,:)
     deallocate(my_result)

  endsubroutine LA_Matrix_QR_Solve_Matrix

  subroutine LA_Matrix_QR_Solve_Vector(this,vector,result,error)
     type(LA_Matrix), intent(inout) :: this
     real(qp), dimension(:), intent(in) :: vector
     real(qp), dimension(:), intent(out) :: result
     integer, intent(out), optional :: error

     real(qp), dimension(:,:), allocatable :: my_result
     integer :: n, m

!     INIT_ERROR(error)

     n = size(vector)
     m = size(result)

     allocate(my_result(m,1))

     call LA_Matrix_QR_Solve_Matrix(this,reshape(vector,(/n,1/)),my_result,error=error)
     result = my_result(:,1)

     deallocate(my_result)

  endsubroutine LA_Matrix_QR_Solve_Vector
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X  Array size and shape checking routines
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


  ! The one dimensional checks seem unnecessarily complicated, but they can be
  ! easily generalised to higher dimensions (an extra : in intarray declaration
  ! for example)

  subroutine check_size_int_dim1(arrayname,intarray,n,caller,error)

    character(*),           intent(in) :: arrayname ! The name of the array to be checked
    integer, dimension(:),  intent(in) :: intarray  ! The array to be tested
    integer, dimension(:),  intent(in) :: n         ! The size that intarray should be
    character(*),           intent(in) :: caller    ! The name of the calling routine
    integer, intent(out), optional :: error

    integer, dimension(:), allocatable :: actual_size ! container for the actual size of intarray
    logical                            :: failed      ! Keeps track of any failures 
    integer                            :: i           ! dimension iterator

!    INIT_ERROR(error)

    failed = .false.
    allocate( actual_size( size(shape(intarray)) ) )
    actual_size = shape(intarray)

    if (size(actual_size) /= size(n)) then
       write(line,'(a,i0,a,i0,a)') caller//': '//arrayname//' is ',size(actual_size), &
            ' dimensional and not ',size(n),' dimensional as expected'
!       call print(line)
       failed = .true.
    else
       do i = 1, size(actual_size)
          if (actual_size(i) /= n(i)) then
             write(line,'(3(a,i0),a)') caller//': The size of dimension ',i,' of '//arrayname//' is ', &
                  actual_size(i),' and not ',n(i),' as expected'
!             call print(line)
             failed = .true.
          end if
       end do
    end if

    if (failed) then
!       RAISE_ERROR(trim(caller)//': Size checking failed', error)
    end if

  end subroutine check_size_int_dim1

  !overloaded subroutine which allows a scalar 'n' in the one dimensional case
  subroutine check_size_int_dim1_s(arrayname,intarray,n,caller,error)
    character(*),           intent(in) :: arrayname 
    integer, dimension(:),  intent(in) :: intarray  
    integer,                intent(in) :: n         
    character(*),           intent(in) :: caller
    integer, intent(out), optional :: error
    
!    INIT_ERROR(error)
    call check_size(arrayname,intarray,(/n/),caller,error)
!    PASS_ERROR(error)

  end subroutine check_size_int_dim1_s

  subroutine check_size_int_dim2(arrayname,intarray,n,caller,error)

    character(*),             intent(in) :: arrayname 
    integer, dimension(:,:),  intent(in) :: intarray  
    integer, dimension(:),    intent(in) :: n         
    character(*),             intent(in) :: caller
    integer, intent(out), optional :: error

    integer, dimension(:), allocatable :: actual_size 
    logical                            :: failed       
    integer                            :: i          

!    INIT_ERROR(error)
    failed = .false.
    allocate( actual_size( size(shape(intarray)) ) )
    actual_size = shape(intarray)

    if (size(actual_size) /= size(n)) then
       write(line,'(a,i0,a,i0,a)') caller//': '//arrayname//' is ',size(actual_size), &
            ' dimensional and not ',size(n),' dimensional as expected'
!       call print(line)
       failed = .true.
    else
       do i = 1, size(actual_size)
          if (actual_size(i) /= n(i)) then
             write(line,'(3(a,i0),a)') caller//': The size of dimension ',i,' of '//arrayname//' is ', &
                  actual_size(i),' and not ',n(i),' as expected'
!             call print(line)
             failed = .true.
          end if
       end do
    end if

    if (failed) then
!       RAISE_ERROR(trim(caller)//': Size checking failed', error)
    end if

  end subroutine check_size_int_dim2

  subroutine check_size_real_dim1(arrayname,realarray,n,caller,error)

    character(*),            intent(in) :: arrayname
    real(dp), dimension(:),  intent(in) :: realarray 
    integer,  dimension(:),  intent(in) :: n         
    character(*),            intent(in) :: caller
    integer, intent(out), optional :: error     

    integer, dimension(:), allocatable :: actual_size 
    logical                            :: failed      
    integer                            :: i           

!    INIT_ERROR(error)
    failed = .false.
    allocate( actual_size( size(shape(realarray)) ) )
    actual_size = shape(realarray)

    if (size(actual_size) /= size(n)) then
       write(line,'(a,i0,a,i0,a)') caller//': '//arrayname//' is ',size(actual_size), &
            ' dimensional and not ',size(n),' dimensional as expected'
!       call print(line)
       failed = .true.
    else
       do i = 1, size(actual_size)
          if (actual_size(i) /= n(i)) then
             write(line,'(3(a,i0),a)') caller//': The size of dimension ',i,' of '//arrayname//' is ', &
                  actual_size(i),' and not ',n(i),' as expected'
!             call print(line)
             failed = .true.
          end if
       end do
    end if

    if (failed) then
!       RAISE_ERROR(trim(caller)//': Size checking failed', error)
    end if

  end subroutine check_size_real_dim1

  !overloaded subroutine which allows a scalar 'n' in the one dimensional case
  subroutine check_size_real_dim1_s(arrayname,realarray,n,caller,error)
    character(*),           intent(in) :: arrayname 
    real(dp), dimension(:), intent(in) :: realarray  
    integer,                intent(in) :: n         
    character(*),           intent(in) :: caller
    integer, intent(out), optional :: error     
    
!    INIT_ERROR(error)
    call check_size(arrayname,realarray,(/n/),caller,error)
!    PASS_ERROR(error)

  end subroutine check_size_real_dim1_s

  subroutine check_size_real_dim2(arrayname,realarray,n,caller, error)

    character(*),              intent(in) :: arrayname 
    real(dp), dimension(:,:),  intent(in) :: realarray 
    integer,  dimension(:),    intent(in) :: n        
    character(*),              intent(in) :: caller
    integer, intent(out), optional :: error  

    integer, dimension(:), allocatable :: actual_size
    logical                            :: failed      
    integer                            :: i          

!    INIT_ERROR(error)
    failed = .false.
    allocate( actual_size( size(shape(realarray)) ) )
    actual_size = shape(realarray)

    if (size(actual_size) /= size(n)) then
       write(line,'(a,i0,a,i0,a)') caller//': '//arrayname//' is ',size(actual_size), &
            ' dimensional and not ',size(n),' dimensional as expected'
!       call print(line)
       failed = .true.
    else
       do i = 1, size(actual_size)
          if (actual_size(i) /= n(i)) then
             write(line,'(3(a,i0),a)') caller//': The size of dimension ',i,' of '//arrayname//' is ', &
                  actual_size(i),' and not ',n(i),' as expected'
!             call print(line)
             failed = .true.
          end if
       end do
    end if

    if (failed) then
!       RAISE_ERROR(trim(caller) //': Size checking failed. Expected: ' // n // ', got: ' // actual_size, error)
    end if

  end subroutine check_size_real_dim2

  subroutine check_size_real_dim3(arrayname,realarray,n,caller, error)

    character(*),              intent(in) :: arrayname 
    real(dp), dimension(:,:,:),intent(in) :: realarray 
    integer,  dimension(:),    intent(in) :: n        
    character(*),              intent(in) :: caller
    integer, intent(out), optional :: error  

    integer, dimension(:), allocatable :: actual_size
    logical                            :: failed      
    integer                            :: i          

!    INIT_ERROR(error)
    failed = .false.
    allocate( actual_size( size(shape(realarray)) ) )
    actual_size = shape(realarray)

    if (size(actual_size) /= size(n)) then
       write(line,'(a,i0,a,i0,a)') caller//': '//arrayname//' is ',size(actual_size), &
            ' dimensional and not ',size(n),' dimensional as expected'
!       call print(line)
       failed = .true.
    else
       do i = 1, size(actual_size)
          if (actual_size(i) /= n(i)) then
             write(line,'(3(a,i0),a)') caller//': The size of dimension ',i,' of '//arrayname//' is ', &
                  actual_size(i),' and not ',n(i),' as expected'
!             call print(line)
             failed = .true.
          end if
       end do
    end if

    if (failed) then
!       RAISE_ERROR(trim(caller) //': Size checking failed. Expected: ' // n // ', got: ' // actual_size, error)
    end if

  end subroutine check_size_real_dim3

  subroutine LA_Matrix_SVD_Allocate(this,s,u,v,error)

    type(LA_Matrix), intent(in) :: this
    real(dp), dimension(:), allocatable, intent(inout), optional :: s
    real(dp), dimension(:,:), allocatable, intent(inout), optional :: u, v
    integer, optional, intent(out) :: error

    !INIT_ERROR(error)

    if(.not.this%initialised) then
!       RAISE_ERROR('LA_Matrix_SVD: not initialised',error)
    endif

    if(present(s)) then
       if(allocated(s)) deallocate(s)
       allocate(s(min(this%n,this%m)))
    endif

    if(present(u)) then
       if(allocated(u)) deallocate(u)
       if(this%n <= this%m) then
          allocate(u(this%n,this%n))
       else
          allocate(u(this%n,this%m))
       endif
    endif

    if(present(v)) then
       if(allocated(v)) deallocate(v)
       if(this%n <= this%m) then
          allocate(v(this%m,this%n))
       else
          allocate(v(this%m,this%m))
       endif
    endif

  endsubroutine LA_Matrix_SVD_Allocate

  subroutine LA_Matrix_SVD(this,s,u,v,error)

    type(LA_Matrix), intent(in) :: this
    real(dp), dimension(:), intent(out), target, optional :: s
    real(dp), dimension(:,:), intent(out), target, optional :: u, v
    integer, optional, intent(out) :: error

    real(dp), dimension(:,:), allocatable :: a
    real(dp), dimension(:), allocatable :: work
    real(dp), dimension(:), pointer :: my_s
    real(dp), dimension(:,:), pointer :: my_u, my_vt
    real(dp), dimension(1,1), target :: dummy_u, dummy_vt
    real(dp) :: tmp
    character(len=1) :: jobu, jobvt
    integer :: lwork, info, i, j
    integer(kind=8) :: lwork_min

    !INIT_ERROR(error)

    if(.not.this%initialised) then
!       RAISE_ERROR('LA_Matrix_SVD: not initialised',error)
    endif

    allocate(a(this%n,this%m))
    a = this%matrix

    if(present(s)) then
       call check_size('s',s,min(this%n,this%m),'LA_Matrix_SVD',error=error)
       my_s => s
    else
       allocate(my_s(min(this%n,this%m)))
    endif

    if(present(u)) then
       if(this%n <= this%m) then
          call check_size('u',u,(/this%n,this%n/),'LA_Matrix_SVD',error=error)
          jobu = "A"
          my_u => u
       else
          call check_size('u',u,(/this%n,this%m/),'LA_Matrix_SVD',error=error)
          jobu = "S"
          my_u => u
       endif
    else
       jobu = "N"
       my_u => dummy_u
    endif

    if(present(v)) then
       if(this%n <= this%m) then
          call check_size('v',v,(/this%m,this%n/),'LA_Matrix_SVD',error=error)
          jobvt = "O"
          my_vt => dummy_vt
       else
          call check_size('v',v,(/this%m,this%m/),'LA_Matrix_SVD',error=error)
          jobvt = "A"
          my_vt => v
       endif
    else
       jobvt = "N"
       my_vt => dummy_vt
    endif

    allocate(work(1))
    lwork = -1
    call dgesvd(jobu, jobvt, this%n, this%m, a, this%n, my_s, my_u, this%n, my_vt, this%m, work, lwork, info)
    if( work(1) > huge(lwork) ) then ! optimal size of work is a bit too large.
       ! that's the minimum size of the work array
       lwork_min = max( 3*min(this%m, this%n) + max(this%n,this%m), 5*min(this%m,this%n) )
       if( lwork_min > huge(lwork) ) then
!          RAISE_ERROR("LA_Matrix_SVD: temporary array for SVD would be too large.", error)
       else
          ! max out the work array
          lwork = huge(lwork)
       endif
    else
       ! otherwise just go with the optimum
       lwork = ceiling(work(1))
    endif
       
    deallocate(work)

    allocate(work(lwork))
    call dgesvd(jobu, jobvt, this%n, this%m, a, this%n, my_s, my_u, this%n, my_vt, this%m, work, lwork, info)
    deallocate(work)

    if( this%n <= this%m ) then
       if(present(v)) then
          do i = 1, this%n
             do j = 1, this%m
                v(j,i) = a(i,j)
             enddo
          enddo
       endif
    else
       if(present(v)) then
          do i = 1, this%m
             do j = i+1, this%m
                tmp = v(i,j)
                v(i,j) = v(j,i)
                v(j,i) = tmp
             enddo
          enddo
       endif
    endif

    if(allocated(a)) deallocate(a)
    my_u => null()
    my_vt => null()

    if(present(s)) then
       my_s => null()
    else
       deallocate(my_s)
    endif

    if(info < 0) then
!       RAISE_ERROR('LA_Matrix_SVD: '//(-info)//'-th parameter had an illegal value.',error)
    elseif( info > 0) then
!       RAISE_ERROR('LA_Matrix_SVD: singular value decomposition of the bidiagonal matrix did not converge, '//info//' superdiagonals did not converge to zero.',error)
    endif
  endsubroutine LA_Matrix_SVD

  subroutine pseudo_inverse(this,inverse,error)
    real(dp),intent(in), dimension(:,:) :: this
    real(dp),intent(out), dimension(:,:) :: inverse
    integer, optional, intent(out) :: error

    type(LA_Matrix) :: LA_this

    !INIT_ERROR(error)

    call LA_Matrix_initialise(LA_this,this)
    call LA_Matrix_PseudoInverse(LA_this,inverse,error=error)
    call LA_Matrix_finalise(LA_this)

  endsubroutine pseudo_inverse

  subroutine LA_Matrix_PseudoInverse(this, inverse, error)
    type(LA_Matrix), intent(in) :: this
    real(dp), dimension(:,:), intent(out) :: inverse
    integer, optional, intent(out) :: error

    real(dp), dimension(:), allocatable :: s
    real(dp), dimension(:,:), allocatable :: u, v, inverse_tmp

  !  INIT_ERROR(error)

    if(.not.this%initialised) then
   !    RAISE_ERROR('LA_Matrix_SVD: not initialised',error)
    endif

    call check_size('inverse',inverse,(/this%m,this%n/),'LA_Matrix_PseudoInverse',error=error)

    call LA_Matrix_SVD_Allocate(this,s=s,u=u,v=v,error=error)
    call LA_Matrix_SVD(this,s=s,u=u,v=v,error=error)

    where(abs(s) > TOL_SVD)
       s = 1.0_dp / s
    elsewhere
       s = 0.0_dp
    endwhere

    allocate(inverse_tmp(size(v,1),size(v,2)))
    inverse_tmp = v .multd. s
    inverse = inverse_tmp .mult. transpose(u)
    deallocate(inverse_tmp)

    if(allocated(s)) deallocate(s)
    if(allocated(u)) deallocate(u)
    if(allocated(v)) deallocate(v)

  endsubroutine LA_Matrix_PseudoInverse

   subroutine heap_sort_i(array, r_data)
     integer,            dimension(:), intent(inout)  :: array
     real(dp), optional, dimension(:), intent(inout)  :: r_data

     ! ---

     integer   :: N, i, j, root, tmpi
     real(DP)  :: tmpr

     ! ---

     N = size(array)

     do i = N/2, 1, -1

        j = i
        call siftdown(j, N)

     enddo

     do i = N, 2, -1
        
        ! Swap
        tmpi       = array(1)
        array(1)   = array(i)
        array(i)   = tmpi
        if (present(r_data)) then
           tmpr       = r_data(1)
           r_data(1)  = r_data(i)
           r_data(i)  = tmpr
        endif

        root = 1
        j = i -1
        call siftdown(root, j)
    
     enddo

   contains

     subroutine siftdown(root, bottom)
       integer, intent(inout)  :: root
       integer, intent(in)     :: bottom

       ! ---

       logical  :: done
       integer  :: maxchild

       ! ---

       done = .false.

       do while ((root*2 <= bottom) .and. .not. done)

          if (root*2 == bottom) then
             maxchild = root * 2
          else if (array(root*2) > array(root*2+1)) then
             maxchild = root * 2
          else
             maxchild = root*2 + 1
          endif

          if (array(root) < array(maxchild)) then

             ! Swap
             tmpi              = array(root)
             array(root)       = array(maxchild)
             array(maxchild)   = tmpi
             if (present(r_data)) then
                tmpr              = r_data(root)
                r_data(root)      = r_data(maxchild)
                r_data(maxchild)  = tmpr
             endif

             root = maxchild
          else
             done = .true.
          endif

       enddo

     endsubroutine siftdown

   endsubroutine heap_sort_i


   !% Sort an array of integers into ascending order.
   !% The function uses heapsort, which always scales as N log N.
   !% i_data and r_data are a accompanying arrays of integers and reals
   !% on which the same reordering is performed
   !% (Initial implementation by Andreas Wonisch)
   subroutine heap_sort_r(array, i_data, r_data)
     real(dp),           dimension(:), intent(inout)  :: array
     integer,  optional, dimension(:), intent(inout)  :: i_data
     real(dp), optional, dimension(:), intent(inout)  :: r_data

     ! ---

     integer   :: N, i, j, root, tmpi
     real(DP)  :: tmpr

     ! ---

     N = size(array)

     do i = N/2, 1, -1

        j = i
        call siftdown(j, N)

     enddo

     do i = N, 2, -1
        
        ! Swap
        tmpr       = array(1)
        array(1)   = array(i)
        array(i)   = tmpr
        if (present(i_data)) then
           tmpi       = i_data(1)
           i_data(1)  = i_data(i)
           i_data(i)  = tmpi
        endif
        if (present(r_data)) then
           tmpr       = r_data(1)
           r_data(1)  = r_data(i)
           r_data(i)  = tmpr
        endif

        root = 1
        j = i -1
        call siftdown(root, j)
    
     enddo

   contains

     subroutine siftdown(root, bottom)
       integer, intent(inout)  :: root
       integer, intent(in)     :: bottom

       ! ---

       logical  :: done
       integer  :: maxchild

       ! ---

       done = .false.

       do while ((root*2 <= bottom) .and. .not. done)

          if (root*2 == bottom) then
             maxchild = root * 2
          else if (array(root*2) > array(root*2+1)) then
             maxchild = root * 2
          else
             maxchild = root*2 + 1
          endif

          if (array(root) < array(maxchild)) then

             ! Swap
             tmpr             = array(root)
             array(root)      = array(maxchild)
             array(maxchild)  = tmpr
             if (present(i_data)) then
                tmpi              = i_data(root)
                i_data(root)      = i_data(maxchild)
                i_data(maxchild)  = tmpi
             endif
             if (present(r_data)) then
                tmpr              = r_data(root)
                r_data(root)      = r_data(maxchild)
                r_data(maxchild)  = tmpr
             endif

             root = maxchild
          else
             done = .true.
          endif

       enddo

     endsubroutine siftdown

   endsubroutine heap_sort_r

   subroutine heap_sort_r_2dim(array, i_data, r_data)
     real(dp),           dimension(:,:), intent(inout)  :: array
     integer,  optional, dimension(:), intent(inout)    :: i_data
     real(dp), optional, dimension(:), intent(inout)    :: r_data

     ! ---

     integer   :: N, i, j, root, tmpi, d
     real(dp), dimension(:), allocatable  :: tmp_array
     real(dp) :: tmpr

     ! ---

     N = size(array,2)
     d = size(array,1)
     allocate(tmp_array(d))

     do i = N/2, 1, -1

        j = i
        call siftdown(j, N)

     enddo

     do i = N, 2, -1
        
        ! Swap
        tmp_array(:) = array(:,1)
        array(:,1)   = array(:,i)
        array(:,i)   = tmp_array(:)
        if (present(i_data)) then
           tmpi       = i_data(1)
           i_data(1)  = i_data(i)
           i_data(i)  = tmpi
        endif
        if (present(r_data)) then
           tmpr       = r_data(1)
           r_data(1)  = r_data(i)
           r_data(i)  = tmpr
        endif

        root = 1
        j = i -1
        call siftdown(root, j)
    
     enddo

     deallocate(tmp_array)

   contains

     subroutine siftdown(root, bottom)
       integer, intent(inout)  :: root
       integer, intent(in)     :: bottom

       ! ---

       logical  :: done
       integer  :: maxchild

       ! ---

       done = .false.

       do while ((root*2 <= bottom) .and. .not. done)

          if (root*2 == bottom) then
             maxchild = root * 2
          else if ( real_array_gt(array(:,root*2),array(:,root*2+1)) ) then  !array(:,root*2) > array(:,root*2+1)
             maxchild = root * 2
          else
             maxchild = root*2 + 1
          endif

          if ( real_array_lt(array(:,root),array(:,maxchild)) ) then ! array(:,root) < array(:,maxchild)

             ! Swap
             tmp_array(:)       = array(:,root)
             array(:,root)      = array(:,maxchild)
             array(:,maxchild)  = tmp_array(:)
             if (present(i_data)) then
                tmpi              = i_data(root)
                i_data(root)      = i_data(maxchild)
                i_data(maxchild)  = tmpi
             endif
             if (present(r_data)) then
                tmpr              = r_data(root)
                r_data(root)      = r_data(maxchild)
                r_data(maxchild)  = tmpr
             endif

             root = maxchild
          else
             done = .true.
          endif

       enddo

     endsubroutine siftdown

   endsubroutine heap_sort_r_2dim

   subroutine heap_sort_i_2dim(array, i_data, r_data)
     integer,            dimension(:,:), intent(inout)  :: array
     integer,  optional, dimension(:), intent(inout)    :: i_data
     real(dp), optional, dimension(:), intent(inout)    :: r_data

     ! ---

     integer   :: N, i, j, root, tmpi, d
     integer, dimension(:), allocatable  :: tmp_array
     real(dp) :: tmpr

     ! ---

     N = size(array,2)
     d = size(array,1)
     allocate(tmp_array(d))

     do i = N/2, 1, -1

        j = i
        call siftdown(j, N)

     enddo

     do i = N, 2, -1
        
        ! Swap
        tmp_array(:) = array(:,1)
        array(:,1)   = array(:,i)
        array(:,i)   = tmp_array(:)
        if (present(i_data)) then
           tmpi       = i_data(1)
           i_data(1)  = i_data(i)
           i_data(i)  = tmpi
        endif
        if (present(r_data)) then
           tmpr       = r_data(1)
           r_data(1)  = r_data(i)
           r_data(i)  = tmpr
        endif

        root = 1
        j = i -1
        call siftdown(root, j)
    
     enddo

     deallocate(tmp_array)

   contains

     subroutine siftdown(root, bottom)
       integer, intent(inout)  :: root
       integer, intent(in)     :: bottom

       ! ---

       logical  :: done
       integer  :: maxchild

       ! ---

       done = .false.

       do while ((root*2 <= bottom) .and. .not. done)

          if (root*2 == bottom) then
             maxchild = root * 2
          else if ( int_array_gt(array(:,root*2),array(:,root*2+1)) ) then  !array(:,root*2) > array(:,root*2+1)
             maxchild = root * 2
          else
             maxchild = root*2 + 1
          endif

          if ( int_array_lt(array(:,root),array(:,maxchild)) ) then ! array(:,root) < array(:,maxchild)

             ! Swap
             tmp_array(:)       = array(:,root)
             array(:,root)      = array(:,maxchild)
             array(:,maxchild)  = tmp_array(:)
             if (present(i_data)) then
                tmpi              = i_data(root)
                i_data(root)      = i_data(maxchild)
                i_data(maxchild)  = tmpi
             endif
             if (present(r_data)) then
                tmpr              = r_data(root)
                r_data(root)      = r_data(maxchild)
                r_data(maxchild)  = tmpr
             endif

             root = maxchild
          else
             done = .true.
          endif

       enddo

     endsubroutine siftdown

   endsubroutine heap_sort_i_2dim

  pure function real_array_gt(array1,array2)
     real(dp), dimension(:), intent(in) :: array1
     real(dp), dimension(size(array1)), intent(in) :: array2
     logical :: real_array_gt

     integer :: i

     do i = 1, size(array1)
        if(array1(i) .feq. array2(i)) then
           cycle
        elseif(array1(i) > array2(i)) then
           real_array_gt = .true.
           return
        else !array1(i) < array2(i)
           real_array_gt = .false.
           return
        endif
     enddo

     real_array_gt = .false.

  endfunction real_array_gt

  pure function real_array_ge(array1,array2)
     real(dp), dimension(:), intent(in) :: array1
     real(dp), dimension(size(array1)), intent(in) :: array2
     logical :: real_array_ge

     integer :: i

     do i = 1, size(array1)
        if(array1(i) .feq. array2(i)) then
           cycle
        elseif(array1(i) > array2(i)) then
           real_array_ge = .true.
           return
        else !array1(i) < array2(i)
           real_array_ge = .false.
           return
        endif
     enddo

     real_array_ge = .true.

  endfunction real_array_ge

  pure function real_array_lt(array1,array2)
     real(dp), dimension(:), intent(in) :: array1
     real(dp), dimension(size(array1)), intent(in) :: array2
     logical :: real_array_lt

     real_array_lt = .not.(real_array_ge(array1,array2))
  endfunction real_array_lt

   pure function complex_feq(x,y) result(feq)

     complex(dp), intent(in) :: x, y
     logical              :: feq

     if ( (abs(real(x-y)) > NUMERICAL_ZERO * abs(real(x))) .or. &
          (abs(aimag(x-y)) > NUMERICAL_ZERO * abs(aimag(x))) ) then
        feq = .false.
     else
        feq = .true.
     end if
     
   end function complex_feq

   pure function real_feq(x,y) result(feq)

     real(dp), intent(in) :: x, y
     logical              :: feq

     feq = (abs(x-y) <= NUMERICAL_ZERO * (abs(x)+abs(y))/2.0_dp) .or. (abs(x-y) <= NUMERICAL_ZERO)
     
   end function real_feq

  ! function matrix_feq(matrix1,matrix2) result (feq)
  !   real(dp),intent(in), dimension(:,:) :: matrix1
  !   real(dp),intent(in), dimension(:,:) :: matrix2

  !   integer::i,j
  !   logical::feq
  !   
  !   call check_size('Matrix2',matrix2,shape(matrix1),'Matrix_FEQ')
  !   
  !   feq =.true.
  !   do j=1,size(matrix1,2)
  !      do i=1,size(matrix1,1)
  !         if (matrix1(i,j).fne.matrix2(i,j)) then
  !            feq=.false.
  !            return
  !         end if
  !      end do
  !   end do
  !   
  ! end function matrix_feq

  !function vector_feq(vector1,vector2) 
  !  real(dp),intent(in), dimension(:) ::vector1,vector2
  !  logical::vector_feq
  !  integer::i
  ! 
  !  if(size(vector1) /= size(vector2)) &
  !       call check_size('Vector2',vector2,size(vector1),'Vector_FEQ')
  !
  !  vector_feq=.true.
  !  do i=1,size(vector1)
  !     if (vector1(i).fne.vector2(i)) then 
  !        vector_feq=.false.
  !        return
  !     end if
  !  end do

  !end function vector_feq
   function matrix_product_vect_asdiagonal_dd(matrix,vect) result (prodmatrix)
     real(dp),intent(in),dimension(:,:)::matrix
     real(dp),intent(in), dimension(:) :: vect
     real(dp),dimension(size(matrix,1),size(matrix,2))::prodmatrix

     call matrix_product_vect_asdiagonal_sub(prodmatrix, matrix, vect)

   end function matrix_product_vect_asdiagonal_dd

   ! m(:,:) .multd. v(:)
   !
   ! returns product of matrix and vector as diagonal of another matrix
   function matrix_product_vect_asdiagonal_zz(matrix,vect) result (prodmatrix)
     complex(dp),intent(in), dimension(:) :: vect
     complex(dp),dimension(size(vect),size(vect))::prodmatrix
     complex(dp),intent(in),dimension(:,:)::matrix

     call matrix_product_vect_asdiagonal_sub(prodmatrix, matrix, vect)

   end function matrix_product_vect_asdiagonal_zz

   subroutine matrix_product_vect_asdiagonal_sub_ddd(lhs, matrix, vect) 
    real(dp), dimension(:,:), intent(out) :: lhs
    real(dp), dimension(:,:), intent(in) :: matrix
    real(dp), dimension(:), intent(in) :: vect

    integer  :: i
     
!$omp parallel do private(i) shared(lhs,vect,matrix)
    do i = 1, size(vect)
       lhs(:,i) = vect(i) * matrix(:,i)
    enddo
     
   endsubroutine matrix_product_vect_asdiagonal_sub_ddd

   ! subroutine form of m(:,:) .multd. v(:) for complex = real * complex
   !
   ! first argument set to product of matrix and a vector as a diagonal of another matrix (no temp on stack)
   subroutine matrix_product_vect_asdiagonal_sub_zdz(lhs, matrix, vect) 
    complex(dp), dimension(:,:), intent(out) :: lhs
    real(dp), dimension(:,:), intent(in) :: matrix
    complex(dp), dimension(:), intent(in) :: vect

    integer  :: i
     
!$omp parallel do private(i) shared(lhs,vect,matrix)
    do i = 1, size(vect)
       lhs(:,i) = vect(i) * matrix(:,i)
    enddo
     
   endsubroutine matrix_product_vect_asdiagonal_sub_zdz

   ! subroutine form of m(:,:) .multd. v(:) for complex = complex * real
   !
   ! first argument set to product of matrix and a vector as a diagonal of another matrix (no temp on stack)
   subroutine matrix_product_vect_asdiagonal_sub_zzd(lhs, matrix, vect) 
    complex(dp), dimension(:,:), intent(out) :: lhs
    complex(dp), dimension(:,:), intent(in) :: matrix
    real(dp), dimension(:), intent(in) :: vect

    integer :: i
     
!$omp parallel do private(i) shared(lhs,vect,matrix)
    do i = 1, size(vect)
       lhs(:,i) = vect(i) * matrix(:,i)
    enddo
     
   endsubroutine matrix_product_vect_asdiagonal_sub_zzd

   ! subroutine form of m(:,:) .multd. v(:) for complex = complex * complex
   !
   ! first argument set to product of matrix and a vector as a diagonal of another matrix (no temp on stack)
   subroutine matrix_product_vect_asdiagonal_sub_zzz(lhs, matrix, vect) 
    complex(dp), dimension(:,:), intent(out) :: lhs
    complex(dp), dimension(:,:), intent(in) :: matrix
    complex(dp), dimension(:), intent(in) :: vect

    integer :: i
     
!$omp parallel do private(i) shared(lhs,vect,matrix)
    do i = 1, size(vect)
       lhs(:,i)=vect(i)*matrix(:,i)
    end do
     
   endsubroutine matrix_product_vect_asdiagonal_sub_zzz

   ! subroutine form of v(:) .multd. m(:,:) for real = real * real
   !
   ! first argument set to product of matrix and a vector as a diagonal of another matrix (no temp on stack)
   subroutine vect_asdiagonal_product_matrix_sub_ddd(lhs, vectL, matrix) 
     real(dp), dimension(:,:), intent(out) :: lhs
     real(dp), dimension(:), intent(in) :: vectL
     real(dp), dimension(:,:), intent(in) :: matrix

     integer :: i
     
     do i = 1, size(matrix,2)
       lhs(:,i)=vectL(:)*matrix(:,i)
     end do
     
   end subroutine vect_asdiagonal_product_matrix_sub_ddd

   ! subroutine form of v(:) .multd. m(:,:) for complex = complex * real
   !
   ! first argument set to product of matrix and a vector as a diagonal of another matrix (no temp on stack)
   subroutine vect_asdiagonal_product_matrix_sub_zzd(lhs, vectL, matrix) 
     complex(dp), dimension(:,:), intent(out) :: lhs
     complex(dp), dimension(:), intent(in) :: vectL
     real(dp), dimension(:,:), intent(in) :: matrix

     integer :: i
     
     do i = 1, size(matrix,2)
       lhs(:,i)=vectL(:)*matrix(:,i)
     end do
     
   end subroutine vect_asdiagonal_product_matrix_sub_zzd

   ! subroutine form of v(:) .multd. m(:,:) for complex = real * complex
   !
   ! first argument set to product of matrix and a vector as a diagonal of another matrix (no temp on stack)
   subroutine vect_asdiagonal_product_matrix_sub_zdz(lhs, vectL, matrix) 
     complex(dp), dimension(:,:), intent(out) :: lhs
     real(dp), dimension(:), intent(in) :: vectL
     complex(dp), dimension(:,:), intent(in) :: matrix

     integer :: i
     
     do i = 1, size(matrix,2)
       lhs(:,i)=vectL(:)*matrix(:,i)
     end do
     
   end subroutine vect_asdiagonal_product_matrix_sub_zdz

   ! subroutine form of v(:) .multd. m(:,:) for complex = complex * complex
   !
   ! first argument set to product of matrix and a vector as a diagonal of another matrix (no temp on stack)
   subroutine vect_asdiagonal_product_matrix_sub_zzz(lhs, vectL, matrix) 
     complex(dp), dimension(:,:), intent(out) :: lhs
     complex(dp), dimension(:), intent(in) :: vectL
     complex(dp), dimension(:,:), intent(in) :: matrix

     integer :: i
     
     do i = 1, size(matrix,2)
       lhs(:,i)=vectL(:)*matrix(:,i)
     end do
   end subroutine vect_asdiagonal_product_matrix_sub_zzz

  function matrix_product_vect(matrix,vect) result (prodvect)
    real(dp),intent(in), dimension(:)   :: vect
    real(dp),intent(in), dimension(:,:) :: matrix
    real(dp), dimension(size(matrix,1)) ::prodvect
    integer::N,M,i,j

    N=size(matrix,1)
    M=size(matrix,2)

    if (M /= size(vect)) call check_size('Vector',vect,M,'Matrix_Product_Vect')
    prodvect = 0.0_dp
    do j=1,M
       forall (i=1:N) prodvect(i) = prodvect(i) + matrix(i,j)*vect(j)
    end do

  end function matrix_product_vect

   function matrix_product_ddd(matrix1,matrix2) result (prodmatrix)
     real(dp),intent(in), dimension(:,:) :: matrix1
     real(dp),intent(in), dimension(:,:) :: matrix2
     real(dp), dimension(size(matrix1,1),size(matrix2,2)) ::prodmatrix

     call matrix_product_sub(prodmatrix, matrix1, matrix2)
   end function matrix_product_ddd
  function matrix_product_int_vect(matrix,intvect) result(prodvect)

     real(dp), intent(in), dimension(:,:) :: matrix
     integer,  intent(in), dimension(:)   :: intvect
     real(dp)                             :: prodvect(size(matrix,1))
     integer                              :: M,N,i,j

     N=size(matrix,1)
     M=size(matrix,2)
     if(M /= size(intvect)) &
          call check_size('Integer Vector',intvect,M,'Matrix_product_int_vect')
     prodvect = 0.0_dp
     do j=1,M
        forall(i=1:N) prodvect(i) = prodvect(i)+matrix(i,j)*intvect(j)
     end do

   end function matrix_product_int_vect

   function matrix_product_zzz(matrix1,matrix2) result (prodmatrix)
     complex(dp),intent(in), dimension(:,:) :: matrix1
     complex(dp),intent(in), dimension(:,:) :: matrix2
     complex(dp), dimension(size(matrix1,1),size(matrix2,2)) ::prodmatrix

     call matrix_product_sub(prodmatrix, matrix1, matrix2)
   end function matrix_product_zzz

  pure function int_array_ge(array1,array2) result(ge)

    integer, intent(in) :: array1(:)
    integer, intent(in) :: array2(size(array1))
    logical             :: ge
    integer             :: i

    ge = .true.
    i = 1
    do while(i <= size(array1))
       if (array1(i) < array2(i)) then
          ge =.false.
          return
       else if (array1(i) > array2(i)) then
          return
       end if
       i = i + 1
    end do

  end function int_array_ge

  pure function int_array_gt(array1,array2) result(gt)

    integer, intent(in) :: array1(:)
    integer, intent(in) :: array2(size(array1))
    logical             :: gt
    integer             :: i

    gt = .true.
    i = 1
    do while(i <= size(array1))
       if (array1(i) < array2(i)) then
          gt =.false.
          return
       else if (array1(i) > array2(i)) then
          return
       end if
       i = i + 1
    end do
    gt = .false.

  end function int_array_gt

  pure function int_array_lt(array1,array2) result(lt)

    integer, intent(in) :: array1(:)
    integer, intent(in) :: array2(size(array1))
    logical             :: lt

    lt = .not.int_array_ge(array1,array2)

  end function int_array_lt

  pure function int_array_le(array1,array2) result(le)

    integer, intent(in) :: array1(:)
    integer, intent(in) :: array2(size(array1))
    logical             :: le

    le = .not.int_array_gt(array1,array2)

  end function int_array_le

  pure function optional_default_l(def, opt_val)
    logical, intent(in) :: def
    logical, intent(in), optional :: opt_val
    logical :: optional_default_l

    if (present(opt_val)) then
      optional_default_l = opt_val
    else
      optional_default_l = def
    endif

  end function optional_default_l

  pure function optional_default_i(def, opt_val)
    integer, intent(in) :: def
    integer, intent(in), optional :: opt_val
    integer :: optional_default_i

    if (present(opt_val)) then
      optional_default_i = opt_val
    else
      optional_default_i = def
    endif

  end function optional_default_i

  pure function optional_default_ia(def, opt_val)
    integer, intent(in) :: def(:)
    integer, intent(in), optional :: opt_val(size(def))
    integer :: optional_default_ia(size(def))

    if (present(opt_val)) then
      optional_default_ia = opt_val
    else
      optional_default_ia = def
    endif

  end function optional_default_ia


  pure function optional_default_r(def, opt_val)
    real(dp), intent(in) :: def
    real(dp), intent(in), optional :: opt_val
    real(dp) :: optional_default_r

    if (present(opt_val)) then
      optional_default_r = opt_val
    else
      optional_default_r = def
    endif

  end function optional_default_r

  pure function optional_default_ra(def, opt_val)
    real(dp), intent(in) :: def(:)
    real(dp), intent(in), optional :: opt_val(size(def))
    real(dp) :: optional_default_ra(size(def))

    if (present(opt_val)) then
      optional_default_ra = opt_val
    else
      optional_default_ra = def
    endif

  end function optional_default_ra


  pure function optional_default_z(def, opt_val)
    complex(dp), intent(in) :: def
    complex(dp), intent(in), optional :: opt_val
    complex(dp) :: optional_default_z

    if (present(opt_val)) then
      optional_default_z = opt_val
    else
      optional_default_z = def
    endif

  end function optional_default_z

  pure function optional_default_c(def, opt_val)
    character(len=*), intent(in) :: def
    character(len=*), intent(in), optional :: opt_val
    character(SYSTEM_STRING_LENGTH) :: optional_default_c

    if (present(opt_val)) then
      optional_default_c = opt_val
    else
      optional_default_c = def
    endif

  end function optional_default_c

  pure function optional_default_ca(def, opt_val)
    character(len=*), dimension(:), intent(in) :: def
    character(len=*), dimension(:), intent(in), optional :: opt_val
    character(SYSTEM_STRING_LENGTH), dimension(size(def)) :: optional_default_ca

    if (present(opt_val)) then
      optional_default_ca = opt_val
    else
      optional_default_ca = def
    endif

  end function optional_default_ca

  function ran_uniform()
    real(dp)::ran_uniform
    logical :: system_use_fortran_random
    system_use_fortran_random = .true.

       call random_number(ran_uniform)
  end function ran_uniform

   subroutine matrix_product_sub_ddd(lhs, matrix1, matrix2, m1_transpose, m2_transpose, &
    lhs_factor, rhs_factor)
      real(dp), intent(out) :: lhs(:,:)
      real(dp), intent(in) :: matrix1(:,:), matrix2(:,:)
      logical, intent(in), optional :: m1_transpose, m2_transpose
      real(dp), intent(in), optional :: lhs_factor, rhs_factor

     integer::M,N,maxd,K
     character(len=1) :: m1_transp, m2_transp
     integer :: m1_r, m1_c, m2_r, m2_c
     real(dp) :: my_lhs_factor, my_rhs_factor

     m1_transp = 'N'
     m2_transp = 'N'
     if (present(m1_transpose)) then
      if (m1_transpose) m1_transp = 'T'
     endif
     if (present(m2_transpose)) then
      if (m2_transpose) m2_transp = 'T'
     endif

     if (m1_transp == 'N') then
       m1_r = 1; m1_c = 2
     else
       m1_r = 2; m1_c = 1
     endif
     if (m2_transp == 'N') then
       m2_r = 1; m2_c = 2
     else
       m2_r = 2; m2_c = 1
     endif

     my_lhs_factor = optional_default(0.0_dp, lhs_factor)
     my_rhs_factor = optional_default(1.0_dp, rhs_factor)

     call check_size('lhs',lhs,(/size(matrix1,m1_r),size(matrix2,m2_c)/),'Matrix_Product')
     if (m2_transp == 'N') then
       call check_size('Matrix2',matrix2,(/size(matrix1,m1_c),size(matrix2,m2_c)/),'Matrix_Product')
     else
       call check_size('Matrix2',matrix2,(/size(matrix2,m2_c),size(matrix1,m1_c)/),'Matrix_Product')
     endif

     N = size(lhs,1) !# of rows of lhs
     M = size(lhs,2) !# of columns of rhs
     K = size(matrix1,m1_c) !!shared dimension
     maxd=max(N,M,K) 

     if (use_intrinsic_blas) then
        if (m1_transp == 'T') then
          if (m2_transp == 'T') then
            lhs=lhs*my_lhs_factor + my_rhs_factor*matmul(transpose(matrix1),transpose(matrix2))
          else
            lhs=lhs*my_lhs_factor + my_rhs_factor*matmul(transpose(matrix1),matrix2)
          endif
        else
          if (m2_transp == 'T') then
            lhs=lhs*my_lhs_factor + my_rhs_factor*matmul(matrix1,transpose(matrix2))
          else
            lhs=lhs*my_lhs_factor + my_rhs_factor*matmul(matrix1,matrix2)
          endif
        endif
     else
        call DGEMM(m1_transp,m2_transp, N,M,K,my_rhs_factor,matrix1,size(matrix1,1),matrix2,size(matrix2,1),&
	  my_lhs_factor,lhs,size(lhs,1))
     endif

   end subroutine matrix_product_sub_ddd

   ! subroutine form of m1(:,:) .mult. m2(:,:)
   !
   ! set first argument to matrix product (no temporary on the stack)
   subroutine matrix_product_sub_zzz(lhs, matrix1, matrix2, m1_transpose, m1_conjugate, &
    m2_transpose, m2_conjugate, lhs_factor, rhs_factor)
      complex(dp), intent(out) :: lhs(:,:)
      complex(dp), intent(in) :: matrix1(:,:), matrix2(:,:)
      logical, intent(in), optional :: m1_transpose, m1_conjugate, m2_transpose, m2_conjugate
      complex(dp), intent(in), optional :: lhs_factor, rhs_factor

     integer::M,N,maxd,K
     logical :: m1_transp, m2_transp, m1_conjg, m2_conjg
     character(len=1) :: m1_op, m2_op
     integer :: m1_r, m1_c, m2_r, m2_c
     complex(dp) :: my_lhs_factor, my_rhs_factor

     m1_transp = .false.
     m2_transp = .false.
     m1_conjg = .false.
     m2_conjg = .false.
     if (present(m1_transpose)) m1_transp = m1_transpose
     if (present(m2_transpose)) m2_transp = m2_transpose
     if (present(m1_conjugate)) m1_conjg = m1_conjugate
     if (present(m2_conjugate)) m2_conjg = m2_conjugate

     !if (m1_conjg .and. m1_transp) call system_abort("Called matrix_product_sub_zzz with m1_transp and m1_conjg true")
     !if (m2_conjg .and. m2_transp) call system_abort("Called matrix_product_sub_zzz with m2_transp and m2_conjg true")

     if (m1_transp) then
       m1_op = 'T'
     else if (m1_conjg) then
       m1_op = 'C'
     else
       m1_op = 'N'
     end if

     if (m2_transp) then
       m2_op = 'T'
     else if (m2_conjg) then
       m2_op = 'C'
     else
       m2_op = 'N'
     end if

     if (m1_op == 'N') then
       m1_r = 1; m1_c = 2
     else
       m1_r = 2; m1_c = 1
     endif
     if (m2_op == 'N') then
       m2_r = 1; m2_c = 2
     else
       m2_r = 2; m2_c = 1
     endif

     my_lhs_factor = optional_default(cmplx(0.0_dp, 0.0_dp, dp), lhs_factor)
     my_rhs_factor = optional_default(cmplx(1.0_dp, 0.0_dp, dp), rhs_factor)

     call check_size('lhs',lhs,(/size(matrix1,m1_r),size(matrix2,m2_c)/),'Matrix_Product')
     if (m2_op == 'N') then
       call check_size('Matrix2',matrix2,(/size(matrix1,m1_c),size(matrix2,m2_c)/),'Matrix_Product')
     else
       call check_size('Matrix2',matrix2,(/size(matrix2,m2_c),size(matrix1,m1_c)/),'Matrix_Product')
     endif

     N = size(lhs,1) !# of rows of lhs
     M = size(lhs,2) !# of columns of rhs
     K = size(matrix1,m1_c) !!shared dimension
     maxd=max(N,M,K) 

     if (use_intrinsic_blas) then
	if (m1_transp) then
	  if (m2_transp) then
	    lhs = lhs*my_lhs_factor + my_rhs_factor*matmul(transpose(matrix1),transpose(matrix2))
	  else if (m2_conjg) then
	    lhs = lhs*my_lhs_factor + my_rhs_factor*matmul(transpose(matrix1),conjg(transpose(matrix2)))
	  else
	    lhs = lhs*my_lhs_factor + my_rhs_factor*matmul(transpose(matrix1),matrix2)
	  endif
	else if (m1_conjg) then
	  if (m2_transp) then
	    lhs = lhs*my_lhs_factor + my_rhs_factor*matmul(conjg(transpose(matrix1)),transpose(matrix2))
	  else if (m2_conjg) then
	    lhs = lhs*my_lhs_factor + my_rhs_factor*matmul(conjg(transpose(matrix1)),conjg(transpose(matrix2)))
	  else
	    lhs = lhs*my_lhs_factor + my_rhs_factor*matmul(conjg(transpose(matrix1)),matrix2)
	  endif
	else
	  if (m2_transp) then
	    lhs = lhs*my_lhs_factor + my_rhs_factor*matmul(matrix1,transpose(matrix2))
	  else if (m2_conjg) then
	    lhs = lhs*my_lhs_factor + my_rhs_factor*matmul(matrix1,conjg(transpose(matrix2)))
	  else
	    lhs = lhs*my_lhs_factor + my_rhs_factor*matmul(matrix1,matrix2)
	  endif
	endif
     else
        call ZGEMM(m1_op,m2_op, N,M,K,my_rhs_factor,matrix1,size(matrix1,1),matrix2,size(matrix2,1),my_lhs_factor,lhs,size(lhs,1))
     endif

   end subroutine matrix_product_sub_zzz

   subroutine matrix_vector_product_sub_ddd(lhs, matrix, vector, m_transpose, &
    lhs_factor, rhs_factor)
      real(dp), intent(out) :: lhs(:)
      real(dp), intent(in) :: matrix(:,:), vector(:)
      logical, intent(in), optional :: m_transpose
      real(dp), intent(in), optional :: lhs_factor, rhs_factor

     integer::M,N,maxd
     character(len=1) :: m_transp
     integer :: m_r, m_c
     real(dp) :: my_lhs_factor, my_rhs_factor

     m_transp = 'N'
     if (present(m_transpose)) then
      if (m_transpose) m_transp = 'T'
     endif

     if (m_transp == 'N') then
       m_r = 1; m_c = 2
     else
       m_r = 2; m_c = 1
     endif

     my_lhs_factor = optional_default(0.0_dp, lhs_factor)
     my_rhs_factor = optional_default(1.0_dp, rhs_factor)

     call check_size('lhs',lhs,(/size(matrix,m_r)/),'Matrix_Product')

     N = size(lhs) !# of rows of lhs
     M = size(matrix,m_c) !!shared dimension
     maxd=max(N,M) 

     if (use_intrinsic_blas) then
        if (m_transp == 'T') then
          lhs=lhs*my_lhs_factor + my_rhs_factor*matmul(transpose(matrix),vector)
        else
          lhs=lhs*my_lhs_factor + my_rhs_factor*matmul(matrix,vector)
        endif
     else
        call DGEMV(m_transp,N,M,my_rhs_factor,matrix,size(matrix,1),vector,1,&
	  my_lhs_factor,lhs,1)
     endif

   end subroutine matrix_vector_product_sub_ddd

  subroutine check_size_complex_dim1(arrayname,realarray,n,caller,error)

    character(*),            intent(in) :: arrayname
    complex(dp), dimension(:),  intent(in) :: realarray 
    integer,  dimension(:),  intent(in) :: n         
    character(*),            intent(in) :: caller
    integer, intent(out), optional :: error     
   
    integer, dimension(:), allocatable :: actual_size 
    logical                            :: failed      
    integer                            :: i           

    !INIT_ERROR(error)
    failed = .false.
    allocate( actual_size( size(shape(realarray)) ) )
    actual_size = shape(realarray)

    if (size(actual_size) /= size(n)) then
       write(line,'(a,i0,a,i0,a)') caller//': '//arrayname//' is ',size(actual_size), &
            ' dimensional and not ',size(n),' dimensional as expected'
!       call print(line)
       failed = .true.
    else
       do i = 1, size(actual_size)
          if (actual_size(i) /= n(i)) then
             write(line,'(3(a,i0),a)') caller//': The size of dimension ',i,' of '//arrayname//' is ', &
                  actual_size(i),' and not ',n(i),' as expected'
!             call print(line)
             failed = .true.
          end if
       end do
    end if

    if (failed) then
    !   RAISE_ERROR(trim(caller)//': Size checking failed', error)
    end if

  end subroutine check_size_complex_dim1

  !overloaded subroutine which allows a scalar 'n' in the one dimensional case
  subroutine check_size_complex_dim1_s(arrayname,realarray,n,caller,error)
    character(*),           intent(in) :: arrayname 
    complex(dp), dimension(:), intent(in) :: realarray  
    integer,                intent(in) :: n         
    character(*),           intent(in) :: caller
    integer, intent(out), optional :: error     
    
    !INIT_ERROR(error)
    call check_size(arrayname,realarray,(/n/),caller)
    !PASS_ERROR(error)

  end subroutine check_size_complex_dim1_s

  subroutine check_size_complex_dim2(arrayname,realarray,n,caller, error)

    character(*),              intent(in) :: arrayname 
    complex(dp), dimension(:,:),  intent(in) :: realarray 
    integer,  dimension(:),    intent(in) :: n        
    character(*),              intent(in) :: caller 
    integer, intent(out), optional :: error     

    integer, dimension(:), allocatable :: actual_size
    logical                            :: failed      
    integer                            :: i          

    !INIT_ERROR(error)
    failed = .false.
    allocate( actual_size( size(shape(realarray)) ) )
    actual_size = shape(realarray)

    if (size(actual_size) /= size(n)) then
       write(line,'(a,i0,a,i0,a)') caller//': '//arrayname//' is ',size(actual_size), &
            ' dimensional and not ',size(n),' dimensional as expected'
!       call print(line)
       failed = .true.
    else
       do i = 1, size(actual_size)
          if (actual_size(i) /= n(i)) then
             write(line,'(3(a,i0),a)') caller//': The size of dimension ',i,' of '//arrayname//' is ', &
                  actual_size(i),' and not ',n(i),' as expected'
!             call print(line)
             failed = .true.
          end if
       end do
    end if

    if (failed) then
     !  RAISE_ERROR(trim(caller)//': Size checking failed', error)
    end if

  end subroutine check_size_complex_dim2

  subroutine check_size_log_dim1(arrayname,logarray,n,caller,error)

    character(*),           intent(in) :: arrayname 
    logical, dimension(:),  intent(in) :: logarray  
    integer, dimension(:),  intent(in) :: n         
    character(*),           intent(in) :: caller
    integer, intent(out), optional :: error     
    
    integer, dimension(:), allocatable :: actual_size
    logical                            :: failed     
    integer                            :: i          

    !INIT_ERROR(error)
    failed = .false.
    allocate( actual_size( size(shape(logarray)) ) )
    actual_size = shape(logarray)

    if (size(actual_size) /= size(n)) then
       write(line,'(a,i0,a,i0,a)') caller//': '//arrayname//' is ',size(actual_size), &
            ' dimensional and not ',size(n),' dimensional as expected'
!       call print(line)
       failed = .true.
    else
       do i = 1, size(actual_size)
          if (actual_size(i) /= n(i)) then
             write(line,'(3(a,i0),a)') caller//': The size of dimension ',i,' of '//arrayname//' is ', &
                  actual_size(i),' and not ',n(i),' as expected'
!             call print(line)
             failed = .true.
          end if
       end do
    end if

    if (failed) then
     !  RAISE_ERROR(trim(caller)//': Size checking failed', error)
    end if

  end subroutine check_size_log_dim1

  !overloaded subroutine which allows a scalar 'n' in the one dimensional case
  subroutine check_size_log_dim1_s(arrayname,logarray,n,caller,error)
    character(*),           intent(in) :: arrayname 
    logical, dimension(:),  intent(in) :: logarray  
    integer,                intent(in) :: n         
    character(*),           intent(in) :: caller
    integer, intent(out), optional :: error     
    
    !INIT_ERROR(error)
    call check_size(arrayname,logarray,(/n/),caller)
    !PASS_ERROR(error)

  end subroutine check_size_log_dim1_s

  subroutine check_size_log_dim2(arrayname,logarray,n,caller,error)

    character(*),             intent(in) :: arrayname 
    logical, dimension(:,:),  intent(in) :: logarray  
    integer, dimension(:),    intent(in) :: n         
    character(*),             intent(in) :: caller    
    integer, intent(out), optional :: error     
    
    integer, dimension(:), allocatable :: actual_size
    logical                            :: failed     
    integer                            :: i          

    !INIT_ERROR(error)
    failed = .false.
    allocate( actual_size( size(shape(logarray)) ) )
    actual_size = shape(logarray)

    if (size(actual_size) /= size(n)) then
       write(line,'(a,i0,a,i0,a)') caller//': '//arrayname//' is ',size(actual_size), &
            ' dimensional and not ',size(n),' dimensional as expected'
!       call print(line)
       failed = .true.
    else
       do i = 1, size(actual_size)
          if (actual_size(i) /= n(i)) then
             write(line,'(3(a,i0),a)') caller//': The size of dimension ',i,' of '//arrayname//' is ', &
                  actual_size(i),' and not ',n(i),' as expected'
!             call print(line)
             failed = .true.
          end if
       end do
    end if

    if (failed) then
     !  RAISE_ERROR(trim(caller)//': Size checking failed', error)
    end if

  end subroutine check_size_log_dim2  
end module

