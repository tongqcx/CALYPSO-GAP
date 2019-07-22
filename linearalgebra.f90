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
  interface check_size
     module procedure check_size_int_dim1, check_size_int_dim1_s, check_size_int_dim2  
     module procedure check_size_real_dim1, check_size_real_dim1_s, check_size_real_dim2, check_size_real_dim3
  end interface check_size
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

end module

