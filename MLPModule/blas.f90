  function dasum(n,dx,incx)
!
!     takes the sum of the absolute values.
!     jack dongarra, linpack, 3/11/78.
!
  !use datatypes
  implicit none
  integer, parameter :: dp  = kind(1.0d0)
  integer, parameter :: i4  = selected_int_kind(5)
  real(dp) dx(*),dtemp,dasum
  integer(i4) i,incx,m,mp1,n,nincx
!
  dasum = 0.0_dp
  dtemp = 0.0_dp
  if (n.le.0) return
  if (incx.eq.1) go to 20
!
!        code for increment not equal to 1
!
  nincx = n*incx
  do 10 i = 1,nincx,incx
    dtemp = dtemp + dabs(dx(i))
10 continue
  dasum = dtemp
  return
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
20 m = mod(n,6_i4)
  if( m .eq. 0 ) go to 40
  do 30 i = 1,m
    dtemp = dtemp + dabs(dx(i))
30 continue
  if( n .lt. 6 ) go to 60
40 mp1 = m + 1
  do 50 i = mp1,n,6
    dtemp = dtemp + dabs(dx(i)) + dabs(dx(i + 1)) + dabs(dx(i + 2)) &
    + dabs(dx(i + 3)) + dabs(dx(i + 4)) + dabs(dx(i + 5))
50 continue
60 dasum = dtemp
  return
  end
  subroutine daxpy(n,da,dx,incx,dy,incy)
!
!     constant times a vector plus a vector.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!
  !use datatypes
  implicit none
  integer, parameter :: dp  = kind(1.0d0)
  integer, parameter :: i4  = selected_int_kind(5)
  real(dp) dx(*),dy(*),da
  integer(i4) i,incx,incy,ix,iy,m,mp1,n
!
  if (n.le.0) return
  if (da .eq. 0.0_dp) return
  if (incx.eq.1.and.incy.eq.1) go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
  ix = 1
  iy = 1
  if(incx.lt.0)ix = (-n+1)*incx + 1
  if(incy.lt.0)iy = (-n+1)*incy + 1
  do 10 i = 1,n
    dy(iy) = dy(iy) + da*dx(ix)
    ix = ix + incx
    iy = iy + incy
10 continue
  return
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
20 m = mod(n,4_i4)
  if( m .eq. 0 ) go to 40
  do 30 i = 1,m
    dy(i) = dy(i) + da*dx(i)
30 continue
  if( n .lt. 4 ) return
40 mp1 = m + 1
  do 50 i = mp1,n,4
    dy(i) = dy(i) + da*dx(i)
    dy(i + 1) = dy(i + 1) + da*dx(i + 1)
    dy(i + 2) = dy(i + 2) + da*dx(i + 2)
    dy(i + 3) = dy(i + 3) + da*dx(i + 3)
50 continue
  return
  end
  subroutine dcopy(n,dx,incx,dy,incy)
!
!     copies a vector, x, to a vector, y.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!
  !use datatypes
  implicit none
  integer, parameter :: dp  = kind(1.0d0)
  integer, parameter :: i4  = selected_int_kind(5)
  real(dp) dx(*),dy(*)
  integer(i4) i,incx,incy,ix,iy,m,mp1,n
!
  if(n.le.0)return
  if(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
  ix = 1
  iy = 1
  if(incx.lt.0)ix = (-n+1)*incx + 1
  if(incy.lt.0)iy = (-n+1)*incy + 1
  do 10 i = 1,n
    dy(iy) = dx(ix)
    ix = ix + incx
    iy = iy + incy
10 continue
  return
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
20 m = mod(n,7_i4)
  if( m .eq. 0 ) go to 40
  do 30 i = 1,m
    dy(i) = dx(i)
30 continue
  if( n .lt. 7 ) return
40 mp1 = m + 1
  do 50 i = mp1,n,7
    dy(i) = dx(i)
    dy(i + 1) = dx(i + 1)
    dy(i + 2) = dx(i + 2)
    dy(i + 3) = dx(i + 3)
    dy(i + 4) = dx(i + 4)
    dy(i + 5) = dx(i + 5)
    dy(i + 6) = dx(i + 6)
50 continue
  return
  end
  function ddot(n,dx,incx,dy,incy)
!
!     forms the dot product of two vectors.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!
  !use datatypes
  implicit none
  integer, parameter :: dp  = kind(1.0d0)
  integer, parameter :: i4  = selected_int_kind(5)
  real(dp) dx(*),dy(*),dtemp,ddot
  integer(i4) i,incx,incy,ix,iy,m,mp1,n
!
  ddot = 0.0_dp
  dtemp = 0.0_dp
  if(n.le.0)return
  if(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
  ix = 1
  iy = 1
  if(incx.lt.0)ix = (-n+1)*incx + 1
  if(incy.lt.0)iy = (-n+1)*incy + 1
  do 10 i = 1,n
    dtemp = dtemp + dx(ix)*dy(iy)
    ix = ix + incx
    iy = iy + incy
10 continue
  ddot = dtemp
  return
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
20 m = mod(n,5_i4)
  if( m .eq. 0 ) go to 40
  do 30 i = 1,m
    dtemp = dtemp + dx(i)*dy(i)
30 continue
  if( n .lt. 5 ) go to 60
40 mp1 = m + 1
  do 50 i = mp1,n,5
    dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) + &
     dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
50 continue
60 ddot = dtemp
  return
  end
  subroutine  dscal(n,da,dx,incx)
!
!     scales a vector by a constant.
!     uses unrolled loops for increment equal to one.
!     jack dongarra, linpack, 3/11/78.
!
  !use datatypes
  implicit none
  integer, parameter :: dp  = kind(1.0d0)
  integer, parameter :: i4  = selected_int_kind(5)
  real(dp) da,dx(*)
  integer(i4) i,incx,m,mp1,n,nincx
!
  if(n.le.0)return
  if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
  nincx = n*incx
  do 10 i = 1,nincx,incx
    dx(i) = da*dx(i)
10 continue
  return
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
20 m = mod(n,5_i4)
  if( m .eq. 0 ) go to 40
  do 30 i = 1,m
    dx(i) = da*dx(i)
30 continue
  if( n .lt. 5 ) return
40 mp1 = m + 1
  do 50 i = mp1,n,5
    dx(i) = da*dx(i)
    dx(i + 1) = da*dx(i + 1)
    dx(i + 2) = da*dx(i + 2)
    dx(i + 3) = da*dx(i + 3)
    dx(i + 4) = da*dx(i + 4)
50 continue
  return
  end
  subroutine  dswap (n,dx,incx,dy,incy)
!
!     interchanges two vectors.
!     uses unrolled loops for increments equal one.
!     jack dongarra, linpack, 3/11/78.
!
  !use datatypes
  implicit none
  integer, parameter :: dp  = kind(1.0d0)
  integer, parameter :: i4  = selected_int_kind(5)
  real(dp) dx(*),dy(*),dtemp
  integer(i4) i,incx,incy,ix,iy,m,mp1,n
!
  if(n.le.0)return
  if(incx.eq.1.and.incy.eq.1)go to 20
!
!       code for unequal increments or equal increments not equal
!         to 1
!
  ix = 1
  iy = 1
  if(incx.lt.0)ix = (-n+1)*incx + 1
  if(incy.lt.0)iy = (-n+1)*incy + 1
  do 10 i = 1,n
    dtemp = dx(ix)
    dx(ix) = dy(iy)
    dy(iy) = dtemp
    ix = ix + incx
    iy = iy + incy
10 continue
  return
!
!       code for both increments equal to 1
!
!
!       clean-up loop
!
20 m = mod(n,3_i4)
  if( m .eq. 0 ) go to 40
  do 30 i = 1,m
    dtemp = dx(i)
    dx(i) = dy(i)
    dy(i) = dtemp
30 continue
  if( n .lt. 3 ) return
40 mp1 = m + 1
  do 50 i = mp1,n,3
    dtemp = dx(i)
    dx(i) = dy(i)
    dy(i) = dtemp
    dtemp = dx(i + 1)
    dx(i + 1) = dy(i + 1)
    dy(i + 1) = dtemp
    dtemp = dx(i + 2)
    dx(i + 2) = dy(i + 2)
    dy(i + 2) = dtemp
50 continue
  return
  end
  subroutine zaxpy(n,za,zx,incx,zy,incy)
!
!     constant times a vector plus a vector.
!     jack dongarra, 3/11/78.
!
  !use datatypes
  implicit none
  integer, parameter :: dp  = kind(1.0d0)
  integer, parameter :: i4  = selected_int_kind(5)
  integer, parameter :: dpc = kind((1.0d0,1.0d0))
  complex(dpc) :: zx(*),zy(*),za
  real(dp)     :: dcabs1
  integer(i4)  :: n,incx,incy,ix,iy,i
!
  if (n.le.0)return
  if (dcabs1(za) .eq. 0.0_dp) return
  if (incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
  ix = 1
  iy = 1
  if(incx.lt.0)ix = (-n+1)*incx + 1
  if(incy.lt.0)iy = (-n+1)*incy + 1
  do 10 i = 1,n
    zy(iy) = zy(iy) + za*zx(ix)
    ix = ix + incx
    iy = iy + incy
10 continue
  return
!
!        code for both increments equal to 1
!
20 do 30 i = 1,n
    zy(i) = zy(i) + za*zx(i)
30 continue
  return
  end
  subroutine zcopy(n,zx,incx,zy,incy)
!
!     copies a vector, x, to a vector, y.
!     jack dongarra, linpack, 4/11/78.
!
  !use datatypes
  implicit none
  integer, parameter :: dp  = kind(1.0d0)
  integer, parameter :: i4  = selected_int_kind(5)
  integer, parameter :: dpc = kind((1.0d0,1.0d0))

  integer(i4) incx,incy,ix,iy,n,i
  complex(dpc) zx(*),zy(*)
!
  if (n.le.0) return
  if (incx.eq.1.and.incy.eq.1) go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
  ix = 1
  iy = 1
  if (incx.lt.0) ix = (-n+1)*incx + 1
  if (incy.lt.0) iy = (-n+1)*incy + 1
  do 10 i = 1,n
    zy(iy) = zx(ix)
    ix = ix + incx
    iy = iy + incy
10 continue
  return
!
!        code for both increments equal to 1
!
20 do 30 i = 1,n
    zy(i) = zx(i)
30 continue
  return
  end
  function zdotc(n,zx,incx,zy,incy)
!
!     forms the dot product of a vector.
!     jack dongarra, 3/11/78.
!
  !use datatypes
  implicit none
  integer, parameter :: dp  = kind(1.0d0)
  integer, parameter :: i4  = selected_int_kind(5)
  integer, parameter :: dpc = kind((1.0d0,1.0d0))

  complex(dpc) zx(*),zy(*),ztemp,zdotc
  real(dp)     conjg
  integer(i4)  n,incx,incy,ix,iy,i
!
  ztemp = (0.0_dp,0.0_dp)
  zdotc = (0.0_dp,0.0_dp)
  if(n.le.0)return
  if(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
  ix = 1
  iy = 1
  if(incx.lt.0)ix = (-n+1)*incx + 1
  if(incy.lt.0)iy = (-n+1)*incy + 1
  do 10 i = 1,n
    ztemp = ztemp + conjg(zx(ix))*zy(iy)
    ix = ix + incx
    iy = iy + incy
10 continue
  zdotc = ztemp
  return
!
!        code for both increments equal to 1
!
20 do 30 i = 1,n
    ztemp = ztemp + conjg(zx(i))*zy(i)
30 continue
  zdotc = ztemp
  return
  end
  subroutine  zswap(n,zx,incx,zy,incy)
!
!     interchanges two vectors.
!     jack dongarra, 3/11/78.
!
  !use datatypes
  implicit none
  integer, parameter :: dp  = kind(1.0d0)
  integer, parameter :: i4  = selected_int_kind(5)
  integer, parameter :: dpc = kind((1.0d0,1.0d0))

  complex(dpc) zx(*),zy(*),ztemp
  integer(i4)  incx,incy,n,ix,iy,i
!
  if(n.le.0)return
  if(incx.eq.1.and.incy.eq.1)go to 20
!
!       code for unequal increments or equal increments not equal
!         to 1
!
  ix = 1
  iy = 1
  if(incx.lt.0)ix = (-n+1)*incx + 1
  if(incy.lt.0)iy = (-n+1)*incy + 1
  do 10 i = 1,n
    ztemp = zx(ix)
    zx(ix) = zy(iy)
    zy(iy) = ztemp
    ix = ix + incx
    iy = iy + incy
10 continue
  return
!
!       code for both increments equal to 1
20 do 30 i = 1,n
    ztemp = zx(i)
    zx(i) = zy(i)
    zy(i) = ztemp
30 continue
  return
  end
  function idamax(n,dx,incx)
!
!     finds the index of element having max. absolute value.
!     jack dongarra, linpack, 3/11/78.
!
  !use datatypes
  implicit none
  integer, parameter :: dp  = kind(1.0d0)
  integer, parameter :: i4  = selected_int_kind(5)
  real(dp) dx(*),dmax
  integer(i4) i,incx,ix,n,idamax
!
  idamax = 0
  if( n .lt. 1 ) return
  idamax = 1
  if(n.eq.1)return
  if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
  ix = 1
  dmax = dabs(dx(1))
  ix = ix + incx
  do 10 i = 2,n
     if(dabs(dx(ix)).le.dmax) go to 5
     idamax = i
     dmax = dabs(dx(ix))
5    ix = ix + incx
10 continue
  return
!
!        code for increment equal to 1
!
20 dmax = dabs(dx(1))
  do 30 i = 2,n
     if(dabs(dx(i)).le.dmax) go to 30
     idamax = i
     dmax = dabs(dx(i))
30 continue
  return
  end
  function dcabs1(z)
  !use datatypes
  implicit none
  integer, parameter :: dp  = kind(1.0d0)
  integer, parameter :: i4  = selected_int_kind(5)
  integer, parameter :: dpc = kind((1.0d0,1.0d0))

  complex(dpc) z,zz
  real(dp) t(2),dcabs1
  equivalence (zz,t(1))
  zz = z
  dcabs1 = dabs(t(1)) + dabs(t(2))
  return
  end
  function izamax(n,zx,incx)
!
!     finds the index of element having max. absolute value.
!     jack dongarra, 1/15/85.
!
  !use datatypes
  implicit none
  integer, parameter :: dp  = kind(1.0d0)
  integer, parameter :: i4  = selected_int_kind(5)
  integer, parameter :: dpc = kind((1.0d0,1.0d0))

  complex(dpc) zx(*)
  real(dp) smax
  integer(i4) i,incx,ix,n,izamax
  real(dp) dcabs1
!
  izamax = 0
  if(n.lt.1)return
  izamax = 1
  if(n.eq.1)return
  if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
  ix = 1
  smax = dcabs1(zx(1))
  ix = ix + incx
  do 10 i = 2,n
     if(dcabs1(zx(ix)).le.smax) go to 5
     izamax = i
     smax = dcabs1(zx(ix))
5    ix = ix + incx
10 continue
  return
!
!        code for increment equal to 1
!
20 smax = dcabs1(zx(1))
  do 30 i = 2,n
     if(dcabs1(zx(i)).le.smax) go to 30
     izamax = i
     smax = dcabs1(zx(i))
30 continue
  return
  end
  subroutine matrix_inversion(a,ia,n,wrk,ifail)
!
!  Matrix inverter
!
!  On entry :
!
!  a     = matrix to be inverted
!  ia    = lower dimension of a
!  n     = actual size of matrix to be inverted
!  wrk   = workspace array of length 2*n
!
!  On exit :
!
!  a     = inverse matrix
!  ifail = 0, if OK
!
!   3/14 Renamed from matinv for benefit of ChemShell
!
!  use times
  implicit none
  integer, parameter :: dp  = kind(1.0d0)
  integer, parameter :: i4  = selected_int_kind(5)
!
!  Passed variables
!
  integer(i4)    :: ia
  integer(i4)    :: ifail
  integer(i4)    :: n
  real(dp)       :: a(ia,*)
  real(dp)       :: wrk(*)
!
!  Local variables
!
  integer(i4)    :: i
  integer(i4)    :: i1
  integer(i4)    :: ins
  integer(i4)    :: j
  integer(i4)    :: jj1
  integer(i4)    :: k
  integer(i4)    :: kk
  integer(i4)    :: l
  integer(i4)    :: l1
  integer(i4)    :: m
  integer(i4)    :: nupper
  real(dp)       :: acc
  real(dp)       :: aloc
  real(dp)       :: g_cpu_time
  real(dp)       :: t
  real(dp)       :: t1
  real(dp)       :: t2
!
!  t1 = g_cpu_time()
  acc = 1.0d-8
  ifail = 0
  do j = 1,n
    if (j.ne.1) then
      call mxm1(a(j,1),ia,a(1,j),ia,wrk,ia,n-j+1_i4,j-1_i4)
      nupper = n - j + 1
      do l = 1,nupper
        a(j+l-1,j) = a(j+l-1,j) + wrk(l)
      enddo
    endif
    t = abs(a(j,j))
    k = j
    if (j.ne.n) then
      do i = j+1,n
        if (abs(a(i,j)).gt.t) then
          t = abs(a(i,j))
          k = i
        endif
      enddo
    endif
    wrk(j+n) = dble(k)
    if (t.le.acc) then
      ifail = 1
!      t2 = g_cpu_time()
!      tmati = tmati + t2 - t1
      return
    endif
    if (k.ne.j) then
      do m = 1,n
        t = a(j,m)
        a(j,m) = a(k,m)
        a(k,m) = t
      enddo
    endif
    a(j,j) = 1.0_dp/a(j,j)
    if (j.ne.n) then
      if (j.ne.1) then
        call mxm2(a(1,j+1),ia,a(j,1),ia,wrk,n-j,j-1_i4)
        nupper = n - j
        do l1 = 1,nupper
          a(j,j+l1) = a(j,j+l1) + wrk(l1)
        enddo
      endif
      t = - a(j,j)
      nupper = n - (j+1)
      do i1 = j+1,n
        a(j,i1) = t*a(j,i1)
      enddo
    endif
  enddo
!
!  Use cminv method to solve for a**-1
!
  do k = 2,n
    nupper = k - 1
    do m = 1,nupper
      wrk(m) = 0.0_dp
    enddo
    do j = 1,k-1
      aloc = a(k,j)
      do m = 1,j
        wrk(m) = wrk(m)-aloc*a(j,m)
      enddo
    enddo
    aloc = a(k,k)
    nupper = k - 1
    do m = 1,nupper
      a(k,m) = wrk(m)*aloc
    enddo
  enddo
!
!  Now back substitution
!
  k = n
  do kk = 2,n
    k = k - 1
    jj1 = kk - 1
    call mxm2(a(k+1,1),ia,a(k,k+1),ia,wrk,n,jj1)
    do l = 1,k
      wrk(l) = wrk(l)+a(k,l)
    enddo
    do j = 1,n
      a(k,j) = wrk(j)
    enddo
  enddo
!
!  Multiply solution by inverse of permutation matrix
!
  k = n + 1
  do i = 1,n
    k = k-1
    ins = int(wrk(k+n))
    if (ins.ne.k) then
      do j = 1,n
        wrk(j) = a(j,ins)
        a(j,ins) = a(j,k)
        a(j,k) = wrk(j)
      enddo
    endif
  enddo
!  t2 = g_cpu_time()
!  tmati = tmati + t2 - t1
  return
  end
!******************************
!  Ancillary matrix routines  *
!******************************
  subroutine mxm1(rarr,ir,sarr,jr,tarr,kr,id1,id2)
!
!  Matrix multiplier
!
!  use datatypes
  implicit none
  integer, parameter :: dp  = kind(1.0d0)
  integer, parameter :: i4  = selected_int_kind(5)
!
!  Passed variables
!
  integer(i4)    :: id1
  integer(i4)    :: id2
  integer(i4)    :: ir
  integer(i4)    :: jr
  integer(i4)    :: kr
  real(dp)       :: rarr(*)
  real(dp)       :: sarr(*)
  real(dp)       :: tarr(*)
!
!  Local variables
!
  integer(i4)    :: i
  integer(i4)    :: ia
  integer(i4)    :: j
  real(dp)       :: sum
!
  do i = 1,id1
    ia = i - ir
    sum = 0.0_dp
    do j = 1,id2
      ia = ia + ir
      sum = sum + rarr(ia)*sarr(j)
    enddo
    tarr(i) = sum
  enddo
  return
  end
  subroutine mxm2(rarr,ic,sarr,jc,tarr,id1,id2)
!
!  Matrix multiplier
!
!  use datatypes
  implicit none
  integer, parameter :: dp  = kind(1.0d0)
  integer, parameter :: i4  = selected_int_kind(5)
!
!  Passed variables
!
  integer(i4)    :: id1
  integer(i4)    :: id2
  integer(i4)    :: ic
  integer(i4)    :: jc
  real(dp)       :: rarr(*)
  real(dp)       :: sarr(*)
  real(dp)       :: tarr(*)
!
!  Local variables
!
  integer(i4)    :: i
  integer(i4)    :: ia
  integer(i4)    :: ira
  integer(i4)    :: j
  integer(i4)    :: ja
  real(dp)       :: sum
!
!
  ira = 1 - ic
  do i = 1,id1
    ira = ira + ic
    ia = ira - 1
    ja = 1 - jc
    sum = 0.0_dp
    do j = 1,id2
      ja = ja + jc
      sum = sum + rarr(ia+j)*sarr(ja)
    enddo
    tarr(i) = sum
  enddo
  return
  end
