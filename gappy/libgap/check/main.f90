program test
integer,parameter  :: na = 64
integer,parameter  :: nf = 66
real(8)   :: lat(3,3)
real(8)   :: pos(na,3)
INTEGER   :: elements(na)
REAL(8)   :: xx(nf, na), dxdy(nf, na, na, 3) , strs(3,3, na, nf)
REAL(8)   :: rcut
logical   :: lgrad
open(111,file='POSCAR')
read(111, *)
read(111, *)
do i= 1,3
    read(111,*) lat(i,:)
enddo
read(111,*)
read(111,*)
read(111,*)
do i= 1,na
    read(111,*) pos(i,:)
enddo
rcut = 6.0
lgrad= .false.
do i= 1, 8
    elements(i) = 5
enddo
do i= 9, 64
    elements(i) =6
enddo
call fcar2wacsf(na, nf, lat, elements, pos, xx, dxdy, strs, rcut, lgrad)
end
