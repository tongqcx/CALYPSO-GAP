MODULE GAP_INIT
integer, parameter                :: natoms_max = 300
integer, parameter                :: nsf_max = 100
integer, parameter                :: nsparseX_max = 1000
double precision,parameter        :: delta = 1.0
double precision,parameter        :: GPa2eVPang =6.24219D-3
type Structure
    integer                       :: na
    integer                       :: ntype
    double precision              :: types(natoms_max)
    double precision              :: lat(3,3)
    double precision              :: volume
    double precision              :: xyz(natoms_max,3)
    double precision              :: xx(nsf_max,natoms_max)
    double precision              :: kk(natoms_max,nsf_max)
    double precision              :: ckm(natoms_max,nsparseX_max)  ! k in the number of atoms
    double precision              :: e(natoms_max)
    double precision              :: covf(natoms_max)
    double precision              :: energy_ref
    double precision              :: energy_cal
    double precision              :: stress_ref(6)
    double precision              :: stress_cal(6)
    double precision              :: force_ref(natoms_max,3)
    double precision              :: force_cal(natoms_max,3)
    double precision              :: dxdy(nsf_max,natoms_max,natoms_max,3)
    double precision              :: strs(3,3,nsf_max,natoms_max)
    double precision              :: dedg(natoms_max,nsf_max)
    double precision              :: dkdg(natoms_max,nsf_max,nsparseX_max)
end type

contains
SUBROUTINE get_ele_weights(cc,nw)!{{{
implicit none
integer           ::  cc
double precision  ::  nw
select case(cc)
case (1)
    nw = -1.0
case (3)
    nw = 1.0
case (5)
    nw = -1.0
case (6)
    nw = 1.0
case (8)
    nw = 2.0
case (12)
    nw = -1.0
case (13)
    nw = 1.0
case (15)
    nw = 2.0
case (16)
    nw = -1.0
case (14)
    nw = -2.0
case (20)
    nw = -1.0
case (28)
    nw = -1.0
case (55)
    nw = 1.0
end select
END SUBROUTINE!}}}
SUBROUTINE GET_COV(n,m,xx,yy,theta, cov)!{{{
    integer i,j,k
    integer,intent(in)       ::  n,m
    double precision, intent(in)       ::  xx(:,:)
    double precision, intent(in)       ::  yy(:,:)
    double precision, intent(in)       ::  theta(:)
    double precision, intent(out)      ::  cov(n,m)
    double precision                   ::  temp
    integer                  :: nf
    nf = size(xx,2)
    temp = 0.d0
    do i = 1,n
        do  j=1,m 
            temp = 0.d0
            do k = 1,nf
                temp = temp+((xx(i,k)-yy(j,k))/theta(k))**2
            enddo
            cov(i,j) = delta*exp(-0.5d0*temp)
        enddo
    enddo
END SUBROUTINE!}}}
FUNCTION Det(matrix)!{{{
   double precision,  intent(in) :: Matrix(3,3)
   double precision              :: Det


   Det = Matrix(1,1)*(Matrix(2,2)*Matrix(3,3)-Matrix(2,3)*Matrix(3,2))&
      -Matrix(1,2)*(Matrix(2,1)*Matrix(3,3)-Matrix(2,3)*Matrix(3,1))&
      +Matrix(1,3)*(Matrix(2,1)*Matrix(3,2)-Matrix(3,1)*Matrix(2,2))

END FUNCTION!}}}
END MODULE

