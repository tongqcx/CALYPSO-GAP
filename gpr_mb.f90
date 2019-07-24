MODULE GPR_MB
use constants
use io
use math
use linearalgebra
use struct

CONTAINS
SUBROUTINE INI_GAP_MB(GAP, ACSF, nsparse, nobf)
type(GAP_type),intent(inout)           :: GAP
type(ACSF_type),intent(in)             :: acsf
integer,intent(in)                     :: nsparse, nobf

!local

call READ_ACSF('neural.in', ACSF)
GAP%nsparse = nsparse
GAP%dd = ACSF%NSF * 2     ! D_tot = D_topology + D_species
GAP%nglobalY = nobf

allocate(GAP%cmm(GAP%nsparse, GAP%nsparse))
allocate(GAP%cmo(GAP%nsparse, GAP%nglobalY, 1))
allocate(GAP%sparseX(GAP%nsparse, GAP%dd))
allocate(GAP%obe(GAP%nglobalY))
allocate(GAP%coeff(GAP%nsparse, 1))
allocate(GAP%lamda(GAP%nglobalY))
allocate(GAP%lamdaobe(GAP%nglobalY))

END SUBROUTINE INI_GAP_MB

SUBROUTINE CAR2ACSF(at, GAP, ACSF)
type(Structure),intent(inout)          :: at
type(GAP_type),intent(in)              :: GAP
type(ACSF_type),intent(in)             :: acsf

!local
INTEGER                                :: i_nsf, i_atom

allocate(at%xx(GAP%dd, at%natoms))
allocate(at%dxdy(GAP%dd, at%natoms, at%natoms, 3))
allocate(at%strs(3, 3, GAP%dd, at%natoms))

nnn = ACSF%nsf
do ii = 1, nnn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!G1 = SUM_j{exp(-alpha*rij**2)fc(rij)}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (ACSF%sf(i)%types.eq.1) then
        do i = 1, at%natoms
            do i_type = 1, nspecies
                do i_neighbor = 1, at%atom(i)%count(i_type)
                    rij = at%atom(i)%neighbor(i_type,i_neighbor,4)
                    xyz = at%atom(i)%neighbor(i_type,i_neighbor,1:3)
                    n = int(at%atom(i)%neighbor(i_type,i_neighbor,5))
                    cutoff = ACSF%sf(ii)%cutoff
                    alpha = ACSF%sf(ii)%alpha
                    weights = at%mlp_weights(n)
                    if (rij.le.cutoff) then
                        deltaxj = -1.d0*(at%atom(i)%pos(1) - xyz(1))
                        deltayj = -1.d0*(at%atom(i)%pos(2) - xyz(2))
                        deltazj = -1.d0*(at%atom(i)%pos(3) - xyz(3))
                        drijdxi = -1.d0*deltaxj/rij
                        drijdyi = -1.d0*deltayj/rij
                        drijdzi = -1.d0*deltazj/rij
                        drijdxj = -1.d0*drijdxi
                        drijdyj = -1.d0*drijdyi
                        drijdzj = -1.d0*drijdzi
                        fcutij = 0.5d0 * (dcos(pi*rij/cutoff) + 1.d0)
                        temp1=0.5d0*(-dsin(pi*rij/cutoff))*(pi/cutoff)
                        dfcutijdxi=temp1*drijdxi
                        dfcutijdyi=temp1*drijdyi
                        dfcutijdzi=temp1*drijdzi
                        dfcutijdxj=-1.d0*dfcutijdxi
                        dfcutijdyj=-1.d0*dfcutijdyi
                        dfcutijdzj=-1.d0*dfcutijdzi
                        at%xx(ii,i) = at%xx(ii,i) + dexp(-1.d0*alpha*rij**2)*fcutij
                        at%xx(ii + nnn, i) = at%xx(ii + nnn, i) + dexp(-1.d0*alpha*rij**2)*fcutij * weights !!!!!!! 

                        !dxx/dx
                        temp1=-2.d0*alpha*rij*dexp(-1.d0*alpha*rij**2)*fcutij
                        temp2= dexp(-1.d0*alpha*rij**2)

                        dxdy(ii,i,i,1)=dxdy(ii,i,i,1)+(drijdxi*temp1 + temp2*dfcutijdxi)
                        dxdy(ii+nnn,i,i,1)=dxdy(ii+nnn,i,i,1) + &
                        (drijdxi*temp1+ temp2*dfcutijdxi) * weights
                
                        temp3=drijdxj*temp1 + temp2*dfcutijdxj
                        dxdy(ii,i,n,1)=dxdy(ii,i,n,1)+temp3
                
                        temp4=temp3*weights
                        dxdy(ii + nnn,i,n,1)=dxdy(ii + nnn,i,n,1)+temp4
                
                        strs(1,1,ii,i)=strs(1,1,ii,i)+deltaxj*temp3
                        strs(2,1,ii,i)=strs(2,1,ii,i)+deltayj*temp3
                        strs(3,1,ii,i)=strs(3,1,ii,i)+deltazj*temp3
                
                        strs(1,1,ii+nnn,i)=strs(1,1,ii+nnn,i)+deltaxj*temp4
                        strs(2,1,ii+nnn,i)=strs(2,1,ii+nnn,i)+deltayj*temp4
                        strs(3,1,ii+nnn,i)=strs(3,1,ii+nnn,i)+deltazj*temp4
                        !dxx/dy
                        dxdy(ii,i,i,2)=dxdy(ii,i,i,2)+(drijdyi*temp1+temp2*dfcutijdyi)
                        dxdy(ii+nnn,i,i,2)=dxdy(ii+nnn,i,i,2)+(drijdyi*temp1+temp2*dfcutijdyi)*weights
                        temp3= drijdyj*temp1 + temp2*dfcutijdyj
                        dxdy(ii,i,n,2)=dxdy(ii,i,n,2)+temp3
                        temp4=temp3 * weights
                        dxdy(ii + nnn,i,n,2)=dxdy(ii + nnn ,i,n,2)+temp4
                
                        strs(1,2,ii,i)=strs(1,2,ii,i)+deltaxj*temp3
                        strs(2,2,ii,i)=strs(2,2,ii,i)+deltayj*temp3
                        strs(3,2,ii,i)=strs(3,2,ii,i)+deltazj*temp3
                
                        strs(1,2,ii + nnn,i)=strs(1,2,ii + nnn,i)+deltaxj*temp4
                        strs(2,2,ii + nnn,i)=strs(2,2,ii + nnn,i)+deltayj*temp4
                        strs(3,2,ii + nnn,i)=strs(3,2,ii + nnn,i)+deltazj*temp4
                        !dxx/dz
                        dxdy(ii,i,i,3)=dxdy(ii,i,i,3)+&
                               (drijdzi*temp1&
                              + temp2*dfcutijdzi)
                        dxdy(ii + nnn,i,i,3)=dxdy(ii + nnn,i,i,3)+&
                               (drijdzi*temp1&
                              + temp2*dfcutijdzi)*types(n)
                        temp3=drijdzj*temp1 + temp2*dfcutijdzj
                        dxdy(ii,i,n,3)=dxdy(ii,i,n,3)+temp3
                        temp4=temp3*types(n)
                        dxdy(ii + nnn,i,n,3)=dxdy(ii + nnn,i,n,3)+temp4
                
                        strs(1,3,ii,i)=strs(1,3,ii,i)+deltaxj*temp3
                        strs(2,3,ii,i)=strs(2,3,ii,i)+deltayj*temp3
                        strs(3,3,ii,i)=strs(3,3,ii,i)+deltazj*temp3
                
                        strs(1,3,ii + nnn,i)=strs(1,3,ii + nnn,i)+deltaxj*temp4
                        strs(2,3,ii + nnn,i)=strs(2,3,ii + nnn,i)+deltayj*temp4
                        strs(3,3,ii + nnn,i)=strs(3,3,ii + nnn,i)+deltazj*temp4
                    endif ! rij .le. cutoff
                enddo ! i_neighbor
            enddo ! i_type
        enddo ! i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!G1 = SUM_j{exp(-alpha*rij**2)fc(rij)}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    elseif (ACSF%sf(i)%types.eq.2) then
        do i = 1, at%natoms
            do j_type = 1, nspecies
                do j_neighbor = 1, at%atom(i)%count(j_type)
                    rij = at%atom(i)%neighbor(j_type,j_neighbor,4)
                    xyz_j = at%atom(i)%neighbor(j_type,j_neighbor,1:3)
                    n = int(at%atom(i)%neighbor(j_type,j_neighbor,5))
                    cutoff = ACSF%sf(ii)%cutoff
                    alpha = ACSF%sf(ii)%alpha
                    weights = at%mlp_weights(n)
                    if (rij.le.cutoff) then
                        deltaxj = -1.d0*(at%atom(i)%pos(1) - xyz_j(1))
                        deltayj = -1.d0*(at%atom(i)%pos(2) - xyz_j(2))
                        deltazj = -1.d0*(at%atom(i)%pos(3) - xyz_j(3))
                        drijdxi = -1.d0*deltaxj/rij
                        drijdyi = -1.d0*deltayj/rij
                        drijdzi = -1.d0*deltazj/rij
                        drijdxj = -1.d0*drijdxi
                        drijdyj = -1.d0*drijdyi
                        drijdzj = -1.d0*drijdzi
                        drijdxk = 0.d0
                        drijdyk = 0.d0
                        drijdzk = 0.d0
                        fcutij=0.5d0*(dcos(pi*rij/cutoff)+1.d0)
                        temp1=0.5d0*(-dsin(pi*rij/cutoff))*(pi/cutoff)
                        dfcutijdxi=temp1*drijdxi
                        dfcutijdyi=temp1*drijdyi
                        dfcutijdzi=temp1*drijdzi
                        dfcutijdxj=-1.d0*dfcutijdxi
                        dfcutijdyj=-1.d0*dfcutijdyi
                        dfcutijdzj=-1.d0*dfcutijdzi
                        dfcutijdxk=0.0d0
                        dfcutijdyk=0.0d0
                        dfcutijdzk=0.0d0
                        do k_type = 1, nspecies
                            do k_neighbor = 1, at%atom(i)%count(k_type)
                                if ((k_type == j_type) .and. (k_neighbor <= j_neighbor)) cycle
                                rik = at%atom(i)%neighbor(j_type,j_neighbor,4)
                                xyz_k = at%atom(i)%neighbor(j_type,j_neighbor,1:3)
                                m = int(at%atom(i)%neighbor(j_type,j_neighbor,5))
                                deltaxk = -1.d0*(at%atom(i)%pos(1) - xyz_k(1))
                                deltayk = -1.d0*(at%atom(i)%pos(2) - xyz_k(2))
                                deltazk = -1.d0*(at%atom(i)%pos(3) - xyz_k(3))
                                drikdxi = -deltaxk/rik
                                drikdyi = -deltayk/rik
                                drikdzi = -deltazk/rik
                                drikdxk = -1.d0*drikdxi
                                drikdyk = -1.d0*drikdyi
                                drikdzk = -1.d0*drikdzi
                                drikdxj = 0.d0
                                drikdyj = 0.d0
                                drikdzj = 0.d0
                                fcutik=0.5d0*(dcos(pi*rik/cutoff)+1.d0)
                                temp1=0.5d0*(-dsin(pi*rik/cutoff))*(pi/cutoff)
                                dfcutikdxi=temp1*drikdxi
                                dfcutikdyi=temp1*drikdyi
                                dfcutikdzi=temp1*drikdzi
                                dfcutikdxj=0.0d0
                                dfcutikdyj=0.0d0
                                dfcutikdzj=0.0d0
                                dfcutikdxk=-1.d0*dfcutikdxi
                                dfcutikdyk=-1.d0*dfcutikdyi
                                dfcutikdzk=-1.d0*dfcutikdzi
                                rjk = (xyz_j(1) - xyz_k(1))**2 + (xyz_j(2) - xyz_k(2))**2 + (xyz_j(3) - xyz_k(3))**2
                                rjk = dsqrt(rjk)

                                if (rij < Rmin) then
                                    print*, ''
                                    stop
                                endif
                                drjkdxj = (xyz_j(1) - xyz_k(1))/rjk
                                drjkdyj = (xyz_j(2) - xyz_k(2))/rjk
                                drjkdzj = (xyz_j(3) - xyz_k(3))/rjk
                                drjkdxk = -1.d0*drjkdxj
                                drjkdyk = -1.d0*drjkdyj
                                drjkdzk = -1.d0*drjkdzj
                                drjkdxi = 0.d0
                                drjkdyi = 0.d0
                                drjkdzi = 0.d0
                                fcutjk=0.5d0*(dcos(pi*rjk/cutoff)+1.d0)
                                temp1=0.5d0*(-dsin(pi*rjk/cutoff))*(pi/cutoff)
                                dfcutjkdxj=temp1*drjkdxj
                                dfcutjkdyj=temp1*drjkdyj
                                dfcutjkdzj=temp1*drjkdzj
                                dfcutjkdxk=-1.d0*dfcutjkdxj
                                dfcutjkdyk=-1.d0*dfcutjkdyj
                                dfcutjkdzk=-1.d0*dfcutjkdzj
                                dfcutjkdxi=0.0d0
                                dfcutjkdyi=0.0d0
                                dfcutjkdzi=0.0d0
! costheta=(rjk**2 - rij**2 -rik**2)/(-2.d0*rij*rik)
                                f=rjk**2 - rij**2 -rik**2
                                g=-2.d0*rij*rik
                                costheta=f/g
                                costheta=costheta+1.d0 ! avoid negative values
                                dfdxi=-2.d0*rij*drijdxi - 2.d0*rik*drikdxi
                                dfdyi=-2.d0*rij*drijdyi - 2.d0*rik*drikdyi
                                dfdzi=-2.d0*rij*drijdzi - 2.d0*rik*drikdzi

                                dfdxj=2.d0*rjk*drjkdxj - 2.d0*rij*drijdxj
                                dfdyj=2.d0*rjk*drjkdyj - 2.d0*rij*drijdyj
                                dfdzj=2.d0*rjk*drjkdzj - 2.d0*rij*drijdzj

                                dfdxk=2.d0*rjk*drjkdxk - 2.d0*rik*drikdxk
                                dfdyk=2.d0*rjk*drjkdyk - 2.d0*rik*drikdyk
                                dfdzk=2.d0*rjk*drjkdzk - 2.d0*rik*drikdzk

                                dgdxi=-2.d0*(drijdxi*rik + rij*drikdxi)
                                dgdyi=-2.d0*(drijdyi*rik + rij*drikdyi)
                                dgdzi=-2.d0*(drijdzi*rik + rij*drikdzi)

                                dgdxj=-2.d0*drijdxj*rik
                                dgdyj=-2.d0*drijdyj*rik
                                dgdzj=-2.d0*drijdzj*rik

                                dgdxk=-2.d0*rij*drikdxk
                                dgdyk=-2.d0*rij*drikdyk
                                dgdzk=-2.d0*rij*drikdzk

                                temp1=1.d0/g**2
                                dcosthetadxi=(dfdxi*g - f*dgdxi)*temp1
                                dcosthetadyi=(dfdyi*g - f*dgdyi)*temp1
                                dcosthetadzi=(dfdzi*g - f*dgdzi)*temp1
                                dcosthetadxj=(dfdxj*g - f*dgdxj)*temp1
                                dcosthetadyj=(dfdyj*g - f*dgdyj)*temp1
                                dcosthetadzj=(dfdzj*g - f*dgdzj)*temp1
                                dcosthetadxk=(dfdxk*g - f*dgdxk)*temp1
                                dcosthetadyk=(dfdyk*g - f*dgdyk)*temp1
                                dcosthetadzk=(dfdzk*g - f*dgdzk)*temp1



    



END MODULE
