subroutine roechar(rho1,u1,p1,rho4,u4,p4,gam,frho,frhou,frhoe,s)
!--------------------------------------------------------------
! Solutore approssimato di Roe del problema di Riemann con uso variabili caratteristiche
!--------------------------------------------------------------
implicit none
real(kind=8) rho1,u1,p1,rho4,u4,p4,gam,frho,frhou,frhoe,s
real(kind=8), dimension(3):: lambda,flux,dw,du
real(kind=8), dimension(3,3):: R
real(kind=8) rhos,rhod,u,h,a,rho,h1,h4
integer j,k

! Calcolo medie di Roe
rhos=sqrt(rho1)
rhod=sqrt(rho4)
h1=(gam*p1/rho1)/(gam-1) + .5*u1**2 
h4=(gam*p4/rho4)/(gam-1) + .5*u4**2

rho=rhos*rhod
u=(rhod*u4+rhos*u1)/(rhod+rhos)
h=(rhos*h1+rhod*h4)/(rhod+rhos)
a=sqrt((gam-1)*(h-.5*u**2))

! Valori assoluti degli autovalori
lambda(1)=abs(u-a)
lambda(2)=abs(u)
lambda(3)=abs(u+a)
s=maxval(lambda)

! Differenza delle varibili primitive
du(1)=rho4-rho1
!du(2)=rho4*u4-rho1*u1
!du(3)=(p4-p1)/(gam-1) + .5*(rho4*u4**2-rho1*u1**2)
du(2)=u4-u1
du(3)=p4-p1

! (Differenza delle) Variabili caratteristiche
!dw(1)=du(3)-(u+a/(gam-1))*du(2)+(.5*u**2+a*u/(gam-1))*du(1)
!dw(2)=du(3)-u*du(2)+(u**2-h)*du(1)
!dw(3)=du(3)-(u-a/(gam-1))*du(2)+(.5*u**2-a*u/(gam-1))*du(1)
dw(1)=.5*(du(3)-rho*a*du(2))/(a**2)
dw(2)=du(1)-du(3)/(a**2)
dw(3)=.5*(du(3)+rho*a*du(2))/(a**2)


! Autovettori destri
R(1,1) = 1
R(2,1) = u - a
R(3,1) = h - u*a

R(1,2) = 1
R(2,2) = u
R(3,2) = .5*u**2

R(1,3) = 1
R(2,3) = u + a
R(3,3) = h + u*a


!Flusso medio
   flux(1)=(rho1*u1 + rho4*u4)/2
   flux(2)=(rho1*u1**2 + p1 + rho4*u4**2 + p4)/2
   flux(3)=(rho1*u1*h1 + rho4*u4*h4)/2

!Flusso di Roe
  do j = 1, 3
   do k = 1, 3
    flux(j) = flux(j) -.5*R(j,k) *lambda(k)*dw(k)
   end do
  end do

frho=flux(1)
frhou=flux(2)
frhoe=flux(3)
end subroutine roechar

