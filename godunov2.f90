program godunov2
!-------------------------------------------------------------
! Metodo di Godunov al secondo ordine con condizione TVD 
! per le equazioni di Eulero unidimensionali
! Condizioni al bordo trasmissive
!-------------------------------------------------------------
implicit none
real(kind=8) a, b, dx, dt, d, del, Tf, t, s, smax, c
real(kind=8), parameter:: gam=1.4 !1.330  ! rapporto tra calori specifici: acqua a 20°
real(kind=8), dimension(:), allocatable:: nodi,r,u,p,rn,un,pn,fr,fu,fp,rho,rhou,rhoe,frho,frhou,frhoe,rhon,rhoun,rhoen
integer n, np, nm, j, lim, sol, i
real(kind=8), dimension(3):: rjn,rjnp,der,derp
real(kind=8) rl,ul,pl,rr,ur,pr
real(kind=8) slope,limiter,lmax,tau
!-------------------------------------------------------
! Lettura dati e inizializzazione
!-------------------------------------------------------
open(unit=1,file='iniz.txt')
read(1,*) a, b
read(1,*) n
np=n+1
read(1,*) dx
allocate(nodi(0:n))
read(1,*) nodi
allocate(r(0:np),u(0:np),p(0:np))
read(1,*) r(0:np)
read(1,*) u(0:np)
read(1,*) p(0:np)
close(1)


allocate(frho(0:n),frhou(0:n),frhoe(0:n))  ! flussi delle variabili conservate
allocate(rho(0:np),rhou(0:np),rhoe(0:np))  ! variabili conservate
allocate(rhon(0:np),rhoun(0:np),rhoen(0:np))
allocate(rn(0:np),un(0:np),pn(0:np))       ! variabili primitive al passo successivo

!-------------------------------------------------------
! Scelta Limitatore e Solutore
!-------------------------------------------------------
print *, 'Scegli il limitatore'
print *, '1= SUPERBEE'
print *, '2= MINMOD'
print *, '3= MC'
print *, '4= LAX-WENDROFF'
print *, '5= NULLO (1° ORDINE)'
read *, lim

print *, 'Scegli il solutore'
print *, '1= Esatto'
print *, '2= Roe'
read *, sol


print *, 'tempo finale'
read *, Tf
t=0

! passaggio alle variabili conservate
do j=0,np
	rho(j)=r(j)
	rhou(j)=r(j)*u(j)
	rhoe(j)=p(j)/(gam-1)+.5*rhou(j)*u(j)
end do

!-------------------------------------------------------
! Ciclo principale del metodo di Godunov
!-------------------------------------------------------
do while (t.lt.Tf)	
	! queste celle fittizie forniscono le condizioni al bordo trasmissive
	r(0)=r(1)
	r(n+1)=r(n)
	u(0)=u(1)
	u(n+1)=u(n)
	p(0)=p(1)
	p(n+1)=p(n)
	
	! Passo temporale
	lmax=0
	do i=0,np
		lmax=max(lmax,abs(u(i))+sqrt(gam*p(i)/r(i)))
	end do
	tau=.6*dx/lmax
	if ((t+tau).gt.Tf) then
		tau=(Tf-t)
	end if
	t=t+tau
	
	
	do j=0,n
		call stati(r,u,p,j,np,tau,dx,rl,ul,pl,rr,ur,pr,gam,lim)
		if (sol==1) then
			call riemann(rl,ul,pl,rr,ur,pr,gam,frho(j),frhou(j),frhoe(j),s)
		else
			call roechar(rl,ul,pl,rr,ur,pr,gam,frho(j),frhou(j),frhoe(j),s)
		end if
	end do

	write (6,*) 't = ', t
	!write (6,*) 'smax =', smax
	c=tau/dx
	do j=1,n
		rhon(j)=rho(j) + c*(frho(j-1)-frho(j))
		rhoun(j)=rhou(j) + c*(frhou(j-1)-frhou(j))
		rhoen(j)=rhoe(j) + c*(frhoe(j-1)-frhoe(j))
	end do
	
	do j=1,n
		! ritorno alle variabili primitive
		r(j)=rhon(j)
		u(j)=rhoun(j)/r(j)
		p(j)=(gam-1)*(rhoen(j)-.5*rhoun(j)*u(j))
	end do
	rho=rhon
	rhou=rhoun
	rhoe=rhoen
end do

!-------------------------------------------------------
! Scrittura su file
!-------------------------------------------------------
nm=n-1
open(unit=2,file='rho.dat')
do j=0,n
	write(2,*) nodi(j), r(j+1)
end do
close(2)
open(unit=3,file='u.dat')
do j=0,n
	write(3,*) nodi(j), u(j+1)
end do
close(3)
open(unit=4,file='p.dat')
do j=0,n
	write(4,*) nodi(j), p(j+1)
end do
close(4)

end program godunov2












subroutine stati(r,u,p,j,np,dt,dx,rl,ul,pl,rr,ur,pr,gam,lim)
!-------------------------------------------------------
! Calcolo stati al tempo n+1/2
!-------------------------------------------------------
implicit none
integer np, j, lim, jm1, jp2
real(kind=8) dx,rl,ul,pl,rr,ur,pr,gam
real(kind=8), dimension(3):: der,derp
real(kind=8), dimension(0:np):: r,u,p
real(kind=8) dt,slope,limiter

if (j==0) then
	jm1=j
else
	jm1=j-1
end if
if (j==np-1) then
	jp2=j+1
else
	jp2=j+2
end if

der(1)=limiter(slope(r(jm1),r(j),r(j+1)),lim)*(r(j+1)-r(j))
der(2)=limiter(slope(u(jm1),u(j),u(j+1)),lim)*(u(j+1)-u(j))
der(3)=limiter(slope(p(jm1),p(j),p(j+1)),lim)*(p(j+1)-p(j))

derp(1)=limiter(slope(r(j),r(j+1),r(jp2)),lim)*(r(jp2)-r(j+1))
derp(2)=limiter(slope(u(j),u(j+1),u(jp2)),lim)*(u(jp2)-u(j+1))
derp(3)=limiter(slope(p(j),p(j+1),p(jp2)),lim)*(p(jp2)-p(j+1))

rl= r(j) - .5*dt/dx*(u(j)*der(1)+r(j)*der(2)) +.5*der(1)
rr= r(j+1) - .5*dt/dx*(u(j+1)*derp(1)+r(j+1)*derp(2)) -.5*derp(1)
ul= u(j) - .5*dt/dx*(u(j)*der(2)+der(3)/r(j)) +.5*der(2)
ur= u(j+1) - .5*dt/dx*(u(j+1)*derp(2)+derp(3)/r(j+1)) -.5*derp(2)
pl= p(j) - .5*dt/dx*(u(j)*der(3)+gam*p(j)*der(2)) +.5*der(3)
pr= p(j+1) - .5*dt/dx*(u(j+1)*derp(3)+gam*p(j+1)*derp(2)) -.5*derp(3)
end subroutine stati



real(kind=8) function slope(um1,u,up1)
implicit none
real(kind=8) um1,u,up1
if (up1==u) then
	if ((u-um1).ge.0) then
		slope=huge(1.0)
	else
		slope=-huge(1.0)
	end if
else
	slope=(u-um1)/(up1-u)
end if
end


real(kind=8) function limiter(r,lim)
real(kind=8) r
integer lim
if (r.lt.0) then
	if (lim==4) then
		limiter=1
	else
		limiter=0
	end if
else
    select case (lim)
		case (1)
			if (r .lt. .5) then
				limiter=2*r
			else if (r.lt.1) then
				limiter=1
			else if (r.lt.2) then
				limiter=r
			else
				limiter=2
			end if
		case (2)
			if (r .lt. 1) then
				limiter=r
			else
				limiter=1
			endif
		case (3)
			if (r .lt. 1/3) then
				limiter=2*r
			else if (r.lt.3) then
				limiter=(1+r)/2
			else
				limiter=2
			end if
		case (4)
			limiter=1
		case (5)
			limiter=0
	end select
end if
end
