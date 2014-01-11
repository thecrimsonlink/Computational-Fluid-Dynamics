program godunov
!-------------------------------------------------------------
! Metodo di Godunov per le equazioni di Eulero unidimensionali
! Condizioni al bordo trasmissive
!-------------------------------------------------------------
implicit none
real(kind=8) a, b, dx, dt, d, del, Tf, t, s, smax, c
real(kind=8), parameter:: gam=1.4 !1.330  ! rapporto tra calori specifici: acqua a 20°
real(kind=8), dimension(:), allocatable:: nodi,r,u,p,rn,un,pn,fr,fu,fp,rho,rhou,rhoe,frho,frhou,frhoe,rhon,rhoun,rhoen
integer n, np, nm , j, sol
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
! Ciclo principale
!-------------------------------------------------------
do while (t.lt.Tf)
	! queste celle fittizie forniscono le condizioni al bordo trasmissive
	r(0)=r(1)
	r(n+1)=r(n)
	u(0)=u(1)
	u(n+1)=u(n)
	p(0)=p(1)
	p(n+1)=p(n)
	
	smax=0
	do j=0,n
		!write(6,*) j
		if (sol==1) then
			call riemann(r(j),u(j),p(j),r(j+1),u(j+1),p(j+1),gam,frho(j),frhou(j),frhoe(j),s)
		else
			call roechar(r(j),u(j),p(j),r(j+1),u(j+1),p(j+1),gam,frho(j),frhou(j),frhoe(j),s)
		end if
		smax=max(s,smax)
	end do
	!print *, smax
	c=.9/smax
	dt=c*dx
	if ((t+dt).gt.Tf) then
		dt=Tf-t
	end if
	t=t+dt
	write (6,*) 't = ', t
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
end program godunov
