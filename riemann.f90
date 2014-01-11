subroutine riemann(rho1,u1,p1,rho4,u4,p4,gam,frho,frhou,frhoe,s)
!--------------------------------------------------------------
! Solutore esatto del problema di Riemann con determinazione dei flussi alle interfacce
!--------------------------------------------------------------
implicit none
real(kind=8) gam,rho1,u1,p1,rho4,u4,p4,a1,a4,del
real(kind=8) rho2,rho3,p2,p3,u,z,s
real(kind=8) af,pf,rhof,uf,frho,frhou,frhoe
real(kind=8) v1,v2,v3,v4,a2,a3,w1,w2,dp2,dp3,u2,u3
real(kind=8), parameter:: eps= 1.E-5
integer iter, flag1, flag2

del=(gam-1)/2
!velocità del suono
a1=sqrt(gam*p1/rho1)
a4=sqrt(gam*p4/rho4)

!--------------------------------------------------------------
! Controllo ammissibilità dati iniziali
!--------------------------------------------------------------
if ((u4-u1).ge.(a1+a4)/del) then
	write(6,*) 'i dati iniziali generano una regione di vuoto'
	write(6,*) 'il programma e'' stato fermato'
	stop
end if

!--------------------------------------------------------------
! Calcolo u di prima approssimazione e Newton-Raphson
!--------------------------------------------------------------
z=(p1/p4)**(del/gam)*a4/a1
u=(z*(a1+del*u1)-(a4-del*u4))/(del*(1+z))

iter=1
do 
	if (iter.gt.100) then
		write(6,*) 'non converge'
		write(6,*) flag1, flag2
		stop
	end if
	!print *, 'u=', u
	if (u.lt.u1) then
		!urto
		call leftshock(rho1,u1,p1,a1,u,p2,a2,gam,dp2,w1)
		!call shockleft(u1,a1,p1,u,a2,p2,dp2,w1,gam)
		flag1=2
		!print *, 'w1=', w1
	else
		!espansione
		call leftrar(rho1,u1,p1,a1,u,p2,a2,gam,del,dp2,v1,v2)
		flag1=1
		!print *, 'v1=', v1, '    v2=', v2
	end if
	
	if (u.gt.u4) then
		!urto
		call rightshock(rho4,u4,p4,a4,u,p3,a3,gam,dp3,w2)
		!call shockright(u4,a4,p4,u,a3,p3,dp3,w2,gam)
		flag2=2
		!print *, 'w2=', w2
	else
		!espansione
		call rightrar(rho4,u4,p4,a4,u,p3,a3,gam,del,dp3,v3,v4)
		flag2=1
		!print *, 'v3=', v3, '    v4=', v4
	end if

	if (abs(1-p2/p3).le.eps) then
		exit
	else
		u=u-(p2-p3)/(dp2-dp3)
	end if
	iter=iter+1
end do

rho2 = gam*p2/(a2**2)
rho3 = gam*p3/(a3**2)
rho1 = gam*p1/(a1**2)
rho4 = gam*p4/(a4**2)
u2 = u
u3 = u

!------------------------------------------------------------------
! Determinazione del flusso all'interfaccia
!------------------------------------------------------------------
if (u.gt.0) then
	af=a2
	pf=p2
	uf=u2
	rhof=gam*pf/(af**2)
	if (u.le.u1) then
		if (w1.gt.0) then
			af=a1
			uf=u1
			pf=p1
			rhof=gam*pf/(af**2)
		end if
	else
		if (v1.gt.0) then
			af=a1
			uf=u1
			pf=p1
			rhof=gam*pf/(af**2)
		else
			if (v2.gt.0) then
				af=(u1+a1/del)/(1+1/del)
				uf=af
				pf=p1*(af/a1)**(gam/del)
				rhof=gam*pf/(af**2)
			end if
		end if
	end if
else
	af=a3
	pf=p3
	uf=u3
	rhof=gam*pf/(af**2)
	if (u.ge.u4) then
		if (w2.lt.0) then
			af=a4
			uf=u4
			pf=p4
			rhof=gam*pf/(af**2)
		end if
	else
		if (v4.lt.0) then
			af=a1
			uf=u1
			pf=p1
			rhof=gam*pf/(af**2)
		else
			if (v3.lt.0) then
				af=-(u4-a4/del)/(1+1/del)
				uf=-af
				pf=p4*(af/a4)**(gam/del)
				rhof=gam*pf/(af**2)
			end if
		end if
	end if
end if
frho=rhof*uf
frhou=pf+frho*uf
frhoe=gam/(gam-1)*pf*uf+(rhof*uf**3)/2

!------------------------------------------------------------------
! Stima massima velocità di propagazione
!------------------------------------------------------------------

if (u.ge.u1) then
	!espansione sx
	s=abs(v1)
else
	!urto sx
	s=abs(w1)
end if

if (u.le.u4) then
	!espansione sx
	s=max(s,abs(v4))
else
	!urto sx
	s=max(s,abs(w2))
end if
end subroutine riemann




subroutine leftrar(rho1,u1,p1,a1,u,p,a,gam,del,dp,v1,v2)
!----------------------------------------------------------------------------------------------
! Determina lo stato a destra di una espansione, noto lo stato a sinistra e la velocità a destra
!----------------------------------------------------------------------------------------------
implicit none
real(kind=8) rho1,u1,p1,a1,u,p,a,gam,del,dp,v1,v2

a=a1+del*(u1-u)
p=p1*(a/a1)**(gam/del)
dp=-gam*p/a
v1=u1-a1
v2=u-a

end subroutine leftrar


subroutine leftshock(rho1,u1,p1,a1,u,p,a,gam,dp,w)
!----------------------------------------------------------------------------------------------
! Determina lo stato a destra di un urto sulla u-a, noto lo stato a sinistra e la velocità a destra
!----------------------------------------------------------------------------------------------
implicit none
real(kind=8) rho1,u1,p1,a1,u,p,a,gam,del,w,dp
real(kind=8) x,v,M1,rp

x=(gam+1)*(u-u1)/(4*a1)
M1=x-sqrt(1+x**2)    
w=a1*M1 + u1
rp=(2*gam*M1**2-gam+1)/(1+gam)
p=rp*p1
a=a1*sqrt((gam+1+(gam-1)*rp)/(gam+1+(gam-1)/rp))
dp=-gam*p/a

!print *, 'rho =', gam*p/(a**2)
!print *, 'u =', u
!print *, 'p =', p
!print *, 'a =', a
end subroutine leftshock


subroutine rightrar(rho4,u4,p4,a4,u,p,a,gam,del,dp,v3,v4)
!----------------------------------------------------------------------------------------------
! Determina lo stato a sinistra di una espansione, noto lo stato a destra e la velocità a sinistra
!----------------------------------------------------------------------------------------------
implicit none
real(kind=8) rho4,u4,p4,a4,u,p,a,gam,del,dp,v3,v4

a=a4-del*(u4-u)
p=p4*(a/a4)**(gam/del)
dp=gam*p/a
v3=u+a
v4=u4+a4

end subroutine rightrar


subroutine rightshock(rho4,u4,p4,a4,u,p,a,gam,dp,w)
!----------------------------------------------------------------------------------------------
! Determina lo stato a sinistra di un urto sulla u+a, noto lo stato a destra e la velocità a sinistra
!----------------------------------------------------------------------------------------------
implicit none
real(kind=8) rho4,u4,p4,a4,u,p,a,gam,del,w,dp
real(kind=8) x,v,M4,rp

x=(gam+1.)*(u-u4)/(4.*a4)
M4=x+sqrt(1+x**2)   
w=u4 + a4*M4
rp=(2*gam*M4**2-gam+1)/(1+gam)
p=rp*p4
a=a4*sqrt((gam+1+(gam-1)*rp)/(gam+1+(gam-1)/rp))
dp=gam*p/a

end subroutine rightshock
