 module lie_m
  use global_m
  use transform_m
  use output_m
  use cutoffmerge_m

  implicit none

public::Int_LieSW  !Lie Integrator with adaptive stepsize control
public::Int_LieSY  !Symmetric Lie Series Integrator with fixed stepsize
public::Int_LieST  !Lie Integrator with adaptive stepsize control and adaptive series truncation

contains

!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!--------------Lie Integrator with adaptive stepsize control----------
!
!                                         by Siegfried Eggl  20080110 
!--------------------------------------------------------------------------------

subroutine Int_LieSW(p,rv, mass,add)

use transform_m
use output_m

integer(kind=intdble)::lterms,allocstat
integer(kind=intdble)::i,count,m0body,body,dim,bodyp1
integer(kind=intdble)::cutbody,m0cutbody,totbody,mbody
                     
type(Problem),intent(inout)::p 
type(additional),intent(in)::add

integer(kind=intdble)::order(1:p%NrKoerper)
real(kind=double),dimension(1:p%NrKoerper),intent(inout)::mass
real(kind=double),dimension(1:p%NrKoerper,1:6)::rv,rvini
real(kind=double)::fac(0:add%lterms+1),lan(1:add%lterms,1:add%lterms),tmod,dtmin
real(kind=double)::lbn(0:add%lterms,0:add%lterms),nonu(1:add%lterms,0:add%lterms),dtk,t,dtnew,dtkmax,dtkmin
real(kind=double),allocatable,dimension(:,:)::rho
real(kind=double),allocatable,dimension(:,:,:)::alloctest

real(kind=double)::kahan_y,kahan_t,kahan_c

logical::mflag,m0flag
logical,dimension(:,:),allocatable::merger


lterms=add%lterms
t=0._double
kahan_y=0._double
kahan_t=0._double
kahan_c=0._double

count=1
dtmin=add%lstepmin
dtk=dtmin
dtkmin=1.d12
dtkmax=0._double
body=p%NrKoerper-m0count
m0body=m0count
totbody=p%NrKoerper
cutbody=0
m0cutbody=0


dim=3

bodyp1=body+1
mbody=body



!for cutoff subroutine: are there massless particles?
if(m0body.gt.0) then
   m0flag=.true.
end if   


do i=1,totbody
   order(i)=i
end do

!-----------------------------------------------
! Benoetigte Koeffizienten fuer xi,lambda,phi, in Subroutine liestep
!-----------------------------------------------
call koeffan(lterms,lan)
call koeffbn(lterms,lbn)
call novernu(lterms,nonu)

!----------------------------------------------
! Tabelle fuer Fakultaeten
!----------------------------------------------
do i=0,lterms+1
   call faculty(i,fac(i))
end do


!alloctest: see if enough memory for lietimestep is available
    allocate(rho(1:totbody,1:body),alloctest(1:totbody,1:body,0:(2*lterms+3)),merger(1:totbody,1:body), stat=allocstat)
     if(allocstat.ne.0) then
      write(*,*)'Not enough memory to allocate mutual distance tensor. ', &
                'Decrease number of particles! Terminating program.'
      STOP
     else
      deallocate(alloctest)
     end if
merger(:,:)=.false.

!---------------------------------------------
!Programmschleife
!---------------------------------------------


rvini=rv
call lietimestep(rv,mass,totbody,body,cutbody,m0cutbody,dtk,lterms,lan,lbn,nonu,fac,p%eps,dtnew,rho)
rv=rvini


!--------------------------------------------------------------
prgs: do while(t.lt.p%tEnde)

  !Schrittweiten - Ausgabecheck
   if(p%tAusg.ne.0._double) then

      tmod=abs(p%tAusg-modulo(t,p%tAusg))*k  !Der nächste Schritt würde über den gewünschten Ausgabezeitpunkt hinausgehen, deswegen wird er auf die benötigte Länge zusammengestutzt         

      if(dtk>tmod) then
         dtk=tmod
      end if

      if(dtk.ge.(p%tAusg*k-dtmin).and.dtk.ne.p%tAusg*k) then !Damit wird sichergestellt,dass die Minimale Schrittweite nicht Aufgrund der fixen Ausgabezeiten unterschritten wird
         dtk=dtk-dtmin
      end if
   
      if(dtk<dtmin) then
   !      write(unit=*,fmt=*) 'Stepsize becoming too small!',dtk,'< ',dtmin
         dtk=dtmin
      end if
   end if 
 
!do a step
call lietimestep(rv,mass,totbody,body,cutbody,m0cutbody,dtk,lterms,lan,lbn,nonu,fac,p%eps,dtnew,rho)


   !t=t+dtk/k
   kahan_y=dtk/k-kahan_c
   kahan_t=t+kahan_y
   kahan_c=(kahan_t-t)-kahan_y
   t=kahan_t


   if(dtk.gt.dtkmax) then
      dtkmax=dtk
    end if
    if (dtk.lt.dtkmin) then
       dtkmin=dtk
   end if

   dtk=dtnew

   !Ausgabe 
   if (mod(real(t/p%tAusg,KIND=double),1._double) .lt.1.d-14.and.t.ge.p%tAusg.or. dble(p%tAusg).eq.0.d0)then
           call Out(t,rv(order(:),:),mass(order(:)),p)  

           !cutoff active?
           if(p%cutoff.gt.0._double) then
              call cutoffr(m0flag,dim,t,totbody,body,m0body,p%cutoff,cutbody,m0cutbody,rv,mass,order)
           end if
       end if
 
if(p%merge) then
   call detectmerge(totbody-m0cutbody,body,bodyp1,p%mergermin,p%mergermax,rho,mass,mflag,merger)

   if(mflag) then
      call scattermerge(p%mergermin,t,dim,totbody,body,m0body,cutbody,m0cutbody,rv,rho,mass,merger,order)
   end if
end if
if(p%cen) then
  call closeenc(t,totbody-m0cutbody,body,bodyp1,p%cend,p%cemlim,rho,rv,mass,order)
end if  


 count=count+1
end do prgs

write(*,*)'time accomplished:',t,'with an average stepsize of:',t/dble(count)
write(*,*)'minimum stepsize',dtkmin/k
write(*,*)'maximum stepsize',dtkmax/k


write(*,*)'number of massive particles at the beginning / end of integration:',body+cutbody,' / ',body
write(*,*)'number of massless particles at the beginning / end of integration:',m0body+m0cutbody,' / ',m0body
if(p%outegy) then
  
   write(*,*)'error in total energy (sum of dE)',dEtot
   write(*,*)'error in total angular momentum (sum of dL)',dLtot
end if


end subroutine Int_LieSW

!****************************************************************
subroutine lietimestep(rv,mass,totbody,body,cutbody,m0cutbody,dt,lt,lan,lbn,nonu,fac,eps,dtnew,rho)

implicit none
integer(kind=intdble)::i,j,body,dim,lt,l,n,nu,bodyp1,bodyges,cutbody,m0cutbody,totbody
real(kind=double),dimension(1:totbody,1:6),intent(inout)::rv
real(kind=double),dimension(1:totbody),intent(in)::mass
real(kind=double),dimension(1:totbody,1:body+cutbody)::rho,rho2
real(kind=double)::eps,dt,sum1,sumnu,dtnew
real(kind=double),intent(in)::lan(1:lt,1:lt),lbn(0:lt,0:lt),nonu(1:lt,0:lt)
real(kind=double)::rlk(1:totbody,1:body+cutbody,1:3),dxi(1:totbody,1:3,0:lt)
real(kind=double)::dphi(1:totbody,1:body+cutbody,0:lt),dlambda(1:totbody,1:body+cutbody,0:lt),fac(0:lt+1)




dim=3

bodyp1=body+cutbody+1
bodyges=totbody-m0cutbody

!-------------------------------------------------
! initialconditions :
! D°xi=xi
! D¹xi=eta
!------------------------------------------------

!$omp parallel default(NONE)   &
!$omp shared(rlk,rv,rho,rho2,body,bodyp1,bodyges,dim,dlambda,dxi,dphi,mass,lan,lbn,lt,nonu,fac,dt) &
!$omp private(l,j,i,sum1,nu,sumnu,n)
!$omp do
do j=1,body
 do i = 1,dim
   dxi(j,i,0)=rv(j,i)
   dxi(j,i,1)=rv(j,i+dim)
 end do 
end do  
!$omp end do nowait

!$omp do
do j=bodyp1,bodyges
 do i = 1,dim
   dxi(j,i,0)=rv(j,i)
   dxi(j,i,1)=rv(j,i+dim)
 end do  
end do
!$omp end do

!-------------------------------------------------
! massive particles: r_lj
! 
! r12=r2-r1
!-------------------------------------------------

!$omp do
do l=1,body-1
   do j=l+1,body
         rlk(l,j,:)=rv(j,1:dim)-rv(l,1:dim)
         rlk(j,l,:)=-rlk(l,j,:)
    end do
end do
!$omp end do nowait

!$omp do
do l=1,body
   rlk(l,l,:)=0._double
end do
!$omp end do nowait

!-------------------------------------------------
!  massless particles: r_lk
!-------------------------------------------------
!$omp do 
do l=bodyp1,bodyges
   do j=1,body
      rlk(l,j,:)=rv(j,1:dim)-rv(l,1:dim)
   end do
end do
!$omp end do

!-------------------------------------------------
!  massive particles: respective distances rho_lj
!------------------------------------------------
!$omp do
do l=1,body-1
   do j=l+1,body
      rho(l,j)=Sqrt(Dot_product(rlk(l,j,1:dim),rlk(l,j,1:dim)))
      rho(j,l)=rho(l,j)
   end do
end do
!$omp end do nowait

!$omp do
do l=1,body
   rho(l,l)=0._double
end do
!$omp end do nowait

!------------------------------------------------------------------------------
!   massless particles: respective distances to massive particles rho_lj , j>l
!------------------------------------------------------------------------------


!$omp do 
do l=bodyp1,bodyges
   do j=1,body
     rho(l,j)=Sqrt(Dot_product(rlk(l,j,1:dim),rlk(l,j,1:dim)))
    end do
end do
!$omp end do    

!$omp do
do l=1,body
   rho2(l,:)=rho(l,:)*rho(l,:)
end do
!$omp end do nowait

!$omp do
do l=bodyp1,bodyges
   rho2(l,:)=rho(l,:)*rho(l,:)
end do
!$omp end do

!-------------------------------------------------
! D°lambda=xi.eta=Sum(D°xi.D¹xi)
!------------------------------------------------
!
!
!  massive
!
!$omp do 
   do l=1,body
      do j=1,body
         dlambda(l,j,0)=Dot_Product ((dxi(j,:,0)-dxi(l,:,0)),(dxi(j,:,1)-dxi(l,:,1)))   
       end do
    end do 
!$omp end do nowait

! massless

!$omp do 
   do l=bodyp1,bodyges
      do j=1,body    
          dlambda(l,j,0)=Dot_Product ((dxi(j,:,0)-dxi(l,:,0)),(dxi(j,:,1)-dxi(l,:,1)))   
       end do
   end do
!$omp end do 

!------------------------------------
! D°phi=phi
!------------------------------------
!
!   massive 
!

!$omp do 
do l=1,body
   do j=1,body  
     if(l.eq.j) then
     else 
      dphi(l,j,0)=1._double/(rho(l,j)*rho2(l,j))
     end if
   end do
end do
!$omp end do nowait

! massless

!$omp do 
do l=bodyp1,bodyges
   do j=1,body  
      dphi(l,j,0)=1._double/(rho(l,j)*rho2(l,j))
   end do 
end do
!$omp end do nowait


!$omp do
do l=1,body
   dphi(l,l,:)=0._double
end do
!$omp end do 
!--------------------------------------------------
! D²xi=Sum(mass(l)*(D°xi(l)-D°xi(j))*D°phi(l,j)) equals m*r_lj/|r_lj|³
!------------------------------------------------
!
!   massive
!

!$omp do
do l =1,body       !bodies
   do i = 1,dim      !coordinates
      sum1=0._double
      do j=1,body       !force on body l = sum over all influencing bodies j
        if (l.eq.k)then
        else
             sum1=sum1-mass(j)*rlk(j,l,i)*dphi(j,l,0)
         end if
      end do
      do j=bodyp1,bodyges     !massless on massive
          if (mass(j).eq.0._double) then
          else
             sum1=sum1-mass(j)*rlk(j,l,i)*dphi(j,l,0)
          end if
      end do
      dxi(l,i,2)=sum1
   end do
end do
!$omp end do 

!
!   massless
!

!$omp do
do l =bodyp1,bodyges     !bodies
   do i = 1,dim      !coordinates
      sum1=0._double
      do j=1,body       !force on body l = sum over all influencing bodies j
             sum1=sum1+mass(j)*rlk(l,j,i)*dphi(l,j,0)
      end do
      dxi(l,i,2)=sum1
   end do
end do
!$omp end do 

!-------------------------------------------------------------
! CALCULATION OF HIGHER LIE SERIES TERMS VIA RECURSION FORMULAE
!--------------------------------------------------------------

do n=1,lt-2     !SEEMS TO WORK FINE WITH GFORTRAN4.4 OPENMP (ORDER IS OK FOR EACH THREAD)
!Lie-Terms: lt-2 because calculation starts at D^(n+2)xi. 

!----------------------------------
!phi_lj
!----------------------------------
!massive

!$omp do
    do l=1,body                                                               
       do j=1,body
          if(l.eq.j) then
          else
            sum1=0._double
           do nu=0,n-1
              sum1=sum1+lan(n,nu+1)*dphi(j,l,n-1-nu)*dlambda(j,l,nu)
           end do           
               dphi(j,l,n)=sum1*1._double/(rho2(j,l))
          end if
        end do
!massless on massive
        do j=bodyp1,bodyges
          if(l.eq.j) then
          else
            sum1=0._double
           do nu=0,n-1
              sum1=sum1+lan(n,nu+1)*dphi(j,l,n-1-nu)*dlambda(j,l,nu)
           end do           
               dphi(j,l,n)=sum1*1._double/(rho2(j,l))
          end if
        end do
     end do
!$omp end do nowait

!   massless
!$omp do
    do l=bodyp1,bodyges                                                             
        do j=1,body
            sum1=0._double
           do nu=0,n-1
              sum1=sum1+lan(n,nu+1)*dphi(l,j,n-1-nu)*dlambda(l,j,nu)
           end do           
               dphi(l,j,n)=sum1*1._double/(rho2(l,j))
         end do
    end do
!$omp end do 


!----------------------------------------------------
! D^(n+2)xi_lj
!---------------------------------------------------
!massive

!$omp do
  do l=1,body
     do i=1,dim
        sum1=0._double
        do j=1,body
           sumnu=0._double
           if (l.eq.j) then
                sumnu=0._double
           else
                do nu=0,n            
                   sumnu=sumnu+nonu(n,nu)*dphi(j,l,nu)*(dxi(j,i,n-nu)-dxi(l,i,n-nu))   
                end do
                sumnu=sumnu*mass(j)
                sum1=sum1+sumnu                
           end if
        end do
!massless on massive
        do j=bodyp1,bodyges
           sumnu=0._double
           if (mass(j).eq.0._double) then
                sumnu=0._double
           else
                do nu=0,n            
                   sumnu=sumnu+nonu(n,nu)*dphi(j,l,nu)*(dxi(j,i,n-nu)-dxi(l,i,n-nu))   
                end do
                sumnu=sumnu*mass(j)
                sum1=sum1+sumnu                
           end if
        end do

        dxi(l,i,n+2)=sum1
     end do
  end do
!$omp end do nowait

!massless

!$omp do 
  do l=bodyp1,bodyges
     do i=1,dim
        sum1=0._double
        do j=1,body
           sumnu=0._double
           if (l.ne.j.and.mass(j).ne.0._double) then
                do nu=0,n            
                   sumnu=sumnu+nonu(n,nu)*dphi(l,j,nu)*(dxi(j,i,n-nu)-dxi(l,i,n-nu))   
                end do
                sumnu=sumnu*mass(j)
                sum1=sum1+sumnu
             else 
                sumnu=0._double
             end if
        end do
        dxi(l,i,n+2)=sum1
     end do
  end do
!$omp end do 

!----------------------------------------------------------
!Lambda_lj 
!----------------------------------------------------------
!massive
!$omp do 
    do l=1,body                                                       
        do j=1,body
           if(l.eq.j) then
           else 
              sum1=0._double                                     
              do nu=0,nint(real(n)/2._double)
                     sum1=sum1+lbn(n,nu)*Dot_Product((dxi(j,:,nu)-dxi(l,:,nu)),(dxi(j,:,n+1-nu)-dxi(l,:,n+1-nu)))             
              end do            
            dlambda(j,l,n)=sum1
           end if
       end do
       
        do j=bodyp1,bodyges
              sum1=0._double                                     
              do nu=0,nint(real(n)/2._double)
                     sum1=sum1+lbn(n,nu)*Dot_Product((dxi(j,:,nu)-dxi(l,:,nu)),(dxi(j,:,n+1-nu)-dxi(l,:,n+1-nu)))             
              end do            
            dlambda(j,l,n)=sum1
       end do
     end do  
!$omp end do nowait
  
!massless
!$omp do
    do l=bodyp1,bodyges                                                    
        do j=1,body 
              sum1=0._double                                     
              do nu=0,nint(real(n)/2._double)
                     sum1=sum1+lbn(n,nu)*Dot_Product((dxi(j,:,nu)-dxi(l,:,nu)),(dxi(j,:,n+1-nu)-dxi(l,:,n+1-nu)))             
              end do
            dlambda(l,j,n)=sum1
       end do
     end do   
!$omp end do 
end do !lterms: end do
!-----------------------------------------------------
! rvout(1:3)=exp(dt*D)xi
! rvout(4:6)=exp(dt*D)eta       
!
! with eta=D xi, -> D eta=D² xi  / in the equation for velocities dxi(:,:,n+1)
!----------------------------------------------------

!$omp do
do l=1,body
 do i=1,dim
       rv(l,i)=0._double 
       rv(l,i+3)=0._double
    do j=0,lt-1                                   ! lt-1 because 0 is counted as first lie-series term                                                    
       rv(l,i)=rv(l,i)+dxi(l,i,j)/fac(j)*dt**j
       rv(l,i+3)=rv(l,i+3)+dxi(l,i,j+1)/fac(j)*dt**j
     end do
  end do
end do
!$omp end do nowait

!$omp do
do l=bodyp1,bodyges
 do i=1,dim
       rv(l,i)=0._double
       rv(l,i+3)=0._double
    do j=0,lt-1 
       rv(l,i)=rv(l,i)+dxi(l,i,j)/fac(j)*dt**j
       rv(l,i+3)=rv(l,i+3)+dxi(l,i,j+1)/fac(j)*dt**j
    end do
 end do
end do
!$omp end do
!$omp end parallel

!-------------------------------------------------------------------------------
! new stepsize via
! |R_(lt)|< 2*|(dt*D)|^(lt)/(lt)! < epsilon
! consequently
! dtnew=(eps*(lt)!/2)/|maximum(D^(lt))|)^(1/t)
!--------------------------------------------------------------------------------

dtnew = (0.5_double*eps*fac(lt)/max(maxval(abs(dxi(1:body,:,lt))),maxval(abs(dxi(bodyp1:bodyges,:,lt)))))**(1._double/lt)

return
end subroutine lietimestep
!/////////////////////////////////////////////////////////////////////////////////////////////////
subroutine koeffan(lt,lan)
!-------------------------------
! berechnet die Koeffizienten
! die bei der Ableitung von 
! phi = rho³ entstehen
!
!paper: Hanselmeier and Dvorak
! Numerical integration with Lie-Series
!
! dependencies: none
! written by: Siegfried Eggl 20070118
! -----------------------------

implicit none
integer(kind=intdble)::i,j,lt
real(kind=double),dimension(1:lt,1:lt)::lan

lan(:,:)=0._double

do i=1,lt
   lan(i,i)=-3._double
end do

do i=2,lt
   lan(i,1)=lan(i-1,1)-2._double
   do j=2,lt-1
      lan(i,j)=lan(i-1,j-1)+lan(i-1,j)
   end do
end do

return

end subroutine koeffan

!****************************
subroutine novernu(n,nonu)
!----------------------------------------------------
! berechnet die Binomialkoeffizienten
! n ueber nu fuer die Liereihenmethode
! zur Loesung von Differentialgleichungen
! Da die 1. und 2.  Ableitung nicht im Rekursions-
! schema vorhanden sind, beginnt die Tafel bei
! den Koeffizienten fuer  D³ xi 
!
! i:       1...n
! nu : 0...i
!
! dependencies: function faculty
! written by: Siegfried Eggl  20070118
!---------------------------------------------------- 

implicit none
integer(kind=intdble)::i,nu,n
real(kind=double)::fac1,fac2,fac3
real(kind=double)::nonu(1:n,0:n)

nonu(1:n,0:n)=0._double


do i=1,n
   do nu=0,i
      if (i-nu<0)then 
         write(*,*)'subroutine novernu: argument error, i-nu: ',i-nu
      end if
        call faculty(i,fac1)
        call faculty(nu,fac2)
        call faculty(i-nu,fac3)
         nonu(i,nu)=fac1/(fac2*fac3)
         
   end do
end do

return
end subroutine novernu

!******************************
subroutine koeffbn(lt,lbn)
!-------------------------------
! berechnet die Koeffizienten
! die bei der Ableitung von 
! lambda= xi . eta entstehen
!
!paper: Hanselmeier and Dvorak
! Numerical integration with Lie-Series
!
! dependencies: none
! written by: Siegfried Eggl 20070118
! -----------------------------

implicit none
integer(kind=intdble)::i,j,m,lt
real(kind=double),dimension(0:lt,0:lt)::lbn



lbn(:,:)=0._double

lbn(:,0)=1._double


do i=1,lt
   m=nint(real(i)/2._double)

      do j=1,m
         lbn(i,j)=lbn(i-1,j-1)+lbn(i-1,j)
      end do
      if (mod(int(i),2).eq.0)then 
        lbn(i,m)=2._double*lbn(i-1,m)+lbn(i-1,m-1)
      else
         lbn(i,m)=lbn(i-1,m-1)     
 
      end if    
 
end do



return

end subroutine koeffbn
!*********************************************
subroutine dist2(rij,n,d)
!----------------------------------
! problemspezifisches Programm
! berechnet die Paarabstände zwischen den Körpern
! dependencies: none
!----------------------------------

implicit none
integer(kind=intdble)::i,j,n
real(kind=double)::rij(1:n,1:n,1:3),d(1:n,1:n)

do i=1,n-1
  do j=i+1,n
   d(i,j) = SQRT(rij(i,j,1)*rij(i,j,1)+rij(i,j,2)*rij(i,j,2 )+rij(i,j,3)*rij(i,j,3))
   d(j,i)=d(i,j)
  end do
end do

do i=1,n
   d(i,i)=0._double
end do

return

end subroutine dist2
!************************************************
subroutine faculty(a,value)
implicit none
integer(kind=intdble)::i
integer(kind=intdble),intent(in)::a
real(kind=double),intent(out)::value

value =1._double

if (a<0)then
   write(*,*)'function faculty: argument error',a
end if

if (a.eq.0)then
   value=1._double
else

   do i=1,a
      value=value*real(i)
   end do
end if

end subroutine faculty


!///////////////////////////////////////////////////////////////////////////////////////////////////////////////
 !--------------Symmetric Lie Series Integrator with fixed stepsize and fixed Series truncation-----------------
!
!                                         by Siegfried Eggl  20080110 
!--------------------------------------------------------------------------------
 subroutine Int_LieSY(p, rv, mass,add)
              use transform_m
              use output_m
        integer(kind=intdble)::i,count,ltnew,ltsum,ltmax,ltmin,lterms

        type(additional),intent(in)::add
        type(Problem),intent(in)::p
        real(kind=double),dimension(:),intent(in)::mass
        real(kind=double),dimension(1:p%NrKoerper,1:6)::rv,rvout,rvold
        real(kind=double)::fac(0:add%lterms),lan(1:add%lterms,1:add%lterms)
        real(kind=double)::lbn(0:add%lterms,0:add%lterms),nonu(1:add%lterms,0:add%lterms)
        real(kind=double)::dt,dtk,t,tt

dtk=p%Ns*k
dt=1._double !in order to get rid of roundoff errors due to readin 
t=0._double
lterms=add%lterms
count=0
ltsum=0
ltmin=lterms
ltmax=0

 !-----------------------------------------------
! Benoetigte Koeffizienten fuer xi,lambda,phi, in Subroutine liestep
!-----------------------------------------------
call koeffan(lterms,lan)
call koeffbn(lterms,lbn)
call novernu(lterms,nonu)

!----------------------------------------------
! Tabelle fuer Fakultaeten
!--------------------------------------------
do i=0,lterms
   call faculty(i,fac(i))
end do

ltnew=lterms

 !  rvold=rv
 ! rvold(:,4:6)=-rvold(:,4:6)

!   call lietermstep1(rv,rvout,mass,p%NrKoerper,dtk,lterms,lan,lbn,nonu,fac,p%eps,ltnew)


! 'LIE BACKWARDS ' 
   call lietermstep1(rv,rvout,mass,p%NrKoerper,-dtk,lterms,lan,lbn,nonu,fac,p%eps,ltnew)


! write(*,*)'LIE BACKWARDS ' 
! do i=1,2
!    write(*,*)rvout(i,:)
! end do
! 
! 
!    call lietermstep1(rv,rvout,mass,p%NrKoerper,dtk,lterms,lan,lbn,nonu,fac,p%eps,ltnew)
! 
! 
! write(*,*)'LIE FORWARDS ' 
! do i=1,2
!    write(*,*)rvout(i,:)
! end do
!  
! STOP

rvold=rvout
tt=t

prgs: do while(tt.lt.p%tEnde)

   call lietermstep(rv,rvout,rvold,mass,p%NrKoerper,dtk,lterms,lan,lbn,nonu,fac,ltnew)
   
   rvold(:,:)=rv(:,:)
   rv(:,:)=rvout(:,:)
  
   t=t+dt
   tt=t*p%Ns
      if (mod(real(t/p%tAusg,KIND=double),1._double) .lt.1.d-14 .and.t.ge.p%tAusg .or. dble(p%tAusg).eq.0.d0)then
       call Out(tt,rv,mass,p)
      end if
! !   
! !     if (modulo(tt/p%tAusg,1._double) .le. 1.d-4) then
! !       write(*,*)modulo(tt/p%tAusg,1._double),tt,p%tAusg,dt
! !     end if

      ltsum=ltsum+ltnew
      count=count+1
      if(ltnew.gt.ltmax) then
         ltmax=ltnew
      end if
      if (ltnew.lt.ltmin) then
         ltmin=ltnew
      end if

end do prgs

  write(*,*) 'average number of Lieterms:', real(ltsum)/real(count)
  write(*,*)'minimum number of Lieterms:',ltmin
  write(*,*)'maximum number of Lieterms:',ltmax
if(p%outegy) then
  
  write(*,*)'error in total energy (sum of dE)',dEtot
  write(*,*)'error in total angular momentum (sum of dL)',dLtot
end if
end subroutine Int_LieSY

!****************************************************************
subroutine lietermstep(rv,rvout,rvold,mass,body,dt,lt,lan,lbn,nonu,fac,ltnew)

implicit none
integer(kind=intdble)::i,body,dim,lt,l,n,nu,k,ltnew
parameter(dim=3)
real(kind=double),dimension(1:body,1:2*dim)::rv,rvout,rvold
real(kind=double),dimension(1:body)::mass
real(kind=double),dimension(1:body,1:body)::rho,rho2
real(kind=double)::dt,sum,sumnu,error
real(kind=double)::lan(1:lt,1:lt),lbn(0:lt,0:lt),nonu(1:lt,0:lt),rlk(1:body,1:body,1:3),dxi(1:body,1:dim,0:lt)
real(kind=double)::dphi(1:body,1:body,0:lt),dlambda(1:body,1:body,0:lt),fac(0:lt)

ltnew=lt

!-------------------------------------------------
! Anfangsbedingungen :
! D°xi=xi
! D¹xi=eta
!------------------------------------------------

do i = 1,dim
   dxi(:,i,0)=rv(:,i)
   dxi(:,i,1)=rv(:,i+3)
end do


!-------------------------------------------------
! Paarabstaende r_lk und deren Betraege rho_lk
! VORZEICHEN KONVENTION WIE UEBLICH!!!
! r12=r2-r1
!-------------------------------------------------
do l=1,body-1
   do k=l+1,body
         rlk(l,k,:)=rv(k,1:3)-rv(l,1:3)
         rlk(k,l,:)=-rlk(l,k,:)
    end do
end do

do l=1,body
   rlk(l,l,:)=0._double
end do

call dist2(rlk,body,rho)

rho2(:,:)=rho(:,:)*rho(:,:)

!-------------------------------------------------
! D°lambda=xi.eta=Sum(D°xi.D¹xi)
!------------------------------------------------

!   do l=1,body-1
!      do k=l+1,body
!         dlambda(l,k,0)=Dot_Product ((dxi(k,:,0)-dxi(l,:,0)),(dxi(k,:,1)-dxi(l,:,1)))   
!       end do
!    end do

   do l=1,body-1
      do k=l+1,body
         dlambda(l,k,0)=Dot_Product (rlk(l,k,:),(dxi(k,:,1)-dxi(l,:,1)))   
       end do
    end do

!------------------------------------
! D°phi=phi
!------------------------------------
do l=1,body-1
   do k=l+1,body   
     dphi(l,k,0)=1._double/(rho(l,k)*rho2(l,k))
     dphi(k,l,0)=dphi(l,k,0)
   end do
end do

do l=1,body
   dphi(l,l,:)=0._double
end do

!--------------------------------------------------
! D²xi=Sum(mass(l)*(D°xi(l)-D°xi(k))*D°phi(l,k)) enspricht m*r_lk/|r_lk|³
!-------------------------------------------------

do l =1,body       !Koerper
   do i = 1,dim      !Koordinaten
      sum=0._double
      do k=1,body       !Kraft auf Koerper l = Summe über alle Beiträge von Köpern k
         sum=sum+mass(k)*rlk(l,k,i)*dphi(l,k,0)
      end do
      dxi(l,i,2)=sum
   end do
end do


!-------------------------------------------------
! Brechnung der hoeheren Lie-Terme durch
! Rekursionen
!------------------------------------------------

do n=1,lt-2                                                                    !Lie-Terme: lt-2 da wir bei D^(n+2)xi starten
    

  !phi_lk
    do l=1,body-1                                                                
        do k=l+1,body
            sum=0._double
           do nu=0,n-1
              sum=sum+lan(n,nu+1)*dphi(l,k,n-1-nu)*dlambda(l,k,nu)
           end do           
               dphi(l,k,n)=sum*1._double/(rho2(l,k))
               dphi(k,l,n)=dphi(l,k,n)
         end do
      end do



! D^(n+2)xi_lk
  do l=1,body
     do i=1,dim
        sum=0._double
        do k=1,body
           sumnu=0._double
           if (l.ne.k) then
              do nu=0,n            
                 sumnu=sumnu+nonu(n,nu)*dphi(l,k,nu)*(dxi(k,i,n-nu)-dxi(l,i,n-nu))   
              end do
              sumnu=sumnu*mass(k)
              sum=sum+sumnu
              else 
                 sumnu=0._double
           end if
        end do
        dxi(l,i,n+2)=sum
     end do
  end do


  !Lambda_lk 

    do l=1,body-1                                                       
        do k=l+1,body 
              sum=0._double                                     
              do nu=0,nint(real(n)/2._double)
                     sum=sum+lbn(n,nu)*Dot_Product((dxi(k,:,nu)-dxi(l,:,nu)),(dxi(k,:,n+1-nu)-dxi(l,:,n+1-nu)))             
              end do
            dlambda(l,k,n)=sum
       end do
     end do  

!ist die Genauigkeitsgrenze schon erfüllt?
!-------------------------------------------------------------------------------
! Fehler der Ordnung lt
!nach der Formel zur Abschaetzung des Restgliedes der Exponentialreihe
! |R_(lt+1)|< 2*|(dt*D)|^(lt+1)/(lt+1)! < epsilon
!--------------------------------------------------------------------------------

  error=2._double*(dt**(n)*maxval(abs(dxi(:,:,(n)))))/fac(n)

! if(error.lt.eps) then
!      ltnew=n+2
!      exit
! end if



end do
  
!-----------------------------------------------------
! rvout(1:3)=exp(dt*D)xi
! rvout(4:6)=exp(dt*D)eta       
!
! mit eta=D xi, sodaß D eta=D² xi  / somit in der Gleichngung für die Geschwindigkeiten: dxi(:,:,n+1)
!----------------------------------------------------
rvout(:,:)=0._double

do i=1,dim
   
    do n=0,lt-2,2                                   ! lt-1 da 0 als erster Lie-Term mitgezählt wird
       rvout(:,i)=rvout(:,i)+2._double*dxi(:,i,n)/fac(n)*dt**n
       rvout(:,i+3)=rvout(:,i+3)+2._double*dxi(:,i,n+1)/fac(n)*dt**n
   end do
end do

       rvout(:,:)=rvout(:,:)-rvold(:,:)

!if(error .gt. eps) then
!   write(*,*) 'liestep error: failed to reach desired truncation error epsilon:',  eps, 'current error:',error
!end if
ltnew=lt

return
end subroutine lietermstep

!**************************************************************************************************

subroutine lietermstep1(rv,rvout,mass,body,dt,lt,lan,lbn,nonu,fac,eps,ltnew)

implicit none
integer(kind=intdble)::i,body,dim,lt,l,n,nu,k,ltnew
parameter(dim=3)
real(kind=double),dimension(1:body,1:2*dim)::rv,rvout
real(kind=double),dimension(1:body)::mass
real(kind=double),dimension(1:body,1:body)::rho,rho2
real(kind=double)::eps,dt,sum,sumnu,error
real(kind=double)::lan(1:lt,1:lt),lbn(0:lt,0:lt),nonu(1:lt,0:lt),rlk(1:body,1:body,1:3),dxi(1:body,1:dim,0:lt)
real(kind=double)::dphi(1:body,1:body,0:lt),dlambda(1:body,1:body,0:lt),fac(0:lt)

ltnew=lt
!-------------------------------------------------
! Anfangsbedingungen :
! D°xi=xi
! D¹xi=eta
!------------------------------------------------

do i = 1,dim
   dxi(:,i,0)=rv(:,i)
   dxi(:,i,1)=rv(:,i+3)
end do


!-------------------------------------------------
! Paarabstaende r_lk und deren Betraege rho_lk
! VORZEICHEN KONVENTION WIE UEBLICH!!!
! r12=r2-r1
!-------------------------------------------------
do l=1,body-1
   do k=l+1,body
         rlk(l,k,:)=rv(k,1:3)-rv(l,1:3)
         rlk(k,l,:)=-rlk(l,k,:)
    end do
end do

do l=1,body
   rlk(l,l,:)=0._double
end do

call dist2(rlk,body,rho)

rho2(:,:)=rho(:,:)*rho(:,:)

!-------------------------------------------------
! D°lambda=xi.eta=Sum(D°xi.D¹xi)
!------------------------------------------------

!   do l=1,body-1
!      do k=l+1,body
!         dlambda(l,k,0)=Dot_Product ((dxi(k,:,0)-dxi(l,:,0)),(dxi(k,:,1)-dxi(l,:,1)))   
!       end do
!    end do

   do l=1,body-1
      do k=l+1,body
         dlambda(l,k,0)=Dot_Product (rlk(l,k,:),(dxi(k,:,1)-dxi(l,:,1)))   
       end do
    end do

!------------------------------------
! D°phi=phi
!------------------------------------
do l=1,body-1
   do k=l+1,body   
     dphi(l,k,0)=1._double/(rho(l,k)*rho2(l,k))
     dphi(k,l,0)=dphi(l,k,0)
   end do
end do

do l=1,body
   dphi(l,l,:)=0._double
end do

!--------------------------------------------------
! D²xi=Sum(mass(l)*(D°xi(l)-D°xi(k))*D°phi(l,k)) enspricht m*r_lk/|r_lk|³
!-------------------------------------------------

do l =1,body       !Koerper
   do i = 1,dim      !Koordinaten
      sum=0._double
      do k=1,body       !Kraft auf Koerper l = Summe über alle Beiträge von Köpern k
         sum=sum+mass(k)*rlk(l,k,i)*dphi(l,k,0)
      end do
      dxi(l,i,2)=sum
   end do
end do


!-------------------------------------------------
! Brechnung der hoeheren Lie-Terme durch
! Rekursionen
!------------------------------------------------

do n=1,lt-2                                                                    !Lie-Terme: lt-2 da wir bei D^(n+2)xi starten
    

  !phi_lk
    do l=1,body-1                                                                
        do k=l+1,body
            sum=0._double
           do nu=0,n-1
              sum=sum+lan(n,nu+1)*dphi(l,k,n-1-nu)*dlambda(l,k,nu)
           end do           
               dphi(l,k,n)=sum*1._double/(rho2(l,k))
               dphi(k,l,n)=dphi(l,k,n)
         end do
      end do



! D^(n+2)xi_lk
  do l=1,body
     do i=1,dim
        sum=0._double
        do k=1,body
           sumnu=0._double
           if (l.ne.k) then
              do nu=0,n            
                 sumnu=sumnu+nonu(n,nu)*dphi(l,k,nu)*(dxi(k,i,n-nu)-dxi(l,i,n-nu))   
              end do
              sumnu=sumnu*mass(k)
              sum=sum+sumnu
              else 
                 sumnu=0._double
           end if
        end do
        dxi(l,i,n+2)=sum
     end do
  end do


  !Lambda_lk 

    do l=1,body-1                                                       
        do k=l+1,body 
              sum=0._double                                     
              do nu=0,nint(real(n)/2._double)
                     sum=sum+lbn(n,nu)*Dot_Product((dxi(k,:,nu)-dxi(l,:,nu)),(dxi(k,:,n+1-nu)-dxi(l,:,n+1-nu)))             
              end do
            dlambda(l,k,n)=sum
       end do
     end do  

!ist die Genauigkeitsgrenze schon erfüllt?
!-------------------------------------------------------------------------------
! Fehler der Ordnung lt
!nach der Formel zur Abschaetzung des Restgliedes der Exponentialreihe
! |R_(lt+1)|< 2*|(dt*D)|^(lt+1)/(lt+1)! < epsilon
!--------------------------------------------------------------------------------

  error=2._double*(dt**(n)*maxval(abs(dxi(:,:,(n)))))/fac(n)

 ! if(error.lt.eps) then
 !      ltnew=n+2
 !      exit
 ! end if

end do
  
!-----------------------------------------------------
! rvout(1:3)=exp(dt*D)xi
! rvout(4:6)=exp(dt*D)eta       
!
! mit eta=D xi, sodaß D eta=D² xi  / somit in der Gleichngung für die Geschwindigkeiten: dxi(:,:,n+1)
!----------------------------------------------------
rvout(:,:)=0._double

do i=1,dim
    do n=0,ltnew-1                                   ! lt-1 da 0 als erster Lie-Term mitgezählt wird
       rvout(:,i)=rvout(:,i)+dxi(:,i,n)/fac(n)*dt**n
       rvout(:,i+3)=rvout(:,i+3)+dxi(:,i,n+1)/fac(n)*dt**n
   end do
end do


if(error .gt. eps) then
   write(*,*) 'liestep error: failed to reach desired truncation error epsilon:',  eps, 'current error:',error
end if

return
end subroutine lietermstep1



!///////////////////////////////////////////////////////////////////////////////////////////////////////////////
!--------------Lie Integrator with adaptive stepsize and adaptive series truncation-----------------
!
!                                         by Siegfried Eggl  20080110 
!--------------------------------------------------------------------------------
 
subroutine Int_LieST(p, rv, mass,add)
 
  use transform_m
  use output_m
  integer(kind=intdble)::i,count,ltsum,ltmax,ltmin,tltmax,tltmin,cycnum,ltdid,tltmid
  type(Problem),intent(in)::p
  type(additional),intent(in)::add
!  integer,parameter::lterms=170 !subroutine 'faculty' cannot process a higher number of lie-series terms
  real(kind=double),dimension(:),intent(in)::mass
  real(kind=double),dimension(1:p%NrKoerper,1:6)::rv,rvout
  real(kind=double)::fac(0:add%ltmax),lan(1:add%ltmax,1:add%ltmax),dtkmax,dtkmin,dtmin,dtnew
  real(kind=double)::lbn(0:add%ltmax,0:add%ltmax),nonu(1:add%ltmax,0:add%ltmax),dt,dtk,t,tmod
  logical::repstep

dtk=k
dtnew=dtk
dt=p%Ns
t=0
count=0
ltsum=0
ltmin=add%ltmax
tltmax=add%ltmax
tltmin=add%ltmin
ltmax=0
dtmin=add%lstepmin
dtkmin=1.d12
dtkmax=0._double
tltmid=nint(real(tltmax+tltmin,kind=double)/2._double)

cycnum=0

 !-----------------------------------------------
! Benoetigte Koeffizienten fuer xi,lambda,phi, in Subroutine liestep
!-----------------------------------------------
call koeffan(tltmax,lan)
call koeffbn(tltmax,lbn)
call novernu(tltmax,nonu)

!----------------------------------------------
! Tabelle fuer Fakultaeten
!--------------------------------------------
do i=0,tltmax
   call faculty(i,fac(i))
end do

prgs: do while(t.lt.p%tEnde)

 
   call lieststep(rv,rvout,mass,p%NrKoerper,dtk,lan,lbn,nonu,fac,p%eps, &
                  tltmin,tltmax,tltmid,ltdid,repstep,dtnew)

!-------------------------------------------------------------
!    repeat a step due to not having reached desired precision
    if (repstep) then
      dtk=dtnew
      if (dtk<dtmin) then
         dtk=dtmin
       end if 
        cycnum=cycnum+1
        repstep=.false.
       cycle prgs
    end if
!-------------------------------------------------------------
 
     rv(:,:)=rvout(:,:)
      t=t+dtk/k


        dtk=dtnew
        

    if(dtk.gt.dtkmax) then
      dtkmax=dtk
    end if
    if (dtk.lt.dtkmin) then
       dtkmin=dtk
    end if

  
   !Heliozentrisches Koordinatensystem?
   if(p%center.eq.'h') then
                          call hekoo(rv,p%NrKoerper)
   end if   
  
   
    !Schrittweiten - Ausgabecheck
   if(p%tAusg.ne.0._double) then

      tmod=abs(dble(p%tAusg)-modulo(dble(t),dble(p%tAusg)))*k  !Der nächste Schritt würde über den gewünschten Ausgabezeitpunkt hinausgehen, deswegen wird er auf die benötigte Länge zusammengestutzt         
  
      if(dtk>tmod) then
         dtk=tmod
      end if

      if(dtk.ge.(p%tAusg*k-dtmin).and.dtk.ne.p%tAusg*k) then !Damit wird sichergestellt,dass die Minimale Schrittweite nicht Aufgrund der fixen Ausgabezeiten unterschritten wird
         dtk=abs(dtk-dtmin)
      end if
   
      if(dtk<dtmin) then
         write(unit=*,fmt=*) 'Stepsize becoming too small!',dtk,'< ',dtmin
      end if

 end if 

if(count.gt.0) then
   !Ausgabe
      if (mod(real(t/p%tAusg,KIND=double),1._double) .lt.1.d-15 .and.t.ge.p%tAusg .or. dble(p%tAusg).eq.0.d0)then
      call Out(t,rv,mass,p)    
    end if
end if
 
      ltsum=ltsum+ltdid
      count=count+1
      if(ltdid.gt.ltmax) then
         ltmax=ltdid
      end if
      if (ltdid.lt.ltmin) then
         ltmin=ltdid
      end if

end do prgs

  write(*,*) 'average number of Lieterms:', dble(ltsum)/dble(count)
  write(*,*)'minimum number of Lieterms:',ltmin
  write(*,*)'maximum number of Lieterms:',ltmax
write(*,*)'time accomplished:',t,'with an average stepsize of:',t/dble(count)
write(*,*)'minimum stepsize',dtkmin/k
write(*,*)'maximum stepsize',dtkmax/k
write(*,*)'number of repeated steps',cycnum

if(p%outegy) then
  
   write(*,*)'error in total energy (sum of dE)',dEtot
   write(*,*)'error in total angular momentum (sum of dL)',dLtot
   write(*,*)'accepted/rejected steps:',dble(count)/dble(cycnum)
end if
end subroutine Int_LieST
!****************************************************************
subroutine lieststep(rv,rvout,mass,body,dt,lan,lbn,nonu,fac,eps,tltmin,tltmax,tltmid,lt,repeatstep,dtnew)

implicit none
integer(kind=intdble)::i,body,dim,l,n,nu,k,tltmin,tltmax,tltmid
integer(kind=intdble),intent(out)::lt
integer(kind=intdble)::np2
parameter(dim=3)
real(kind=double),dimension(1:body,1:2*dim)::rv,rvout
real(kind=double),dimension(1:body)::mass
real(kind=double),dimension(1:body,1:body)::rho,rho2
real(kind=double)::eps,dt,sum,sumnu,error
real(kind=double)::lan(1:tltmax,1:tltmax),lbn(0:tltmax,0:tltmax),nonu(1:tltmax,0:tltmax)
real(kind=double)::rlk(1:body,1:body,1:3),dxi(1:body,1:dim,0:tltmax)
real(kind=double)::dphi(1:body,1:body,0:tltmax),dlambda(1:body,1:body,0:tltmax),fac(0:tltmax)
real(kind=double),intent(out)::dtnew
logical::repeatstep,step2big,step2small

dxi(:,:,:)=0._double
step2small=.false.
step2big=.false.
repeatstep=.false.

!----------------------------------
! receive: lt, tltmin, tltmax
! intent out: ltnew
!----------------------------------

!-------------------------------------------------
! Initial conditions :
! D°xi=xi
! D¹xi=eta
!------------------------------------------------


do i = 1,dim
   dxi(:,i,0)=rv(:,i)
   dxi(:,i,1)=rv(:,i+3)
end do


!-------------------------------------------------
! Pairdistances r_lk and their norms rho_lk
! BEWARE SIGNUM!!
! r12=r2-r1
!-------------------------------------------------
do l=1,body-1
   do k=l+1,body
         rlk(l,k,:)=rv(k,1:3)-rv(l,1:3)
         rlk(k,l,:)=-rlk(l,k,:)
    end do
end do

do l=1,body
   rlk(l,l,:)=0._double
end do

call dist2(rlk,body,rho)

rho2(:,:)=rho(:,:)*rho(:,:)

!-------------------------------------------------
! D°lambda=xi.eta=Sum(D°xi.D¹xi)
!------------------------------------------------

!   do l=1,body-1
!      do k=l+1,body
!         dlambda(l,k,0)=Dot_Product ((dxi(k,:,0)-dxi(l,:,0)),(dxi(k,:,1)-dxi(l,:,1)))   
!       end do
!    end do

   do l=1,body-1
      do k=l+1,body
         dlambda(l,k,0)=Dot_Product (rlk(l,k,:),(dxi(k,:,1)-dxi(l,:,1)))   
       end do
    end do
!------------------------------------
! D°phi=phi
!------------------------------------
do l=1,body-1
   do k=l+1,body   
     dphi(l,k,0)=1._double/(rho(l,k)*rho2(l,k))
     dphi(k,l,0)=dphi(l,k,0)
   end do
end do

do l=1,body
   dphi(l,l,:)=0._double
end do

!--------------------------------------------------
! D²xi=Sum(mass(l)*(D°xi(l)-D°xi(k))*D°phi(l,k)) enspricht m*r_lk/|r_lk|³
!-------------------------------------------------

do l =1,body       !Bodies
   do i = 1,dim      !Coordinates
      sum=0._double
      do k=1,body       !Force to body l = Sum over all bodies k
         sum=sum+mass(k)*rlk(l,k,i)*dphi(l,k,0)
      end do
      dxi(l,i,2)=sum
   end do
end do


!-------------------------------------------------
! Calculation of higher terms by recursive relations
!------------------------------------------------


do n=1,tltmax-2                                                                    !Lie-Terms: tltmax-2 because we start with D^(n+2)xi
    

  !phi_lk
    do l=1,body-1                                                                
        do k=l+1,body
            sum=0._double
           do nu=0,n-1
              sum=sum+lan(n,nu+1)*dphi(l,k,n-1-nu)*dlambda(l,k,nu)
           end do           
               dphi(l,k,n)=sum*1._double/(rho2(l,k))
               dphi(k,l,n)=dphi(l,k,n)
         end do
      end do



! D^(n+2)xi_lk
  do l=1,body
     do i=1,dim
        sum=0._double
        do k=1,body
           sumnu=0._double
           if (l.ne.k) then
              do nu=0,n            
                 sumnu=sumnu+nonu(n,nu)*dphi(l,k,nu)*(dxi(k,i,n-nu)-dxi(l,i,n-nu))   
              end do
              sumnu=sumnu*mass(k)
              sum=sum+sumnu
              else 
                 sumnu=0._double
           end if
        end do
        dxi(l,i,n+2)=sum
     end do
  end do


  !Lambda_lk 

    do l=1,body-1                                                       
        do k=l+1,body 
              sum=0._double                                     
              do nu=0,nint(real(n)/2._double)
                     sum=sum+lbn(n,nu)*Dot_Product((dxi(k,:,nu)-dxi(l,:,nu)),(dxi(k,:,n+1-nu)-dxi(l,:,n+1-nu)))             
              end do
            dlambda(l,k,n)=sum
       end do
     end do  


!precision reached?
!-------------------------------------------------------------------------------
!Estimate error of order lt 
!by using the remainder of exponential series
! 
! |R_(lt)|< 2*|(dt*D)|^(lt)/(lt)! < epsilon
!
!--------------------------------------------------------------------------------
np2=n+2

  if(np2.ge.tltmin) then
    error=2._double*(dt**(np2)*maxval(abs(dxi(:,:,np2))))/fac(np2)

    if(error.lt.eps) then
         if(np2.eq.tltmin) then
           step2small=.true.  !step dt is probably too small, so go to higher order in Lie series (tltmid) and adapt stepsize
           lt=np2
          end if
          
          if(np2.gt.tltmin.and..not.step2small) then !everything is fine, precision reached, leave 
            lt=np2
            exit
          end if

    end if
        
    if(step2small.and.np2.eq.tltmid) then     
            exit
    end if

    
    if (np2.ge.tltmax) then
          lt=np2
          step2big=.true.      !step dt is for sure too large: m
          exit
    end if
 
  end if

  
end do
  
!-------------------------------------------------------------------------------
! Calculation of the new stepsize via error in velocities: order lt
! using the remainder estimates of exponential series
!
! |R_(lt)|< 2*|(dt*D)|^(lt)/(lt)! < epsilon
!
! resulting in  
! dtnew=(eps*(lt)!/(2*|maximum(D^(lt))|))^(1/t)
!
! ATTENTION: eps = relative error! 
! dtnew=(eps*(lt)!/(2*maximum(|D^(lt)|))*minimum(||D^(0)|| )^(1/t)
!--------------------------------------------------------------------------------

if (step2big) then
!choose new stepsize for minimal, allowed order 
  dtnew = (eps*fac(tltmin)/(2._double*maxval(abs(dxi(1:body,:,tltmin)))))**(1._double/tltmin)
  repeatstep=.true.

 
!stepsize still too big but wouldn't go further down
  if(dtnew.ge.dt) then
    dtnew=dtnew/10._double
  end if

else

if(step2small) then
    !if possible make stepsize larger
    dtnew = (eps*fac(tltmid)/2._double)**(1._double/tltmid)*maxval(abs(dxi(1:body,:,tltmid)))
    if(dtnew.lt.dt) then
      dtnew=dt
  
    end if
 
    
!everything is ok, proceed with the same stepsize
 else
  dtnew=dt   

end if
 repeatstep=.false.

rvout(:,:)=0._double


!-----------------------------------------------------
! rvout(1:3)=exp(dt*D)xi
! rvout(4:6)=exp(dt*D)eta       
!
! with eta=D xi, sodaß D eta=D² xi  / velocities: dxi(:,:,n+1)
!----------------------------------------------------

do i=1,dim
    do n=0,lt-1                                   ! lt-1 because 0 will be counted as first term in lie-series
       rvout(:,i)=rvout(:,i)+dxi(:,i,n)/fac(n)*dt**n
       rvout(:,i+3)=rvout(:,i+3)+dxi(:,i,n+1)/fac(n)*dt**n
   end do
end do


end if



return
end subroutine lieststep
!/////////////////////////////////////////////////////////////////////////////////////////////


end module lie_m

