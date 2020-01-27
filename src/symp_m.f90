module symp_m
  use global_m
  use transform_m
  use output_m
  use cutoffmerge_m

  implicit none

public::Int_Candy
public::Int_Yoshida
contains

!/////////////////////////////////////////////////////////////////////////////////////////////////////////
!---------------------------Candy-Integration--------------------------
!
!                                               by Siegfried Eggl  20080110 
!------------------------------------------------------------------
Subroutine Int_Candy(p, rv, m)
!--------------------------------------------------------
! Symplektischer Integrator nach Candy
! v=p/m
! a=F/m
!-------------------------------------------------------
type(Problem),intent(in)::p
real(kind=double),dimension(:),intent(inout)::m
real(kind=double),dimension(1:p%NrKoerper,1:6)::rv
real(kind=double),dimension(1:p%NrKoerper,1:3)::v2,v3,v4,q1,q2,q3,a
real(kind=double),dimension(1:4)::coeffa,coeffb
real(kind=double)::dt,t,dtk,tmod
integer(kind=intdble)::i,dim,allocstat
integer(kind=intdble)::outsteps,outcnt
!cutoff and merge variables
integer(kind=intdble)::m0body,mbody,m0cutbody,totbody,cutbody,totm0cut,mbodyp1
integer(kind=intdble)::order(1:p%NrKoerper)
real(kind=double),dimension(:,:),allocatable::rho
real(kind=double),dimension(:,:,:),allocatable::alloctest
logical::mflag,m0flag
logical,dimension(:,:),allocatable::merger


!initialization
dim=3
totbody=p%NrKoerper
mbody=totbody-m0count
mbodyp1=mbody+1 !even though massive bodies may be lost, this will give the starting value for massless body loops
m0body=m0count
m0cutbody=0
totm0cut=totbody
 cutbody=0

mflag=.false.




!for cutoff subroutine: are there massless particles?
if(m0body.gt.0) then
   m0flag=.true.
end if  
 
!store particle ids 
do i=1,totbody
   order(i)=i
end do


!alloctest: see if enough memory is available
allocate(rho(1:totbody,1:mbody),alloctest(1:totbody,1:mbody,1:3),merger(1:totbody,1:mbody),stat=allocstat)
  if(allocstat.ne.0) then
      write(*,*)'Not enough memory to allocate mutual distance tensor. ', &
                'Decrease number of particles! Terminating program.'
      STOP
  else
      deallocate(alloctest)
  end if
merger(:,:)=.false.
rho=0._double  
      !timestep (looks weird, but is necessary in order to get rid of most of the roundoff from input)   

        if (p%NS.lt.1._double) then
         dt=real(nint(1._double/p%Ns),kind=double)*10._double**(2*int(log10(p%Ns)))
        else
         dt=real(nint(p%Ns),kind=double)
        end if
        dtk=dt*k
        t=0._double
        
        outsteps=nint(p%tAusg/dt)-1
        outcnt=0

!Coefficients (multiplied by timestep to save cpu resources)
 coeffa(1)=(2._double+2._double**(1._double/3._double)+2**(-1._double/3._double))/6._double*dtk
 coeffa(4)=coeffa(1)
 coeffa(2)=(1._double-2._double**(1._double/3._double)-2._double**(-1._double/3._double))/6._double*dtk
 coeffa(3)=coeffa(2)

 coeffb(1)=0._double
 coeffb(2)=1._double/(2._double-2._double**(1._double/3._double))*dtk
 coeffb(4)=coeffb(2)
 coeffb(3)=1._double/(1._double-2._double**(2._double/3._double))*dtk



              

! integrate
do while(t.lt.p%tEnde)

                totm0cut=totbody-m0cutbody

! step1
                !massive
                !$omp parallel default(shared)   &
                !$omp private(i)
                !$omp do
                do i=1,mbody
                        q1(i,:)=rv(i,1:3) + coeffa(1)* rv(i,4:6) 
                end do
                !$omp end do

                !massless (mbodyp1 = original number of mbody before merging and scattering + 1)
                
                !$omp do
                do i=mbodyp1,totm0cut
                        q1(i,:)=rv(i,1:3) + coeffa(1)* rv(i,4:6) 
                end do
                !$omp end do
                !$omp end parallel
               
                call force(mbody,mbodyp1,totm0cut,q1,m,a,rho) ! Kraefte berechnen
                        
 !step 2      
                !massive    
                !$omp parallel default(shared)   &
                !$omp private(i)
                !$omp do            
                do i=1,mbody
                  v2(i,:) = rv(i,4:6) + coeffb(2)* a(i,:)
                  q2(i,:) = q1(i,:) +coeffa(2)* v2(i,:)
                end do
                !$omp end do
                
                !massless  
                !$omp do           
                do i=mbodyp1,totm0cut
                  v2(i,:) = rv(i,4:6) + coeffb(2)* a(i,:)
                  q2(i,:) = q1(i,:) +coeffa(2)* v2(i,:)
                end do
                !$omp end do
                !$omp end parallel

                call force(mbody,mbodyp1,totm0cut,q2,m,a,rho)  ! Kraefte berechnen
    !step 3        

                !massive
                !$omp parallel default(shared)   &
                !$omp private(i)
                !$omp do
                do i=1,mbody
                    v3(i,:)= v2(i,:) + coeffb(3) * a(i,:)
                    q3(i,:) = q2(i,:) + coeffa(3) * v3(i,:)
                end do
                !$omp end do nowait
 
                !massless
                !$omp do
                do i=mbodyp1,totm0cut
                    v3(i,:)= v2(i,:) + coeffb(3) * a(i,:)
                    q3(i,:) = q2(i,:) + coeffa(3) * v3(i,:)
                end do
                !$omp end do
                !$omp end parallel

               call force(mbody,mbodyp1,totm0cut,q3,m,a,rho) ! Kraefte berechnen

 !step 4 
                !massive
                 !$omp parallel default(shared)   &
                !$omp private(i)
                !$omp do
                do i=1,mbody   
                        v4(i,:) = v3(i,:) + coeffb(4) * a(i,:)
                        rv(i,1:3)= q3(i,:) + coeffa(4) * v4(i,:)     
                        rv(i,4:6) = v4(i,:)
                end do
                !$omp end do nowait

                !massless
                !$omp do
                do i=mbodyp1,totm0cut   
                        v4(i,:) = v3(i,:) + coeffb(4) * a(i,:)
                        rv(i,1:3)= q3(i,:) + coeffa(4) * v4(i,:)     
                        rv(i,4:6) = v4(i,:)
                end do
                !$omp end do
                !$omp end parallel

  
  t=(real(nint(t/dt,Kind=intdble),kind=double)+1._double)*dt
 ! output
  if (outsteps==outcnt) then
   !tmod=modulo(t/p%tAusg,1._double)
   !if (tmod.lt.1.d-14.or.dble(p%tAusg).eq.0.d0)then
     call Out(t, rv(order(:),:),m(order(:)),p)    
     outcnt=0
  else
     outcnt=outcnt+1
  end if

   
  !cutoff active?
     if(p%cutoff.gt.0._double) then
       call cutoffr(m0flag,dim,t,totbody,mbody,m0body,p%cutoff,cutbody,m0cutbody,rv,m,order)
     end if

 

!Merging?
if(p%merge) then
   call force(mbody,mbodyp1,totm0cut,rv(:,1:3),m,a,rho)
   call detectmerge(totm0cut,mbody,mbodyp1,p%mergermin,p%mergermax,rho,m,mflag,merger)

   if(mflag) then
      call scattermerge(p%mergermin,t,dim,totbody,mbody,m0body,cutbody,m0cutbody,rv,rho,m,merger,order)
      mflag=.false.
   end if
end if

if(p%cen) then
  call force(mbody,mbodyp1,totm0cut,rv(:,1:3),m,a,rho)
  call closeenc(t,totbody-m0cutbody,mbody,mbodyp1,p%cend,p%cemlim,rho,rv,m,order)
end if      

end do

write(*,*)'time accomplished:',t

write(*,*)'number of massive particles at the beginning / end of integration:', &
           mbody+cutbody,' / ',mbody
write(*,*)'number of massless particles at the beginning / end of integration:',& 
           m0body+m0cutbody,' / ',m0body


if(p%outegy) then
   write(*,*)'error in total energy (sum of dE)',dEtot
   write(*,*)'error in total angular momentum (sum of dL)',dLtot
end if

end Subroutine Int_Candy
!/////////////////////////////////////////////////////////////////////////////////////////////////////////
!--------------------------Yoshida-Symmetric-Integration--------------------------
!
!                                               by Siegfried Eggl  20110323 
!------------------------------------------------------------------
Subroutine Int_Yoshida(p, rv, m)
!--------------------------------------------------------
!Yoshida 1990
!-------------------------------------------------------
type(Problem),intent(in)::p
real(kind=double),dimension(:),intent(inout)::m
real(kind=double),dimension(1:p%NrKoerper,1:6)::rv
real(kind=double),dimension(1:p%NrKoerper,1:3)::p1,q1,a
real(kind=double),dimension(1:16)::coeffc
real(kind=double),dimension(1:15)::coeffd
real(kind=double),dimension(0:7)::w
real(kind=double)::dt,t,dtk,tmod
integer(kind=intdble)::outsteps,outcnt
integer(kind=intdble)::i,j,dim,allocstat
!cutoff and merge variables
integer(kind=intdble)::m0body,mbody,m0cutbody,totbody,cutbody,totm0cut,mbodyp1
integer(kind=intdble)::order(1:p%NrKoerper)
real(kind=double),dimension(:,:),allocatable::rho
real(kind=double),dimension(:,:,:),allocatable::alloctest
logical::mflag,m0flag
logical,dimension(:,:),allocatable::merger

!initialization
dim=3
totbody=p%NrKoerper
mbody=totbody-m0count
mbodyp1=mbody+1 !even though massive bodies may be lost, this will give the starting value for massless body loops
m0body=m0count
m0cutbody=0
totm0cut=totbody
 cutbody=0
mflag=.false.



!for cutoff subroutine: are there massless particles?
if(m0body.gt.0) then
   m0flag=.true.
end if  
 
!store particle ids 
do i=1,totbody
   order(i)=i
end do


!alloctest: see if enough memory is available
allocate(rho(1:totbody,1:mbody),alloctest(1:totbody,1:mbody,1:3),merger(1:totbody,1:mbody),stat=allocstat)
  if(allocstat.ne.0) then
      write(*,*)'Not enough memory to allocate mutual distance tensor. ', &
                'Decrease number of particles! Terminating program.'
      STOP
  else
      deallocate(alloctest)
  end if

rho=0._double
merger(:,:)=.false.
!timestep (looks weird, but is necessary in order to get rid of most of the roundoff from input)   
        if (p%NS.lt.1._double) then
         dt=real(nint(1._double/p%Ns),kind=double)*10._double**(2*int(log10(p%Ns)))
        else
         dt=real(nint(p%Ns),kind=double)
        end if
        dtk=dt*k
        t=0._double
             
        outsteps=nint(p%tAusg/dt)-1
        outcnt=0



! w(1)=0.102799849391985_double
! w(2)=-1.96061023297549_double
! w(3)=1.93813913762276_double
! w(4)=-0.158240635368243_double
! w(5)=-1.44485223686048_double
! w(6)=0.253693336566229_double
! w(7)=0.914844246229740_double


w(1)=-0.161582374150097E1
w(2)=-0.244699182370524E1
w(3)=-0.716989419708120E-2
w(4)=0.244002732616735E1
w(5)=0.157739928123617E0
w(6)= 0.182020630970714E1
w(7)=0.104242620869991E1


w(0)=1._double-2._double*sum(w(1:7))

do i=1,7
 coeffd(i)=w(8-i)
 coeffd(16-i)=coeffd(i)

 coeffc(i+1)=0.5_double*(w(8-i)+w(7-i))
 coeffc(16-i)=coeffc(i+1)
end do

 coeffd(8)=w(0)
 coeffc(1)=0.5_double*w(7)
 coeffc(16)=coeffc(1)
    
  !Coefficients (multiplied by timestep to save cpu resources)
        
   coeffc(:)=coeffc(:)*dtk
   coeffd(:)=coeffd(:)*dtk

! integrate
do while(t.lt.p%tEnde)
                totm0cut=totbody-m0cutbody
! step1
                !massive
                !$omp parallel default(shared)   &
                !$omp private(i)
                !$omp do 
                do i=1,mbody
                        q1(i,:)=rv(i,1:3) + coeffc(1)*rv(i,4:6) 
                end do
                !$omp end do nowait

                !massless (mbodyp1 = original number of mbody before merging and scattering + 1)                
                !$omp do
                do i=mbodyp1,totm0cut
                        q1(i,:)=rv(i,1:3) + coeffc(1)* rv(i,4:6) 
                end do
                !$omp end do
                !$omp end parallel
               
                call force(mbody,mbodyp1,totm0cut,q1,m,a,rho) ! Kraefte berechnen

                !massive
                !$omp parallel default(shared)   &
                !$omp private(i)  
                !$omp do 
                 do i=1,mbody
                   p1(i,:)=rv(i,4:6)
                 end do
                 !$omp end do nowait

                !massless  
                !$omp do             
                do i=mbodyp1,totm0cut
                   p1(i,:) = rv(i,4:6)
                end do
                !$omp end do
                !$omp end parallel                        
 !steps 2-15      
do j=2,15
 !step j
                !massive    
                !$omp parallel default(shared)   &
                !$omp private(i)
                !$omp do             
                do i=1,mbody
                  p1(i,:) = p1(i,:) + coeffd(j-1)* a(i,:)
                  q1(i,:) = q1(i,:) +coeffc(j)* p1(i,:)
                end do
                !$omp end do nowait
                
                !massless  
                !$omp do             
                do i=mbodyp1,totm0cut
                  p1(i,:) = p1(i,:) + coeffd(j-1)* a(i,:)
                  q1(i,:) = q1(i,:) + coeffc(j)* p1(i,:)
                end do
                !$omp end do
                !$omp end parallel

                call force(mbody,mbodyp1,totm0cut,q1,m,a,rho)  ! Kraefte berechnen
 end do

!step 16 
                !massive
                !$omp parallel default(shared)   &
                !$omp private(i)
                !$omp do
                do i=1,mbody   
                        rv(i,4:6) = p1(i,:) + coeffd(15) * a(i,:)
                        rv(i,1:3)= q1(i,:) + coeffc(16) * rv(i,4:6)                            
                end do
                !$omp end do nowait

                !massless
                !$omp do 
                do i=mbodyp1,totm0cut   
                       rv(i,4:6)  = p1(i,:) + coeffd(15) * a(i,:)
                       rv(i,1:3)= q1(i,:) + coeffc(16) * rv(i,4:6)     
                end do
                !$omp end do
                !$omp end parallel

!calculate new time (sorry, t=t+dt does not work due to accumulation of roundoff)  
  t=(real(nint(t/dt,Kind=intdble),kind=double)+1._double)*dt

 
! output
!    tmod=modulo(t/p%tAusg,1._double)
!    if (tmod.lt.1.d-14 .or.dble(p%tAusg).eq.0.d0)then
!      call Out(t, rv(order(:),:),m(order(:)),p)    
!    end if
   if (outsteps==outcnt) then
     call Out(t, rv(order(:),:),m(order(:)),p)    
     outcnt=0
   else
     outcnt=outcnt+1
   end if
   
   
  !cutoff active?
     if(p%cutoff.gt.0._double) then
       call cutoffr(m0flag,dim,t,totbody,mbody,m0body,p%cutoff,cutbody,m0cutbody,rv,m,order)
     end if

!Merging?
  if(p%merge) then
   call force(mbody,mbodyp1,totm0cut,rv(:,1:3),m,a,rho)
   call detectmerge(totm0cut,mbody,mbodyp1,p%mergermin,p%mergermax,rho,m,mflag,merger)
   if(mflag) then
      call scattermerge(p%mergermin,t,dim,totbody,mbody,m0body,cutbody,m0cutbody,rv,rho,m,merger,order)
      mflag=.false.
   end if
  end if
 !Close Encounters?
 if(p%cen) then
  call closeenc(t,totbody-m0cutbody,mbody,mbodyp1,p%cend,p%cemlim,rho,rv,m,order)
end if  
end do

write(*,*)'time accomplished:',t

write(*,*)'number of massive particles at the beginning / end of integration:', &
           mbody+cutbody,' / ',mbody
write(*,*)'number of massless particles at the beginning / end of integration:',& 
           m0body+m0cutbody,' / ',m0body


if(p%outegy) then
   write(*,*)'error in total energy (sum of dE)',dEtot
   write(*,*)'error in total angular momentum (sum of dL)',dLtot
end if

end Subroutine Int_Yoshida

!*********************************************************
subroutine force(mbody,mbodyp1,totm0cut,q,m,a,rho)
implicit none
integer(kind=intdble)::i,j
integer(kind=intdble),intent(in)::mbody,mbodyp1,totm0cut
real(kind=double),intent(inout),dimension(:)::m
real(kind=double),intent(in),dimension(:,:)::q
real(kind=double),intent(inout),dimension(:,:)::a,rho
real(kind=double),dimension(totm0cut,mbody,1:3)::r

!real(kind=double),dimension(1:3)::rbc,baryc,rqbc,abc
!real(kind=double)::mbc




!Massive
!$omp parallel default(NONE)   &
!$omp shared(r,q,m,a,rho,mbody,mbodyp1,totm0cut) private(i,j)
!$omp do
do i=1,mbody
 a(i,:)=0._double
end do
!$omp end do nowait
!$omp do
do i=mbodyp1,totm0cut
 a(i,:)=0._double
end do
!$omp end do

!$omp do
do i=1,mbody-1
      do j=i+1,mbody         
                r(i,j,:)=q(j,:)-q(i,:)
                r(j,i,:)=-r(i,j,:)
                rho(i,j)=Sqrt(Dot_Product(r(i,j,:),r(i,j,:)))
                rho(j,i)=rho(i,j)
     end do
end do
!$omp end do nowait

!Massless
!$omp do 
do i=mbodyp1, totm0cut
  do j=1,mbody
                r(i,j,:)=q(j,:)-q(i,:)
                rho(i,j)=Sqrt(Dot_Product(r(i,j,:),r(i,j,:)))
   end do
end do
!$omp end do

!Force on massive
!$omp do 
do i=1,mbody
  do j=1,mbody
          if(i.eq.j) then
           else
              a(i,:)= a(i,:) + m(j)*rho(i,j)**(-3._double)*r(i,j,:)
          end if
   end do  

   do j=mbodyp1,totm0cut
              a(i,:)= a(i,:) - m(j)*rho(j,i)**(-3._double)*r(j,i,:)
   end do               
end do
!$omp end do nowait


!Force on massless 

!calculate small body barycenter
!baryc(:)=0._double

!do i=mbodyp1,totm0cut
!  baryc(1:3)=baryc(1:3)+q(i,1:3)*m(i)
!end do

!$omp do
do i=mbodyp1, totm0cut
  do j=1,mbody
         a(i,:)= a(i,:) + m(j)*rho(i,j)**(-3._double)*r(i,j,:)
  end do   
!  !Force on small body i towards small body barycenter
!   mbc=sum(m(mbodyp1:totm0cut))-m(i)   !mass of barycenter without particle i
!   if(mbc.eq.0._double) then
!   else
!     rbc(:)=(baryc(:)-q(i,:)*m(i))/mbc  !position of barycenter without partilce i
!     rqbc(:)= rbc(:)-q(i,:)         !vector i to barycenter
!     abc(:)=mbc*rqbc(:)*(Dot_Product(rqbc(:),rqbc(:)))**(-1._double) 
!     a(i,:)=a(i,:)+abc(:)
!  end if
end do
!$omp end do
!$omp end parallel

return
end subroutine
!/////////////////////////////////////////////////////////////////////////////////////////////////////////
 end module symp_m
