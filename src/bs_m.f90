 module bs_m
  use global_m
  use transform_m
  use output_m
  use cutoffmerge_m

  implicit none

public::Int_BS

contains
!//////////////////////////////////////////////////////////////////////7
subroutine Int_Bs(p, rv, mass)
!***************************************************************
! Bulirsch Stroer Extrapolation algorithm 
! containing Modified Midpoint Integration 
! and polynomial Extrapolation
! Core taken from Mercury6 (John Chambers) 
! written by 
! Siegfried Eggl   20090227
!***************************************************************
implicit none
        type(Problem),intent(in)::p
!        type(vec),dimension(1:p%NrKoerper)::rn,vn
        real(kind=double),dimension(1:p%NrKoerper,1:6)::rv,rvout
        real(kind=double),dimension(:)::mass
        real(kind=double)::dtk,dt,t,dtkdid,dtknxt,tmod,dtkmin,dtkmax,dtmin
        integer(kind=intdble)::bsrep,body,count,bodyp1,totm0cut
        integer(kind=intdble)::cutbody,m0cutbody,totbody,dim,m0body
        integer(kind=intdble)::order(1:p%NrKoerper),allocstat,i
        logical::mflag,m0flag
        logical,dimension(:,:),allocatable::merger 
        real(kind=double),allocatable,dimension(:,:)::rho
        real(kind=double),allocatable,dimension(:,:,:)::alloctest
        
        real(kind=double)::kahan_y,kahan_t,kahan_c
        
dim=3
dtk=real(p%Ns,kind=double)*K
dt=p%Ns
t=0._double
kahan_y=0._double
kahan_t=0._double
kahan_c=0._double
bsrep=8
dtmin=1.d-12
dtkmin=Huge(dtkmin)
dtkmax=Tiny(dtkmax)
count=0


body=p%NrKoerper-m0count
bodyp1=body+1
m0body=m0count
totbody=p%NrKoerper
cutbody=0
m0cutbody=0


!for cutoff subroutine: are there massless particles?
if(m0body.gt.0) then
   m0flag=.true.
end if   

do i=1,totbody
   order(i)=i
end do

!alloctest: see if enough memory is available
allocate(rho(1:totbody,1:body),alloctest(1:totbody,1:body,1:3),merger(1:totbody,1:body),stat=allocstat)
    if(allocstat.ne.0) then
      write(*,*)'Not enough memory to allocate mutual distance tensor. ', &
                'Decrease number of particles! Terminating program.'
      STOP
    else
      deallocate(alloctest)
    end if

merger(:,:)=.false.

prgs: do while(t.lt.p%tEnde)

!-------------------------------------------
! BS STEP
! der Dimensionsfaktor (Gaussgravitationskonstante) wurde in die Zeit verpackt
!------------------------------------------- 
   totm0cut=totbody-m0cutbody


   call  mdt_bs(dtk,dtkdid,dtknxt,p%eps,bsrep,totbody,body,cutbody,m0cutbody,&
                mass,rv,rvout,rho)
       
   !t=t+dtkdid/k
   kahan_y=dtkdid/k-kahan_c
   kahan_t=t+kahan_y
   kahan_c=(kahan_t-t)-kahan_y
   t=kahan_t

!check for global max and min of stepsize
    dtkmax=max(dtkmax,dtkdid)
    dtkmin=min(dtkmin,dtkdid)
    dtk=dtknxt
   
!massive 
    rv(1:body,:)=rvout(1:body,:)
!massless
    rv(bodyp1:totm0cut,:)=rvout(bodyp1:totm0cut,:)



   !Schrittweiten - Ausgabecheck
   if(p%tAusg.ne.0.d0) then

      tmod=abs(p%tAusg-modulo(t,p%tAusg))*k  !Der nächste Schritt würde über den gewünschten Ausgabezeitpunkt hinausgehen, deswegen wird er auf die benötigte Länge zusammengestutzt        
      if(dtk>tmod) then
         dtk=tmod
      end if

      if(dtk.ge.(p%tAusg*k-dtmin).and.dtk.ne.p%tAusg*k) then !Damit wird sichergestellt,dass die Minimale Schrittweite nicht Aufgrund der fixen Ausgabezeiten unterschritten wird
         dtk=dtk-dtmin
      end if
   
      if(dtk<dtmin) then
    !     write(unit=*,fmt=*) 'Stepsize becoming too small!',dtk,'< ',dtmin
    !     write(unit=*,fmt=*) 'setting stepsize to:',dtmin
         dtk=dtmin
      end if

 end if 

 
   !Ausgabe 
      if (mod(real(t/p%tAusg,KIND=double),1._double) .lt.1.d-15.and.t.ge.p%tAusg )then    
      call Out(t,rv(order(:),:),mass(order(:)),p)    
      !cutoff active?
           if(p%cutoff.gt.0._double) then
              call cutoffr(m0flag,dim,t,totbody,body,m0body,p%cutoff, &
                           cutbody,m0cutbody,rv,mass,order)
           end if
    end if

 

if(p%merge) then
!merging with massive particles only   
!   rhomerge=huge(rhomerge)
!   bcm=body+cutbody+m0body

 !  rhomerge(1:body,1:body)=rho(1:body,1:body)
 !  rhomerge(body+cutbody:bcm,1:body)= rho(body+cutbody:bcm,1:body)
  
 !  deallocate(rho)
 !calculate m0 particles' mutual distances
 !  do i=body+cutbody,bcm-1
 !     do j=i+1,bcm
 !        rhomerge(i,j)= Sqrt(Dot_Product(rv(i,1:dim),rv(j,1:dim)))
 !     end do
 !  end do

   call detectmerge(totm0cut,body,bodyp1,p%mergermin,p%mergermax,rho,& 
                    mass,mflag,merger)

   if(mflag) then
      call scattermerge(p%mergermin,t,dim,totbody,body,m0body,cutbody,&
                        m0cutbody,rv,rho,mass,merger,order)
   end if
end if

if(p%cen) then
  call closeenc(t,totbody-m0cutbody,body,bodyp1,p%cend,p%cemlim,rho,rv,mass,order)
end if  
  count=count+1
end do prgs

write(*,*)'time accomplished:',t,'with an average stepsize of:',t/real(count)
write(*,*)'minimum stepsize',dtkmin/k
write(*,*)'maximum stepsize',dtkmax/k

write(*,*)'number of massive particles at the beginning / end of integration:',body+cutbody,' / ',body
write(*,*)'number of massless particles at the beginning / end of integration:',m0body+m0cutbody,' / ',m0body


return
end subroutine Int_BS
!///////////////////////////////////////////////////////////////////

!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     MDT_BS1.FOR    (ErikSoft   2 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!Author: John E. Chambers
!
!Integrates NBOD bodies (of which NBIG are Big) for one timestep H0
!using the Bulirsch-Stoer method. The accelerations are calculated using the 
!subroutine FORCE. The accuracy of the step is approximately determined 
!by the tolerance parameter TOL.
!
!
!Last Modified by
!
! Siegfried Eggl   20090227
!------------------------------------------------------------------------------
!
subroutine mdt_bs(h0,hdid,hnext,tol,bsrep,totbody,body,cutbody,m0cutbody,mass,rv0,rv,rho)
      implicit none
     
!Input/Output
      integer(kind=intdble):: totbody,body,cutbody,m0cutbody,hred,bsrep,bsrep1
      real(kind=double):: h0,hdid,tol,mass(1:totbody)
      real(kind=double),dimension(1:totbody,1:6),intent(in):: rv0
      real(kind=double),dimension(1:totbody,1:6),intent(out):: rv
   
!Local
      integer(kind=intdble):: j, j1, k, n,tm0c,bpcp1
      real(kind=double):: tmp0,tmp1,tmp2,errmax,tol2,h,hx2,h2(1:bsrep),hnext
      real(kind=double):: a(1:totbody,1:3),a0(1:totbody,1:3),d(1:6,1:totbody,1:bsrep), &
                          xscal(1:totbody),vscal(1:totbody)
      real(kind=double),dimension(1:totbody,1:6)::rvend
      !
! 
      real(kind=double):: SHRINK,GROW
      parameter (SHRINK=.55_double,GROW=1.3_double)

      real(kind=double)::rho(:,:)
      
!------------------------------------------------------------------------------
!
      tol2 = tol*tol
!
      bsrep1=bsrep-1

      tm0c=totbody-m0cutbody
      bpcp1=body+cutbody+1
!massive:      
!Calculate arrays used to scale the relative error (R^2 for position and
!V^2 for velocity).


!$omp parallel default(shared)   &
!$omp private(k,tmp1,tmp2)

!$omp do
      do k = 1, body
         tmp1 = Dot_Product(rv0(k,1:3),rv0(k,1:3))
         tmp2 = Dot_Product(rv0(k,4:6),rv0(k,4:6))
        xscal(k) = 1._double / tmp1
        vscal(k) = 1._double / tmp2
      end do
!$omp end do nowait

!massless:
!$omp do
      do k = bpcp1, tm0c
         tmp1 = Dot_Product(rv0(k,1:3),rv0(k,1:3))
         tmp2 = Dot_Product(rv0(k,4:6),rv0(k,4:6))
        xscal(k) = 1._double / tmp1
        vscal(k) = 1._double / tmp2
      end do
!$omp end do
!$omp end parallel

!Calculate accelerations at the start of the step

      call frho(rv0,totbody,body,cutbody,m0cutbody,mass,a0,rho)
   
      errmax=Huge(errmax)

do while(errmax.gt.tol2)
!
!For each value of N, do a modified-midpoint integration with 2N substeps
      do n = 1, bsrep
        h = h0 / (2._double * dble(n))
        h2(n) = .25_double / dble(n*n)
        hx2 = h * 2._double
!
!massive
!$omp parallel default(shared)   &
!$omp private(k)
!$omp do
        do k = 1, body
           rv(k,1:3)=rv0(k,1:3)+h*rv0(k,4:6)
           rv(k,4:6) = rv0(k,4:6) + h*a0(k,1:3)
        end do
!$omp end do  nowait

!massless
!$omp do 
         do k = bpcp1, tm0c
           rv(k,1:3)=rv0(k,1:3)+h*rv0(k,4:6)
           rv(k,4:6) = rv0(k,4:6) + h*a0(k,1:3)
        end do
!$omp end do
!$omp end parallel

        call frho(rv,totbody,body,cutbody,m0cutbody,mass,a,rho)

!massive   
!$omp parallel default(shared)   &
!$omp private(k)    
!$omp do 
        do k = 1, body
          rvend(k,1:3)=rv0(k,1:3)+hx2*rv(k,4:6)
          rvend(k,4:6) = rv0(k,4:6) + hx2*a(k,1:3)
        end do
!$omp end do nowait

!massless
!$omp do 
        do k = bpcp1, tm0c
          rvend(k,1:3)=rv0(k,1:3)+hx2*rv(k,4:6)
          rvend(k,4:6) = rv0(k,4:6) + hx2*a(k,1:3)
        end do
!$omp end do
!$omp end parallel


        do j = 2, n


        call frho(rvend,totbody,body,cutbody,m0cutbody,mass,a,rho)

!massive        
          !$omp parallel default(shared)   &
          !$omp private(k)
          !$omp do
          do k = 1, body
              rv(k,1:3)=  rv(k,1:3)+hx2*rvend(k,4:6)
              rv(k,4:6) = rv(k,4:6) + hx2*a(k,1:3)
          end do
         !$omp end do nowait
!massless 
         !$omp do
           do k = bpcp1, tm0c
              rv(k,1:3)=  rv(k,1:3)+hx2*rvend(k,4:6)
              rv(k,4:6) = rv(k,4:6) + hx2*a(k,1:3)
          end do
         !$omp end do
         !$omp end parallel
          
        call frho(rv,totbody,body,cutbody,m0cutbody,mass,a,rho)

!massive 
         !$omp parallel default(shared)   &
         !$omp private(k)
          !$omp do
          do k = 1, body
            rvend(k,1:3)=rvend(k,1:3)+hx2*rv(k,4:6)
            rvend(k,4:6) = rvend(k,4:6) + hx2*a(k,1:3)
          end do
         !$omp end do nowait
!massless
        !$omp do
         do k = bpcp1, tm0c
          rvend(k,1:3)=rvend(k,1:3)+hx2*rv(k,4:6)
          rvend(k,4:6) = rvend(k,4:6) + hx2*a(k,1:3)
         end do
        !$omp end do
        !$omp end parallel    

        end do
      
        call frho(rvend,totbody,body,cutbody,m0cutbody,mass,a,rho)
!massive
         !$omp parallel default(shared)   &
         !$omp private(k)
          !$omp do
        do k = 1, body
          d(1:3,k,n) = .5_double*(rvend(k,1:3) + rv(k,1:3) + h*rvend(k,4:6))
          d(4:6,k,n) = .5_double*(rvend(k,4:6) + rv(k,4:6) + h*a(k,1:3))
        end do
        !$omp end do nowait
!massless
        !$omp do
        do k = bpcp1, tm0c
          d(1:3,k,n) = .5_double*(rvend(k,1:3) + rv(k,1:3) + h*rvend(k,4:6))
          d(4:6,k,n) = .5_double*(rvend(k,4:6) + rv(k,4:6) + h*a(k,1:3))
        end do
         !$omp end do
        !$omp end parallel 


!Update the D array, used for polynomial extrapolation
        do j = n - 1, 1, -1
          j1 = j + 1
          tmp0 = 1._double / (h2(j) - h2(n))
          tmp1 = tmp0 * h2(j1)
          tmp2 = tmp0 * h2(n)
!massive
        !$omp parallel default(shared)   &
        !$omp private(k)
        !$omp do
          do k = 1, body
            d(1:6,k,j) = tmp1 * d(1:6,k,j1)  -  tmp2 * d(1:6,k,j)
          end do
         !$omp end do nowait
!massless   
        !$omp do
          do k =  bpcp1, tm0c
            d(1:6,k,j) = tmp1 * d(1:6,k,j1)  -  tmp2 * d(1:6,k,j)
          end do
          !$omp end do
         !$omp end parallel 
       end do
!
!After several integrations, test the relative error on extrapolated values
        if (n.gt.3) then
          errmax = 0._double
!
!Maximum relative position and velocity errors (last D term added)
!massive
          do k = 1, body
            tmp1 = max( d(1,k,1)*d(1,k,1), d(2,k,1)*d(2,k,1), &
                       d(3,k,1)*d(3,k,1))
            tmp2 = max( d(4,k,1)*d(4,k,1), d(5,k,1)*d(5,k,1), &
                       d(6,k,1)*d(6,k,1))
            errmax = max(errmax, tmp1*xscal(k), tmp2*vscal(k))
          end do
!massless
          do k =  bpcp1, tm0c
            tmp1 = max( d(1,k,1)*d(1,k,1), d(2,k,1)*d(2,k,1), &
                       d(3,k,1)*d(3,k,1))
            tmp2 = max( d(4,k,1)*d(4,k,1), d(5,k,1)*d(5,k,1), &
                       d(6,k,1)*d(6,k,1))
            errmax = max(errmax, tmp1*xscal(k), tmp2*vscal(k))
          end do

!If error is smaller than TOL, update position and velocity arrays, and exit
          if (errmax.le.tol2) then
!massive
            !$omp parallel default(shared)   &
            !$omp private(k)
            !$omp do
            do k = 1, body
              rv(k,1:6) = d(1:6,k,1)
            end do
            !$omp end do nowait
!massless
           !$omp do
            do k =  bpcp1, tm0c
              rv(k,1:6) = d(1:6,k,1)
            end do
            !$omp end do
            !$omp end parallel 

         
           do j = 2, n
!massive
            !$omp parallel default(shared)   &
            !$omp private(k)
            !$omp do 
            do k = 1, body
                rv(k,1:6) = rv(k,1:6) + d(1:6,k,j)
            end do
            !$omp end do nowait
!massless
             !$omp do
            do k =  bpcp1, tm0c
                rv(k,1:6) = rv(k,1:6) + d(1:6,k,j)
            end do
            !$omp end do 
            !$omp end parallel 
          end do
!
        
!Save the actual stepsize used
            hdid = h0

         
!Recommend a new stepsize for the next call to this subroutine
            if (n.eq.bsrep) then 
              hnext = h0 * SHRINK
              ! write(*,*)'helloexit',n,h0
            end if

            if (n.lt.bsrep1) then
              hnext = h0 * GROW
            end if
!
 
            exit

         end if
      end if
!
   end do
  
          if (n.ge.8) then
!If errors were too large, redo the step with half the previous step size.
              h0 = h0 * 0.5_double
              hred=hred+1
           end if

 end do
 
            return
            end subroutine mdt_bs
!//////////////////////////////////////////////////////////////////////


!************************************************************************
subroutine frho(rv,totbody,body,cutbody,m0cutbody,mass,force,d) 
implicit none
integer(kind=intdble)::totbody,body,m0cutbody,i,j,l,tm0c,cutbody,bpcp1
real(kind=double)::rv(1:totbody,1:6),force(1:totbody,1:3),mass(1:totbody)
real(kind=double)::rij(1:totbody,1:body,1:3),d(:,:)

force(:,:)=0._double

tm0c=totbody-m0cutbody
bpcp1=body+cutbody+1

!massive pairdistances
!$omp parallel default(NONE)   &
!$omp shared(rij,rv,d,body,bpcp1,tm0c,force,mass) &
!$omp private(i,j,l)

!$omp do 
do i=1,body-1
   do j=i+1,body
      do l=1,3
         rij(i,j,l)=rv(j,l)-rv(i,l)
      end do
         rij(j,i,:)=-rij(i,j,:)
!norm of pair distances
        d(i,j)=Sqrt(Dot_Product(rij(i,j,:),rij(i,j,:)))
        d(j,i)=d(i,j)
    end do
end do
!$omp end do nowait

!$omp do
do i=1,body
   rij(i,i,:)=0._double
end do
!$omp end do nowait

!massless distances to massive
!$omp do 
do i=bpcp1,tm0c
   do j=1,body
      do l=1,3
         rij(i,j,l)=rv(j,l)-rv(i,l)
      end do
!norm of pair distances
        d(i,j)=Sqrt(Dot_Product(rij(i,j,:),rij(i,j,:)))
    end do
end do
!$omp end do

!massive: accelerations

!$omp do
   do i=1,body
      do j=1,body!massive on massive
         if(i.eq.j) then
            else
          force(i,:)=force(i,:)-mass(j)*rij(j,i,:)/d(j,i)**3._double
         end if
      end do

      do j=bpcp1,tm0c !massive on massless
          force(i,:)=force(i,:)-mass(j)*rij(j,i,:)/d(j,i)**3._double
      end do
      
   end do
!$omp end do nowait


!massless: accelerations
!$omp do 
   do i=bpcp1,tm0c
      do j=1,body
          force(i,:)=force(i,:)+mass(j)*rij(i,j,:)/d(i,j)**3._double
      end do
   end do
!$omp end do
!$omp end parallel
return
end subroutine frho


!************************************************************
end  module bs_m
