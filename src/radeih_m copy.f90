 module radeih_m
!contains Gauss Radau quadrature 15th order 5 (Everhart 1974, 1985)
! with Einstein Inffeld Hoffman PPN N-Body metric  

  use global_m
  use transform_m
  use output_m
  use cutoffmerge_m

  implicit none

public::Int_Radau_eih
public::Int_Radau_yeih
 contains

!***************************************************************************
!###########################################################
subroutine Int_Radau_eih(p,rv,mass,add)
implicit none
integer(kind=intdble)::allocstat
integer(kind=intdble)::cnt,m0body,body,dim,bodyp1
integer(kind=intdble)::cutbody,m0cutbody,totbody,mbody


type(Problem),intent(inout)::p 
type(additional),intent(in)::add

integer(kind=intdble)::order(1:p%NrKoerper)
real(kind=double),dimension(1:p%NrKoerper),intent(inout)::mass
real(kind=double),dimension(1:p%NrKoerper,1:6)::rv
real(kind=double)::tmod,dtmin
real(kind=double)::dtk,t,dtkmax,dtkmin,dtdid
real(kind=double),allocatable,dimension(:,:,:)::c_b,c_e
real(kind=double),allocatable,dimension(:,:)::rho
real(kind=double),allocatable,dimension(:)::c_h,c_vc,c_xc,c_c,c_d,c_r

real(kind=double)::kahan_y,kahan_t,kahan_c

logical::mflag,m0flag
logical,dimension(:,:),allocatable::merger

integer(kind=intdble)::i,n,j,stiffcnt



t=0._double
kahan_y=0._double
kahan_t=0._double
kahan_c=0._double

 cnt=1
 stiffcnt=0
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

!RA15 constants
allocate(c_b(7,totbody,3),c_e(7,totbody,3),c_h(8),c_xc(8),c_vc(7),c_c(21),c_d(21),c_r(28),rho(totbody,body), & 
         merger(totbody,body),stat=allocstat)

  if (allocstat.ne.0) then
            write(*,*)'error in radau_m: allocation not possible, try reducing the number of particles'
  end if 

merger(:,:)=.false.

! Gauss-Radau spacings for substeps within a sequence, for the 15th order 
! integrator. The sum of the H values should be 3.733333333333333
!
  c_h(:)=(/          0._double,.0562625605269221464656522_double,.1802406917368923649875799_double,  &
       .3526247171131696373739078_double,.5471536263305553830014486_double,.7342101772154105315232106_double,  &
       .8853209468390957680903598_double,.9775206135612875018911745_double/)
!
! Constant coefficients used in series expansions for X and V
!  XC: 1/2,  1/6,  1/12, 1/20, 1/30, 1/42, 1/56, 1/72
!  VC: 1/2,  1/3,  1/4,  1/5,  1/6,  1/7,  1/8
 
!   c_xc(:)=(/.5_double,.1666666666666667_double,.08333333333333333_double,.05_double,  &
!        .03333333333333333_double,.02380952380952381_double,.01785714285714286_double, &
!        .01388888888888889_double/)
! 
!   c_vc(:)=(/.5d0,.3333333333333333_double,.25_double,.2_double, &
!         .1666666666666667_double,.1428571428571429_double,.125_double/)
  c_xc(1)=0.5_double
  c_xc(2)=1._double/6._double
  c_xc(3)=1._double/12._double
  c_xc(4)=0.05_double
  c_xc(5)=1._double/30._double
  c_xc(6)=1._double/42._double
  c_xc(7)=1._double/56._double
  c_xc(8)=1._double/72._double

 c_vc(1)=0.5_double
 c_vc(2)=1._double/3._double
 c_vc(3)=0.25_double
 c_vc(4)=0.2_double
 c_vc(5)=1._double/6._double
 c_vc(6)=1._double/7._double
 c_vc(7)=0.125_double

  c_b(:,:,:)=0._double
  c_e(:,:,:)=0._double
!
! set values of the constant arrays
! (R = R21, R31, R32, R41, R42, R43 in Everhart's paper.)
    
        n = 0
        do j = 2, 8
          do i = 1, j - 1
            n = n + 1
            c_r(n) = 1._double / (c_h(j) - c_h(i))
          end do
        end do
!
! Constants to convert between B and G arrays (C = C21, C31, C32, C41, C42...)
        c_c(1) = - c_h(2)
        c_d(1) =   c_h(2)
        n = 1
        do j = 3, 7
          n = n + 1
          c_c(n) = -c_h(j) * c_c(n-j+2)
          c_d(n) =  c_h(2) * c_d(n-j+2)
          do i = 3, j - 1
            n = n + 1
            c_c(n) = c_c(n-j+1)  -  c_h(j) * c_c(n-j+2)
            c_d(n) = c_d(n-j+1)  +  c_h(i) * c_d(n-j+2)
          end do
          n = n + 1
          c_c(n) = c_c(n-j+1) - c_h(j)
          c_d(n) = c_d(n-j+1) + c_h(j)
        end do



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
! 
call radau_eih_step(rv,mass,totbody,body,cutbody,m0cutbody,dtk,dtdid,p%eps,rho,c_b, &
                c_e,c_h,c_xc,c_vc,c_c,c_d,c_r)
!
   !t=t+dtdid/k

   kahan_y=dtdid/k-kahan_c
   kahan_t=t+kahan_y
   kahan_c=(kahan_t-t)-kahan_y
   t=kahan_t

    if(dtdid.gt.dtkmax) then
      dtkmax=dtdid
    end if
    if (dtdid.le.dtkmin) then
       dtkmin=dtdid
       stiffcnt=stiffcnt+1
      else 
       stiffcnt=0 
   end if

   if (stiffcnt.gt.10000) then
     write(*,*)'possibly stiff region encountered, stopping'
     STOP
   end if 


 

   !Ausgabe 
   if (mod(real(t/p%tAusg,KIND=double),1._double) .lt.1.d-15.and.t.ge.p%tAusg .or. dble(p%tAusg).eq.0.d0)then
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
 cnt=cnt+1
end do prgs

write(*,*)'time accomplished:',t,'with an average stepsize of:',t/dble(cnt)
write(*,*)'minimum stepsize',dtkmin/k
write(*,*)'maximum stepsize',dtkmax/k


write(*,*)'number of massive particles at the beginning / end of integration:',body+cutbody,' / ',body
write(*,*)'number of massless particles at the beginning / end of integration:',m0body+m0cutbody,' / ',m0body
if(p%outegy) then
  
   write(*,*)'error in total energy (sum of dE)',dEtot
   write(*,*)'error in total angular momentum (sum of dL)',dLtot
end if

close(41)


return
end subroutine Int_Radau_eih
!******************************************************************************
!
! Author: John E. Chambers
! Modified: Siegfried Eggl  201104
!
! Integrates bodies for one timestep TDID using
! Everhart's RA15 integrator algorithm. The accelerations are calculated
! using the subroutine FROH. The accuracy of the step is approximately 
! determined by the tolerance parameter TOL.
!
! Based on RADAU by E. Everhart, Physics Department, University of Denver.
! Comments giving equation numbers refer to Everhart (1985) ``An Efficient
! Integrator that Uses Gauss-Radau Spacings'', in The Dynamics of Comets:
! Their Origin and Evolution, p185-202, eds. A. Carusi & G. B. Valsecchi,
! pub Reidel. (A listing of the original subroutine is also given in this 
! paper.)
! 
!------------------------------------------------------------------------------
      subroutine radau_eih_step(rv,mass,totbody,body,cutbody,m0cutbody,& 
                            t,tdid,tol,rho,b,e,h,xc,vc,c,d,r)

      implicit none
!
! Input/Output
      integer(kind=intdble):: totbody,body,cutbody,m0cutbody 
      real(kind=double):: t,tdid,tol
      real(kind=double):: rv(:,:),rho(:,:),mass(:)  
!
! Local
      integer(kind=intdble):: j,k,n,i
      integer(kind=intdble):: bmc,bp1,totm0 
      real(kind=double)::g(7,totbody,3),b(7,totbody,3),e(7,totbody,3)
      real(kind=double)::h(8),xc(8),vc(7),c(21),d(21),r(28),s(9)
      real(kind=double)::q,q2,q3,q4,q5,q6,q7,temp,gk,gkmax
      real(kind=double),dimension(totbody,3)::acc,acc1
      real(kind=double)::rv1(totbody,6),temp3(3)
      logical::iter
!
bmc=body-cutbody
bp1=body+1
totm0=totbody-m0cutbody
!
! Calculate forces at the start of the sequence
      call frho_eih(rv,totbody,bmc,bp1,totm0,mass,acc,rho) 
!
! Find G values from B values predicted at the last call (Eqs. 7 of Everhart)
!
!massive

!$omp parallel default(NONE) &
!$omp shared(g,bmc,bp1,totm0,b,d) &
!$omp private(k)
!$omp do
      do  k= 1, bmc
        g(1,k,:) = b(7,k,:)*d(16) + b(6,k,:)*d(11) + b(5,k,:)*d(7) &
                 + b(4,k,:)*d(4)  + b(3,k,:)*d(2)  + b(2,k,:)*d(1)  + b(1,k,:)
        g(2,k,:) = b(7,k,:)*d(17) + b(6,k,:)*d(12) + b(5,k,:)*d(8) &
                 + b(4,k,:)*d(5)  + b(3,k,:)*d(3)  + b(2,k,:)
        g(3,k,:) = b(7,k,:)*d(18) + b(6,k,:)*d(13) + b(5,k,:)*d(9)  &
                 + b(4,k,:)*d(6)  + b(3,k,:)
        g(4,k,:) = b(7,k,:)*d(19) + b(6,k,:)*d(14) + b(5,k,:)*d(10) + b(4,k,:)
        g(5,k,:) = b(7,k,:)*d(20) + b(6,k,:)*d(15) + b(5,k,:)
        g(6,k,:) = b(7,k,:)*d(21) + b(6,k,:)
        g(7,k,:) = b(7,k,:)
     end do
!$omp end do nowait

!massless
!$omp do 
     do  k=bp1,totm0
        g(1,k,:) = b(7,k,:)*d(16) + b(6,k,:)*d(11) + b(5,k,:)*d(7) &
                 + b(4,k,:)*d(4)  + b(3,k,:)*d(2)  + b(2,k,:)*d(1)  + b(1,k,:)
        g(2,k,:) = b(7,k,:)*d(17) + b(6,k,:)*d(12) + b(5,k,:)*d(8) &
                 + b(4,k,:)*d(5)  + b(3,k,:)*d(3)  + b(2,k,:)
        g(3,k,:) = b(7,k,:)*d(18) + b(6,k,:)*d(13) + b(5,k,:)*d(9)  &
                 + b(4,k,:)*d(6)  + b(3,k,:)
        g(4,k,:) = b(7,k,:)*d(19) + b(6,k,:)*d(14) + b(5,k,:)*d(10) + b(4,k,:)
        g(5,k,:) = b(7,k,:)*d(20) + b(6,k,:)*d(15) + b(5,k,:)
        g(6,k,:) = b(7,k,:)*d(21) + b(6,k,:)
        g(7,k,:) = b(7,k,:)
     end do
!$omp end do
!$omp end parallel
!------------------------------------------------------------------------------
!
!  MAIN  LOOP  STARTS  HERE
!
! For each iteration (six for first call to subroutine, two otherwise)...
 
!   do n = 1, niter
    n=0
    iter=.true.
    do while(iter)
!    do while(n<6)
!    
    n=n+1    
! For each substep within a sequence...
        do j = 2, 8
!
!$omp parallel default(NONE) &
!$omp shared(h,rv1,rv,b,acc,bmc,bp1,totm0,t,tol,iter,j) &
!$omp private(s,k)
! Calculate position predictors using Eqn. 9 of Everhart
          s(1) = t * h(j)
          s(2) = s(1) * s(1) * .5_double
          s(3) = s(2) * h(j) *  1._double/3._double
          s(4) = s(3) * h(j) * .5_double
          s(5) = s(4) * h(j) * .6_double
          s(6) = s(5) * h(j) *  2._double/3._double
          s(7) = s(6) * h(j) * 5._double/7._double
          s(8) = s(7) * h(j) * .75_double
          s(9) = s(8) * h(j) *  7._double/9._double
!
!massive
!$omp do
          do k = 1, bmc           
            rv1(k,1:3) = s(9)*b(7,k,:) + s(8)*b(6,k,:) + s(7)*b(5,k,:) &
                       + s(6)*b(4,k,:) + s(5)*b(3,k,:) + s(4)*b(2,k,:) &
                       + s(3)*b(1,k,:) + s(2)*acc(k,:)  + s(1)*rv(k,4:6) + rv(k,1:3)
          end do
!$omp end do
!massless
!$omp do
         do k = bp1,totm0           
            rv1(k,1:3) = s(9)*b(7,k,:) + s(8)*b(6,k,:) + s(7)*b(5,k,:) &
                       + s(6)*b(4,k,:) + s(5)*b(3,k,:) + s(4)*b(2,k,:) &
                       + s(3)*b(1,k,:) + s(2)*acc(k,:)  + s(1)*rv(k,4:6) + rv(k,1:3)
          end do
!$omp end do
!
! If necessary, calculate velocity predictors too, from Eqn. 10 of Everhart
            s(1) = t * h(j)
            s(2) = s(1) * h(j) * .5_double
            s(3) = s(2) * h(j) * 2._double/3._double
            s(4) = s(3) * h(j) * .75_double
            s(5) = s(4) * h(j) * .8_double
            s(6) = s(5) * h(j) * 5._double/6._double
            s(7) = s(6) * h(j) * 6._double/7._double
            s(8) = s(7) * h(j) * .875_double
!massive
!$omp do
            do k = 1, bmc
              rv1(k,4:6) = s(8)*b(7,k,:) + s(7)*b(6,k,:) + s(6)*b(5,k,:) &
                         + s(5)*b(4,k,:) + s(4)*b(3,k,:) + s(3)*b(2,k,:)    &
                         + s(2)*b(1,k,:) + s(1)*acc(k,:)  + rv(k,4:6)
            end do
!$omp end do nowait
!massless
!$omp do
            do k = bp1,totm0
              rv1(k,4:6) = s(8)*b(7,k,:) + s(7)*b(6,k,:) + s(6)*b(5,k,:) &
                         + s(5)*b(4,k,:) + s(4)*b(3,k,:) + s(3)*b(2,k,:)    &
                         + s(2)*b(1,k,:) + s(1)*acc(k,:)  + rv(k,4:6)
            end do 
!$omp end do
!$omp end parallel

!
! Calculate forces at the current substep
          call frho_eih(rv1,totbody,bmc,bp1,totm0,mass,acc1,rho) 
! Update G values using Eqs. 4 of Everhart, and update B values using Eqs. 5
 
  select case(j)
     case(2)
!$omp parallel default(NONE) &
!$omp shared(b,r,g,acc,acc1,bmc,bp1,totm0,j) &
!$omp private(k,temp3)
!massive
!$omp do
            do k = 1, bmc
              temp3 = g(1,k,:)
              g(1,k,:) = (acc1(k,:) - acc(k,:)) * r(1)
              b(1,k,:) = b(1,k,:) + g(1,k,:) - temp3
            end do
!$omp end do nowait
!massless
!$omp do
            do k = bp1,totm0 
              temp3 = g(1,k,:)
              g(1,k,:) = (acc1(k,:) - acc(k,:)) * r(1)
              b(1,k,:) = b(1,k,:) + g(1,k,:) - temp3
            end do
!$omp end do
!$omp end parallel

      case(3)
!$omp parallel default(NONE) &
!$omp shared(b,c,r,g,acc,acc1,bmc,bp1,totm0,j) &
!$omp private(i,k,temp3,gk,temp)
!massive
!$omp do
            do k = 1,bmc
              do i=1,3             
               temp = g(2,k,i)
               gk = acc1(k,i) - acc(k,i)
               g(2,k,i) = (gk*r(2) - g(1,k,i))*r(3)
               temp = g(2,k,i) - temp
               b(1,k,i) = b(1,k,i)  +  temp * c(1)
               b(2,k,i) = b(2,k,i)  +  temp
              end do
            end do
!$omp end do nowait
!massless
!$omp do
            do k = bp1,totm0
              do i=1,3             
               temp = g(2,k,i)
               gk = acc1(k,i) - acc(k,i)
               g(2,k,i) = (gk*r(2) - g(1,k,i))*r(3)
               temp = g(2,k,i) - temp
               b(1,k,i) = b(1,k,i)  +  temp * c(1)
               b(2,k,i) = b(2,k,i)  +  temp
              end do
            end do
!$omp end do
!$omp end parallel
      case(4)
!$omp parallel default(NONE) &
!$omp shared(b,r,c,g,acc,acc1,bmc,bp1,totm0,j) &
!$omp private(i,k,temp3,gk,temp)
!massive
!$omp do
            do k = 1,bmc
             do i=1,3
              temp = g(3,k,i)
              gk = acc1(k,i) - acc(k,i)
              g(3,k,i) = ((gk*r(4) - g(1,k,i))*r(5) - g(2,k,i))*r(6)
              temp = g(3,k,i) - temp
              b(1,k,i) = b(1,k,i)  +  temp * c(2)
              b(2,k,i) = b(2,k,i)  +  temp * c(3)
              b(3,k,i) = b(3,k,i)  +  temp
             end do
            end do
!$omp end do nowait
!massless
!$omp do
           do k = bp1,totm0
             do i=1,3
              temp = g(3,k,i)
              gk = acc1(k,i) - acc(k,i)
              g(3,k,i) = ((gk*r(4) - g(1,k,i))*r(5) - g(2,k,i))*r(6)
              temp = g(3,k,i) - temp
              b(1,k,i) = b(1,k,i)  +  temp * c(2)
              b(2,k,i) = b(2,k,i)  +  temp * c(3)
              b(3,k,i) = b(3,k,i)  +  temp
             end do
            end do
!$omp end do
!$omp end parallel
       case(5)
!$omp parallel default(NONE) &
!$omp shared(b,c,r,g,acc,acc1,bmc,bp1,totm0,j) &
!$omp private(i,k,temp3,gk,temp)
!massive
!$omp do
            do k = 1, bmc
             do i=1,3
              temp = g(4,k,i)
              gk = acc1(k,i) - acc(k,i)
              g(4,k,i) = (((gk*r(7) - g(1,k,i))*r(8) - g(2,k,i))*r(9) &
                    - g(3,k,i))*r(10)
              temp = g(4,k,i) - temp
              b(1,k,i) = b(1,k,i)  +  temp * c(4)
              b(2,k,i) = b(2,k,i)  +  temp * c(5)
              b(3,k,i) = b(3,k,i)  +  temp * c(6)
              b(4,k,i) = b(4,k,i)  +  temp
             end do
            end do
!$omp end do nowait
!massless
!$omp do
            do k = bp1,totm0
             do i=1,3
              temp = g(4,k,i)
              gk = acc1(k,i) - acc(k,i)
              g(4,k,i) = (((gk*r(7) - g(1,k,i))*r(8) - g(2,k,i))*r(9) &
                    - g(3,k,i))*r(10)
              temp = g(4,k,i) - temp
              b(1,k,i) = b(1,k,i)  +  temp * c(4)
              b(2,k,i) = b(2,k,i)  +  temp * c(5)
              b(3,k,i) = b(3,k,i)  +  temp * c(6)
              b(4,k,i) = b(4,k,i)  +  temp
             end do
            end do
!$omp end do
!$omp end parallel
        case(6)
!$omp parallel default(NONE) &
!$omp shared(b,c,r,g,acc,acc1,bmc,bp1,totm0,j) &
!$omp private(i,k,gk,temp3,temp)
!massive
!$omp do
            do k = 1, bmc
             do i=1,3
              temp = g(5,k,i)
              gk = acc1(k,i) - acc(k,i)
              g(5,k,i) =  ((((gk*r(11) - g(1,k,i))*r(12) - g(2,k,i))*r(13) &
                    - g(3,k,i))*r(14) - g(4,k,i))*r(15)
              temp = g(5,k,i) - temp
              b(1,k,i) = b(1,k,i)  +  temp * c(7)
              b(2,k,i) = b(2,k,i)  +  temp * c(8)
              b(3,k,i) = b(3,k,i)  +  temp * c(9)
              b(4,k,i) = b(4,k,i)  +  temp * c(10)
              b(5,k,i) = b(5,k,i)  +  temp
             end do
            end do
!$omp end do nowait
!massless
!$omp do
           do k = bp1,totm0
             do i=1,3
              temp = g(5,k,i)
              gk = acc1(k,i) - acc(k,i)
              g(5,k,i) =  ((((gk*r(11) - g(1,k,i))*r(12) - g(2,k,i))*r(13) &
                    - g(3,k,i))*r(14) - g(4,k,i))*r(15)
              temp = g(5,k,i) - temp
              b(1,k,i) = b(1,k,i)  +  temp * c(7)
              b(2,k,i) = b(2,k,i)  +  temp * c(8)
              b(3,k,i) = b(3,k,i)  +  temp * c(9)
              b(4,k,i) = b(4,k,i)  +  temp * c(10)
              b(5,k,i) = b(5,k,i)  +  temp
             end do
           end do
!$omp end do
!$omp end parallel
         case(7)
!$omp parallel default(NONE) &
!$omp shared(b,c,r,g,acc,acc1,bmc,bp1,totm0,j) &
!$omp private(i,k,temp3,gk,temp)
!massive
!$omp do
            do k = 1, bmc
             do i=1,3
              temp = g(6,k,i)
              gk = acc1(k,i) - acc(k,i)
              g(6,k,i) = (((((gk*r(16) - g(1,k,i))*r(17) - g(2,k,i))*r(18) &
                     - g(3,k,i))*r(19) - g(4,k,i))*r(20) - g(5,k,i))*r(21)
              temp = g(6,k,i) - temp
              b(1,k,i) = b(1,k,i)  +  temp * c(11)
              b(2,k,i) = b(2,k,i)  +  temp * c(12)
              b(3,k,i) = b(3,k,i)  +  temp * c(13)
              b(4,k,i) = b(4,k,i)  +  temp * c(14)
              b(5,k,i) = b(5,k,i)  +  temp * c(15)
              b(6,k,i) = b(6,k,i)  +  temp
            end do
          end do
!$omp end do nowait
!massless
!$omp do
          do k = bp1,totm0
             do i=1,3
              temp = g(6,k,i)
              gk = acc1(k,i) - acc(k,i)
              g(6,k,i) = (((((gk*r(16) - g(1,k,i))*r(17) - g(2,k,i))*r(18) &
                     - g(3,k,i))*r(19) - g(4,k,i))*r(20) - g(5,k,i))*r(21)
              temp = g(6,k,i) - temp
              b(1,k,i) = b(1,k,i)  +  temp * c(11)
              b(2,k,i) = b(2,k,i)  +  temp * c(12)
              b(3,k,i) = b(3,k,i)  +  temp * c(13)
              b(4,k,i) = b(4,k,i)  +  temp * c(14)
              b(5,k,i) = b(5,k,i)  +  temp * c(15)
              b(6,k,i) = b(6,k,i)  +  temp
            end do
          end do
!$omp end do
!$omp end parallel

         case(8)
           gkmax=0._double
!massive
            do k = 1, bmc
             do i=1,3
              temp = g(7,k,i)
              gk = acc1(k,i) - acc(k,i)
              g(7,k,i) = ((((((gk*r(22) - g(1,k,i))*r(23) - g(2,k,i))*r(24) &
                     - g(3,k,i))*r(25) - g(4,k,i))*r(26) - g(5,k,i))*r(27) &
                     - g(6,k,i))*r(28)
              temp = g(7,k,i) - temp
              b(1,k,i) = b(1,k,i)  +  temp * c(16)
              b(2,k,i) = b(2,k,i)  +  temp * c(17)
              b(3,k,i) = b(3,k,i)  +  temp * c(18)
              b(4,k,i) = b(4,k,i)  +  temp * c(19)
              b(5,k,i) = b(5,k,i)  +  temp * c(20)
              b(6,k,i) = b(6,k,i)  +  temp * c(21)
              b(7,k,i) = b(7,k,i)  +  temp
              gkmax=max(gkmax,abs(temp))             
             end do !i 1..3
            end do  !k 1..nbody
!massless
             do k = bp1,totm0
             do i=1,3
              temp = g(7,k,i)
              gk = acc1(k,i) - acc(k,i)
              g(7,k,i) = ((((((gk*r(22) - g(1,k,i))*r(23) - g(2,k,i))*r(24) &
                     - g(3,k,i))*r(25) - g(4,k,i))*r(26) - g(5,k,i))*r(27) &
                     - g(6,k,i))*r(28)
              temp = g(7,k,i) - temp
              b(1,k,i) = b(1,k,i)  +  temp * c(16)
              b(2,k,i) = b(2,k,i)  +  temp * c(17)
              b(3,k,i) = b(3,k,i)  +  temp * c(18)
              b(4,k,i) = b(4,k,i)  +  temp * c(19)
              b(5,k,i) = b(5,k,i)  +  temp * c(20)
              b(6,k,i) = b(6,k,i)  +  temp * c(21)
              b(7,k,i) = b(7,k,i)  +  temp
              gkmax=max(gkmax,abs(temp))             
             end do !i 1..3
            end do  !k 1..nbody

!more iterations necessary?
                if(gkmax.gt.tol) then
                  iter=.true.
               else
                  iter=.false.
               end if
!              write(*,*)gkmax,n
        end select
      end do      !j 2...8
!     rerun subroutine with smaller stepsize if things do not converge      
      if(n>6) then
         t=t/3.13_double
         tdid=0._double
!         write(*,*)'rvint',rv(2,:)
!         write(*,*)'no convergence, trying smaller initial stepsize!',t*3.13d0,'->',t        
         return
      end if
     end do      !n... while(iter)
!
!------------------------------------------------------------------------------
!
!  END  OF  MAIN  LOOP
!


! Estimate suitable sequence size for the next call to subroutine (Eqs. 15, 16)
      temp = 0._double

!massive
      do k = 1, bmc
       do i=1,3
        temp = max( temp, abs( b(7,k,i) ) )
       end do
      end do
!massless
      do k = bp1,totm0
       do i=1,3
        temp = max( temp, abs( b(7,k,i) ) )
       end do
      end do
!
      temp = temp / (72._double * abs(t)**7)
      tdid = t
      if (temp.eq.0) then
        t = tdid * 1.4_double
      else
        t = sign( (tol/temp)**(1._double/9._double), tdid )
      end if
!
! If new sequence size is much bigger than the current one, reduce it
      if (abs(t/tdid).gt.1.4_double) t = tdid * 1.4_double
!
! Find new position and velocity values at end of the sequence (Eqs. 11, 12)

      temp = tdid * tdid

!$omp parallel default(NONE) &
!$omp shared(rv,xc,vc,b,acc,bmc,bp1,totm0,t,tdid,temp) &
!$omp private(k)
!massive
!$omp do
      do k = 1 , bmc
        rv(k,1:3) = (xc(8)*b(7,k,:) + xc(7)*b(6,k,:) + xc(6)*b(5,k,:) &
                  +  xc(5)*b(4,k,:) + xc(4)*b(3,k,:) + xc(3)*b(2,k,:) &
                  +  xc(2)*b(1,k,:) + xc(1)*acc(k,:))*temp + rv(k,4:6)*tdid + rv(k,1:3)
!
        rv(k,4:6) = (vc(7)*b(7,k,:) + vc(6)*b(6,k,:) + vc(5)*b(5,k,:) &
                  +  vc(4)*b(4,k,:) + vc(3)*b(3,k,:) + vc(2)*b(2,k,:)     &
                  +  vc(1)*b(1,k,:) + acc(k,:))*tdid + rv(k,4:6)
      end do
!$omp end do nowait
!massless
!$omp do
      do k = bp1,totm0
        rv(k,1:3) = (xc(8)*b(7,k,:) + xc(7)*b(6,k,:) + xc(6)*b(5,k,:) &
                  +  xc(5)*b(4,k,:) + xc(4)*b(3,k,:) + xc(3)*b(2,k,:) &
                  +  xc(2)*b(1,k,:) + xc(1)*acc(k,:))*temp + rv(k,4:6)*tdid + rv(k,1:3)
!
        rv(k,4:6) = (vc(7)*b(7,k,:) + vc(6)*b(6,k,:) + vc(5)*b(5,k,:) &
                  +  vc(4)*b(4,k,:) + vc(3)*b(3,k,:) + vc(2)*b(2,k,:)     &
                  +  vc(1)*b(1,k,:) + acc(k,:))*tdid + rv(k,4:6)
      end do
!$omp end do
!$omp end parallel

! Predict new B values to use at the start of the next sequence. The predicted
! values from the last call are saved as E. The correction, BD, between the
! actual and predicted values of B is applied in advance as a correction.
      q = t / tdid
      q2 = q  * q
      q3 = q  * q2
      q4 = q2 * q2
      q5 = q2 * q3
      q6 = q3 * q3
      q7 = q3 * q4
!
!massive
      do k = 1, bmc
       do i =1,3
        s(1) = b(1,k,i) - e(1,k,i)
        s(2) = b(2,k,i) - e(2,k,i)
        s(3) = b(3,k,i) - e(3,k,i)
        s(4) = b(4,k,i) - e(4,k,i)
        s(5) = b(5,k,i) - e(5,k,i)
        s(6) = b(6,k,i) - e(6,k,i)
        s(7) = b(7,k,i) - e(7,k,i)
!
! Estimate B values for the next sequence (Eqs. 13 of Everhart).
        e(1,k,i) = q* (b(7,k,i)* 7._double + b(6,k,i)* 6._double + b(5,k,i)* 5._double &
               +       b(4,k,i)* 4._double + b(3,k,i)* 3._double + b(2,k,i)*2._double + b(1,k,i))
        e(2,k,i) = q2*(b(7,k,i)*21._double + b(6,k,i)*15._double + b(5,k,i)*10._double &
               +       b(4,k,i)* 6._double + b(3,k,i)* 3._double + b(2,k,i))      
        e(3,k,i) = q3*(b(7,k,i)*35._double + b(6,k,i)*20._double + b(5,k,i)*10._double  &
               +       b(4,k,i)*4._double  + b(3,k,i))
        e(4,k,i) = q4*(b(7,k,i)*35._double + b(6,k,i)*15._double + b(5,k,i)*5._double + b(4,k,i))
        e(5,k,i) = q5*(b(7,k,i)*21._double + b(6,k,i)*6._double  + b(5,k,i))
        e(6,k,i) = q6*(b(7,k,i)*7._double  + b(6,k,i))
        e(7,k,i) = q7* b(7,k,i)
!
        b(1,k,i) = e(1,k,i) + s(1)
        b(2,k,i) = e(2,k,i) + s(2)
        b(3,k,i) = e(3,k,i) + s(3)
        b(4,k,i) = e(4,k,i) + s(4)
        b(5,k,i) = e(5,k,i) + s(5)
        b(6,k,i) = e(6,k,i) + s(6)
        b(7,k,i) = e(7,k,i) + s(7)
       end do
      end do
!massless
     do k = bp1,totm0
       do i =1,3
        s(1) = b(1,k,i) - e(1,k,i)
        s(2) = b(2,k,i) - e(2,k,i)
        s(3) = b(3,k,i) - e(3,k,i)
        s(4) = b(4,k,i) - e(4,k,i)
        s(5) = b(5,k,i) - e(5,k,i)
        s(6) = b(6,k,i) - e(6,k,i)
        s(7) = b(7,k,i) - e(7,k,i)
!
! Estimate B values for the next sequence (Eqs. 13 of Everhart).
        e(1,k,i) = q* (b(7,k,i)* 7._double + b(6,k,i)* 6._double + b(5,k,i)* 5._double &
               +       b(4,k,i)* 4._double + b(3,k,i)* 3._double + b(2,k,i)*2._double + b(1,k,i))
        e(2,k,i) = q2*(b(7,k,i)*21._double + b(6,k,i)*15._double + b(5,k,i)*10._double &
               +       b(4,k,i)* 6._double + b(3,k,i)* 3._double + b(2,k,i))      
        e(3,k,i) = q3*(b(7,k,i)*35._double + b(6,k,i)*20._double + b(5,k,i)*10._double  &
               +       b(4,k,i)*4._double  + b(3,k,i))
        e(4,k,i) = q4*(b(7,k,i)*35._double + b(6,k,i)*15._double + b(5,k,i)*5._double + b(4,k,i))
        e(5,k,i) = q5*(b(7,k,i)*21._double + b(6,k,i)*6._double  + b(5,k,i))
        e(6,k,i) = q6*(b(7,k,i)*7._double  + b(6,k,i))
        e(7,k,i) = q7* b(7,k,i)
!
        b(1,k,i) = e(1,k,i) + s(1)
        b(2,k,i) = e(2,k,i) + s(2)
        b(3,k,i) = e(3,k,i) + s(3)
        b(4,k,i) = e(4,k,i) + s(4)
        b(5,k,i) = e(5,k,i) + s(5)
        b(6,k,i) = e(6,k,i) + s(6)
        b(7,k,i) = e(7,k,i) + s(7)
       end do
      end do
      return
      end subroutine 

!******************************************************************************************
subroutine frho_eih(rv,totbody,bmc,bp1,totm0,mass,frel,d) 
!Calculates relativistic forces from Einstein Infeld Hoffman PPN N-Body metric
!NICE MODEL DISABLED
implicit none
integer(kind=intdble)::totbody,bmc,bp1,totm0,i,j,l
real(kind=double)::rv(1:totbody,1:6),mass(1:totbody)
real(kind=double)::frel(1:totbody,1:3),fclas(1:totbody,1:3)
real(kind=double),dimension(1:totbody,1:bmc,1:3)::rij
real(kind=double)::subsum1,subsum2
real(kind=double)::d(:,:)
real(kind=double),parameter::k2=0.00029591220828559115_double !gaussian gravitational constant squared (k^2)
real(kind=double),parameter::km1=58.132440867048956_double    !1/k
real(kind=double),parameter::km2=3379.3806811609434_double    !1/k^2
real(kind=double),parameter::cm2=3.3356611871248456d-5        !1/c^2 vacuum light speed
real(kind=double),parameter::cm2k2=9.870628679746496d-9       !k^2/c^2


 
 
 subsum1=0._double
 subsum2=0._double

!$omp parallel default(shared)   &
!$omp firstprivate(k2,km1,km2,cm2,bmc,bp1,totm0) & 
!$omp private(i,j,l,subsum1,subsum2)
!$omp do 
do i=1,bmc
frel(i,:)=0._double
end do
!$omp end do nowait

!$omp do
do i=bp1,totm0
frel(i,:)=0._double
end do
!$omp end do nowait


!$omp do 
do i=1,bmc
fclas(i,:)=0._double
end do
!$omp end do nowait

!$omp do
do i=bp1,totm0
fclas(i,:)=0._double
end do
!$omp end do nowait

!massive pair distances
!$omp do
do i=1,bmc-1
   do j=i+1,bmc
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
do i=1,bmc
   rij(i,i,:)=0._double
end do
!$omp end do nowait


!massless distances to massive
!$omp do
do i=bp1,totm0
   do j=1,bmc
      do l=1,3
         rij(i,j,l)=rv(j,l)-rv(i,l)
      end do
!norm of pair distances
        d(i,j)=Sqrt(Dot_Product(rij(i,j,:),rij(i,j,:)))
    end do
end do
!$omp end do 


!calculate classical forces

!massive
!$omp do
 do i=1,bmc
      do j=1,bmc!massive on massive
         if(i.eq.j) then
            else
          fclas(i,1:3)=fclas(i,1:3)+mass(j)*rij(i,j,1:3)/d(i,j)**3
         end if
      end do

!      do j=bp1,totm0 !massless on massive
!           if (mass(j).eq.0.d0) then
!           else
!          subsum1(i)=subsum1(i)+mass(j)/d(j,i)
!           end if
!       end do  
  end do
!$omp end do 


!massless
!$omp do
 do i=bp1,totm0
      do j=1,bmc
        fclas(i,1:3)=fclas(i,1:3)+mass(j)*rij(i,j,1:3)/d(i,j)**3
      end do !j
 end do !i
!$omp end do



 !relativistic accelerations
!$omp do
 do i=1,bmc
      do j=1,bmc!massive on massive
         if(i.eq.j) then
            else
            do l=1,bmc
             if(i.eq.l) then
             else
              subsum1=subsum1+mass(l)/d(i,l)
             end if
             if(j.eq.l) then
             else
              subsum2=subsum2+mass(l)/d(j,l)
             end if
            end do !l
         !PLEASE TAKE CARE WRT THE RENORMALIZED TIME (T->t/k!!! so you have to take care of 1/c^2 too!)
   frel(i,1:3)=frel(i,1:3)+ rij(i,j,1:3)*mass(j)/d(i,j)**3 *&
           (1._double+cm2k2*(Dot_Product(rv(i,4:6),rv(i,4:6)) + 2._double*Dot_Product(rv(j,4:6),rv(j,4:6)) &
           - 4._double*Dot_Product(rv(i,4:6),rv(j,4:6)) -1.5_double*Dot_Product(rij(j,i,1:3)/d(i,j),rv(j,4:6))**2 &
           -4._double*subsum1-subsum2+0.5_double*Dot_Product(rij(i,j,1:3),fclas(j,1:3)))) &
           +cm2k2*mass(j)/d(i,j)**2*Dot_Product(rij(j,i,1:3)/d(i,j),(4.d0*rv(i,4:6)-3.d0*rv(j,4:6)))*(rv(i,4:6)-rv(j,4:6)) &
           +cm2k2*3.5_double*mass(j)*fclas(j,1:3)/d(i,j)
           
           subsum1=0._double
           subsum2=0._double
         end if
        end do !j
   
!    do j=bp1,totm0!massless on massive
!      if (mass(j).eq.0.d0) then
!      else
!      frel(i,1:3)=frel(i,1:3)-rij(j,i,1:3)*mass(j)/d(j,i)**3 *&
!              (1._double+cm2*k2*(-4._double*subsum1(i)-subsum2(j,i)-5._double*mass(i)/d(j,i)+ &  !!!I THINK THE PROBLEM IS HERE!!!!!
!              (Dot_Product(rv(i,4:6),rv(i,4:6))-4._double*Dot_Product(rv(i,4:6),rv(j,4:6))+ & 
!              2._double*Dot_Product(rv(j,4:6),rv(j,4:6)))- &
!              1.5_double *k2*(Dot_Product(rv(j,4:6),(-rij(j,i,1:3)))/d(j,i))**2)) +&
!              cm2*k2*(mass(j)*(Dot_Product(-rij(j,i,1:3), &
!              (4._double*rv(i,4:6)-3._double*rv(j,4:6)))/d(j,i)**3*(rv(j,4:6)-rv(i,4:6)))+ &
!              3.5_double *(-subsum3(j,i,1:3)))
!      end if
!    end do !j 
 end do !i
!$omp end do

subsum1=0._double
subsum2=0._double

!massless
!$omp do
 do i=bp1,totm0
      do j=1,bmc
          do l=1,bmc
              subsum1=subsum1+mass(l)/d(i,l)
             if(j.eq.l) then
             else
              subsum2=subsum2+mass(l)/d(j,l)
             end if
            end do !l
   frel(i,1:3)=frel(i,1:3)+ rij(i,j,1:3)*mass(j)/d(i,j)**3 *&
           (1._double+cm2k2*(Dot_Product(rv(i,4:6),rv(i,4:6)) + 2._double*Dot_Product(rv(j,4:6),rv(j,4:6)) &
           - 4._double*Dot_Product(rv(i,4:6),rv(j,4:6)) -1.5_double/d(i,j)*Dot_Product(-rij(i,j,1:3),rv(j,4:6))**2 &
           -4._double*subsum1-subsum2+0.5_double*Dot_Product(rij(i,j,1:3),fclas(j,1:3)))) &
           +cm2k2*mass(j)/d(i,j)**3*Dot_Product(-rij(i,j,1:3),(4.d0*rv(i,4:6)-3.d0*rv(j,4:6)))*(rv(i,4:6)-rv(j,4:6)) &
           +cm2k2*3.5_double*mass(j)*fclas(j,1:3)/d(i,j)
           
           subsum1=0._double
           subsum2=0._double
        end do !j
 end do !i
!$omp end do
!$omp end parallel



return
end subroutine frho_eih
!#################################################################
!          
! subroutine frho_eih(rv,totbody,bmc,bp1,totm0,mass,frel,d) 
! !Calculates relativistic forces from Einstein Infeld Hoffman PPN N-Body metric
! implicit none
! integer(kind=intdble)::totbody,bmc,bp1,totm0,i,j,l
! real(kind=double)::rv(1:totbody,1:6),mass(1:totbody)
! real(kind=double)::frel(1:totbody,1:3),km2,cm2,km1,k2
! real(kind=double),dimension(1:totbody,1:bmc,1:3)::rij,subsum3
! real(kind=double),dimension(1:totbody,1:bmc)::subsum2
! real(kind=double),dimension(1:totbody)::subsum1
! real(kind=double)::d(:,:)
! !
! 
! k2=0.00029591220828559115_double
! km1=58.132440867048956_double
! km2=3379.3806811609434_double
!  cm2=3.3356611871248456d-5
! 
! !$omp parallel default(shared)   &
! !$omp firstprivate(k2,km1,km2,cm2,bmc,bp1,totm0) & 
! !$omp private(i,j,l)
! !$omp do 
! do i=1,bmc
! frel(i,:)=0._double
! end do
! !$omp end do nowait
! 
! !$omp do
! do i=bp1,totm0
! frel(i,:)=0._double
! end do
! !$omp end do nowait
! 
! 
! !$omp do
! do i=1,bmc
! subsum1(i)=0._double
! subsum2(i,:)=0._double
! subsum3(i,:,:)=0._double
! end do
! !$omp end do nowait
! 
! !$omp do
! do i=bp1,totm0
! subsum1(i)=0._double
! subsum2(i,:)=0._double
! subsum3(i,:,:)=0._double
! end do
! !$omp end do
! 
! 
! !massive pair distances
! !$omp do
! do i=1,bmc-1
!    do j=i+1,bmc
!       do l=1,3
!          rij(i,j,l)=rv(j,l)-rv(i,l)
!       end do
!          rij(j,i,:)=-rij(i,j,:)
! !norm of pair distances
!         d(i,j)=Sqrt(Dot_Product(rij(i,j,:),rij(i,j,:)))
!         d(j,i)=d(i,j)
!     end do
! end do
! !$omp end do nowait
! !$omp do
! do i=1,bmc
!    rij(i,i,:)=0._double
! end do
! !$omp end do nowait
! 
! 
! !massless distances to massive
! !$omp do
! do i=bp1,totm0
!    do j=1,bmc
!       do l=1,3
!          rij(i,j,l)=rv(j,l)-rv(i,l)
!       end do
! !norm of pair distances
!         d(i,j)=Sqrt(Dot_Product(rij(i,j,:),rij(i,j,:)))
!     end do
! end do
! !$omp end do 
! 
! 
! !massive
! !$omp do
!  do i=1,bmc
!       do j=1,bmc!massive on massive
!          if(i.eq.j) then
!             else
!           subsum1(i)=subsum1(i)+mass(j)/d(j,i)
!          end if
!       end do
! 
! !      do j=bp1,totm0 !massless on massive
! !           if (mass(j).eq.0.d0) then
! !           else
! !          subsum1(i)=subsum1(i)+mass(j)/d(j,i)
! !           end if
! !       end do  
!   end do
! !$omp end do 
! 
! 
! !massless
! !$omp do
!  do i=bp1,totm0
!       do j=1,bmc
!         subsum1(i)=subsum1(i)+mass(j)/d(i,j)
!       end do !j
!  end do !i
! !$omp end do
! 
! 
! !$omp do
!  do i=1,bmc
!       do j=1,bmc!massive on massive
!          do l=1,bmc
!          if(i.eq.l.or.j.eq.l) then
!             else
!           subsum2(i,j)=subsum2(i,j)+mass(l)*(-1._double/d(j,l)+0.5_double*Dot_Product(rij(i,j,1:3),rij(j,l,1:3))/d(j,l)**3)
!           subsum3(i,j,1:3)=subsum3(i,j,1:3)+mass(j)*mass(l)/d(i,j)*rij(j,l,1:3)/d(j,l)**3
!          end if
!          end do !l
! 
!           !massless on massive
! !          do l=bp1,totm0
! !           if (mass(l).eq.0.d0) then
! !           else
! !           subsum2(i,j)=subsum2(i,j)+mass(l)*(-1._double/d(l,j)+0.5_double*Dot_Product(rij(i,j,1:3),(-rij(l,j,1:3)))/d(l,j)**3)
! !           subsum3(i,j,1:3)=subsum3(i,j,1:3)+mass(j)*mass(l)/d(i,j)*(-rij(l,j,1:3))/d(l,j)**3
! !           end if
! !         end do !l
!      end do !j   
! end do !i
! !$omp end do 
! 
! 
! !massless
! !$omp do
!  do i=bp1,totm0
!       do j=1,bmc
!          do l=1,bmc
!           if(j.eq.l) then
!             else
!           subsum2(i,j)=subsum2(i,j)+mass(l)*(-1._double/d(j,l)+0.5_double*Dot_Product(rij(i,j,1:3),rij(j,l,1:3))/d(j,l)**3)
!           subsum3(i,j,1:3)=subsum3(i,j,1:3)+mass(j)*mass(l)/d(i,j)*rij(j,l,1:3)/d(j,l)**3
!           end if
!         end do !l
!       end do !j
!  end do!i
! !$omp end do 
! 
! 
!  !relativistic accelerations
! !$omp do
!  do i=1,bmc
!       do j=1,bmc!massive on massive
!          if(i.eq.j) then
!             else
!    frel(i,1:3)=frel(i,1:3)+rij(i,j,1:3)*mass(j)/d(i,j)**3 *&
!              (1._double+cm2*k2*(-4._double*subsum1(i)+subsum2(i,j)-5._double*mass(i)/d(i,j)+ &
!              (Dot_Product(rv(i,4:6),rv(i,4:6))-4._double*Dot_Product(rv(i,4:6),rv(j,4:6))+ & 
!              2._double*Dot_Product(rv(j,4:6),rv(j,4:6)))- &
!              1.5_double *k2*(Dot_Product(rv(j,4:6),rij(i,j,1:3))/d(i,j))**2)) +&
!              cm2*k2*(mass(j)*(Dot_Product(rij(i,j,1:3), &
!              (4._double*rv(i,4:6)-3._double*rv(j,4:6)))/d(i,j)**3*(rv(j,4:6)-rv(i,4:6)))+ &
!              3.5_double *subsum3(i,j,1:3))
!          end if
!         end do !j
!    
! !    do j=bp1,totm0!massless on massive
! !      if (mass(j).eq.0.d0) then
! !      else
! !      frel(i,1:3)=frel(i,1:3)-rij(j,i,1:3)*mass(j)/d(j,i)**3 *&
! !              (1._double+cm2*k2*(-4._double*subsum1(i)-subsum2(j,i)-5._double*mass(i)/d(j,i)+ &  !!!I THINK THE PROBLEM IS HERE!!!!!
! !              (Dot_Product(rv(i,4:6),rv(i,4:6))-4._double*Dot_Product(rv(i,4:6),rv(j,4:6))+ & 
! !              2._double*Dot_Product(rv(j,4:6),rv(j,4:6)))- &
! !              1.5_double *k2*(Dot_Product(rv(j,4:6),(-rij(j,i,1:3)))/d(j,i))**2)) +&
! !              cm2*k2*(mass(j)*(Dot_Product(-rij(j,i,1:3), &
! !              (4._double*rv(i,4:6)-3._double*rv(j,4:6)))/d(j,i)**3*(rv(j,4:6)-rv(i,4:6)))+ &
! !              3.5_double *(-subsum3(j,i,1:3)))
! !      end if
! !    end do !j 
!  end do !i
! !$omp end do
! 
! !massless
! !$omp do
!  do i=bp1,totm0
!       do j=1,bmc
!        frel(i,1:3)=frel(i,1:3)+rij(i,j,1:3)*mass(j)/d(i,j)**3 *&
!              (1._double+cm2*k2*(-4._double*subsum1(i)+subsum2(i,j)-5._double*mass(i)/d(i,j)+ &
!              (Dot_Product(rv(i,4:6),rv(i,4:6))-4._double*Dot_Product(rv(i,4:6),rv(j,4:6))+ & 
!              2._double*Dot_Product(rv(j,4:6),rv(j,4:6)))- &
!              1.5_double *k2*(Dot_Product(rv(j,4:6),rij(i,j,1:3))/d(i,j))**2)) +&
!              cm2*k2*(mass(j)*(Dot_Product(rij(i,j,1:3), &
!              (4._double*rv(i,4:6)-3._double*rv(j,4:6)))/d(i,j)**3*(rv(j,4:6)-rv(i,4:6)))+ &
!              3.5_double  *subsum3(i,j,1:3))
!        end do !j
!  end do !i
! !$omp end do
! !$omp end parallel
! 
! 
! return
! end subroutine frho_eih

!######################################################################################
!******************************************************************

subroutine Int_Radau_yeih(p,rv,mass,add)
!Calculates relativistic forces from Einstein Infeld Hoffman PPN N-Body metric including maximum Yarkovsky effect
implicit none
integer(kind=intdble)::allocstat
integer(kind=intdble)::cnt,m0body,body,dim,bodyp1
integer(kind=intdble)::cutbody,m0cutbody,totbody,mbody


type(Problem),intent(inout)::p 
type(additional),intent(in)::add

integer(kind=intdble)::order(1:p%NrKoerper)
real(kind=double),dimension(1:p%NrKoerper),intent(inout)::mass
real(kind=double),dimension(1:p%NrKoerper,1:6)::rv
real(kind=double)::tmod,dtmin
real(kind=double)::dtk,t,dtkmax,dtkmin,dtdid
real(kind=double),allocatable,dimension(:,:,:)::c_b,c_e
real(kind=double),allocatable,dimension(:,:)::rho
real(kind=double),allocatable,dimension(:)::c_h,c_vc,c_xc,c_c,c_d,c_r

real(kind=double)::kahan_y,kahan_t,kahan_c

logical::mflag,m0flag
logical,dimension(:,:),allocatable::merger

integer(kind=intdble)::i,n,j,stiffcnt



t=0._double
kahan_y=0._double
kahan_t=0._double
kahan_c=0._double

 cnt=1
 stiffcnt=0
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

!RA15 constants
allocate(c_b(7,totbody,3),c_e(7,totbody,3),c_h(8),c_xc(8),c_vc(7),c_c(21),c_d(21),c_r(28),rho(totbody,body), & 
         merger(totbody,body),stat=allocstat)

  if (allocstat.ne.0) then
            write(*,*)'error in radau_m: allocation not possible, try reducing the number of particles'
  end if 

merger(:,:)=.false.

! Gauss-Radau spacings for substeps within a sequence, for the 15th order 
! integrator. The sum of the H values should be 3.733333333333333
!
  c_h(:)=(/          0._double,.0562625605269221464656522_double,.1802406917368923649875799_double,  &
       .3526247171131696373739078_double,.5471536263305553830014486_double,.7342101772154105315232106_double,  &
       .8853209468390957680903598_double,.9775206135612875018911745_double/)
!
! Constant coefficients used in series expansions for X and V
!  XC: 1/2,  1/6,  1/12, 1/20, 1/30, 1/42, 1/56, 1/72
!  VC: 1/2,  1/3,  1/4,  1/5,  1/6,  1/7,  1/8
 
!   c_xc(:)=(/.5_double,.1666666666666667_double,.08333333333333333_double,.05_double,  &
!        .03333333333333333_double,.02380952380952381_double,.01785714285714286_double, &
!        .01388888888888889_double/)
! 
!   c_vc(:)=(/.5d0,.3333333333333333_double,.25_double,.2_double, &
!         .1666666666666667_double,.1428571428571429_double,.125_double/)
  c_xc(1)=0.5_double
  c_xc(2)=1._double/6._double
  c_xc(3)=1._double/12._double
  c_xc(4)=0.05_double
  c_xc(5)=1._double/30._double
  c_xc(6)=1._double/42._double
  c_xc(7)=1._double/56._double
  c_xc(8)=1._double/72._double

 c_vc(1)=0.5_double
 c_vc(2)=1._double/3._double
 c_vc(3)=0.25_double
 c_vc(4)=0.2_double
 c_vc(5)=1._double/6._double
 c_vc(6)=1._double/7._double
 c_vc(7)=0.125_double

  c_b(:,:,:)=0._double
  c_e(:,:,:)=0._double
!
! set values of the constant arrays
! (R = R21, R31, R32, R41, R42, R43 in Everhart's paper.)
    
        n = 0
        do j = 2, 8
          do i = 1, j - 1
            n = n + 1
            c_r(n) = 1._double / (c_h(j) - c_h(i))
          end do
        end do
!
! Constants to convert between B and G arrays (C = C21, C31, C32, C41, C42...)
        c_c(1) = - c_h(2)
        c_d(1) =   c_h(2)
        n = 1
        do j = 3, 7
          n = n + 1
          c_c(n) = -c_h(j) * c_c(n-j+2)
          c_d(n) =  c_h(2) * c_d(n-j+2)
          do i = 3, j - 1
            n = n + 1
            c_c(n) = c_c(n-j+1)  -  c_h(j) * c_c(n-j+2)
            c_d(n) = c_d(n-j+1)  +  c_h(i) * c_d(n-j+2)
          end do
          n = n + 1
          c_c(n) = c_c(n-j+1) - c_h(j)
          c_d(n) = c_d(n-j+1) + c_h(j)
        end do



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
! 
call radau_yeih_step(rv,mass,totbody,body,cutbody,m0cutbody,dtk,dtdid,p%eps,rho,c_b, &
                c_e,c_h,c_xc,c_vc,c_c,c_d,c_r)
!
   kahan_y=dtdid/k-kahan_c
   kahan_t=t+kahan_y
   kahan_c=(kahan_t-t)-kahan_y
   t=kahan_t

   if(dtdid.gt.dtkmax) then
      dtkmax=dtdid
    end if
    if (dtdid.le.dtkmin) then
       dtkmin=dtdid
       stiffcnt=stiffcnt+1
      else 
       stiffcnt=0 
   end if

   if (stiffcnt.gt.10000) then
     write(*,*)'possibly stiff region encountered, stopping'
     STOP
   end if 


 

   !Ausgabe 
   if (mod(real(t/p%tAusg,KIND=double),1._double) .lt.1.d-15.and.t.ge.p%tAusg  .or. dble(p%tAusg).eq.0.d0)then
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

 cnt=cnt+1
end do prgs

write(*,*)'time accomplished:',t,'with an average stepsize of:',t/dble(cnt)
write(*,*)'minimum stepsize',dtkmin/k
write(*,*)'maximum stepsize',dtkmax/k


write(*,*)'number of massive particles at the beginning / end of integration:',body+cutbody,' / ',body
write(*,*)'number of massless particles at the beginning / end of integration:',m0body+m0cutbody,' / ',m0body
if(p%outegy) then
  
   write(*,*)'error in total energy (sum of dE)',dEtot
   write(*,*)'error in total angular momentum (sum of dL)',dLtot
end if

close(41)


return
end subroutine Int_Radau_yeih
!******************************************************************************
!
! Author: John E. Chambers
! Modified: Siegfried Eggl  201104
!
! Integrates bodies for one timestep TDID using
! Everhart's RA15 integrator algorithm. The accelerations are calculated
! using the subroutine FROH. The accuracy of the step is approximately 
! determined by the tolerance parameter TOL.
!
! Based on RADAU by E. Everhart, Physics Department, University of Denver.
! Comments giving equation numbers refer to Everhart (1985) ``An Efficient
! Integrator that Uses Gauss-Radau Spacings'', in The Dynamics of Comets:
! Their Origin and Evolution, p185-202, eds. A. Carusi & G. B. Valsecchi,
! pub Reidel. (A listing of the original subroutine is also given in this 
! paper.)
! 
!------------------------------------------------------------------------------
      subroutine radau_yeih_step(rv,mass,totbody,body,cutbody,m0cutbody,& 
                            t,tdid,tol,rho,b,e,h,xc,vc,c,d,r)

      implicit none
!
! Input/Output
      integer(kind=intdble):: totbody,body,cutbody,m0cutbody 
      real(kind=double):: t,tdid,tol
      real(kind=double):: rv(:,:),rho(:,:),mass(:)  
!
! Local
      integer(kind=intdble):: j,k,n,i
      integer(kind=intdble):: bmc,bp1,totm0 
      real(kind=double)::g(7,totbody,3),b(7,totbody,3),e(7,totbody,3)
      real(kind=double)::h(8),xc(8),vc(7),c(21),d(21),r(28),s(9)
      real(kind=double)::q,q2,q3,q4,q5,q6,q7,temp,gk,gkmax
      real(kind=double),dimension(totbody,3)::acc,acc1
      real(kind=double)::rv1(totbody,6),temp3(3)
      logical::iter
!
bmc=body-cutbody
bp1=body+1
totm0=totbody-m0cutbody
!
! Calculate forces at the start of the sequence
      call frho_yeih(rv,totbody,bmc,bp1,totm0,mass,acc,rho) 
!
! Find G values from B values predicted at the last call (Eqs. 7 of Everhart)
!
!massive

!$omp parallel default(NONE) &
!$omp shared(g,bmc,bp1,totm0,b,d) &
!$omp private(k)
!$omp do
      do  k= 1, bmc
        g(1,k,:) = b(7,k,:)*d(16) + b(6,k,:)*d(11) + b(5,k,:)*d(7) &
                 + b(4,k,:)*d(4)  + b(3,k,:)*d(2)  + b(2,k,:)*d(1)  + b(1,k,:)
        g(2,k,:) = b(7,k,:)*d(17) + b(6,k,:)*d(12) + b(5,k,:)*d(8) &
                 + b(4,k,:)*d(5)  + b(3,k,:)*d(3)  + b(2,k,:)
        g(3,k,:) = b(7,k,:)*d(18) + b(6,k,:)*d(13) + b(5,k,:)*d(9)  &
                 + b(4,k,:)*d(6)  + b(3,k,:)
        g(4,k,:) = b(7,k,:)*d(19) + b(6,k,:)*d(14) + b(5,k,:)*d(10) + b(4,k,:)
        g(5,k,:) = b(7,k,:)*d(20) + b(6,k,:)*d(15) + b(5,k,:)
        g(6,k,:) = b(7,k,:)*d(21) + b(6,k,:)
        g(7,k,:) = b(7,k,:)
     end do
!$omp end do nowait

!massless
!$omp do 
     do  k=bp1,totm0
        g(1,k,:) = b(7,k,:)*d(16) + b(6,k,:)*d(11) + b(5,k,:)*d(7) &
                 + b(4,k,:)*d(4)  + b(3,k,:)*d(2)  + b(2,k,:)*d(1)  + b(1,k,:)
        g(2,k,:) = b(7,k,:)*d(17) + b(6,k,:)*d(12) + b(5,k,:)*d(8) &
                 + b(4,k,:)*d(5)  + b(3,k,:)*d(3)  + b(2,k,:)
        g(3,k,:) = b(7,k,:)*d(18) + b(6,k,:)*d(13) + b(5,k,:)*d(9)  &
                 + b(4,k,:)*d(6)  + b(3,k,:)
        g(4,k,:) = b(7,k,:)*d(19) + b(6,k,:)*d(14) + b(5,k,:)*d(10) + b(4,k,:)
        g(5,k,:) = b(7,k,:)*d(20) + b(6,k,:)*d(15) + b(5,k,:)
        g(6,k,:) = b(7,k,:)*d(21) + b(6,k,:)
        g(7,k,:) = b(7,k,:)
     end do
!$omp end do
!$omp end parallel
!------------------------------------------------------------------------------
!
!  MAIN  LOOP  STARTS  HERE
!
! For each iteration (six for first call to subroutine, two otherwise)...
 
!   do n = 1, niter
    n=0
    iter=.true.
    do while(iter)
!    do while(n<6)
!    
    n=n+1    
! For each substep within a sequence...
        do j = 2, 8
!
!$omp parallel default(NONE) &
!$omp shared(h,rv1,rv,b,acc,bmc,bp1,totm0,t,tol,iter,j) &
!$omp private(s,k)
! Calculate position predictors using Eqn. 9 of Everhart
          s(1) = t * h(j)
          s(2) = s(1) * s(1) * .5_double
          s(3) = s(2) * h(j) *  1._double/3._double
          s(4) = s(3) * h(j) * .5_double
          s(5) = s(4) * h(j) * .6_double
          s(6) = s(5) * h(j) *  2._double/3._double
          s(7) = s(6) * h(j) * 5._double/7._double
          s(8) = s(7) * h(j) * .75_double
          s(9) = s(8) * h(j) *  7._double/9._double
!
!massive
!$omp do
          do k = 1, bmc           
            rv1(k,1:3) = s(9)*b(7,k,:) + s(8)*b(6,k,:) + s(7)*b(5,k,:) &
                       + s(6)*b(4,k,:) + s(5)*b(3,k,:) + s(4)*b(2,k,:) &
                       + s(3)*b(1,k,:) + s(2)*acc(k,:)  + s(1)*rv(k,4:6) + rv(k,1:3)
          end do
!$omp end do
!massless
!$omp do
         do k = bp1,totm0           
            rv1(k,1:3) = s(9)*b(7,k,:) + s(8)*b(6,k,:) + s(7)*b(5,k,:) &
                       + s(6)*b(4,k,:) + s(5)*b(3,k,:) + s(4)*b(2,k,:) &
                       + s(3)*b(1,k,:) + s(2)*acc(k,:)  + s(1)*rv(k,4:6) + rv(k,1:3)
          end do
!$omp end do
!
! If necessary, calculate velocity predictors too, from Eqn. 10 of Everhart
            s(1) = t * h(j)
            s(2) = s(1) * h(j) * .5_double
            s(3) = s(2) * h(j) * 2._double/3._double
            s(4) = s(3) * h(j) * .75_double
            s(5) = s(4) * h(j) * .8_double
            s(6) = s(5) * h(j) * 5._double/6._double
            s(7) = s(6) * h(j) * 6._double/7._double
            s(8) = s(7) * h(j) * .875_double
!massive
!$omp do
            do k = 1, bmc
              rv1(k,4:6) = s(8)*b(7,k,:) + s(7)*b(6,k,:) + s(6)*b(5,k,:) &
                         + s(5)*b(4,k,:) + s(4)*b(3,k,:) + s(3)*b(2,k,:)    &
                         + s(2)*b(1,k,:) + s(1)*acc(k,:)  + rv(k,4:6)
            end do
!$omp end do nowait
!massless
!$omp do
            do k = bp1,totm0
              rv1(k,4:6) = s(8)*b(7,k,:) + s(7)*b(6,k,:) + s(6)*b(5,k,:) &
                         + s(5)*b(4,k,:) + s(4)*b(3,k,:) + s(3)*b(2,k,:)    &
                         + s(2)*b(1,k,:) + s(1)*acc(k,:)  + rv(k,4:6)
            end do 
!$omp end do
!$omp end parallel

!
! Calculate forces at the current substep
          call frho_yeih(rv1,totbody,bmc,bp1,totm0,mass,acc1,rho) 
! Update G values using Eqs. 4 of Everhart, and update B values using Eqs. 5
 
  select case(j)
     case(2)
!$omp parallel default(NONE) &
!$omp shared(b,r,g,acc,acc1,bmc,bp1,totm0,j) &
!$omp private(k,temp3)
!massive
!$omp do
            do k = 1, bmc
              temp3 = g(1,k,:)
              g(1,k,:) = (acc1(k,:) - acc(k,:)) * r(1)
              b(1,k,:) = b(1,k,:) + g(1,k,:) - temp3
            end do
!$omp end do nowait
!massless
!$omp do
            do k = bp1,totm0 
              temp3 = g(1,k,:)
              g(1,k,:) = (acc1(k,:) - acc(k,:)) * r(1)
              b(1,k,:) = b(1,k,:) + g(1,k,:) - temp3
            end do
!$omp end do
!$omp end parallel

      case(3)
!$omp parallel default(NONE) &
!$omp shared(b,c,r,g,acc,acc1,bmc,bp1,totm0,j) &
!$omp private(i,k,temp3,gk,temp)
!massive
!$omp do
            do k = 1,bmc
              do i=1,3             
               temp = g(2,k,i)
               gk = acc1(k,i) - acc(k,i)
               g(2,k,i) = (gk*r(2) - g(1,k,i))*r(3)
               temp = g(2,k,i) - temp
               b(1,k,i) = b(1,k,i)  +  temp * c(1)
               b(2,k,i) = b(2,k,i)  +  temp
              end do
            end do
!$omp end do nowait
!massless
!$omp do
            do k = bp1,totm0
              do i=1,3             
               temp = g(2,k,i)
               gk = acc1(k,i) - acc(k,i)
               g(2,k,i) = (gk*r(2) - g(1,k,i))*r(3)
               temp = g(2,k,i) - temp
               b(1,k,i) = b(1,k,i)  +  temp * c(1)
               b(2,k,i) = b(2,k,i)  +  temp
              end do
            end do
!$omp end do
!$omp end parallel
      case(4)
!$omp parallel default(NONE) &
!$omp shared(b,r,c,g,acc,acc1,bmc,bp1,totm0,j) &
!$omp private(i,k,temp3,gk,temp)
!massive
!$omp do
            do k = 1,bmc
             do i=1,3
              temp = g(3,k,i)
              gk = acc1(k,i) - acc(k,i)
              g(3,k,i) = ((gk*r(4) - g(1,k,i))*r(5) - g(2,k,i))*r(6)
              temp = g(3,k,i) - temp
              b(1,k,i) = b(1,k,i)  +  temp * c(2)
              b(2,k,i) = b(2,k,i)  +  temp * c(3)
              b(3,k,i) = b(3,k,i)  +  temp
             end do
            end do
!$omp end do nowait
!massless
!$omp do
           do k = bp1,totm0
             do i=1,3
              temp = g(3,k,i)
              gk = acc1(k,i) - acc(k,i)
              g(3,k,i) = ((gk*r(4) - g(1,k,i))*r(5) - g(2,k,i))*r(6)
              temp = g(3,k,i) - temp
              b(1,k,i) = b(1,k,i)  +  temp * c(2)
              b(2,k,i) = b(2,k,i)  +  temp * c(3)
              b(3,k,i) = b(3,k,i)  +  temp
             end do
            end do
!$omp end do
!$omp end parallel
       case(5)
!$omp parallel default(NONE) &
!$omp shared(b,c,r,g,acc,acc1,bmc,bp1,totm0,j) &
!$omp private(i,k,temp3,gk,temp)
!massive
!$omp do
            do k = 1, bmc
             do i=1,3
              temp = g(4,k,i)
              gk = acc1(k,i) - acc(k,i)
              g(4,k,i) = (((gk*r(7) - g(1,k,i))*r(8) - g(2,k,i))*r(9) &
                    - g(3,k,i))*r(10)
              temp = g(4,k,i) - temp
              b(1,k,i) = b(1,k,i)  +  temp * c(4)
              b(2,k,i) = b(2,k,i)  +  temp * c(5)
              b(3,k,i) = b(3,k,i)  +  temp * c(6)
              b(4,k,i) = b(4,k,i)  +  temp
             end do
            end do
!$omp end do nowait
!massless
!$omp do
            do k = bp1,totm0
             do i=1,3
              temp = g(4,k,i)
              gk = acc1(k,i) - acc(k,i)
              g(4,k,i) = (((gk*r(7) - g(1,k,i))*r(8) - g(2,k,i))*r(9) &
                    - g(3,k,i))*r(10)
              temp = g(4,k,i) - temp
              b(1,k,i) = b(1,k,i)  +  temp * c(4)
              b(2,k,i) = b(2,k,i)  +  temp * c(5)
              b(3,k,i) = b(3,k,i)  +  temp * c(6)
              b(4,k,i) = b(4,k,i)  +  temp
             end do
            end do
!$omp end do
!$omp end parallel
        case(6)
!$omp parallel default(NONE) &
!$omp shared(b,c,r,g,acc,acc1,bmc,bp1,totm0,j) &
!$omp private(i,k,gk,temp3,temp)
!massive
!$omp do
            do k = 1, bmc
             do i=1,3
              temp = g(5,k,i)
              gk = acc1(k,i) - acc(k,i)
              g(5,k,i) =  ((((gk*r(11) - g(1,k,i))*r(12) - g(2,k,i))*r(13) &
                    - g(3,k,i))*r(14) - g(4,k,i))*r(15)
              temp = g(5,k,i) - temp
              b(1,k,i) = b(1,k,i)  +  temp * c(7)
              b(2,k,i) = b(2,k,i)  +  temp * c(8)
              b(3,k,i) = b(3,k,i)  +  temp * c(9)
              b(4,k,i) = b(4,k,i)  +  temp * c(10)
              b(5,k,i) = b(5,k,i)  +  temp
             end do
            end do
!$omp end do nowait
!massless
!$omp do
           do k = bp1,totm0
             do i=1,3
              temp = g(5,k,i)
              gk = acc1(k,i) - acc(k,i)
              g(5,k,i) =  ((((gk*r(11) - g(1,k,i))*r(12) - g(2,k,i))*r(13) &
                    - g(3,k,i))*r(14) - g(4,k,i))*r(15)
              temp = g(5,k,i) - temp
              b(1,k,i) = b(1,k,i)  +  temp * c(7)
              b(2,k,i) = b(2,k,i)  +  temp * c(8)
              b(3,k,i) = b(3,k,i)  +  temp * c(9)
              b(4,k,i) = b(4,k,i)  +  temp * c(10)
              b(5,k,i) = b(5,k,i)  +  temp
             end do
           end do
!$omp end do
!$omp end parallel
         case(7)
!$omp parallel default(NONE) &
!$omp shared(b,c,r,g,acc,acc1,bmc,bp1,totm0,j) &
!$omp private(i,k,temp3,gk,temp)
!massive
!$omp do
            do k = 1, bmc
             do i=1,3
              temp = g(6,k,i)
              gk = acc1(k,i) - acc(k,i)
              g(6,k,i) = (((((gk*r(16) - g(1,k,i))*r(17) - g(2,k,i))*r(18) &
                     - g(3,k,i))*r(19) - g(4,k,i))*r(20) - g(5,k,i))*r(21)
              temp = g(6,k,i) - temp
              b(1,k,i) = b(1,k,i)  +  temp * c(11)
              b(2,k,i) = b(2,k,i)  +  temp * c(12)
              b(3,k,i) = b(3,k,i)  +  temp * c(13)
              b(4,k,i) = b(4,k,i)  +  temp * c(14)
              b(5,k,i) = b(5,k,i)  +  temp * c(15)
              b(6,k,i) = b(6,k,i)  +  temp
            end do
          end do
!$omp end do nowait
!massless
!$omp do
          do k = bp1,totm0
             do i=1,3
              temp = g(6,k,i)
              gk = acc1(k,i) - acc(k,i)
              g(6,k,i) = (((((gk*r(16) - g(1,k,i))*r(17) - g(2,k,i))*r(18) &
                     - g(3,k,i))*r(19) - g(4,k,i))*r(20) - g(5,k,i))*r(21)
              temp = g(6,k,i) - temp
              b(1,k,i) = b(1,k,i)  +  temp * c(11)
              b(2,k,i) = b(2,k,i)  +  temp * c(12)
              b(3,k,i) = b(3,k,i)  +  temp * c(13)
              b(4,k,i) = b(4,k,i)  +  temp * c(14)
              b(5,k,i) = b(5,k,i)  +  temp * c(15)
              b(6,k,i) = b(6,k,i)  +  temp
            end do
          end do
!$omp end do
!$omp end parallel

         case(8)
           gkmax=0._double
!massive
            do k = 1, bmc
             do i=1,3
              temp = g(7,k,i)
              gk = acc1(k,i) - acc(k,i)
              g(7,k,i) = ((((((gk*r(22) - g(1,k,i))*r(23) - g(2,k,i))*r(24) &
                     - g(3,k,i))*r(25) - g(4,k,i))*r(26) - g(5,k,i))*r(27) &
                     - g(6,k,i))*r(28)
              temp = g(7,k,i) - temp
              b(1,k,i) = b(1,k,i)  +  temp * c(16)
              b(2,k,i) = b(2,k,i)  +  temp * c(17)
              b(3,k,i) = b(3,k,i)  +  temp * c(18)
              b(4,k,i) = b(4,k,i)  +  temp * c(19)
              b(5,k,i) = b(5,k,i)  +  temp * c(20)
              b(6,k,i) = b(6,k,i)  +  temp * c(21)
              b(7,k,i) = b(7,k,i)  +  temp
              gkmax=max(gkmax,abs(temp))             
             end do !i 1..3
            end do  !k 1..nbody
!massless
             do k = bp1,totm0
             do i=1,3
              temp = g(7,k,i)
              gk = acc1(k,i) - acc(k,i)
              g(7,k,i) = ((((((gk*r(22) - g(1,k,i))*r(23) - g(2,k,i))*r(24) &
                     - g(3,k,i))*r(25) - g(4,k,i))*r(26) - g(5,k,i))*r(27) &
                     - g(6,k,i))*r(28)
              temp = g(7,k,i) - temp
              b(1,k,i) = b(1,k,i)  +  temp * c(16)
              b(2,k,i) = b(2,k,i)  +  temp * c(17)
              b(3,k,i) = b(3,k,i)  +  temp * c(18)
              b(4,k,i) = b(4,k,i)  +  temp * c(19)
              b(5,k,i) = b(5,k,i)  +  temp * c(20)
              b(6,k,i) = b(6,k,i)  +  temp * c(21)
              b(7,k,i) = b(7,k,i)  +  temp
              gkmax=max(gkmax,abs(temp))             
             end do !i 1..3
            end do  !k 1..nbody

!more iterations necessary?
                if(gkmax.gt.tol) then
                  iter=.true.
               else
                  iter=.false.
               end if
!              write(*,*)gkmax,n
        end select
      end do      !j 2...8
!     rerun subroutine with smaller stepsize if things do not converge      
      if(n>6) then
         t=t/3.13_double
         tdid=0._double
!         write(*,*)'rvint',rv(2,:)
!         write(*,*)'no convergence, trying smaller initial stepsize!',t*3.13d0,'->',t        
         return
      end if
     end do      !n... while(iter)
!
!------------------------------------------------------------------------------
!
!  END  OF  MAIN  LOOP
!


! Estimate suitable sequence size for the next call to subroutine (Eqs. 15, 16)
      temp = 0._double

!massive
      do k = 1, bmc
       do i=1,3
        temp = max( temp, abs( b(7,k,i) ) )
       end do
      end do
!massless
      do k = bp1,totm0
       do i=1,3
        temp = max( temp, abs( b(7,k,i) ) )
       end do
      end do
!
      temp = temp / (72._double * abs(t)**7)
      tdid = t
      if (temp.eq.0) then
        t = tdid * 1.4_double
      else
        t = sign( (tol/temp)**(1._double/9._double), tdid )
      end if
!
! If new sequence size is much bigger than the current one, reduce it
      if (abs(t/tdid).gt.1.4_double) t = tdid * 1.4_double
!
! Find new position and velocity values at end of the sequence (Eqs. 11, 12)

      temp = tdid * tdid

!$omp parallel default(NONE) &
!$omp shared(rv,xc,vc,b,acc,bmc,bp1,totm0,t,tdid,temp) &
!$omp private(k)
!massive
!$omp do
      do k = 1 , bmc
        rv(k,1:3) = (xc(8)*b(7,k,:) + xc(7)*b(6,k,:) + xc(6)*b(5,k,:) &
                  +  xc(5)*b(4,k,:) + xc(4)*b(3,k,:) + xc(3)*b(2,k,:) &
                  +  xc(2)*b(1,k,:) + xc(1)*acc(k,:))*temp + rv(k,4:6)*tdid + rv(k,1:3)
!
        rv(k,4:6) = (vc(7)*b(7,k,:) + vc(6)*b(6,k,:) + vc(5)*b(5,k,:) &
                  +  vc(4)*b(4,k,:) + vc(3)*b(3,k,:) + vc(2)*b(2,k,:)     &
                  +  vc(1)*b(1,k,:) + acc(k,:))*tdid + rv(k,4:6)
      end do
!$omp end do nowait
!massless
!$omp do
      do k = bp1,totm0
        rv(k,1:3) = (xc(8)*b(7,k,:) + xc(7)*b(6,k,:) + xc(6)*b(5,k,:) &
                  +  xc(5)*b(4,k,:) + xc(4)*b(3,k,:) + xc(3)*b(2,k,:) &
                  +  xc(2)*b(1,k,:) + xc(1)*acc(k,:))*temp + rv(k,4:6)*tdid + rv(k,1:3)
!
        rv(k,4:6) = (vc(7)*b(7,k,:) + vc(6)*b(6,k,:) + vc(5)*b(5,k,:) &
                  +  vc(4)*b(4,k,:) + vc(3)*b(3,k,:) + vc(2)*b(2,k,:)     &
                  +  vc(1)*b(1,k,:) + acc(k,:))*tdid + rv(k,4:6)
      end do
!$omp end do
!$omp end parallel

! Predict new B values to use at the start of the next sequence. The predicted
! values from the last call are saved as E. The correction, BD, between the
! actual and predicted values of B is applied in advance as a correction.
      q = t / tdid
      q2 = q  * q
      q3 = q  * q2
      q4 = q2 * q2
      q5 = q2 * q3
      q6 = q3 * q3
      q7 = q3 * q4
!
!massive
      do k = 1, bmc
       do i =1,3
        s(1) = b(1,k,i) - e(1,k,i)
        s(2) = b(2,k,i) - e(2,k,i)
        s(3) = b(3,k,i) - e(3,k,i)
        s(4) = b(4,k,i) - e(4,k,i)
        s(5) = b(5,k,i) - e(5,k,i)
        s(6) = b(6,k,i) - e(6,k,i)
        s(7) = b(7,k,i) - e(7,k,i)
!
! Estimate B values for the next sequence (Eqs. 13 of Everhart).
        e(1,k,i) = q* (b(7,k,i)* 7._double + b(6,k,i)* 6._double + b(5,k,i)* 5._double &
               +       b(4,k,i)* 4._double + b(3,k,i)* 3._double + b(2,k,i)*2._double + b(1,k,i))
        e(2,k,i) = q2*(b(7,k,i)*21._double + b(6,k,i)*15._double + b(5,k,i)*10._double &
               +       b(4,k,i)* 6._double + b(3,k,i)* 3._double + b(2,k,i))      
        e(3,k,i) = q3*(b(7,k,i)*35._double + b(6,k,i)*20._double + b(5,k,i)*10._double  &
               +       b(4,k,i)*4._double  + b(3,k,i))
        e(4,k,i) = q4*(b(7,k,i)*35._double + b(6,k,i)*15._double + b(5,k,i)*5._double + b(4,k,i))
        e(5,k,i) = q5*(b(7,k,i)*21._double + b(6,k,i)*6._double  + b(5,k,i))
        e(6,k,i) = q6*(b(7,k,i)*7._double  + b(6,k,i))
        e(7,k,i) = q7* b(7,k,i)
!
        b(1,k,i) = e(1,k,i) + s(1)
        b(2,k,i) = e(2,k,i) + s(2)
        b(3,k,i) = e(3,k,i) + s(3)
        b(4,k,i) = e(4,k,i) + s(4)
        b(5,k,i) = e(5,k,i) + s(5)
        b(6,k,i) = e(6,k,i) + s(6)
        b(7,k,i) = e(7,k,i) + s(7)
       end do
      end do
!massless
     do k = bp1,totm0
       do i =1,3
        s(1) = b(1,k,i) - e(1,k,i)
        s(2) = b(2,k,i) - e(2,k,i)
        s(3) = b(3,k,i) - e(3,k,i)
        s(4) = b(4,k,i) - e(4,k,i)
        s(5) = b(5,k,i) - e(5,k,i)
        s(6) = b(6,k,i) - e(6,k,i)
        s(7) = b(7,k,i) - e(7,k,i)
!
! Estimate B values for the next sequence (Eqs. 13 of Everhart).
        e(1,k,i) = q* (b(7,k,i)* 7._double + b(6,k,i)* 6._double + b(5,k,i)* 5._double &
               +       b(4,k,i)* 4._double + b(3,k,i)* 3._double + b(2,k,i)*2._double + b(1,k,i))
        e(2,k,i) = q2*(b(7,k,i)*21._double + b(6,k,i)*15._double + b(5,k,i)*10._double &
               +       b(4,k,i)* 6._double + b(3,k,i)* 3._double + b(2,k,i))      
        e(3,k,i) = q3*(b(7,k,i)*35._double + b(6,k,i)*20._double + b(5,k,i)*10._double  &
               +       b(4,k,i)*4._double  + b(3,k,i))
        e(4,k,i) = q4*(b(7,k,i)*35._double + b(6,k,i)*15._double + b(5,k,i)*5._double + b(4,k,i))
        e(5,k,i) = q5*(b(7,k,i)*21._double + b(6,k,i)*6._double  + b(5,k,i))
        e(6,k,i) = q6*(b(7,k,i)*7._double  + b(6,k,i))
        e(7,k,i) = q7* b(7,k,i)
!
        b(1,k,i) = e(1,k,i) + s(1)
        b(2,k,i) = e(2,k,i) + s(2)
        b(3,k,i) = e(3,k,i) + s(3)
        b(4,k,i) = e(4,k,i) + s(4)
        b(5,k,i) = e(5,k,i) + s(5)
        b(6,k,i) = e(6,k,i) + s(6)
        b(7,k,i) = e(7,k,i) + s(7)
       end do
      end do
      return
      end subroutine 

!******************************************************************************************
subroutine frho_yeih(rv,totbody,bmc,bp1,totm0,mass,frel,d) 
use global_m
use yark_m
use transform_m
implicit none
integer(kind=intdble)::totbody,bmc,bp1,totm0,i,j,l
real(kind=double)::rv(1:totbody,1:6),mass(1:totbody)
real(kind=double)::frel(1:totbody,1:3)
real(kind=double),dimension(1:totbody,1:bmc,1:3)::rij,subsum3
real(kind=double),dimension(1:totbody,1:bmc)::subsum2
real(kind=double),dimension(1:totbody)::subsum1
real(kind=double)::d(:,:),dm(1:totbody)
real(kind=double),dimension(1:3)::yadiu !,yasea
real(kind=double),dimension(1:totbody,1:6):: drv!, ele
!real(kind=double),dimension(1:3)::ybxrs,ybxrsxyb,rs
real(kind=double),parameter::k2=0.00029591220828559115_double !gaussian gravitational constant squared (k^2)
real(kind=double),parameter::km1=58.132440867048956_double    !1/k
real(kind=double),parameter::km2=3379.3806811609434_double    !1/k^2
real(kind=double),parameter::cm2=3.3356611871248456d-5        !1/c^2 vacuum light speed
real(kind=double),parameter::cm2k2=9.870628679746496d-9	      !k^2/c^2



!$omp parallel default(shared)   &
!$omp firstprivate(k2,km1,km2,cm2,bmc,bp1,totm0) & 
!$omp private(i,j,l)
!$omp do 
do i=1,bmc
frel(i,:)=0._double
end do
!$omp end do nowait

!$omp do
do i=bp1,totm0
frel(i,:)=0._double
end do
!$omp end do nowait


!$omp do
do i=1,bmc
subsum1(i)=0._double
subsum2(i,:)=0._double
subsum3(i,:,:)=0._double
end do
!$omp end do nowait

!$omp do
do i=bp1,totm0
subsum1(i)=0._double
subsum2(i,:)=0._double
subsum3(i,:,:)=0._double
end do
!$omp end do


!massive pair distances
!$omp do
do i=1,bmc-1
   do j=i+1,bmc
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
do i=1,bmc
   rij(i,i,:)=0._double
end do
!$omp end do nowait


!massless distances to massive
!$omp do
do i=bp1,totm0
   do j=1,bmc
      do l=1,3
         rij(i,j,l)=rv(j,l)-rv(i,l)
      end do
!norm of pair distances
        d(i,j)=Sqrt(Dot_Product(rij(i,j,:),rij(i,j,:)))
    end do
end do
!$omp end do 


!massive
!$omp do
 do i=1,bmc
      do j=1,bmc!massive on massive
         if(i.eq.j) then
            else
          subsum1(i)=subsum1(i)+mass(j)/d(j,i)
         end if
      end do

!      do j=bp1,totm0 !massless on massive
!          subsum1(i)=subsum1(i)+mass(j)/d(j,i)
!       end do  
  end do
!$omp end do 


!massless
!$omp do
 do i=bp1,totm0
      do j=1,bmc
        subsum1(i)=subsum1(i)+mass(j)/d(i,j)
      end do !j
 end do !i
!$omp end do


!$omp do
 do i=1,bmc
      do j=1,bmc!massive on massive
         do l=1,bmc
         if(i.eq.l.or.j.eq.l) then
            else
          subsum2(i,j)=subsum2(i,j)+mass(l)*(-1._double/d(j,l)+0.5_double*Dot_Product(rij(i,j,1:3),rij(j,l,1:3))/d(j,l)**3)
          subsum3(i,j,1:3)=subsum3(i,j,1:3)+mass(j)*mass(l)/d(i,j)*rij(j,l,1:3)/d(j,l)**3
         end if
         end do !l

!           !massless on massive
!          do l=bp1,totm0
!           subsum2(i,j)=subsum2(i,j)+mass(l)*(-1._double/d(l,j)+0.5_double*Dot_Product(rij(i,j,1:3),(-rij(l,j,1:3)))/d(l,j)**3)
!           subsum3(i,j,1:3)=subsum3(i,j,1:3)+mass(j)*mass(l)/d(i,j)*(-rij(l,j,1:3))/d(l,j)**3
!         end do !l
     end do !j   
end do !i
!$omp end do 


!massless
!$omp do
 do i=bp1,totm0
      do j=1,bmc
         do l=1,bmc
          if(j.eq.l) then
            else
          subsum2(i,j)=subsum2(i,j)+mass(l)*(-1._double/d(j,l)+0.5_double*Dot_Product(rij(i,j,1:3),rij(j,l,1:3))/d(j,l)**3)
          subsum3(i,j,1:3)=subsum3(i,j,1:3)+mass(j)*mass(l)/d(i,j)*rij(j,l,1:3)/d(j,l)**3
          end if
        end do !l
      end do !j
 end do!i
!$omp end do 


 !relativistic accelerations
!$omp do
 do i=1,bmc
      do j=1,bmc!massive on massive
         if(i.eq.j) then
            else
   frel(i,1:3)=frel(i,1:3)+rij(i,j,1:3)*mass(j)/d(i,j)**3 *&
             (1._double+cm2k2*(-4._double*subsum1(i)+subsum2(i,j)-5._double*mass(i)/d(i,j)+ &
             (Dot_Product(rv(i,4:6),rv(i,4:6))-4._double*Dot_Product(rv(i,4:6),rv(j,4:6))+ & 
             2._double*Dot_Product(rv(j,4:6),rv(j,4:6)))- &
             1.5_double *(Dot_Product(rv(j,4:6),rij(i,j,1:3))/d(i,j))**2)) +&
             cm2k2*(mass(j)*(Dot_Product(rij(i,j,1:3), &
             (4._double*rv(i,4:6)-3._double*rv(j,4:6)))/d(i,j)**3*(rv(j,4:6)-rv(i,4:6)))+ &
             3.5_double *subsum3(i,j,1:3))
         end if
        end do !j
   
!    do j=bp1,totm0!massless on massive
!      frel(i,1:3)=frel(i,1:3)-rij(j,i,1:3)*mass(j)/d(j,i)**3 *&
!              (1._double+cm2*k2*(-4._double*subsum1(i)-subsum2(j,i)-5._double*mass(i)/d(j,i)+ &
!              (Dot_Product(rv(i,4:6),rv(i,4:6))-4._double*Dot_Product(rv(i,4:6),rv(j,4:6))+ & 
!              2._double*Dot_Product(rv(j,4:6),rv(j,4:6)))- &
!              1.5_double *k2*(Dot_Product(rv(j,4:6),(-rij(j,i,1:3)))/d(j,i))**2)) +&
!              cm2*k2*(mass(j)*(Dot_Product(-rij(j,i,1:3), &
!              (4._double*rv(i,4:6)-3._double*rv(j,4:6)))/d(j,i)**3*(rv(j,4:6)-rv(i,4:6)))+ &
!              3.5_double *(-subsum3(j,i,1:3)))
!    end do !j 
 end do !i
!$omp end do nowait

!massless
!$omp do
 do i=bp1,totm0
      do j=1,bmc
       frel(i,1:3)=frel(i,1:3)+rij(i,j,1:3)*mass(j)/d(i,j)**3 *&
             (1._double+cm2k2*(-4._double*subsum1(i)+subsum2(i,j)-5._double*mass(i)/d(i,j)+ &
             (Dot_Product(rv(i,4:6),rv(i,4:6))-4._double*Dot_Product(rv(i,4:6),rv(j,4:6))+ & 
             2._double*Dot_Product(rv(j,4:6),rv(j,4:6)))- &
             1.5_double *(Dot_Product(rv(j,4:6),rij(i,j,1:3))/d(i,j))**2)) +&
             cm2k2*(mass(j)*(Dot_Product(rij(i,j,1:3), &
             (4._double*rv(i,4:6)-3._double*rv(j,4:6)))/d(i,j)**3*(rv(j,4:6)-rv(i,4:6)))+ &
             3.5_double  *subsum3(i,j,1:3))
       end do !j
 end do !i
!$omp end do
!$omp end parallel
!add Yarkovsky effekt to massless accelerations:

! !$omp do
!  do i=bp1,totm0
!        !calculate vector away from the sun
!        rs(:)=-rij(i,1,1:3)/Sqrt(Dot_product(rij(i,1,:),rij(i,1,:)))
!        !auxiliary vectors yb x rs and yb x (rs x yb)
!        call crossp3d(yb(i,1:3),rs(1:3),ybxrs(1:3))
!        call crossp3d(yb(i,1:3),(-ybxrs(1:3)),ybxrsxyb(1:3))   
!        !make sure vectors are nomralized
!        ybxrs(:)=ybxrs(:)/Sqrt(Dot_product(dble(ybxrs(1:3)),dble(ybxrs(1:3))))
!        ybxrsxyb(:)=ybxrsxyb(:)/Sqrt(Dot_product(ybxrsxyb(1:3),ybxrsxyb(1:3)))
!        !Yarkovsky force
!        frel(i,1:3)=frel(i,1:3)+ & 
!        ycdiu(i)*d(i,1)**(-2._double)*(yse(i)*ybxrs(1:3)+yce(i)*ybxrsxyb(1:3)) + & !diurnal part
!        ycsea(i)*Dot_Product(yb(i,1:3),rs(1:3))* & !seansonal part
!                 d(i,1)**(-2._double) * yb(i,1:3)
!  end do !i
! !$omp end do 

!calculate heliocentric vectors and velocities for massless only
j=1
do i=bp1,totm0
  j=j+1
   drv(j,:)=rv(i,:)-rv(1,:)
   dm(j)=mass(i)
end do  
!add the sun at position 1 for htrnsel 
   drv(1,:)=0.d0
   dm(1)=mass(1)

!calculate Keplerian Orbital Elements ->  not necessary for diurnal effect
!call htrnsel (j,drv(1:j,1:6),dm(1:j),ele(1:j,1:6))

!$omp parallel default(shared)   &
!$omp private(i,l,yadiu)
!$omp do
 do i=2,j
 l=bp1+i-2
 !we only use max diurnal yark effect, spin axis is oriented perpendicular to body's orbital plane 
!$omp critical       
       call yarkdi(drv(i,1:3),drv(i,4:6),dm(1),dm(i),yark_p(l,1:10) ,yadiu(1:3),dadt_yark(l,1:2))
!      call yarkse(ele(i,1:6),yark_p(l,1:10),yasea(1:3))
!$omp end critical
       !Yarkovsky force added
       frel(l,1:3)=frel(l,1:3)+ yadiu(1:3)!+yasea(1:3) 
 end do !i
!$omp end do 
!$omp end parallel


return
end subroutine frho_yeih

end module