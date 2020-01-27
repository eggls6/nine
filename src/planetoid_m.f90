module planetoid_m

use global_m
use transform_m

implicit none

public::ptosinit

contains
!********************************************************

subroutine ptosinit(p,plt,rv,m)
!-------------------------------------------------------------------------------------------------------------------
!    generates positioning of planetoids around all heavy bodies with respect to reference coordinate system
!    returns: ptosrv
!
!    dependencies: equidist, randdist, makeom
!-------------------------------------------------------------------------------------------------------------------

implicit none

type(Problem)::p
type(planetoids)::plt
real(kind=double),intent(inout),dimension(:,:)::rv
real(kind=double),intent(inout),dimension(:)::m

integer(kind=intdble)::i,j,l,allocstat,nhb,n_plt,totnbody
real(kind=double),allocatable,dimension(:,:,:)::ptosrv
real(kind=double),allocatable,dimension(:,:)::dumrv,angle,pltmass2
real(kind=double),dimension(1:3,1:3)::om
real(kind=double),allocatable,dimension(:)::brv,pltmass
real(kind=double)::hbmass


nhb=p%NrKoerper
totnbody=nhb+sum(plt%bodyn(:,2))


allocate(ptosrv(1:plt%n,1:int(maxval(plt%bodyn(:,2))),1:6),angle(1:plt%n,1:2),brv(1:6), &
         pltmass2(1:int(maxval(plt%bodyn(:,2))),1:plt%n), stat=allocstat)
           
         if (allocstat.ne.0) then
            write(*,*)'error in planetoid_m: allocating ptosrv: please reduce number of planetoids'
         end if 

ptosrv(:,:,:)=0.d0
pltmass2(:,:)=0.d0
!------------------------------------------------------------
!             planetoid positioning around respective heavy body
!------------------------------------------------------------
do i=1,plt%n !number of discs

   n_plt=int(plt%bodyn(i,2))!number of planetoids in disc i
   hbmass=m(i)    !mass of heavy body central to disc

   allocate(dumrv(1:n_plt,1:6), pltmass(1:n_plt), stat=allocstat)

      if (allocstat.ne.0) then
            write(*,*)'error in planetoid_m: allocating dumrv'
      end if

      dumrv(:,:)=0.d0
      pltmass(:)=0.d0
        
!-------------------------------------------------------------
!   which kind of planetoid spaceing?
!--------------------------------------------------------------


   call distr(plt,i,n_plt,hbmass,dumrv,pltmass)

    do j=1,n_plt
        pltmass2(j,i)=pltmass(j)
    end do

!--------------------------------------------------------------------
!    use orientation-matrix
!--------------------------------------------------------------------
   angle(i,1)=plt%phi(i)
   angle(i,2)=plt%incl(i)

   call makeom(angle(i,:),om)


!$omp parallel do default(shared)   &
!$omp private(j)
   do j=1,plt%bodyn(i,2)
      dumrv(j,1:3)=matmul(om(:,:),dumrv(j,1:3))
      dumrv(j,4:6)=matmul(om(:,:),dumrv(j,4:6))
    end do
!$omp end parallel do

!-------------------------------------------------------------------
!         add positions and velocities of the heavy bodies
!------------------------------------------------------------------
    if (plt%bodyn(i,1).gt.0) then

       do j=1,plt%bodyn(i,2)     
          ptosrv(i,j,:)=dumrv(j,:)+rv(plt%bodyn(i,1),:)
       end do
    
    else
! calculating barycenter of heavy bodies
       brv(:)=0.d0
       do j=1,nhb
        brv(:)=brv(:)+rv(j,:)*m(j)
       end do
        brv(:)=brv(:)/sum(m(1:nhb))
! adding vector of barycenter to ring planetoids       
       do j=1,plt%bodyn(i,2)  
           ptosrv(i,j,:)=dumrv(j,:)+brv(:)
       end do
     end if
       
   deallocate(dumrv,pltmass)

end do !disc around body i

!rearranging planetoid tensor to fit into rv matrix

l=nhb

do j=1,plt%n
   do i=1,plt%bodyn(j,2)
      l=l+1
      rv(l,:)=ptosrv(j,i,:)
      m(l)=pltmass2(i,j)
   end do
end do


  deallocate(ptosrv,brv)


return   
end subroutine ptosinit



!*****************************************************************************
subroutine makeom(ang,om)
!--------------------------------------------------------------------------------
!  produce orientation matrix for angels: 
!  ang(1)=phi=horizontal angle
!  ang(2)=theta=vertical angle
!  with respect to the reference plane
!--------------------------------------------------------------------------------
implicit none

real(kind=double),intent(in)::ang(1:2)
real(kind=double),intent(out)::om(1:3,1:3) 
real(kind=double)::dang(1:2)

dang(:)=ang(:)/180.d0*pi

om(1,:)=(/cos(dang(2))*cos(dang(1)),-sin(dang(1)),sin(dang(2))*cos(dang(1))/)
om(2,:)=(/cos(dang(2))*sin(dang(1)),cos(dang(1)) ,sin(dang(2))*sin(dang(1))/)
om(3,:)=(/-sin(dang(2))            , 0._double         ,cos(dang(2))/)

return
end subroutine makeom


!**********************************************************************************
subroutine distr(plt,discid,n_plt,hbmass,posvel,pltmass)
!-------------------------------------------------------------------------
! supposed to produce an inital axisymmetric a*r^b and exp(a*r^b) distribution 
! from equidistributed pseudo random numbers
!-------------------------------------------------------------------------
implicit none

type(planetoids)::plt
type(kepele),dimension(1:2)::kele
real(kind=double),intent(in)::hbmass
integer(kind=intdble),intent(in)::n_plt,discid
integer(kind=intdble)::i,j,l,n_p,n_r,cycnt

real(kind=double),dimension(1:6)::discin,dr
real(kind=double),dimension(1:6)::ele,ran
real(kind=double)::rhill,pairdist,pi,pol_r
real(kind=double),dimension(:,:),allocatable::rgrid
real(kind=double),dimension(:,:),intent(inout)::posvel
real(kind=double),dimension(1:2,1:6)::dumposvel
real(kind=double),dimension(1:2)::dummass
real(kind=double),dimension(:),intent(inout)::pltmass
logical::packed_ok

!constant
pi=4.d0*atan(1.d0)


!initialization of global disc properties
  n_p=n_plt
  n_r=plt%gridsize(discid) !default 50
  

allocate(rgrid(1:n_r+1,1:6))

posvel(:,:)=0.d0
ele(:)=0.d0


call probability_distr(plt,n_r,discid,rgrid,dr,discin) 

call makepltmass(plt,n_p,discid,pltmass)

i=1
packed_ok=.false.
cycnt=1

!---------------------------------------------------------------
loop1: do while(i.lt.n_p)
     
   if(packed_ok) then
         i=i+1  !i... id of particle
         cycnt=1
   end if

!generate actual distribution of keplerian elements and planetoid masses
 

  loop2: do l=1,6 !l... 1-6 Keplerian elements, 7 mass
    call random_number(ran(:))
       
    loop3:  do j=1,n_r  !j... disc element

    if(ran(l).ge.sum(rgrid(1:j,l)).and.ran(l).lt.sum(rgrid(1:j+1,l)))then
      ele(l)=discin(l)+(dble(j-1)+ran(2))*dr(l) 
      continue
    end if

   end do loop3  
  end do loop2



  kele(2)%a=ele(1)
  kele(2)%e=ele(2)
  kele(2)%i=ele(3)
  kele(2)%kom=ele(4)
  kele(2)%gom=ele(5)
  kele(2)%man=ele(6)

  dummass(1)=hbmass
  dummass(2)=pltmass(i)

  call htrnsko(int(2,kind=intdble),kele,dummass,dumposvel)
      
  posvel(i,1:3)=dumposvel(2,1:3)
  posvel(i,4:6)=dumposvel(2,4:6)/k


!-----------------------------------------------------------
!make sure planetoids are not packed too closely (hill's radii)
!---------------------------------------------------------- 
 packed_ok=.true.

!$omp parallel do default(private)   &
!$omp shared(posvel,plt,pltmass,hbmass,packed_ok,cycnt,i,discid)
 packed: do j=1,i-1

   dumposvel(1,1:3)=posvel(j,1:3)-posvel(i,1:3) !calculate distance vector r_ij 

   pairdist=Sqrt(Dot_product(dumposvel(1,1:3),dumposvel(1,1:3)))
 
   pol_r=Sqrt(Dot_product(posvel(i,1:3),posvel(i,1:3)))
! maximum Hill's radius of particle i or j, multiplied by xrhill gives minimum allowed distance
! between particles 
   rhill=plt%xrhill(discid)*max(pol_r*(pltmass(i)/(3.d0*hbmass))**(1.d0/3.d0), &
             pol_r*(pltmass(j)/(3.d0*hbmass))**(1.d0/3.d0))


  if(pairdist.lt.rhill) then       
        packed_ok=.false.
        cycnt=cycnt+1
   end if 


   if(cycnt.ge.1d6) then
    write(*,*)'could not place particle',i
    write(*,*)'reduce spaceing restrictions'
    STOP
   end if
 end do packed
!$omp end parallel do

end do loop1

 write(*,*)' '
 write(*,*)'total discmass [M_sun]',sum(pltmass(:))
 write(*,*)'max, min mass of planetoids in  disc [M_sun]',maxval(pltmass(:)),minval(pltmass(:))

return
end subroutine

!*****************************************************
subroutine probability_distr(plt,n_r,discid,rgrid,dr,discin)
type(planetoids)::plt
integer(kind=intdble),intent(in)::n_r,discid
integer(kind=intdble)::i,j
integer(kind=intdble),dimension(1:6)::exp_or_r

real(kind=double),dimension(1:6)::a,b,discout,r0
real(kind=double),dimension(1:6),intent(out)::discin,dr
real(kind=double)::discpos,pi
real(kind=double),dimension(:,:),intent(inout)::rgrid
real(kind=double),dimension(1:n_r,1:6)::rgrid2


pi=4.d0*atan(1.d0)
!calculate probability distributions for each keplerian element
do i=1,6
!initialization of element distribution parameters
  a(i)=plt%a(discid,i)
  b(i)=plt%b(discid,i) 
  r0(i)=plt%r0(discid,i) 
  exp_or_r(i)=plt%exp_or_r(discid,i) !0 or 1 (masses 2,3)

!disc borders and element-spacings
  discin(i)=plt%distinn(discid,i) 
  discout(i)=plt%distout(discid,i) 
  dr(i)=abs(discin(i)-discout(i))/dble(n_r)  
end do


!--------------------------------------------------------
!error messages
if(minval(discin(:)).lt.0.d0) then
    write(*,*)'ERROR: module planetoid_m'
    write(*,*)'Keplerian element boundaries in disc ',discid,'cannot be smaller than 0' 
    STOP
end if

if(discout(2).gt.1.d0) then
        write(*,*)'ERROR: module planetoid_m'
        write(*,*)'eccentricity-boundaries in disc ',discid,' beyond 1' 
         STOP
end if

if(discout(3).gt.180.d0) then
        write(*,*)'ERROR: module planetoid_m'
        write(*,*)'inclination boundaries in disc ',discid,' beyond 180 deg' 
         STOP
end if

if(maxval(discout(4:6)).gt.360.d0) then
        write(*,*)'ERROR: module planetoid_m'
        write(*,*)'Keplerian angle boundaries in disc ',discid,' beyond 360 deg' 
         STOP
end if
!----------------------


!make distribution
!$omp parallel do default(shared)   &
!$omp private(i,j,discpos)
do j=1,6
 do i=1,n_r
   !calculate distribution position
  discpos=discin(j)+dble(i)*dr(j)-r0(j)

  !choose distribution
  if(exp_or_r(j).eq.1) then
     rgrid(i,j)=exp(a(j)*discpos**b(j)) 
  else
     rgrid(i,j)=a(j)*discpos**b(j)      
  end if

  !non negative densities 
  if(rgrid(i,j).lt.0.d0) then
        rgrid(i,j)=0.d0
  end if

  !NAN and INFINITY filter
  if (rgrid(i,j).ne.rgrid(i,j).or.rgrid(i,j).ge.huge(rgrid(i,j)).or. &
           rgrid(i,j).le.-huge(rgrid(i,j))) then
         write(*,*)'ERROR: module planetoid_m'
         call whosfault(j)
         write(*,*)'Signularity in density distribution encountered'

         STOP 
   end if
 end do !i...disc elements


if(sum(rgrid(:,j)).le.0.d0) then
    write(*,*)'ERROR: module planetoid_m'
    call whosfault(j)
    write(*,*)'Distribution completely out of disc boundaries.'
    STOP
end if

 rgrid2(1:n_r,j)=rgrid(1:n_r,j)
 
 rgrid(1,j)=0.d0
 
 do i=1,n_r    
     rgrid(i+1,j)=rgrid2(i,j)
 end do 

!normalize grid so that Sum(f_i(x_i))=1
    rgrid(:,j)=rgrid(:,j)/sum(rgrid(:,j))

end do !j...elements
!$omp end parallel do


return

end subroutine
!******************************************************************************
subroutine makepltmass(plt,n_p,discid,pltmass)
!---------------------------------------------------------------------
! produce and nomralize planetoid masses
!---------------------------------------------------------------------
type(planetoids)::plt
integer(kind=intdble)::i
integer(kind=intdble),intent(in)::n_p,discid
real(kind=double),dimension(:),intent(inout)::pltmass
real(kind=double),dimension(1:2)::ran


!mass distribution logarithmic (to bridge large scales)

if (plt%exp_or_r(discid,7).eq.2) then
  if (plt%distinn(discid,7).le.1.d-35) then
    plt%distinn(discid,7)=1.d-35
  end if

  if (plt%distout(discid,7).le.1.d-35) then
    plt%distout(discid,7)=1.d-35
  end if
end if

if(plt%massmax(discid).le.0.d0.or.plt%distout(discid,7).le.0.d0) then
  pltmass(:)=0.d0
else
  do i=1,n_p
    call random_number(ran(:))
                     !logarithmic equidistribution for masses? 
    if(plt%exp_or_r(discid,7).eq.2) then  !log distribution
      pltmass(i)=10.d0**(log10(plt%distinn(discid,7))+ran(2)*log10(abs(plt%distout(discid,7)- &
                     plt%distinn(discid,7))))
    else !equidistribution
      pltmass(i)=plt%distinn(discid,7)+ran(2)*abs(plt%distout(discid,7)-plt%distinn(discid,7)) 
    end if
  end do
pltmass(:)=pltmass(:)/sum(pltmass(:))*plt%massmax(discid)
end if

return

end subroutine

!******************************************************************

subroutine whosfault(j)
integer(kind=intdble),intent(in)::j

select case(j)
  case(1)
    write(*,*)'distribution parameter: semi-major axis'
  case(2)
    write(*,*)'distribution parameter: eccentricity'
  case(3)
    write(*,*)'distribution parameter: inclination'
  case(4)
    write(*,*)'distribution parameter: argument of pericenter'
  case(5)
    write(*,*)'distribution parameter: arument of the ascending node'
  case(6)
    write(*,*)'distribution parameter: mean anomaly'
  case(7)
    write(*,*)'distribution parameter: mass'
end select

return
end subroutine

!//////////////////////////////////////////////////////////////////////

end module planetoid_m
