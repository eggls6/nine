program nine
!----------------------------------------------------------------------------
! Few-Body integration Package by Siegfried Eggl and Makrus Gyergyovits
! Version 1.8
! last update: 2011 10 26
!---------------------------------------------------------------------------- 

use global_m
use input_m
use transform_m
use output_m
use lie_m
use symp_m
use planetoid_m
use bs_m
use radau_m
use yark_m
use radeih_m


implicit none

integer(kind=intdble)::i,j,bpdum,astat
! Problem-configuration
type(Problem)::vkp
! Planetoid configuration
type(planetoids)::plt
! Orte und Geschwindigkeiten der Koerper  
real(kind=double),dimension(:,:),allocatable::rv,rvinit,dum4,rvd
! Keplerelemente der Koerper
type(kepele),dimension(:),allocatable::ele
!Massen der Koerper
real(kind=double),dimension(:),allocatable::m
!Zusatzinformationen fuer Lie und Hybrid Integratoren
type(additional)::add
!Input  Steuerungsvariable
 character(len=2)::input
 character::yn
!Dummy Variablen 
real(kind=double)::dum1,dum2,dum3(1:3),t0,rvdum(1:6),mdum,mlimit
real(kind=double),dimension(:),allocatable::md

!initialzie varialbes defined in global_m
!Initializing total error of energy and angular momentum 
dEtot=0._double
dLtot=0._double
!Initialization of total energy and angular momentum loss through system escapes
lEtot=0._double
lLtot(:)=0._double
lLsum=0._double
!allocation status
astat=0
toutprev=-1.d99

write(unit=*,fmt='(A)',advance='no')'gathering of initial conditions...'


!read in configuration
 
call ReadConfig(vkp,input,add,plt)

!Planetoids?

 if(plt%pto) then
   
    call ReadPto(plt)

    totnbody=vkp%NrKoerper+sum(plt%bodyn(:,2))

  else

    totnbody=vkp%NrKoerper
  end if

!allocate main variables

allocate(m(totnbody),rv(totnbody,1:6),&
         dum4(totnbody,1:3),rvinit(totnbody,1:6), &
         bodypos(totnbody),rvd(totnbody,1:6),md(totnbody),&
         trojan(totnbody),bname(totnbody),dadt_yark(totnbody,2),stat=astat)

  if(astat.eq.0) then
  else
   write(*,*)'Error: Allocation problem in main.f90'
   STOP
  end if

  rv(:,:)=0._double
  rvinit(:,:)=0._double
  m(:)=0._double
  dadt_yark(:,:)=0
!initial time (necessary for initial condition output)
  t0=0._double

  !read in initial conditions
 select case (input)

  case('he','el')
       allocate(ele(vkp%NrKoerper),stat=astat)      
       call ReadIniEL(ele,m(1:vkp%NrKoerper),vkp%NrKoerper)
       call htrnsko(vkp%NrKoerper,ele,m(1:vkp%NrKoerper),rv(1:vkp%NrKoerper,:))

  case('te','tr')
       allocate(ele(vkp%NrKoerper),stat=astat)            
       call ReadIniELT(ele, m(1:vkp%NrKoerper), vkp%NrKoerper)                                            
       call trojtrnsko(vkp%NrKoerper,ele,m(1:vkp%NrKoerper),rv(1:vkp%NrKoerper,:))

  case('ye')
          allocate(ele(vkp%NrKoerper),yark_p(vkp%NrKoerper,1:10),stat=astat)
          call ReadIniYEL(ele,m(1:vkp%NrKoerper), vkp%NrKoerper)
          call htrnsko(vkp%NrKoerper,ele,m(1:vkp%NrKoerper),rv(1:vkp%NrKoerper,:))

  case('be','bi')
       allocate(ele(vkp%NrKoerper),stat=astat)
       call ReadIniEL(ele,m(1:vkp%NrKoerper),vkp%NrKoerper) 
       !transform elements to rv vectors
       call btrnsko(vkp%NrKoerper,ele,m(1:vkp%NrKoerper),rv(1:vkp%NrKoerper,:))
       !correct binary position
       call binbakoo(rv(1:vkp%NrKoerper,:),m(1:vkp%NrKoerper),vkp%NrKoerper)
       ! call bakoo(rv(1:vkp%NrKoerper,:),m(1:vkp%NrKoerper),vkp%NrKoerper)
                                              
  case('rv','ve','co')
       call ReadIniRV(rv(1:vkp%NrKoerper,:), m(1:vkp%NrKoerper), vkp%NrKoerper)        
       
  case default
     write(unit=*,fmt=*)' '
     write(unit=*,fmt=*)'initial condition error: unknown type of input (elements or rv-coordinates)'
  end select
  
  
  if(astat.eq.0) then
  else
   write(*,*)'Error: Allocation problem in main.f90'
   STOP
  end if
  
!  if(vkp%NrInt.eq.6.and.input.ne.'ye') then
!     write(*,*)'Wrong input format for integrator! Choose "ye" as input elements in config.inn'
!     STOP
!  end if 
!  
!  if(vkp%NrInt.eq.8.and.input.ne.'ye') then
!     write(*,*)'Wrong input format for integrator! Choose "ye" as input elements in config.inn'
!     STOP
!  end if 

  if(plt%pto) then
      
      call ptosinit(vkp,plt,rv,m)

  end if


!Checking for massless particles and generating a linked list for reordering massless particles at the end of the
! rv-array
 
  do i=1,totnbody
     bodypos(i)=i
  end do

!masslimit below which particles are considered to have not mutual interaction (NICE Model)
  mlimit=vkp%mlimit

i=1 
  do while (i<totnbody)
     if(m(i).le.mlimit.and.maxval(m(i:totnbody)).gt.mlimit) then
        rvdum(:)=rv(i,:)
        mdum=m(i)
        bpdum=bodypos(i)
        do j=i+1,totnbody
           rv(j-1,:)=rv(j,:)
           m(j-1)=m(j)
           bodypos(j-1)=bodypos(j)
        end do
        rv(totnbody,:)=rvdum(:)
        m(totnbody)=mdum
        bodypos(totnbody)=bpdum
     else
        i=i+1
     end if
  end do
  

  m0count=0
  do i=1,totnbody
     if(m(i).le.mlimit) then
        m0count=m0count+1
     end if
  end do
  write(unit=*,fmt=*)' '
  write(unit=*,fmt=*)'number of massive bodies:',totnbody-m0count
  write(unit=*,fmt=*)'number of massless bodies:',m0count
  write(unit=*,fmt=*)' '

  if(totnbody-m0count.le.1) then
       write(unit=*,fmt=*)'You are going to calculate one massive body only.' 
       write(unit=*,fmt=*)'Are you sure you want to waste CPU time on a calculation like that? yes/no' 
       read(unit=*,fmt=*)yn
       if(yn.eq.'y') then
                write(unit=*,fmt=*)'So be it!'
       else
                STOP
       end if
  end if                    
        
        
!---------Cutoffradius active?-------------------------------
  if (vkp%cutoff.gt.0._double) then
     write(unit=*,fmt=*)'Cutoff Radius = ',vkp%cutoff
  else
     write(unit=*,fmt=*)'Cutoff Radius deactivated'
  end if
!----------------------------------------------------
  write(unit=*,fmt=*)' '

!--------Merging active?----------------------------
  if (vkp%merge) then
     write(unit=*,fmt=*)'Merging activated'
  else
     write(unit=*,fmt=*)'Merging deactivated'
  end if
!--------------------------------------------------
   write(unit=*,fmt=*)' '

!Transformations 
  if(vkp%center.eq.'b') then
                  call bakoo(rv,m,totnbody)                                             
  elseif(vkp%center.eq.'h') then                                      
                  call hekoo(rv,totnbody)
  end if
 
 !Initial values of Energy and Angular Momentum
if(vkp%outegy) then
       initegy=0._double
       initlp=0._double
       call energy(rv,m,totnbody,initegy,dum1,dum2)   
       call drehimpuls(rv,totnbody,m,dum4,initlp,dum3)
       if(initegy.eq.0._double.or.initlp.eq.0._double) then
          write(unit=*,fmt='(A)')'initial condition error: initial total energy or inital total angular momentum euqal to 0'
          write(unit=*,fmt='(A)')'delta energy or delta angular momentum will produce wrong results!'
       end if
end if


  close(21)
  write(unit=*,fmt=*)'done'
   
  rvinit(:,:)=rv(:,:)

 if(plt%pto) then
!velocities of the bodies from config.inn file have to have velocities adapted
    rv(1:vkp%NrKoerper,4:6)=rv(1:vkp%NrKoerper,4:6)/k
!all bodies (config.inn and planetoid.inn) with and without mass are set to be vkp%NrKoerper
    vkp%NrKoerper=totnbody
 else
  rv(:,4:6)=rv(:,4:6)/k
 end if

 
 !-------------- negative times for backward integration -------------------        
 if(tend.le.0.d0) then
 write(unit=*,fmt=*)'Integration backward in time, sign of velocities switched' 
 do i=1,totnbody  
     rv(:,4:6)=-rv(:,4:6) 
  end do
  tend=abs(tend) 
 end if


! CHOICE OF INTEGRATOR
  select case (vkp%NrInt)
 
  case (0) !Candy
     write(*,*)'Candy 4th Order Symplectic working...'
     call OpenFile(vkp,input,"Candy")
     call Out(t0,rv,m,vkp)
     call Int_Candy(vkp, rv, m)
     
   case (1) !Yoshida
     write(*,*)'Yoshida 8th Order Symplectic working...'
     call OpenFile(vkp,input,"Yoshi")
     call Out(t0,rv,m,vkp)
     call Int_Yoshida(vkp, rv, m)

     
  case(2)!Liereihen Integrator mit Schrittweitensteuerung
     write(*,*)'Lie Integrator with adaptive stepsize  working...'
     call OpenFile(vkp,input,'LieSW')
     call Out(t0,rv,m,vkp)
     call Int_LieSW(vkp, rv, m,add)
     
!      
!   case(3)!Experimenteller Lie Integrator mit Termsteuerung
!      write(*,*)'Experimental Symmetric Lie Series Integrator with fixed stepsize working..'
!      call OpenFile(vkp,input,'LieSY')
!      call Out(t0,rv,m,vkp)
!      call Int_LieSY(vkp, rv, m,add)
! 
! 
!  case(4)!Experimenteller Lie-Integrator mit Schrittweiten- und Termsteuerung
!      write(*,*)'Experimental Lie Integrator with adaptive stepsize and adaptive number of terms  working...'
!      call OpenFile(vkp,input,'LieST')
!      call Out(t0,rv,m,vkp)
!      call Int_LieST(vkp, rv, m,add)           
  

 case(3)!Bulirsch Stoer
     write(*,*)'Bulirsch Stoer working...'
     call OpenFile(vkp,input,'BulSt')
     call Out(t0,rv,m,vkp)
     call Int_BS(vkp, rv, m)

 case(4)!RA-15 (Everhart Gauss Radau)
      write(*,*)'Radau-15 working...'
     call OpenFile(vkp,input,'Radau')
     call Out(t0,rv,m,vkp)
     call Int_Radau(vkp, rv, m,add) 

 case(5)!Relativistic (shperical heliocentric metric) RA-15 (Everhart Gauss Radau)
      write(*,*)'Relativistic Radau-15 working...'
     call OpenFile(vkp,input,'Rarel')
     call Out(t0,rv,m,vkp)
     call Int_Radau_rel(vkp, rv, m,add) 

! case(6)!Relativistic (spherical heliocentric) RA-15 (Everhart Gauss Radau) - max Yarkovsky effect
!      write(*,*)'Relativistic Radau-15 with Yarkovsky working...'
!     call OpenFile(vkp,input,'Rarya')
!     call Out(t0,rv,m,vkp)
!     call Int_Radau_yrel(vkp, rv, m,add) 
     
 case(6)!Relativistic RA-15 (Everhart Gauss Radau) with Einstein Infeld Hoffman barycentric N-Body metric
      write(*,*)'EIH Relativistic Radau-15  working...'
     call OpenFile(vkp,input,'Raeih')
     call Out(t0,rv,m,vkp)
     call Int_Radau_eih(vkp, rv, m,add) 
     
     
! case(8)!Relativistic RA-15 (Everhart Gauss Radau) with Einstein Infeld Hoffman barycentric N-Body metric + max Yarkovsky effect
!      write(*,*)'EIH Relativistic Radau-15 with Yarkovsky working...'
!     call OpenFile(vkp,input,'Reihy')
!     call Out(t0,rv,m,vkp)
!     call Int_Radau_yeih(vkp, rv, m,add) 


 case default 
        write(*,*)'No integration method specified in config.inn -  terminating...'    
  end select


!DEFAULT OUTPUT  
  if (vkp%outdef) then
     call default(vkp,rvinit,initegy,initlp)
  end if
    
  call CloseFile

  deallocate(m,rv, &
             dum4,rvinit, &
             bodypos,rvd,md,&
             trojan,bname)
  if(input.eq.'ye') then
   deallocate(dadt_yark,ele,yark_p)
  end if

  write(unit=*,fmt=*)'...done'

end program 

