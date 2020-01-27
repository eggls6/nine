module global_m
  implicit none

!/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!                                        global variables definition
!
!                                         by Siegfried Eggl  20080110 
!------------------------------------------------------------------------------
!********************************************************************

! Definition of Precision

  integer,parameter::double=selected_real_kind(18,60)
  integer,parameter::intdble=selected_int_kind(16)
! Position Array
  integer(kind=intdble), dimension(:),allocatable::bodypos
! Number of Massless Bodies  
  integer(kind=intdble)::m0count
! Total number of Integration Objects
  integer(kind=intdble)::totnbody
  ! PI
  real(kind=double),parameter::PI=3.14159265358979323846264338328_double
  ! Gaussian Gravitational Constant
  real(kind=double),parameter::K=0.01720209895_double
  ! Degree to Rad conversion
  real(kind=double),parameter::grad=0.0174532925199432957692369076849_double
  !rad to degree conversion
  real(kind=double),parameter::invgrad=57.2957795130823208767981548141_double
  

  !Initial energy and angular momentum (lp)
  real(kind=double)::initegy,initlp

 !Total deviation of Energy and Angular Momentum
  real(kind=double)::dEtot,dLtot
 
 !energy and angular momentum lost through particle escape (cutoff)
  real(kind=double)::lEtot,lLsum,lLtot(1:3)

 !print progress to stdout, stop calculation when particle goes beyond cutoff
  logical::showprog,cutoffstop(2)

 !time output style variables
  character(len=2)::tdigc,tendc

 !time output dummy variable
  real(kind=double)::toutprev
  
 ! Auxilliary data for N-Body problem
  type Problem
        ! Id of integration algorithm /eg. 2...Lie Series
        integer(kind=intdble)::NrInt
        ! Number of Bodies
        integer(kind=intdble)::NrKoerper
        ! start time
        real(kind=double)::tstart
        ! End Time
        real(kind=double)::tEnde
        ! Output interval
        real(kind=double)::tAusg
        ! stepsize for fixed step integrators
        real(kind=double)::Ns
        ! one step precision
        real(kind=double)::eps
        !Cutoffradius for exclusion of particle form the system
        real(kind=double) ::cutoff 
      
        !Coordinate center (not used)
        character::center 
        !output
        logical::outegy,outbco,outhco,outbele,outhele,outdef
        logical::outbinbco,outbinhele,outbinhco
        logical::outtroj
        !mergeing
        logical::merge
        real(kind=double)::mergermin,mergermax
        !NICE MODEL minimum mass
        real(kind=double)::mlimit
          !close encouter output
        logical::cen
        !close encounter distance for output
        real(kind=double)::cend,cemlim
        
 end type Problem

!Kepler Element-Type
  type kepele
                       real(kind=double)::a,e,i,kom,gom,man
  end type kepele

!Additional Variables
  type additional
    integer(kind=intdble)::lterms,ltmin,ltmax
    real(kind=double)::hyrlimit,lstepmin
    logical::autodt
  end type additional

!planetoids
  type planetoids
     logical::pto
     integer(kind=intdble)::n
     integer,allocatable,dimension(:,:)::bodyn,exp_or_r
     integer,allocatable,dimension(:)::box_flare,gridsize
     real(kind=double),allocatable,dimension(:)::xrhill, &
             incl,phi,dtheta_z,massmax
     real(kind=double),allocatable,dimension(:,:):: &
             distout,distinn,a,b,r0
   
  end type planetoids
!trojan motion

  integer(kind=intdble),allocatable,dimension(:)::trojan
  !start and end times
    real(kind=double)::tstart,tend 
    
!Yarkovsky (see Bottke)
   !solar constant 1367 W/m^2
   real(kind=double),parameter::solconst=1367._double   
    !C_diurnal, C_seasonal, Sin(eps_lag), Cos(eps_lag)
     real(kind=double),allocatable,dimension(:,:)::yark_p
!    for output drift in semimajor axis due to yarkovsky effect [AU/Myrs] 
!     dadt_yark(1): diurnal effect,     dadt_yark(2): seasonal effect.
     real(kind=double),allocatable,dimension(:,:)::dadt_yark

!     !spin direction (stays constant)
!     real(kind=double),allocatable,dimension(:,:)::yb
!     real(kind=double)::yeps


!Name for bodies

    character(len=15),dimension(:),allocatable::bname

!integration backwards?
   logical::intback
    
end module global_m
