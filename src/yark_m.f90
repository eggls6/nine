! ==========MODULE yark_m============
! CONTAINS
!                 yarkdi   yarkovsky diurnal                            
!                 yarkse   yarkovsky seasonal                           
!               yarkinit   initialisation of yarkovsky 
! MODULES AND HEADERS
! yark_pert.o: \
! ../include/sysdep.h90 \
! ../suit/FUND_CONST.mod 

MODULE yark_m
 !USE fund_const
 !USE YARKO, only : yarkp
 !USE CONST
 USE global_m
 
 
! USE Tools, only : scal, vsize

 IMPLICIT NONE
 !PRIVATE

 !PUBLIC yarkdi, yarkse, yarkinit, secular_nongrav
 PUBLIC yarkdi, yarkse, init_yark !,  secular_nongrav

 ! private common data
 ! former yarkom.h
 ! controls and directory for Yarkovski force
 integer(kind=intdble) iyark,iyarpt
! CHARACTER(len=80) yardir

 ! fromer yarkov.h 
 ! common containing all physical information
 ! on the current asteroid needed to compute Yarkovsky acceleration
 !      real(kind=double) yarkp(10),beya(7),alya(7)
 real(kind=double) beya(7),alya(7)

 real(kind=double) spya,sqya,etaya75,fmeaya,thfacya,radfluya
 ! logical flggs: availability of physical data, 
 !  has yarkinit routine been called
 LOGICAL yarfil,yarini
 ! secular perturbation in semimajor axis 
! real(kind=double) dadt(2)
! PUBLIC iyark,iyarpt,yardir,yarfil,yarini,dadt

CONTAINS 


 SUBROUTINE init_yark(yarkp)

   real(kind=double), INTENT(INOUT), DIMENSION(3) :: yarkp
   real(kind=double),                DIMENSION(2) :: eclip_axis

   eclip_axis(1) = yarkp(1)*grad !! longitude
   eclip_axis(2) = yarkp(2)*grad !! latitude


   yarkp(1) = cos(eclip_axis(2))*cos(eclip_axis(1))
   yarkp(2) = cos(eclip_axis(2))*sin(eclip_axis(1))
   yarkp(3) = sin(eclip_axis(2))

 END SUBROUTINE init_yark


 ! =========================================================
 ! SECACC   added 17/9/2008
 ! secular perturbation on semimajor axis, presumably due to
 ! non gravitational perturbations (including Yarkovsky)
 ! implemented as acceleration along the velocity
!  SUBROUTINE secular_nongrav(x,v,secacc)
!    !  interface: INPUT
!    real(kind=double), INTENT(IN), DIMENSION(3) :: x,v ! position and velocity, heliocentric
!    !  interface: OUTPUT
!    real(kind=double), INTENT(OUT), DIMENSION(3) :: secacc
!    ! end interface
!    !        real(kind=double) :: vvec(3),vv2, rr, vsize, prscal, factor, conv
!    real(kind=double) :: rr,rr2,vv2,rv,angm_x,angm_y,angm_z,angm2
!    real(kind=double) :: factor,conv,yark_acc,trans_size,trans(3)
!    real(kind=double) :: vsize, prscal
!    real(kind=double) :: gms = const_gauss**2
!    ! =========================================================
!    !! first attempt (David V.)
!    !         rr=vsize(x)
!    !         vv2=prscal(v,v)
!    !         vvec=v*(1/vv2)
!    !         factor=(2*gms/rr-vv2)**2/gms
!    !! warning: dadt is in AU/My, thus it has to be converted in AU/d
!    !         conv=1.d6*365.25d0 
!    !         factor=factor/conv
!    !         secacc=0.5*factor*vvec
!    ! Steve Chesley's implementation
!    rr=sqrt(vsize(x))
!    rr2=scal(x,x,3)
!    vv2=scal(v,v,3)
!    rv =scal(x,v,3)
!    angm_x=x(2)*v(3)-x(3)*v(2)
!    angm_y=x(3)*v(1)-x(1)*v(3)
!    angm_z=x(1)*v(2)-x(2)*v(1)
!    angm2=angm_x*angm_x+angm_y*angm_y+angm_z*angm_z
!    factor=angm2*sqrt(2*gms/rr-vv2)/gms
!    ! warning: dadt is in AU/My, thus it has to be converted in AU/d
!    conv=1.d6*365.25d0 
!    yark_acc=0.5d0*factor/conv/rr2
!    trans=v-(x*rv/rr2)
!    trans_size=sqrt(vsize(trans))
!    secacc=yark_acc*trans/trans_size
! 
! 
!  END SUBROUTINE secular_nongrav


 ! ******************************************************************    
 SUBROUTINE yarkdi(xast,vast,msun,mass,yarkp,a,dadt) 
! SUBROUTINE yarkdi(xast,elkep,a,dadt_yark,iparti) 
   ! ******************************************************************    
   ! eggl 04.2012                                                                     
   ! This subroutine computes the heliocentric components of the           
   ! Yarkovsky thermal acceleration -- the diurnal variant only.           
                          
   !                                                                       
   ! Input parameters:                                                     
   ! -----------------                                                     
   !                                                                       
   ! - via header   xast(1-3) ... heliocentric coordinates of the body (in AU)
   !                vast(1-3) ... heliocentric velocities of the body (in AU/D^2/gaussk)
   !                msun      ... mass of central body (in M_sol)
   !                mass      ... mass of asteroid (in M_sol)
   !                
   !                yarkp(2) ...  (0 or 180) prograd or retrograde orientation of          
   !                               body's spin axis)      
   !                yarkp(3) ...   thermal capacity [J/kg/K] 
   !    
   !                yarkp(4-5) ... k_0 and k_1 parameters of the           
   !                               surface thermal conductivity            
   !                               [K(T) = k_0 + k_1 T_av^3]               
   !                yarkp(6) ... density of the surface layer              
   !                yarkp(7) ... radius of the body                        
   !                yarkp(8) ... rotation frequency                        
   !                yarkp(9) ... surface absorptivity   
   !                yarkp(10)... bulk density            
   !                                                                       
   !                                                                       
   ! Output parameters: a(1-3) ... diurnal acceleration      
   !              
   ! --NOT WORKING:--   a(4-6) ... partials wrt the radius of the body     
   !                    a(7-9) ... partials wrt the thermal conductivity   
   !                               parameter k_0                           
   !                    a(10-12) ... partials wrt the thermal conductivity 
   !                                 parameter k_1                         
   !                    a(13-15) ... partials wrt the x-component of the   
   !                                 spin axis unit vector                 
   !                    a(16-18) ... partials wrt the y-component of the   
   !                                 spin axis unit vector                 
   !                    a(19-21) ... partials wrt the z-component of the   
   !                                 spin axis unit vector                 
   !                                                                       
   ! SI units are assumed internally in the subroutine, but the results    
   ! (e.g. accelerations) are given in AU and days.                        
   !                                                                       
   ! Written by: D. Vokrouhlicky, Oct 99 , modified by S. Eggl 04.2012                                  
   ! (queries to vokrouhl@mbox.cesnet.cz)                                  
   ! ..................................................................    
   implicit real(kind=double) (a-h,o-z) 
   ! here we specify two more parameters that eventually might be changed: 
   ! -- the average bulk density of the body (densityb) which is now set   
   !    to 2 g/cm^3                                                        
   ! -- the heat capacity of the surface layer (capacity) which is now set 
   !    to 680 J/kg/K     
                                           
   real(kind=double), parameter :: speed_light=2.99792458d8
   real(kind=double), parameter :: stefboltz=5.66962d-8,clight3=8.99377374d8 
 !  real(kind=double), parameter :: dsqrt2=1.414213562373d0,dsqrt23=1414.213562373d0,emiss=0.9d0 
   real(kind=double):: dsqrt2,dsqrt23
   real(kind=double), parameter :: aceuni=0.049900176d0,ua=1.49597870691d11,gaussk=0.01720209895_double
   real(kind=double) :: densityb 
   ! input: asteroid position, flag for partials                           
   real(kind=double) xast(3),vast(3),msun,mass
 !  real(kind=double),intent(in):: elkep(6)
  ! integer(kind=intdble) iparti 
   ! output: acceleration and partials                                     
!   real(kind=double) a(21) 
     real(kind=double) a(3),da
   ! internal variables                                                    
   real(kind=double) aster_obliq,semia,rbar,tbar,esinf,ecosfp1,mu,ecc,hh,p,lhxr,cobl
   real(kind=double),dimension(1:3):: vprod1,vprod2,h,hxr
   real(kind=double) rau2,rau,xn,yn,zn,tstar, tav1000,surcon,bgama,theta,diudepth,rp,phi_a,rdotv,rdot
   ! physical data on the current asteroid                   
   real(kind=double),intent(in)::yarkp(1:10)        
   real(kind=double),intent(out)::dadt(1:2)       
   
              
   ! ----------------------------------------------------------------------
!eggl 04.2012  Kirchhoff: emissivity=absorbtivity
   emiss=yarkp(9)

   dsqrt2=sqrt(2._double)
   dsqrt23=sqrt(2.d6)
   
   capacity=yarkp(3)
 !  aster_obliq = acos(yarkp(3)/sqrt(yarkp(1)**2+yarkp(2)**2+yarkp(3)**2))
!eggl 04.2012 We use maximum possible Yark -> obliquity is always = 0 because the spin axis is perpendicular to the orbital plane
    
    aster_obliq=yarkp(2)
    cobl=cos(aster_obliq) 
   
   !yarkp(1:3) = spin_temp

   densityb=yarkp(10)
   rau2=xast(1)*xast(1)+xast(2)*xast(2)+xast(3)*xast(3) 
   rau=sqrt(rau2) 
   xn=xast(1)/rau 
   yn=xast(2)/rau 
   zn=xast(3)/rau 
   ! initializations & constants                                           
   radflu=solconst/rau2 
   
  !eggl 04.2012 calculate semimajor axis from vis viva equation
   vnorm=sqrt(dot_product(vast,vast))
    
   semia=1._double/(2._double/rau-vnorm*vnorm/(msun+mass))
 
 
  
   ! - mean motion & solar radiation flux at r=a                           
   fmeaya=(1.9909837d-7)/semia/sqrt(semia)  !! mean motion 1.99d-7 is gauss constant converted in AU^(3/2)/s
   radfluya=solconst/semia/semia  !! flux at r=a 1371d0 is the solar constant
   radfluya = radfluya/(UA**2)  !! flux in W/m² -> no its not  

   phi_a = 3._double*radfluya/(4._double*yarkp(7)*yarkp(10)*speed_light) !! factor phi at distance a in m/s²

   ! - subsolar temperature                                                
   tstar=(yarkp(9)*radflu/emiss/stefboltz)**0.25d0 
   tav1000=tstar/dsqrt23 
   ! - surface conductivity                                                
   surcon=yarkp(4)+yarkp(5)*(tav1000**3) 
   ! - thermal inertia & diurnal thermal parameter                         
   bgama=sqrt(surcon*yarkp(6)*capacity) 
   theta=bgama*sqrt(yarkp(8))/emiss/stefboltz/(tstar**3) 
   diudepth=sqrt(surcon/yarkp(6)/capacity/yarkp(8)) 
   ! - radius of the body scaled by the depth of the diurnal wave          
   rp=yarkp(7)/diudepth 
   al=dsqrt2*rp 

   tau=theta/al 
   tau1=1._double+tau 
   ! - the auxiliary functions A-D, a,b                                    
   cal=cos(al) 
   sal=sin(al) 
   if (al.lt.90._double) then 
      ealm=exp(-al) 
   else 
      ealm=0._double 
   endif
   af=3._double*(al+2._double)*ealm+(3._double*(al-2._double)*cal+al*(al-3._double)*sal) 
   bf=al*(al+3._double)*ealm+(-al*(al-3._double)*cal+3._double*(al-2._double)*sal) 
   caf=-(al+2._double)*ealm+(-(al-2._double)*cal+al*sal) 
   cbf=-al*ealm-(al*cal+(al-2._double)*sal) 
   ccf=caf+tau*af/tau1 
   cdf=cbf+tau*bf/tau1 
   ! - G exp(i delta) & amplitude computed                                 
   facp=aceuni*yarkp(9)*radflu/yarkp(7)/densityb/clight3 
   deno=ccf*ccf+cdf*cdf 
   deno1=deno*tau1 
   gcosd=(caf*ccf+cbf*cdf)/deno1 
   gsind=(cbf*ccf-caf*cdf)/deno1 

!********************OUTDATED SEE FURTHER DOWN******************************************
   !!    Computation of dadt, diurnal(1) and seasonal(2), for output -> via global_m

  ! dadt(1) = -((8._double*yarkp(9)*phi_a*gsind)/(9._double*fmeaya*tau1))*cos(aster_obliq) !! diurnal effect in m/s
  ! dadt(2) = -((4._double*yarkp(9)*phi_a*gsind)/(9._double*fmeaya*tau1))*(sin(aster_obliq))**2 !! seasonal effect in m/s
 !eggl 04 2012
 !Vokrouhlicky 1999 A&A 344,362-366
 !...aster_obliq -> gamma
 !...g -> E_R'
 !...phi_a -> PHI
 !...tau -> CHI, tau1 = 1+CHI
 !...yarkp(9) -> alpha


  !for output 
 !  dadt = dadt*UA*86400._double*1.d6*365.25d0 !! in UA/Myr

  ! write(*,*)'dadt diurnal AU/Myr', dadt_yark(1),'dadt seasonal AU/Myr',dadt_yark(2)
!*************************************************************

   ! geometric products                                                    
   ! - r x s                                                               
!    vprod1(1)=yn*yarkp(3)-zn*yarkp(2) 
!    vprod1(2)=zn*yarkp(1)-xn*yarkp(3) 
!    vprod1(3)=xn*yarkp(2)-yn*yarkp(1) 

   ! - s x (r x s) = r - (r.s) s   !this is wrong (-)                                         
!    scalar=xn*yarkp(1)+yn*yarkp(2)+zn*yarkp(3) 
!    vprod2(1)=xn-scalar*yarkp(1) 
!    vprod2(2)=yn-scalar*yarkp(2) 
!    vprod2(3)=zn-scalar*yarkp(3) 
!  


!eggl 04.2012 spin axis = orbital angular momentum axis ->
   !  r x s = r x h = -h x r =~ - v   
 
  mu=(msun+mass)  
    
  h(1)=xast(2)*vast(3)-xast(3)*vast(2)
  h(2)=xast(3)*vast(1)-xast(1)*vast(3)
  h(3)=xast(1)*vast(2)-xast(2)*vast(1)
  hh=Dot_Product(h,h)
 
  hxr(1)=h(2)*xast(3)-h(3)*xast(2)
  hxr(2)=h(3)*xast(1)-h(1)*xast(3)
  hxr(3)=h(1)*xast(2)-h(2)*xast(1)
  lhxr=sqrt(Dot_Product(hxr,hxr))
  
   vprod1(1)=-hxr(1)/lhxr*cobl
   vprod1(2)=-hxr(2)/lhxr*cobl
   vprod1(3)=-hxr(3)/lhxr*cobl
   

!eggl 04.2012 spin axis = orbital angular momentum axis ->   
    ! s x (r x s) = s x (-v) =  h x (-v) = (r x v) x (-v) = v x (r x v) = r - v(r.v)     (only true for unit vectors)
   rdotv=xn*vprod1(1)+yn*vprod1(2)+zn*vprod1(3)
   
   vprod2(1)=xn-vprod1(1)*rdotv
   vprod2(2)=yn-vprod1(2)*rdotv
   vprod2(3)=zn-vprod1(3)*rdotv
   
   ! diurnal acceleration                                                  
   a(1)=facp*(gsind*vprod1(1)+gcosd*vprod2(1)) 
   a(2)=facp*(gsind*vprod1(2)+gcosd*vprod2(2)) 
   a(3)=facp*(gsind*vprod1(3)+gcosd*vprod2(3)) 
 
  
 !eggl 04.2012 try to compute da/dt diurnal
  !da/dt Murray & Dermott Solar System Dynamics p.54
  
  rbar=Dot_Product(a,xast/rau)
  tbar=Dot_Product(a,hxr/lhxr)
  
  ecc=sqrt(1._double-hh/(mu*semia))
  rdot=sign(sqrt(vnorm*vnorm-hh/rau2),Dot_Product(xast,vast))
  
  p=semia*(1._double-ecc*ecc)
  esinf=p/sqrt(hh)*rdot
  ecosfp1=p/rau  
  
   da=2._double*sqrt(semia)/(sqrt(mu)*gaussk*p)*(rbar*esinf+tbar*ecosfp1)

 !for output 
   dadt(1) = da*1.d6*365.25d0 !! in AU/Myr
 
 !*******************************************************************
   !!    Computation of mean diurnal da/dt -> dadt(2) according to Vokrouhlicky 1999 A&A 344,362-366 
   !!     
    !Vokrouhlicky 1999 A&A 344,362-366
    !...aster_obliq -> gamma
    !...g -> E_R'
    !...phi_a -> PHI
    !...tau -> CHI, tau1 = 1+CHI
    !...yarkp(9) -> alpha_1
 
   !calculate average insolation onto asteroid (Eggl et al 2012, ApJ), equivalence r = a (1-e^2)
   radflu=solconst/(p*p)
  ! - subsolar temperature                                                
   tstar=(yarkp(9)*radflu/emiss/stefboltz)**0.25d0 
   tav1000=tstar/dsqrt23 
   ! - surface conductivity                                                
   surcon=yarkp(4)+yarkp(5)*(tav1000**3) 
   ! - thermal inertia & diurnal thermal parameter                         
   bgama=sqrt(surcon*yarkp(6)*capacity) 
   theta=bgama*sqrt(yarkp(8))/emiss/stefboltz/(tstar**3) 
   diudepth=sqrt(surcon/yarkp(6)/capacity/yarkp(8)) 
   ! - radius of the body scaled by the depth of the diurnal wave          
   rp=yarkp(7)/diudepth 
   al=dsqrt2*rp 

   tau=theta/al 
   tau1=1._double+tau 
   ! - the auxiliary functions A-D, a,b                                    
   cal=cos(al) 
   sal=sin(al) 
   if (al.lt.90._double) then 
      ealm=exp(-al) 
   else 
      ealm=0._double 
   endif
   af=3._double*(al+2._double)*ealm+(3._double*(al-2._double)*cal+al*(al-3._double)*sal) 
   bf=al*(al+3._double)*ealm+(-al*(al-3._double)*cal+3._double*(al-2._double)*sal) 
   caf=-(al+2._double)*ealm+(-(al-2._double)*cal+al*sal) 
   cbf=-al*ealm-(al*cal+(al-2._double)*sal) 
   ccf=caf+tau*af/tau1 
   cdf=cbf+tau*bf/tau1 
   ! - G exp(i delta) & amplitude computed                                 
   facp=aceuni*yarkp(9)*radflu/yarkp(7)/densityb/clight3 
   deno=ccf*ccf+cdf*cdf 
   deno1=deno*tau1 
   gcosd=(caf*ccf+cbf*cdf)/deno1 
   gsind=(cbf*ccf-caf*cdf)/deno1 

   dadt(2) = -((8._double*yarkp(9)*phi_a*gsind)/(9._double*fmeaya*tau1))*cos(aster_obliq) !! diurnal effect in m/s
   dadt(2)=dadt(2)*UA*86400._double*1.d6*365.25d0 !! in UA/Myr
   
  ! dadt(2) = -((4._double*yarkp(9)*phi_a*gsind)/(9._double*fmeaya*tau1))*(sin(aster_obliq))**2 !! seasonal effect in m/s
  ! write(*,*)'dadt diurnal AU/Myr', dadt_yark(1),'dadt seasonal AU/Myr',dadt_yark(2)
!************************************************************
 
 
 
   ! Partials?                                                             
!    if (iparti.eq.0) return 
!    ! - general                                                             
!    cafp=-ealm+cal+(2._double*al-1._double)*sal 
!    cbfp=-ealm-(2._double*al-1._double)*cal+sal 
!    afp=3._double*ealm+(al*al-3._double)*cal+(al*(al-4._double)+3._double)*sal 
!    bfp=(2._double*al+3._double)*ealm-(al*(al-4._double)+3._double)*cal                   &
!         &     +(al*al-3._double)*sal                                            
!    ! - thermal conductivity parameters (k_0,k_1)                           
!    xi1r=caf*ccf-cbf*cdf 
!    xi1i=cbf*ccf+caf*cdf 
!    xi2r=cafp*af-cbfp*bf 
!    xi2i=cbfp*af+cafp*bf 
!    xi2r=xi2r-caf*afp+cbf*bfp 
!    xi2i=xi2i-cbf*afp-caf*bfp 
!    deno=xi1r*xi1r+xi1i*xi1i 
!    facr=1._double+0.5d0*al*(xi2r*xi1r+xi2i*xi1i)/deno 
!    faci=     0.5d0*al*(xi2i*xi1r-xi2r*xi1i)/deno 
!    derikr=-tau*(gcosd*facr-gsind*faci)/tau1 
!    deriki=-tau*(gsind*facr+gcosd*faci)/tau1 
!    a(7)=facp*(deriki*vprod1(1)+derikr*vprod2(1)) 
!    a(8)=facp*(deriki*vprod1(2)+derikr*vprod2(2)) 
!    a(9)=facp*(deriki*vprod1(3)+derikr*vprod2(3)) 
!    a(10)=a(7)*(tav1000**3) 
!    a(11)=a(8)*(tav1000**3) 
!    a(12)=a(9)*(tav1000**3) 
!    ! - radius of the body                                                  
!    rfac=(tau+tau1)/tau1 
!    a(4)=-a(1)*rfac-2._double*a(7) 
!    a(5)=-a(2)*rfac-2._double*a(8) 
!    a(6)=-a(3)*rfac-2._double*a(9) 
!    ! - partials d_K (a), d_R (a) ...                                       
!    a(4)=a(4)/yarkp(7) 
!    a(5)=a(5)/yarkp(7) 
!    a(6)=a(6)/yarkp(7) 
!    a(7)=a(7)/surcon 
!    a(8)=a(8)/surcon 
!    a(9)=a(9)/surcon 
!    a(10)=a(10)/surcon 
!    a(11)=a(11)/surcon 
!    a(12)=a(12)/surcon 
!    ! - spin axis components                                                
!    ! ... sx                                                                
!    a(13)=-facp*gcosd*(xn*yarkp(1)+scalar) 
!    a(14)=facp*(gsind*zn-gcosd*xn*yarkp(2)) 
!    a(15)=-facp*(gsind*yn+gcosd*xn*yarkp(3)) 
!    ! ... sy                                                                
!    a(16)=-facp*(gsind*zn+gcosd*yn*yarkp(1)) 
!    a(17)=-facp*gcosd*(yn*yarkp(2)+scalar) 
!    a(18)=facp*(gsind*xn-gcosd*yn*yarkp(3)) 
!    ! ... sz                                                                
!    a(19)=facp*(gsind*yn-gcosd*zn*yarkp(1)) 
!    a(20)=-facp*(gsind*xn+gcosd*zn*yarkp(2)) 
!    a(21)=-facp*gcosd*(zn*yarkp(3)+scalar) 
   return 
 END SUBROUTINE yarkdi
 ! ******************************************************************    
SUBROUTINE yarkse(kepele,yarkp,a) 
! SUBROUTINE yarkse(xast,vast,elkep,yarkp,a) 
   ! ******************************************************************    
   !                                                                       
   ! This subroutine computes the heliocentric components of the           
   ! Yarkovsky thermal acceleration -- the seasonal variant only.          
   ! If the flag (iparti) is set to 1, one gets at the output also         
   ! partials wrt to some parameters of the thermal model (if these        
   ! are desired to be adjusted).                                          
   !                                                                       
   ! Input parameters:                                                     
   ! -----------------                                                     
   !                                                                       
   ! - via header                       
   !                (xast,vast) ... state vector of the asteroid           
   ! - via common   yarkp(1-3) ... sx, sy, sz (unit vector of the          
   !                               body's spin axis orientation)           
   !                yarkp(4-5) ... k_0 and k_1 parameters of the           
   !                               surface thermal conductivity            
   !                               [K(T) = k_0 + k_1 T_av^3]               
   !                yarkp(6) ... density of the surface layer [kg/m^3]              
   !                yarkp(7) ... radius of the body                        
   !                yarkp(8) ... rotation frequency  [s]                      
   !                yarkp(9) ...  surface absorptivity  = 1-albedo  
   !                yarkp(10)... bulk density   
   !                + some more precomputed useful variables               
   !                                                                       
   ! Output parameters: a(1-3) ... seasonal acceleration                   
   ! ------------------ a(4-6) ... partials wrt the radius of the body     
   !                    a(7-9) ... partials wrt the thermal conductivity   
   !                                                                       
   ! REM. PARTIALS ARE DISABLED AT THIS MOMENT                             
   !                                                                       
   ! SI units are assumed throughout the subroutine, but the results       
   ! (e.g. accelerations) are given in AU and days.                        
   !                                                                       
   ! Written by: D. Vokrouhlicky, Oct 99                                   
   ! (queries to vokrouhl@mbox.cesnet.cz)                                  
   ! ..................................................................    
   implicit real(kind=double) (a-h,o-z) 
   ! modified 17/4/2007, to comply with the Golevka code at JPL.
   integer(kind=intdble), parameter:: napprox=7
   !      integer(kind=intdble), parameter:: napprox=12
   parameter (capacity=680._double,stefboltz=5.66962d-8,speed_light=2.99792458d8) 
   parameter (clight3=8.99377374d8,aceuni=0.049900176d0,ua=1.49597870691d11) !emiss=0.9d0
  ! dimension xast(3),vast(3) 
   real(kind=double) elkep(6),pvya(3),qvya(3), nvya(3)
   dimension brac(napprox),bras(napprox),gcosd(napprox),gsind(napprox),a(3)
   real(kind=double) densityb, dsqrt2 
!   integer(kind=intdble) iparti
   integer(kind=intdble) i
!   real(kind=double) ::  =k**2
   real(kind=double),intent(in)::yarkp(1:10),kepele(6)
   ! ----------------------------------------------------------------------
!eggl 04.2012  Kirchhoff: emissivity=absorbtivity
   emiss=yarkp(9) 
 
 !eggl 04.2012 reordering keplerian elements (switch of omega <-> Omega for this routine + degree->rad)
   elkep(1:2)=kepele(1:2)
   elkep(3)=kepele(3)*grad
   elkep(4)=kepele(5)*grad
   elkep(5)=kepele(4)*grad
   elkep(6)=kepele(6)*grad
 
 
   densityb=yarkp(10)

   dsqrt2=sqrt(2._double)
   ! - mean motion & solar radiation flux at r=a                           
   fmeaya=(1.9909837d-7)/elkep(1)/sqrt(elkep(1)) 
   radfluya=solconst/elkep(1)/elkep(1) 


   ! - subsolar temperature                                                
   tstarya=(yarkp(9)*radfluya/emiss/stefboltz)**0.25d0 
   thfacya=sqrt(fmeaya)/stefboltz/(tstarya**3) 
   ! - projections s_P and s_Q of the spin axis computed                   
   pvya(1)=cos(elkep(4))*cos(elkep(5))-                         &
        &           cos(elkep(3))*sin(elkep(4))*sin(elkep(5))           
   pvya(2)=sin(elkep(4))*cos(elkep(5))+                         &
        &           cos(elkep(3))*cos(elkep(4))*sin(elkep(5))           
   pvya(3)=sin(elkep(3))*sin(elkep(5)) 
   qvya(1)=-cos(elkep(4))*sin(elkep(5))-                        &
        &           cos(elkep(3))*sin(elkep(4))*cos(elkep(5))           
   qvya(2)=-sin(elkep(4))*sin(elkep(5))+                        &
        &         cos(elkep(3))*cos(elkep(4))*cos(elkep(5))             
   qvya(3)=sin(elkep(3))*cos(elkep(5)) 
   nvya(1)=sin(elkep(3))*sin(elkep(4)) 
   nvya(2)=-sin(elkep(3))*cos(elkep(4)) 
   nvya(3)=cos(elkep(3)) 


   spya=yarkp(1)*pvya(1)+yarkp(2)*pvya(2)+yarkp(3)*pvya(3) 
   sqya=yarkp(1)*qvya(1)+yarkp(2)*qvya(2)+yarkp(3)*qvya(3) 
   cgam=yarkp(1)*nvya(1)+yarkp(2)*nvya(2)+yarkp(3)*nvya(3) 
   obli=acos(cgam)*invgrad

  ! print *, obli

! 
!         write(*,*)'Incl, Omega, omega (deg) ',elkep(3)/grad,elkep(4)/grad,elkep(5)/grad
!         write(*,*)' Obliquity of the spin axis; Yarkovsky: ',obli 
   ! - compute the \alpha(k) and \beta(k) coefficients                     
   eta=sqrt(1._double-elkep(2)*elkep(2)) 
   etaya75=eta**0.75d0 
   ! -- \beta_1(x) ... \beta_7(x) functions                                
   argu=elkep(2) 
   argu2=argu*argu 
   beya(1)=eta*(1._double+argu2*(-1152._double+argu2*(48._double-argu2))/9216._double) 
   argu=2._double*elkep(2) 
   argu2=argu*argu 
   beya(2)=eta*argu*(1._double+argu2*(-1920._double+argu2*(60._double-argu2))    &
        &           /23040._double)                                             
   argu=3._double*elkep(2) 
   argu2=argu*argu 
   beya(3)=3._double*eta*argu2*(1._double+argu2*(-40._double+argu2)/640._double)/8._double 
   argu=4._double*elkep(2) 
   argu2=argu*argu 
   beya(4)=eta*argu2*argu*(1._double+argu2*(-48._double+argu2)/960._double)/12._double 
   argu=5._double*elkep(2) 
   argu2=argu*argu 
   beya(5)=5._double*eta*argu2*argu2*(1._double-argu2/24._double)/384._double 
   argu=6._double*elkep(2) 
   argu2=argu*argu 
   beya(6)=eta*argu2*argu2*argu*(1._double-argu2/28._double)/640._double 
   argu=7._double*elkep(2) 
   argu2=argu*argu 
   beya(7)=7._double*eta*argu2*argu2*argu2/46080._double 
   ! -- \alpha_1(x) ... \alpha_7(x) functions                              
   argu=elkep(2) 
   argu2=argu*argu 
   alya(1)=1._double+argu2*(-3456._double+argu2*(240._double-7._double*argu2))/9216._double 
   argu=2._double*elkep(2) 
   argu2=argu*argu 
   alya(2)=argu*(1._double+argu2*(-960._double+argu2*(45._double-argu2))/5760._double) 
   argu=3._double*elkep(2) 
   argu2=argu*argu 
   alya(3)=3._double*argu2*(1._double+argu2*(-200._double+7._double*argu2)/1920._double)/  &
        &           8._double                                                   
   argu=4._double*elkep(2) 
   argu2=argu*argu 
   alya(4)=argu*argu2*(1._double+argu2*(-36._double+argu2)/480._double)/12._double 
   argu=5._double*elkep(2) 
   argu2=argu*argu 
   alya(5)=argu2*argu2*(1._double-7._double*argu2/120._double)/76.8d0 
   argu=6._double*elkep(2) 
   argu2=argu*argu 
   alya(6)=argu*argu2*argu2*(1._double-argu2/21._double)/640._double 
   argu=7._double*elkep(2) 
   argu2=argu*argu 
   alya(7)=7._double*argu2*argu2*argu2/46080._double 


   ! - thermal inertia & seasonal thermal parameter                        
   bgama=sqrt(yarkp(4)*yarkp(6)*capacity) 
   theta=bgama*thfacya/emiss 
   seadepth=sqrt(yarkp(4)/yarkp(6)/capacity/fmeaya) 

   ! - radius of the body scaled by the depth of the seasonal wave         
   rp=yarkp(7)/seadepth 
   rp2=dsqrt2*rp 
   tau=theta*etaya75/rp2 
   tau1=1._double+tau 

   ! - amplitude of the effect                                             
   fac=aceuni*yarkp(9)*radfluya/yarkp(7)/densityb/clight3/tau1 
   ! - G_k cos(d_k) & G_K sin(d_k) functions computed                      
   do 10 i=1,napprox 
      fk=i 
      alk=sqrt(fk)*rp2 
      ! - the auxiliary functions A-D, a,b                                    
      cal=cos(alk) 
      sal=sin(alk) 
      if (alk.lt.90._double) then 
         ealm=exp(-alk) 
      else 
         ealm=0._double 
      endif
      af=3._double*(alk+2._double)*ealm+(3._double*(alk-2._double)*cal+alk*(alk-3._double)*sal) 
      bf=alk*(alk+3._double)*ealm+(-alk*(alk-3._double)*cal+3._double*(alk-2._double)*sal) 
      caf=-(alk+2._double)*ealm+(-(alk-2._double)*cal+alk*sal) 
      cbf=-alk*ealm-(alk*cal+(alk-2._double)*sal) 
      ccf=caf+tau*af/tau1 
      cdf=cbf+tau*bf/tau1 
      ! - G exp(i delta)                                                      
      deno=ccf*ccf+cdf*cdf 
      gcosd(i)=(caf*ccf+cbf*cdf)/deno 
      gsind(i)=(cbf*ccf-caf*cdf)/deno 
      ! compute cos- & sin-related brackets                                   
      brac(i)=spya*alya(i)*gcosd(i)+sqya*beya(i)*gsind(i) 
      bras(i)=sqya*beya(i)*gcosd(i)-spya*alya(i)*gsind(i) 
10     continue 


!       ! mean anomaly determined                                               
!       ! - computed from the state vector                                      
!        r2=xast(1)*xast(1)+xast(2)*xast(2)+xast(3)*xast(3) 
!        v2=vast(1)*vast(1)+vast(2)*vast(2)+vast(3)*vast(3) 
!        rdot=xast(1)*vast(1)+xast(2)*vast(2)+xast(3)*vast(3) 
!        r=sqrt(r2) 
!        aaxi=1._double/(2._double/r-(v2/gms)) 
!        esinu=rdot/sqrt(aaxi*gms) 
!        ecosu=(r*v2/gms)-1._double 
!        uano=atan2(esinu,ecosu) 
!        anomaly=uano-esinu 
! !   !   if (anomaly.lt.0._double) anomaly=anomaly+dpig
!       write(*,*) 'M',anomaly,elekep(6)

      !eggl 04.2012 mean anomaly computed already for keplerian elements
       anomaly=elkep(6)
      ! compute the sum...                                                    
      fact=0._double 
      do 100 i=napprox,1,-1 
         fk=i 
         canomaly=cos(fk*anomaly) 
         sanomaly=sin(fk*anomaly) 
         fact=fact+(brac(i)*canomaly+bras(i)*sanomaly) 
100       continue 
         fact=fact*fac 


         ! seasonal acceleration (~ factor * {\bf s})                            
         a(1)=fact*yarkp(1) 
         a(2)=fact*yarkp(2) 
         a(3)=fact*yarkp(3) 
 
         
         ! Partials? -- DISABLED AT THE MOMENT                                   
         !      if (iparti.eq.0) return                                          
         ! - general                                                             
         !      cafp=-ealm+cal+(2._double*al-1._double)*sal                                
         !      cbfp=-ealm-(2._double*al-1._double)*cal+sal                                
         !      afp=3._double*ealm+(al*al-3._double)*cal+(al*(al-4._double)+3._double)*sal           
         !      bfp=(2._double*al+3._double)*ealm-(al*(al-4._double)+3._double)*cal                  
         !     .     +(al*al-3._double)*sal                                           
         ! - thermal conductivity parameters (k_0,k_1)                           
         !      xi1r=caf*ccf-cbf*cdf                                             
         !      xi1i=cbf*ccf+caf*cdf                                             
         !      xi2r=cafp*af-cbfp*bf                                             
         !      xi2i=cbfp*af+cafp*bf                                             
         !      xi2r=xi2r-caf*afp+cbf*bfp                                        
         !      xi2i=xi2i-cbf*afp-caf*bfp                                        
         !      deno=xi1r*xi1r+xi1i*xi1i                                         
         !      facr=1._double+0.5d0*al*(xi2r*xi1r+xi2i*xi1i)/deno                    
         !      faci=     0.5d0*al*(xi2i*xi1r-xi2r*xi1i)/deno                    
         !      derikr=-tau*(gcosd*facr-gsind*faci)/tau1                         
         !      deriki=-tau*(gsind*facr+gcosd*faci)/tau1                         
         !      a(7)=fac*(deriki*vprod1(1)+derikr*vprod2(1))                     
         !      a(8)=fac*(deriki*vprod1(2)+derikr*vprod2(2))                     
         !      a(9)=fac*(deriki*vprod1(3)+derikr*vprod2(3))                     
         !      a(10)=a(7)*(tav1000**3)                                          
         !      a(11)=a(8)*(tav1000**3)                                          
         !      a(12)=a(9)*(tav1000**3)                                          
         ! - radius of the body                                                  
         !      rfac=(tau+tau1)/tau1                                             
         !      a(4)=-a(1)*rfac-2._double*a(7)                                        
         !      a(5)=-a(2)*rfac-2._double*a(8)                                        
         !      a(6)=-a(3)*rfac-2._double*a(9)                                        
         ! - partials d_K (a), d_R (a) ...                                       
         !      a(4)=a(4)/yarkp(7)                                               
         !      a(5)=a(5)/yarkp(7)                                               
         !      a(6)=a(6)/yarkp(7)                                               
         !      a(7)=a(7)/surcon!!!!!! --> yarkp(4)                              
         !      a(8)=a(8)/surcon                                                 
         !      a(9)=a(9)/surcon                                                 
         !      a(10)=a(10)/surcon                                               
         !      a(11)=a(11)/surcon                                               
         !      a(12)=a(12)/surcon                                               
         ! - spin axis components                                                
         return 
END SUBROUTINE yarkse                                          
       ! ==================================================================   



       ! yarkinit: initialisation of the yarkovsky force model for a given aste
       !           written by A. Milani & D. Vokrouhlicky, Oct 99              
!!$SUBROUTINE yarkinit(elem)
!!$  !USE orbit_elements 
!!$  USE TOOLS, only : ORB2CAR
!!$
!!$  IMPLICIT NONE 
!!$  CHARACTER*(*), INTENT(IN) :: astnam 
!!$  TYPE(orbit_elem), INTENT(IN) :: elem
!!$  CHARACTER*80 file 
!!$  real(kind=double) lat,long,emiss,stefboltz,argu,argu2,tstarya,eta 
!!$  TYPE(orbit_elem) :: elekep
!!$  real(kind=double) elkep(6),pvya(3),qvya(3),nvya(3),enne,cgam,obli 
!!$  integer(kind=intdble) unit,le, fail_flag
!!$  INCLUDE 'sysdep.h90' 
!!$! yar is the logical flag for the existence of the physical data        
!!$! allowing computation of Yarkovsky; otherwise, the non gravitational   
!!$! force is set to zero                                                  
!!$!                                                                       
!!$  IF(iyark.eq.0.or.iyark.ge.3)RETURN 
!!$  yarini=.true. 
!!$! convert elements to keplerian                                         
!!$  !call coo_cha(elem,'KEP',elekep,fail_flag)
!!$  !CALL ORB2CAR
!!$  IF(fail_flag.ge.4)THEN
!!$     WRITE(*,*)' yarkinit: not possible with comet ', fail_flag, elekep
!!$     STOP
!!$  ELSE
!!$     elkep=elekep%coord
!!$  ENDIF 
!!$! compute the name of the file which could contain the yarkovsky data   
!!$  CALL filnam(yardir,astnam,'yar',file,le) 
!!$  INQUIRE(file=file(1:le),exist=yarfil) 
!!$  IF(yarfil)THEN 
!!$     call filopn(unit,file(1:le),'old') 
!!$     read(unit,*,end=111) 
!!$! equatorial RA and DEC of the spin axis; beware of old names
!!$     read(unit,*,end=111)long 
!!$     read(unit,*,end=111)lat 
!!$! - via common   yarkp(1-3) ... sx, sy, sz (unit vector of the          
!!$!                               body's spin axis orientation)           
!!$!                yarkp(4-5) ... k_0 and k_1 parameters of the           
!!$!                               surface thermal conductivity            
!!$!                               [K(T) = k_0 + k_1 T_av^3]               
!!$!                yarkp(6) ... density of the surface layer              
!!$!                yarkp(7) ... radius of the body                        
!!$!                yarkp(8) ... rotation frequency                        
!!$!                yarkp(9) ... surface absorptivity
!!$!                yarkp(10)... bulk density                      
!!$     yarkp(1)=cos(lat*radeg)*cos(long*radeg) 
!!$     yarkp(2)=cos(lat*radeg)*sin(long*radeg) 
!!$     yarkp(3)=sin(lat*radeg) 
!!$! rotate to ecliptic
!!$     yarkp(1:3)=MATMUL(roteqec,yarkp(1:3))
!!$     read(unit,*,end=111)yarkp(4) 
!!$     read(unit,*,end=111)yarkp(5) 
!!$     read(unit,*,end=111)yarkp(6) 
!!$     read(unit,*,end=111)yarkp(7) 
!!$     read(unit,*,end=111)yarkp(8) 
!!$     read(unit,*,end=111)yarkp(9) 
!!$     read(unit,*,end=111)yarkp(10) 
!!$! precompute some variables for the seasonal variant of the Yarkovsky   
!!$! effect:                                                               
!!$! - constants                                                           
!!$     emiss=0.9d0 
!!$     stefboltz=5.66962d-8 
!!$! - mean motion & solar radiation flux at r=a                           
!!$     fmeaya=(1.9909837d-7)/elkep(1)/sqrt(elkep(1)) 
!!$     radfluya=1371._double/elkep(1)/elkep(1) 
!!$! - subsolar temperature                                                
!!$     tstarya=(yarkp(9)*radfluya/emiss/stefboltz)**0.25d0 
!!$     thfacya=sqrt(fmeaya)/stefboltz/(tstarya**3) 
!!$! - projections s_P and s_Q of the spin axis computed                   
!!$     pvya(1)=cos(elkep(4))*cos(elkep(5))-                         &
!!$     &           cos(elkep(3))*sin(elkep(4))*sin(elkep(5))           
!!$     pvya(2)=sin(elkep(4))*cos(elkep(5))+                         &
!!$     &           cos(elkep(3))*cos(elkep(4))*sin(elkep(5))           
!!$     pvya(3)=sin(elkep(3))*sin(elkep(5)) 
!!$     qvya(1)=-cos(elkep(4))*sin(elkep(5))-                        &
!!$     &           cos(elkep(3))*sin(elkep(4))*cos(elkep(5))           
!!$     qvya(2)=-sin(elkep(4))*sin(elkep(5))+                        &
!!$     &         cos(elkep(3))*cos(elkep(4))*cos(elkep(5))             
!!$     qvya(3)=sin(elkep(3))*cos(elkep(5)) 
!!$     nvya(1)=sin(elkep(3))*sin(elkep(4)) 
!!$     nvya(2)=-sin(elkep(3))*cos(elkep(4)) 
!!$     nvya(3)=cos(elkep(3)) 
!!$     spya=yarkp(1)*pvya(1)+yarkp(2)*pvya(2)+yarkp(3)*pvya(3) 
!!$     sqya=yarkp(1)*qvya(1)+yarkp(2)*qvya(2)+yarkp(3)*qvya(3) 
!!$     cgam=yarkp(1)*nvya(1)+yarkp(2)*nvya(2)+yarkp(3)*nvya(3) 
!!$     obli=acos(cgam)/radeg 
!!$     write(*,*)'Incl, Omega, omega (deg) ',elkep(3)/radeg,elkep(4)/radeg,elkep(5)/radeg
!!$     write(*,*)' Obliquity of the spin axis; Yarkovsky: ',obli 
!!$! - compute the \alpha(k) and \beta(k) coefficients                     
!!$     eta=sqrt(1._double-elkep(2)*elkep(2)) 
!!$     etaya75=eta**0.75d0 
!!$! -- \beta_1(x) ... \beta_7(x) functions                                
!!$     argu=elkep(2) 
!!$     argu2=argu*argu 
!!$     beya(1)=eta*(1._double+argu2*(-1152._double+argu2*(48._double-argu2))/9216._double) 
!!$     argu=2._double*elkep(2) 
!!$     argu2=argu*argu 
!!$     beya(2)=eta*argu*(1._double+argu2*(-1920._double+argu2*(60._double-argu2))    &
!!$     &           /23040._double)                                             
!!$     argu=3._double*elkep(2) 
!!$     argu2=argu*argu 
!!$     beya(3)=3._double*eta*argu2*(1._double+argu2*(-40._double+argu2)/640._double)/8._double 
!!$     argu=4._double*elkep(2) 
!!$     argu2=argu*argu 
!!$     beya(4)=eta*argu2*argu*(1._double+argu2*(-48._double+argu2)/960._double)/12._double 
!!$     argu=5._double*elkep(2) 
!!$     argu2=argu*argu 
!!$     beya(5)=5._double*eta*argu2*argu2*(1._double-argu2/24._double)/384._double 
!!$     argu=6._double*elkep(2) 
!!$     argu2=argu*argu 
!!$     beya(6)=eta*argu2*argu2*argu*(1._double-argu2/28._double)/640._double 
!!$     argu=7._double*elkep(2) 
!!$     argu2=argu*argu 
!!$     beya(7)=7._double*eta*argu2*argu2*argu2/46080._double 
!!$! -- \alpha_1(x) ... \alpha_7(x) functions                              
!!$     argu=elkep(2) 
!!$     argu2=argu*argu 
!!$     alya(1)=1._double+argu2*(-3456._double+argu2*(240._double-7._double*argu2))/9216._double 
!!$     argu=2._double*elkep(2) 
!!$     argu2=argu*argu 
!!$     alya(2)=argu*(1._double+argu2*(-960._double+argu2*(45._double-argu2))/5760._double) 
!!$     argu=3._double*elkep(2) 
!!$     argu2=argu*argu 
!!$     alya(3)=3._double*argu2*(1._double+argu2*(-200._double+7._double*argu2)/1920._double)/  &
!!$     &           8._double                                                   
!!$     argu=4._double*elkep(2) 
!!$     argu2=argu*argu 
!!$     alya(4)=argu*argu2*(1._double+argu2*(-36._double+argu2)/480._double)/12._double 
!!$     argu=5._double*elkep(2) 
!!$     argu2=argu*argu 
!!$     alya(5)=argu2*argu2*(1._double-7._double*argu2/120._double)/76.8d0 
!!$     argu=6._double*elkep(2) 
!!$     argu2=argu*argu 
!!$     alya(6)=argu*argu2*argu2*(1._double-argu2/21._double)/640._double 
!!$     argu=7._double*elkep(2) 
!!$     argu2=argu*argu 
!!$     alya(7)=7._double*argu2*argu2*argu2/46080._double 
!!$! close the input file                                                  
!!$     call filclo(unit,' ') 
!!$  ELSE 
!!$     WRITE(*,*)' Yarkovsky datafile not found:',file(1:le) 
!!$     stop 
!!$  ENDIF
!!$  WRITE(*,*)' Yarkovsky data loaded for asteroid ', astnam 
!!$  RETURN 
!!$111 yarfil=.false. 
!!$  WRITE(*,*)' incomplete yarkovsky file for asteroid ', astnam 
!!$END SUBROUTINE yarkinit

END MODULE yark_m

!Tutorial for yarkovsky subroutine


! N.B. Some other external modules are needed (see at the top of the program) like- YARKO which contains the yarkp parameters and an other (that you don't need) that is the aster_obliq (obliquity of the spin)
! - CONST which contains some basic constant like invgrad = conversion from radian to degre, WP = kind(0d0), speed_light =  2.99792458d8 (m/s),UA = 1.49597870691d11 m 
! 
! The yarkovsky parameters are global parameters named yarkp(dim=10):

! *yarkp(1-3) = (sx,sy,sz) cartesian coordinates of the spin vectors 
! *yarkp(4-5) = k_0 and k_1 parameters of the surface thermal conductivity [K(T) = k_0 + k_1 T^3] 
! *yarkp(6) = density of the surface layer
! *yarkp(7) = radius of the body                        
! *yarkp(8) = rotation frequency                        
! *yarkp(9) = surface absorptivity   
! *yarkp(10)=  bulk density 


 
