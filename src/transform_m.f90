module transform_m
  use global_m
  implicit none

  public::htrnsko
  public::htrnsel
  public::btrnsko
 
  public::trojtrnsko
  public::trojtrnsel
  public::hekoo
  public::bakoo

  public::binbakoo
  public::binhekoo
  contains

!/////////////////////////////////////////////////////////////////////////////////
!                              Transformationen fuer Eingabe und Ausgabe 
!
!                                         by Siegfried Eggl  20080110 
!--------------------------------------------------------------------------------
!*****************************************************************

!....BERECHNUNG DER HELIOZENTRISCHEN KOORDINATEN UND GESCHWINDIGKEITEN
!....AUS DEN BAHNELEMENTEN


        SUBROUTINE htrnsko(n,ele,mass,rv)
         use global_m
         implicit none   
         integer(kind=intdble)::i,n
         real(kind=double),dimension(1:n,1:6),intent(out)::rv
          type(kepele),dimension(1:n)::ele
          real(kind=double)::xh1,xh2,xh3,vh1,vh2,vh3,a,e,incl,om,gom,m       
          real(kind=double)::cosgom,singom,cosi,sini,px,py,pz,qx,qy,qz,ah,exan,cappaq
          real(kind=double)::cosex,sinex,rh,bh,cosom,sinom,ea,mass(1:n),grd,p,kk
  
         p=pi
         grd=grad
         kk=k

!$omp parallel do default(private)   &
!$omp firstprivate(kk,grd,p,n) &
!$omp shared(rv,mass,ele) 
    do  i=2,n
          a=ele(i)%a
          e=ele(i)%e
          incl=ele(i)%i*grd
          om=ele(i)%kom*grd
          gom=ele(i)%gom*grd
          m=ele(i)%man*grd

          cappaq=mass(1)+mass(i)

          call nr(m,e,ea)
 !        cosea=cos(ea)
 !         cosphi=(cosea-e)/(1._double-e*cosea)
         
          IF (E .EQ.0.D0)  then
             OM  = 0.D0
          end if
          IF (INCL .EQ.0.D0) then
             GOM  = 0.D0
          end if
          COSOM = COS(OM )
          SINOM = SIN(OM )
          COSGOM = COS(GOM)
          SINGOM = SIN(GOM )
          COSI = COS(INCL )
          SINI = SIN(INCL ) 
          PX= COSOM*COSGOM-SINOM*SINGOM*COSI
          PY= COSOM*SINGOM+SINOM*COSGOM*COSI
          PZ= SINOM*SINI
          QX=-SINOM*COSGOM-COSOM*SINGOM*COSI
          QY=-SINOM*SINGOM+COSOM*COSGOM*COSI
          QZ= COSOM*SINI
!          AH = SQRT((1-E )/(1+E ))*SQRT((1._double-cosphi(n))/(1._double+cosphi(n)))
!frueher... *DTAN(V /2._double)
!          EXAN = 2.D0*DATAN(AH)
          EXAN=ea
 !         ta=2._double*DATAN(SQRT((1._double+e)/(1._double-e))*DTAN(ea/2._double))
 
          COSEX = COS(EXAN)
          SINEX =SIN(EXAN)

          AH = A *SQRT(1._double-E *E )
          XH1  = A *PX*(COSEX-E )+AH*QX*SINEX
          XH2  = A *PY*(COSEX-E )+AH*QY*SINEX
          XH3  = A *PZ*(COSEX-E )+AH*QZ*SINEX
          RH = SQRT(XH1 *XH1 +XH2 *XH2 +XH3 *XH3 )

          if(A.eq.0._double.or.RH.eq.0._double) then
             AH=1._double
          else
             AH = SQRT(cappaq )/(SQRT(A )*RH)
          end if
  
          BH = A *SQRT(1._double-E *E )*COSEX
          VH1  = AH*(-A *PX*SINEX+BH*QX)
          VH2  = AH*(-A *PY*SINEX+BH*QY)
          VH3  = AH*(-A *PZ*SINEX+BH*QZ)

   
          rv(i,1)=  XH1
          rv(i,2)=  XH2
          rv(i,3)=  XH3
          rv(i,4)=VH1*kk
          rv(i,5)=VH2*kk
          rv(i,6)=VH3*kk                   
end do
  !$omp end parallel do  

        RETURN
      END  subroutine htrnsko


!*****************************************
!...Calculates Keplerian Orbital Elements from heliocentric position and velocity vectors

SUBROUTINE htrnsel (n,rv,mass,ele)
use global_m
implicit none
integer(kind=intdble)::i,j,n,dim      
real(kind=double):: ele(1:n,1:6),rv(1:n,1:6),mass(1:n),e2,eanom,Lprojxy,U(1:3,1:3)
real(kind=double):: tanom,Hnorm,sinE,cosE
real(kind=double),dimension(1:3)::LRL,L,r,v,H,node,rU
real(kind=double)::cappaq,a,e,incl,om,gom,mmm,L2,rnorm,nnorm,igrad,p
 
p=pi
igrad=invgrad
       
dim=3
!$omp parallel do default(private)   &
!$omp firstprivate (igrad,p,n,dim)  &
!$omp shared(rv,mass,ele)
 do i=2,n    
       cappaq=mass(1)+mass(i)

!*******************************************************
! LETS BE SERIOUS HERE

          r=rv(i,1:3)
          v=rv(i,4:6)
!calculation of unit mass angular momentum vector L and its projection onto xy-plane
        call crossp3d(r,v,L)
        L2=Dot_Product(L,L)
        Lprojxy=SQRT( L(1)*L(1)+L(2)*L(2) )

 ! inclination "incl"
       incl=ATAN2(Lprojxy,L(3))*igrad

!argument of the ascending node "gom"
       gom=Atan2(L(1),-L(2))*igrad
       if (gom.lt.0._double) then
           gom = gom+360._double
        end if

! Laplace-Runge-Lenz vector
        call crossp3d(v,L,LRL)
        rnorm=Sqrt(Dot_Product(r(:),r(:)))

        LRL(:)=LRL(:)/cappaq-r(:)/rnorm
         
! eccentricity "e" = |LRL|
        e=Sqrt(Dot_Product(LRL,LRL))

        if(e.lt.1.d-13)then 
           e=0._double
           e2=0._double
           om=0._double
         else
           e2=e*e
         end if
         
! semi-major axis "a"
     ! a=L2/(cappaq*(1._double+e2))
       a=L2/(cappaq*(1._double-e2))    
   
   !calculation of vector "node" pointing to ascending node    
    ! attention: will change node for retrograde orbits!!!!
         node(1)=-L(2)*L(3)
         node(2)=L(1)*L(3)
         node(3)=0._double

         nnorm=Sqrt(Dot_Product(node,node))
 
       if(incl.lt.1.d-13) then
          gom=0._double
         ! argument of pericenter "om" = angle of LRL to x-axis
          om=Acos(LRL(1)/Sqrt(Dot_Product(LRL,LRL)))
              if(LRL(2).lt.0._double) then
                 om=360._double-om*igrad
              else
                 om=om*igrad
              end if
  
        else
         ! argument of pericenter "om" = angle of LRL to node
    !      om=Acos(Dot_Product(LRL(1:2),node(1:2))/ &
    !         (Sqrt(Dot_Product(LRL,LRL))*nnorm))

         !calculation of H in order to fix quadrant of LRL 
    !      call crossp3d(L,LRL,H)
    !          if((node(1)*H(1)+node(2)*H(2)).gt.0._double) then
    !             om=360._double-om*grad
    !          else
    !             om=om*grad
    !          end if 

          call crossp3d(L,node,H)
          Hnorm=Sqrt(Dot_Product(H,H))
          om=Atan2(Dot_Product(LRL,H)*nnorm,Dot_Product(LRL,node)*Hnorm)*igrad
              if(om.lt.0._double) om=om+360._double 
  
        end if  

          if(e.lt.1.2d-13.and.incl.le.1.d-13) then
              gom=0._double
              om=0._double

              !mean anomaly 'mmm'
              mmm=Atan2(r(2),r(1))*igrad
             if(mmm.lt.0._double) mmm=mmm+360._double

        

           elseif (e.lt.1.2d-13.and.incl.gt.1.d-13) then
             
              !calculation of true anomaly via vector "H" at right angle to LRL and L (alternative)
              !        call crossp3d(L,LRL,H)
              !        Hnorm=Sqrt(Dot_Product(H,H))
              !        tanom= Atan2(Dot_Product(H,r)*e,Dot_Product(LRL,r)*Hnorm)   


              !calculation of true anomaly via vector "node" pointing to ascending node DOES NOT WORK! 
            !  help1=r(1)*node(1)+r(2)*node(2)
             ! tanom=Acos(help1/(rnorm*nnorm))
             
          !    tanom=Acos(Dot_Product(r,node)/(rnorm*nnorm))
          !    if(node(1)*v(1)+node(2)*v(2).gt.0._double) then
          !       tanom=2._double*PI-tanom
          !    end if

              !transformation of positionvector r into coordinate-system U spanned by Angular Momentum vector, vector Sun-ascending Node, 
              !and their crossproduct H, corresponding to motion of rU in the orbital plane. 
              

              call crossp3d(L,node,H)
              Hnorm=Sqrt(Dot_Product(H,H))
              ! calculation of the base - transformation Matrix
              U(:,1)=node(:)/nnorm
              U(:,2)=H/Hnorm
              U(:,3)=L(:)/Sqrt(L2)
              
              !coordinate transformation
              call gauss(dim,U,r,rU)

              !true anomaly in orbital plane
              tanom=Atan2(rU(2),rU(1))

              !mean anomaly 'mmm'
              mmm=tanom*igrad     
              if(mmm.lt.0._double) mmm=mmm+360._double
      
   !  ??       om=0._double

            elseif (incl.lt.1.d-13.and.e.gt.1.d-13) then
              !calculation of true anomaly via vector "H" at right angle to LRL and L (alternative)
                      call crossp3d(L,LRL,H)
                      Hnorm=Sqrt(Dot_Product(H,H))
                      tanom= Atan2(Dot_Product(H,r)*e,Dot_Product(LRL,r)*Hnorm) 
  

!!$             call crossp3d(L,node,H)
!!$              Hnorm=Sqrt(Dot_Product(H,H))
!!$              ! calculation of the base - transformation Matrix
!!$              U(:,1)=node(:)/nnorm
!!$              U(:,2)=H/Hnorm
!!$              U(:,3)=L(:)/Sqrt(L2)
!!$              
!!$              !coordinate transformation
!!$              call gauss(3,U,r,rU)
!!$
!!$              !true anomaly in orbital plane
!!$              tanom=Atan2(rU(2),rU(1))
           

              !calculation of eccentric anomaly eanom (=E)
               cosE=(e+Cos(tanom))/(1._double+e*Cos(tanom))
               sinE=Sqrt(1._double-e2)*Sin(tanom)/(1._double+e*Cos(tanom))
               eanom=ATAN2(sinE,cosE)  

              !mean anomaly 'mmm' via Kepler's equation
              mmm=(eanom-e*sinE)*igrad              
 
!                write(*,*)'hey'

                 if(mmm.lt.0._double) mmm=mmm+360._double

              else
              !calculation of true anomaly via vector "node" pointing to ascending node  NOT WORKING!
           !   help1=r(1)*node(1)+r(2)*node(2)
           !   tanom=Acos(help1/(rnorm*nnorm))
             
           !   if(node(1)*v(1)+node(2)*v(2).gt.0._double) then
           !      tanom=2._double*PI-tanom
           !   end if
           !  

              !calculation of true anomaly via vector "H" at right angle to LRL and L (alternative)
          !            call crossp3d(L,LRL,H)
          !            Hnorm=Sqrt(Dot_Product(H,H))
          !            tanom= Atan2(Dot_Product(H,r)*e,Dot_Product(LRL,r)*Hnorm) 

 !   write(*,*)'whathe'
  
             call crossp3d(L,LRL,H)
              Hnorm=Sqrt(Dot_Product(H,H))
              ! calculation of the base - transformation Matrix
              U(:,1)=LRL(:)/e
              U(:,2)=H/Hnorm
              U(:,3)=L(:)/Sqrt(L2)
              
              !coordinate transformation
              call gauss(dim,U,r,rU)

              !true anomaly in orbital plane
              tanom=Atan2(rU(2),rU(1))


              !calculation of eccentric anomaly eanom (=E)
               cosE=(e+Cos(tanom))/(1._double+e*Cos(tanom))
               sinE=Sqrt(1._double-e2)*Sin(tanom)/(1._double+e*Cos(tanom))
               eanom=Atan2(sinE,cosE)  
               
                
               
              !mean anomaly 'mmm' via Kepler's equation
             mmm=(eanom-e*sinE)*igrad
             if(mmm.lt.0._double) mmm=mmm+360._double
           end if
      
          if(e.lt.0._double.or.e.gt.1._double) e=1._double
          if(a.lt.0._double.or.a.gt.huge(a)) a=0._double 
          if(incl.gt.180._double.or.incl.lt.0._double) om=-100._double

          ele(i,1)=a
          ele(i,2)=e
          ele(i,3)=incl
          ele(i,4)=om
          ele(i,5)=gom
          ele(i,6)=mmm

          do j=4,6
            if(ele(i,j).gt.360._double.or.ele(i,j).lt.0._double.or.ele(i,j).ne.ele(i,j)) then
              ele(i,j)=-100._double
            end if
          end do

       end do
!$omp end parallel do
     

     
       RETURN
     END subroutine htrnsel
!*************************************************************************************************
!...CALCULATION OF BINARY ORBITAL ELEMENTS FROM POSITION- AND VELOCITY-VEKTORS
!....INCLUDING QUADRANT CORRECTION 
SUBROUTINE btrnsel (n,rv,mass,ele)
use global_m
implicit none
integer(kind=intdble)::i,j,n,dim      
real(kind=double):: ele(1:n,1:6),rv(1:n,1:6),mass(1:n),e2,eanom,Lprojxy,U(1:3,1:3)
real(kind=double):: tanom,Hnorm,sinE,cosE
real(kind=double),dimension(1:3)::LRL,L,r,v,H,node,rU
real(kind=double)::cappaq,a,e,incl,om,gom,mmm,L2,igrad,rnorm,nnorm
 
igrad=invgrad       
dim=3
!$omp parallel do default(private)   &
!$omp firstprivate(n,igrad,dim) &
!$omp shared(rv,mass,ele)
 do i=2,n    

       if(i.eq.2) then
       cappaq=mass(1)+mass(i)
         else
       cappaq=mass(1)+mass(2)+mass(i)
       end if
!*******************************************************
! LETS BE SERIOUS HERE

          r=rv(i,1:3)
          v=rv(i,4:6)
!calculation of unit mass angular momentum vector L and its projection onto xy-plane
        call crossp3d(r,v,L)
        L2=Dot_Product(L,L)
        Lprojxy=SQRT( L(1)*L(1)+L(2)*L(2) )

 ! inclination "incl"
       incl=ATAN2(Lprojxy,L(3))*igrad

!argument of the ascending node "gom"
       gom=Atan2(L(1),-L(2))*igrad
       if (gom.lt.0._double) then
           gom = gom+360._double
        end if

! Laplace-Runge-Lenz vector
        call crossp3d(v,L,LRL)
        rnorm=Sqrt(Dot_Product(r(:),r(:)))

        LRL(:)=LRL(:)/cappaq-r(:)/rnorm
         
! eccentricity "e" = |LRL|
        e=Sqrt(Dot_Product(LRL,LRL))

        if(e.lt.1.d-13)then 
           e=0._double
           e2=0._double
           om=0._double
         else
           e2=e*e
         end if
         
! semi-major axis "a"
     ! a=L2/(cappaq*(1._double+e2))
       a=L2/(cappaq*(1._double-e2))    
   
   !calculation of vector "node" pointing to ascending node    
    ! attention: will change node for retrograde orbits!!!!
         node(1)=-L(2)*L(3)
         node(2)=L(1)*L(3)
         node(3)=0._double

         nnorm=Sqrt(Dot_Product(node,node))
 
       if(incl.lt.1.d-13) then
          gom=0._double
         ! argument of pericenter "om" = angle of LRL to x-axis
          om=Acos(LRL(1)/Sqrt(Dot_Product(LRL,LRL)))
              if(LRL(2).lt.0._double) then
                 om=360._double-om*igrad
              else
                 om=om*igrad
              end if
  
        else
          call crossp3d(L,node,H)
          Hnorm=Sqrt(Dot_Product(H,H))
          om=Atan2(Dot_Product(LRL,H)*nnorm,Dot_Product(LRL,node)*Hnorm)*igrad
              if(om.lt.0._double) om=om+360._double 
  
        end if  

          if(e.lt.1.2d-13.and.incl.le.1.d-13) then
              gom=0._double
              om=0._double

              !mean anomaly 'mmm'
              mmm=Atan2(r(2),r(1))*igrad
             if(mmm.lt.0._double) mmm=mmm+360._double

           elseif (e.lt.1.2d-13.and.incl.gt.1.d-13) then
             

              call crossp3d(L,node,H)
              Hnorm=Sqrt(Dot_Product(H,H))
              ! calculation of the base - transformation Matrix
              U(:,1)=node(:)/nnorm
              U(:,2)=H/Hnorm
              U(:,3)=L(:)/Sqrt(L2)
              
              !coordinate transformation
              call gauss(dim,U,r,rU)

              !true anomaly in orbital plane
              tanom=Atan2(rU(2),rU(1))

              !mean anomaly 'mmm'
              mmm=tanom*igrad     
              if(mmm.lt.0._double)  mmm=mmm+360._double

             elseif (incl.lt.1.d-13.and.e.gt.1.d-13) then
              !calculation of true anomaly via vector "H" at right angle to LRL and L (alternative)
                      call crossp3d(L,LRL,H)
                      Hnorm=Sqrt(Dot_Product(H,H))
                      tanom= Atan2(Dot_Product(H,r)*e,Dot_Product(LRL,r)*Hnorm) 
 
              !calculation of eccentric anomaly eanom (=E)
               cosE=(e+Cos(tanom))/(1._double+e*Cos(tanom))
               sinE=Sqrt(1._double-e2)*Sin(tanom)/(1._double+e*Cos(tanom))
               eanom=ATAN2(sinE,cosE)  

              !mean anomaly 'mmm' via Kepler's equation
              mmm=(eanom-e*sinE)*igrad              
                 if(mmm.lt.0._double)  mmm=mmm+360._double

             else
   
             call crossp3d(L,LRL,H)
              Hnorm=Sqrt(Dot_Product(H,H))
              ! calculation of the base - transformation Matrix
              U(:,1)=LRL(:)/e
              U(:,2)=H/Hnorm
              U(:,3)=L(:)/Sqrt(L2)
              
              !coordinate transformation
              call gauss(dim,U,r,rU)

              !true anomaly in orbital plane
              tanom=Atan2(rU(2),rU(1))

              !calculation of eccentric anomaly eanom (=E)
               cosE=(e+Cos(tanom))/(1._double+e*Cos(tanom))
               sinE=Sqrt(1._double-e2)*Sin(tanom)/(1._double+e*Cos(tanom))
               eanom=Atan2(sinE,cosE)  
                            
              !mean anomaly 'mmm' via Kepler's equation
             mmm=(eanom-e*sinE)*igrad
             if(mmm.lt.0._double) mmm=mmm+360._double
           end if
      
          if(e.lt.0._double.or.e.gt.1._double) e=1._double
          if(a.lt.0._double.or.a.gt.huge(a)) a=0._double 
          if(incl.gt.180._double.or.incl.lt.0._double) om=-100._double

          ele(i,1)=a
          ele(i,2)=e
          ele(i,3)=incl
          ele(i,4)=om
          ele(i,5)=gom
          ele(i,6)=mmm

          do j=4,6
            if(ele(i,j).gt.360._double.or.ele(i,j).lt.0._double.or.ele(i,j).ne.ele(i,j)) then
              ele(i,j)=-100._double
            end if
          end do

       end do
!$omp end parallel do
     
       RETURN
     END subroutine btrnsel
!*********************************************************
subroutine gauss(dim,a,b,c)
!----------------------------------------------------------
! Gauß'sches Eliminationsverfahren zur Lösung eines 
! Systems Linearer Gleichungen der Form 
!                                    A.c=b
! mit Teilpivotisierung
!
! dim....[Int] Dimension des Gleichungssystems
! a...[Real] dim x dim Koeffizientenmatrix des Gleichungssystems
! b...[Real] dim Vektor der rechten Seite der Gleichung
! c...[Real] dim  Vektor der Unbekannten
! dependencies: none
!----------------------------------------------------------
implicit none
  integer(kind=intdble),intent(in)::dim
  integer(kind=intdble)::i,j,k,pivot(1:dim),amax(1:1),idum
  real(kind=double)::a(1:dim,1:dim+1),b(1:dim),dum(1:dim+1)
  real(kind=double)::c(1:dim)

  do i=1,dim
     pivot(i)=i
  end do

  dum(:)=0._double
  c(:)=0._double
  a(:,dim+1)=b(:)
 
 do i=1,dim
!-----------------------------------------------------------------------
!  Partial Pivoting (nur Zeilen werden vertauscht)
!-----------------------------------------------------------------------
    amax(:)=maxloc(abs(a(:,i))) !Welche Zeile hat den größten Wert?
    idum=amax(1)
   if (idum.gt.i) then
      pivot(i)=idum
       dum(:)=a(i,:)                                !Sichere die Zeile die ersetzt wird
      a(i,:)=a(pivot(i),:)                       !Ersetze Zeile i mit jener mit größtem Wert 
       a(pivot(i),:)=dum(:)                 !Beende die Vertauschung, der gesicherte Wert ersetzt die ersetzende Zeile
     end if
!------------------------------------------------------
!   Gauss Elimination
!-----------------------------------------------------     
      a(i,:)=a(i,:)/a(i,i)
      if(i+1.le.dim) then
         do j=i+1,dim
            a(j,:)=a(j,:)-a(j,i)*a(i,:)         
        end do  
      end if
  end do
  
!-----------------------------------------------
!   Rückwärtssubstitution
!-----------------------------------------------

do i=0,dim-1
   j=dim-i  
   dum(1)=0._double
   if (j+1.le.dim) then
      do k=j+1,dim
         dum(1)=dum(1)+a(j,k)*c(k)
      end do
   end if
   c(j)=(a(j,dim+1)-dum(1))
end do

!-----------------------------------------------------------------------
!  Back Pivoting enfällt, da nur Zeilen vertauscht wurden, und sich
!  somit zwar die Reihenfolge der Gleichungen, nicht aber die Koeffizienten
!  der Variablen geändert haben.
!-----------------------------------------------------------------------

return
end subroutine gauss

!***********************************************************************

!....CALCULATES POSITION AND VELOCITIES FROM BINARY ORBITAL ELEMENTS
!....ORBITAL ELEMENTS OF SECONDARY ARE WITH RESPECT TO PRIMARY
!....THE REST WITH RESPECT TO THE BINARIES' BARYCENTER

        SUBROUTINE btrnsko(n,ele,mass,rv)
          use global_m
         implicit none   
         integer(kind=intdble)::i,n
         real(kind=double),dimension(1:n,1:6),intent(out)::rv
          type(kepele),dimension(1:n)::ele
          real(kind=double)::xh1,xh2,xh3,vh1,vh2,vh3,a,e,incl,om,gom,m       
          real(kind=double)::cosgom,singom,cosi,sini,px,py,pz,qx,qy,qz,ah,exan,cappaq
          real(kind=double)::cosex,sinex,rh,bh,cosom,sinom,ea,mass(1:n),summass,grd,p

           grd=grad
            p=pi
           summass=mass(1)+mass(2)
        
        
!$omp parallel default(private)   &
!$omp firstprivate(grd,p,summass,n) &
!$omp shared(rv,mass,ele)
!$omp do
    do  i=2,n

          a=ele(i)%a
          e=ele(i)%e
          incl=ele(i)%i*grd
          om=ele(i)%kom*grd
          gom=ele(i)%gom*grd
          m=ele(i)%man*grd

          if(i.eq.2) then
              cappaq=mass(1)+mass(i)
          else
              cappaq=summass+mass(i)
          end if
          call nr(m,e,ea)
 !        cosea=cos(ea)
 !         cosphi=(cosea-e)/(1._double-e*cosea)
         
          IF (E .EQ.0.D0)  then
             OM  = 0.D0
          end if
          IF (INCL .EQ.0.D0) then
             GOM  = 0.D0
          end if
          COSOM = COS(OM )
          SINOM = SIN(OM )
          COSGOM = COS(GOM )
          SINGOM = SIN(GOM )
          COSI = COS(INCL )
          SINI = SIN(INCL ) 
          PX= COSOM*COSGOM-SINOM*SINGOM*COSI
          PY= COSOM*SINGOM+SINOM*COSGOM*COSI
          PZ= SINOM*SINI
          QX=-SINOM*COSGOM-COSOM*SINGOM*COSI
          QY=-SINOM*SINGOM+COSOM*COSGOM*COSI
          QZ= COSOM*SINI
!          AH = SQRT((1-E )/(1+E ))*SQRT((1._double-cosphi(n))/(1._double+cosphi(n)))
!frueher... *DTAN(V /2._double)
!          EXAN = 2.D0*DATAN(AH)
          EXAN=ea
 !         ta=2._double*DATAN(SQRT((1._double+e)/(1._double-e))*DTAN(ea/2._double))
 
          COSEX = COS(EXAN)
          SINEX =SIN(EXAN)

          AH = A *SQRT(1._double-E *E )
          XH1  = A *PX*(COSEX-E )+AH*QX*SINEX
          XH2  = A *PY*(COSEX-E )+AH*QY*SINEX
          XH3  = A *PZ*(COSEX-E )+AH*QZ*SINEX
          RH = SQRT(XH1 *XH1 +XH2 *XH2 +XH3 *XH3 )

          if(A.eq.0._double.or.RH.eq.0._double) then
             AH=1._double
          else
             AH = SQRT(cappaq )/(SQRT(A )*RH)
          end if
  
          BH = A *SQRT(1._double-E *E )*COSEX
          VH1  = AH*(-A *PX*SINEX+BH*QX)
          VH2  = AH*(-A *PY*SINEX+BH*QY)
          VH3  = AH*(-A *PZ*SINEX+BH*QZ)

   
          rv(i,1)=  XH1
          rv(i,2)=  XH2
          rv(i,3)=  XH3
          rv(i,4)=VH1*k
          rv(i,5)=VH2*k
          rv(i,6)=VH3*k                   
end do
!$omp end do
!$omp end parallel
   

        RETURN
      END  subroutine btrnsko



!******************************************************************************************
!....BERECHNUNG DER HELIOZENTRISCHEN KOORDINATEN UND GESCHWINDIGKEITEN
!....AUS DEN BAHNELEMENTEN FUER TROJANER
        SUBROUTINE trojtrnsko(n,ele,mass,rv)
        use global_m
         implicit none   
         integer(kind=intdble)::i,n
         real(kind=double),dimension(1:n,1:6),intent(out)::rv
          type(kepele),dimension(1:n)::ele
 
          real(kind=double)::xh1,xh2,xh3,vh1,vh2,vh3,a,e,incl,om,gom,m       
          real(kind=double)::cosgom,singom,cosi,sini,px,py,pz,qx,qy,qz,ah,exan,cappaq
          real(kind=double)::cosex,sinex,rh,bh,cosom,sinom,ea,mass(1:n),p
          real(kind=double)::avv(3)
    
     p=pi  

    do  i=2,n

          a=ele(i)%a
          e=ele(i)%e
          incl=ele(i)%i/180._double*p
          om=ele(i)%kom/180._double*p
          gom=ele(i)%gom/180._double*p
          m=ele(i)%man/180._double*p

!trojan defined in global_m

          if (trojan(i).le.1) then
            cappaq=mass(1)+mass(i)
          else
            cappaq=mass(1)+mass(trojan(i))
          end if

          call nr(m,e,ea)
 !        cosea=cos(ea)
 !         cosphi=(cosea-e)/(1._double-e*cosea)
         
          IF (E .EQ.0.D0)  then
             OM  = 0.D0
          end if
          IF (INCL .EQ.0.D0) then
             GOM  = 0.D0
          end if
          COSOM = COS(OM )
          SINOM = SIN(OM )
          COSGOM = COS(GOM )
          SINGOM = SIN(GOM )
          COSI = COS(INCL )
          SINI = SIN(INCL ) 
          PX= COSOM*COSGOM-SINOM*SINGOM*COSI
          PY= COSOM*SINGOM+SINOM*COSGOM*COSI
          PZ= SINOM*SINI
          QX=-SINOM*COSGOM-COSOM*SINGOM*COSI
          QY=-SINOM*SINGOM+COSOM*COSGOM*COSI
          QZ= COSOM*SINI
!          AH = SQRT((1-E )/(1+E ))*SQRT((1._double-cosphi(n))/(1._double+cosphi(n)))
!frueher... *DTAN(V /2._double)
!          EXAN = 2.D0*DATAN(AH)
          EXAN=ea
 !         ta=2._double*DATAN(SQRT((1._double+e)/(1._double-e))*DTAN(ea/2._double))
 
          COSEX = COS(EXAN)
          SINEX =SIN(EXAN)

          AH = A *SQRT(1._double-E *E )
          XH1  = A *PX*(COSEX-E )+AH*QX*SINEX
          XH2  = A *PY*(COSEX-E )+AH*QY*SINEX
          XH3  = A *PZ*(COSEX-E )+AH*QZ*SINEX
          RH = SQRT(XH1 *XH1 +XH2 *XH2 +XH3 *XH3 )
    
          
          if(A.eq.0._double.or.RH.eq.0._double) then
             AH=1._double
          else
             AH = SQRT(cappaq )/(SQRT(A )*RH)
          end if
           
  
          BH = A *SQRT(1._double-E *E )*COSEX
          VH1  = AH*(-A *PX*SINEX+BH*QX)
          VH2  = AH*(-A *PY*SINEX+BH*QY)
          VH3  = AH*(-A *PZ*SINEX+BH*QZ)

   
          rv(i,1)=  XH1
          rv(i,2)=  XH2
          rv(i,3)=  XH3

          rv(i,4)=VH1*k
          rv(i,5)=VH2*k
          rv(i,6)=VH3*k                   


!         now rescale velocities to make them have the same angular velocity as the trojan host
          ! dr/dt = avv x r  so it should rest in a corotating frame
          
      if(trojan(i).le.1) then
      else    
          if(ele(trojan(i))%a.eq.a) then
          else
            call crossp3d(rv(trojan(i),1:3),rv(trojan(i),4:6),avv)
            avv=avv/dot_product(rv(trojan(i),1:3),rv(trojan(i),1:3))
  
           call crossp3d(avv,rv(i,1:3),rv(i,4:6))
          end if
     end if
end do
        RETURN
      END  subroutine trojtrnsko

!*****************************************
!...BERECHNUNG DER BAHNELEMENTE AUS ORTS- UND GESCHWINDIGKEITS-VEKTOREN
!....INCLUSIVE QUADRANTENKORREKTUR 

SUBROUTINE trojtrnsel (n,rv,mass,ele)
use global_m
implicit none
integer(kind=intdble)::i,j,n,dim  
real(kind=double):: ele(1:n,1:6),rv(1:n,1:6),mass(1:n),e2,eanom,Lprojxy,U(1:3,1:3)
real(kind=double):: tanom,Hnorm,sinE,cosE
real(kind=double),dimension(1:3)::LRL,L,r,v,H,node,rU
real(kind=double)::cappaq,a,e,incl,om,gom,mmm,L2,igrad,rnorm,nnorm,p
!real(kind=double)::avv
       
p=pi
igrad=invgrad 
dim=3

!$omp parallel default(private)   &
!$omp firstprivate(p,n,igrad) &
!$omp shared(rv,mass,ele)
!$omp do
 do i=2,n    

          r=rv(i,1:3)
         
         if (trojan(i).le.1) then
            cappaq=mass(1)+mass(i)
            
          else
            cappaq=mass(1)+mass(trojan(i))
          end if

          v=rv(i,4:6)
          
          
!calculation of unit mass angular momentum vector L and its projection onto xy-plane
        call crossp3d(r,v,L)
        L2=Dot_Product(L,L)
       Lprojxy=SQRT( L(1)*L(1)+L(2)*L(2) )

 ! inclination "incl"
       incl=ATAN2(Lprojxy,L(3))*igrad

!argument of the ascending node "gom"
       gom=Atan2(L(1),-L(2))*igrad
       if (gom.lt.0._double) then
           gom = gom+360._double
        end if

! Laplace-Runge-Lenz vector
        call crossp3d(v,L,LRL)
        rnorm=Sqrt(Dot_Product(r(:),r(:)))

        LRL(:)=LRL(:)/cappaq-r(:)/rnorm
         
! eccentricity "e" = |LRL|
        e=Sqrt(Dot_Product(LRL,LRL))

        if(e.lt.1.d-13)then 
           e=0._double
           e2=0._double
           om=0._double
         else
           e2=e*e
         end if
         
! semi-major axis "a"
     ! a=L2/(cappaq*(1._double+e2))
       a=L2/(cappaq*(1._double-e2))    
   
   !calculation of vector "node" pointing to ascending node    
    ! attention: will change node for retrograde orbits!!!!
         node(1)=-L(2)*L(3)
         node(2)=L(1)*L(3)
         node(3)=0._double

         nnorm=Sqrt(Dot_Product(node,node))
 
       if(incl.lt.1.d-13) then
          gom=0._double
         ! argument of pericenter "om" = angle of LRL to x-axis
          om=Acos(LRL(1)/Sqrt(Dot_Product(LRL,LRL)))
              if(LRL(2).lt.0._double) then
                 om=360._double-om*igrad
              else
                 om=om*igrad
              end if
  
        else
   

          call crossp3d(L,node,H)
          Hnorm=Sqrt(Dot_Product(H,H))
          om=Atan2(Dot_Product(LRL,H)*nnorm,Dot_Product(LRL,node)*Hnorm)*igrad
              if(om.lt.0._double) om=om+360._double 
  
        end if  

          if(e.lt.1.2d-13.and.incl.le.1.d-13) then
              gom=0._double
              om=0._double

              !mean anomaly 'mmm'
              mmm=Atan2(r(2),r(1))*igrad
             if(mmm.lt.0._double) mmm=mmm+360._double

        

           elseif (e.lt.1.2d-13.and.incl.gt.1.d-13) then
          

              !transformation of positionvector r into coordinate-system U spanned by Angular Momentum vector, vector Sun-ascending Node, 
              !and their crossproduct H, corresponding to motion of rU in the orbital plane. 
              

              call crossp3d(L,node,H)
              Hnorm=Sqrt(Dot_Product(H,H))
              ! calculation of the base - transformation Matrix
              U(:,1)=node(:)/nnorm
              U(:,2)=H/Hnorm
              U(:,3)=L(:)/Sqrt(L2)
              
              !coordinate transformation
              call gauss(dim,U,r,rU)

              !true anomaly in orbital plane
              tanom=Atan2(rU(2),rU(1))

              !mean anomaly 'mmm'
              mmm=tanom*igrad     
              if(mmm.lt.0._double) mmm=mmm+360._double
      

            elseif (incl.lt.1.d-13.and.e.gt.1.d-13) then
              !calculation of true anomaly via vector "H" at right angle to LRL and L (alternative)
                      call crossp3d(L,LRL,H)
                      Hnorm=Sqrt(Dot_Product(H,H))
                      tanom= Atan2(Dot_Product(H,r)*e,Dot_Product(LRL,r)*Hnorm) 
 
              !calculation of eccentric anomaly eanom (=E)
               cosE=(e+Cos(tanom))/(1._double+e*Cos(tanom))
               sinE=Sqrt(1._double-e2)*Sin(tanom)/(1._double+e*Cos(tanom))
               eanom=ATAN2(sinE,cosE)  

              !mean anomaly 'mmm' via Kepler's equation
              mmm=(eanom-e*sinE)*igrad              
 
!                write(*,*)'hey'

                 if(mmm.lt.0._double) mmm=mmm+360._double

              else
        
  
             call crossp3d(L,LRL,H)
              Hnorm=Sqrt(Dot_Product(H,H))
              ! calculation of the base - transformation Matrix
              U(:,1)=LRL(:)/e
              U(:,2)=H/Hnorm
              U(:,3)=L(:)/Sqrt(L2)
              
              !coordinate transformation
              call gauss(dim,U,r,rU)

              !true anomaly in orbital plane
              tanom=Atan2(rU(2),rU(1))


              !calculation of eccentric anomaly eanom (=E)
               cosE=(e+Cos(tanom))/(1._double+e*Cos(tanom))
               sinE=Sqrt(1._double-e2)*Sin(tanom)/(1._double+e*Cos(tanom))
               eanom=Atan2(sinE,cosE)  
               
                
               
              !mean anomaly 'mmm' via Kepler's equation
             mmm=(eanom-e*sinE)*igrad
             if(mmm.lt.0._double) mmm=mmm+360._double
           end if
      
          if(e.lt.0._double.or.e.gt.1._double) e=1._double
          if(a.lt.0._double.or.a.gt.huge(a)) a=0._double 
          if(incl.gt.180._double.or.incl.lt.0._double) om=-100._double

          ele(i,1)=a
          ele(i,2)=e
          ele(i,3)=incl
          ele(i,4)=om
          ele(i,5)=gom
          ele(i,6)=mmm

          do j=4,6
            if(ele(i,j).gt.360._double.or.ele(i,j).lt.0._double.or.ele(i,j).ne.ele(i,j)) then
              ele(i,j)=-100._double
            end if
          end do

       end do
!$omp end do       
!$omp end parallel
     
       RETURN
     END subroutine trojtrnsel

!   *****************************************************************
!...BERECHNUNG DER BAHNELEMENTE AUS DEN KARTHESISCHEN KOORDINATEN
!...UND GESCHWINDIGKEITEN 

  subroutine xtoel(x,ele,mu)

        IMPLICIT none
        real(kind=double)::mu,rsquare,dotrv,vsquare,crossi,crossj,crossk
        real(kind=double)::rvet,rsobrea,semia,ecosu,esinu,esquare
        real(kind=double)::excen,anomu,anomalia,nodo,sinimod,incli
        real(kind=double)::psubz,qsubz,argumento
        real(kind=double):: X(1:6),ELE(1:6),V(1:3)
!	
        V(1)=X(4)/SQRT(mu)
        V(2)=X(5)/SQRT(mu)
        V(3)=X(6)/SQRT(mu)
        RSQUARE=X(1)**2+X(2)**2+X(3)**2
        DOTRV=X(1)*V(1)+X(2)*V(2)+X(3)*V(3)
        VSQUARE=V(1)**2+V(2)**2+V(3)**2
        CROSSI=X(2)*V(3)-X(3)*V(2)
        CROSSJ=X(3)*V(1)-X(1)*V(3)
        CROSSK=X(1)*V(2)-X(2)*V(1)
        RVET=SQRT(RSQUARE)
        RSOBREA=2.D0-RVET*VSQUARE
        SEMIA=RVET/RSOBREA
!  write(*,*)'xtoel',semia,rvet,vsquare,rsobrea
        ECOSU=1.D0-RSOBREA
        ESINU=DOTRV/SQRT(SEMIA)
        ESQUARE=ECOSU**2._double+ESINU**2._double
        EXCEN=SQRT(ESQUARE)
        ANOMU=ATAN2(ESINU,ECOSU)
!		A operacao ATAN2 da resultado entre -PI e +PI.
!     write(*,*)'esinu,ecosu,anomu',esinu,ecosu,anomu
        ANOMALIA=ANOMU-EXCEN*SIN(ANOMU)
!	MODCROSS=SQRT(CROSSI**2+CROSSJ**2+CROSSK**2)
        NODO=ATAN2(CROSSI,-CROSSJ)
        SINIMOD=SQRT(CROSSI**2+CROSSJ**2)      
        INCLI=ATAN2(SINIMOD,CROSSK)
        PSUBZ=X(3)/RVET*COS(ANOMU)-V(3)*SQRT(SEMIA)*SIN(ANOMU)
        QSUBZ=X(3)/RVET*SIN(ANOMU)+V(3)*SQRT(SEMIA)* (COS(ANOMU)-EXCEN)
        QSUBZ=QSUBZ/SQRT(1-ESQUARE)
        ARGUMENTO=ATAN2(PSUBZ,QSUBZ)

 !      write(*,*)'ANOMU',ANOMU
        ELE(1)=SEMIA
        ELE(2)=EXCEN
        ELE(3)=INCLI
        ELE(5)=NODO
        ELE(4)=ARGUMENTO
        ELE(6)=ANOMALIA
 
 
        RETURN
END subroutine xtoel
!**************************************************************************************


subroutine nr(m,ecc,ea)
!***********************************************
! solves kepler's equation by applying  Newton Raphson Method
!with an accuracy limit of E-14
!
! m[real]...............Mean Anomaly (radian!)
!ecc[real]..............Eccentricity (numerical Eccentricity <1!)
!ea[real]................Eccentric Anomaly (radian!)
!step[int]...............Iteration steps
!
!written by Siegfried Eggl  20061222
!modified                   20111026
!dependencies: none
!**********************************************
implicit none
        integer(kind=intdble)::i
        real(kind=double)::m,ecc,ea,ea0,dea,deat


ea=0._double
ea0=1.3421d0
deat=1.E-12
dea=ABS(ea-ea0)


do i=1,100

   if (dea>deat)then

      ea=ea0-(ea0-ecc*SIN(ea0)-m)/(1._double-ecc*COS(ea0))

      dea=ABS(ea-ea0)

      ea0=ea
   else 
      continue

   end if

end do

if (ea.le.1.E-14) then
   ea=0._double
end if
!if precision is not achieved try with different initial condition
if(dea>deat) then

ea=0._double
ea0=0.3421d0
deat=1.E-12
dea=ABS(ea-ea0)
 do i=1,100

   if (dea>deat)then

      ea=ea0-(ea0-ecc*SIN(ea0)-m)/(1._double-ecc*COS(ea0))

      dea=ABS(ea-ea0)

      ea0=ea
   else 
      continue

   end if

 end do
end if


if(dea>deat) then

ea=0._double
ea0=2.3421d0
deat=1.E-12
dea=ABS(ea-ea0)
 do i=1,100

   if (dea>deat)then

      ea=ea0-(ea0-ecc*SIN(ea0)-m)/(1._double-ecc*COS(ea0))

      dea=ABS(ea-ea0)

      ea0=ea
   else 
      continue

   end if

 end do
end if



if(dea>deat) then
   write(unit=*,fmt=*)'convergence-error in subroutine nr'
   write(unit=*,fmt=*)'target precision:',deat
   write(unit=*,fmt=*)'achieved precision:',dea
end if
return

end subroutine nr
!*****************************************************************
SUBROUTINE hekoo(rv,body)
!------------------------------------------------------
! Transformation ins Heliozentirsche Koordinatensystem
! Die Koordinaten und Geschwindigkeiten  der Sonne werden von denen der massiven Körper
! abgezogen. Anschließend wird die Sonne in den Koordinaten Mittelpunkt gesetzt.
! Würde man die Schleife j bei 1 starten, hätte die Sonne schon nach dem ersten
! Schritt den Nullvektor als Koordinaten, und man bekäme kein brauchbares Ergebnis
!-------------------------------------------------------
implicit none
integer(kind=intdble)::body,j
real(kind=double),intent(inout)::rv(1:body,1:6)

!$omp parallel do default(shared)   &
!$omp private(j)
do j=2,body
   rv(j,1:6)=rv(j,1:6)-rv(1,1:6)
end do
!$omp end parallel do

   rv(1,1:6)=0._double

return
end subroutine hekoo
!*******************************************************
SUBROUTINE bakoo(rv,mass,body)
!-----------------------------------------------------
! Transformation ins Baryzentrische Koordinatensystem
!-----------------------------------------------------
implicit none
integer(kind=intdble)::body,j
real(kind=double),intent(inout)::rv(1:body,1:6)
real(kind=double),intent(in)::mass(1:body)
real(kind=double)::rb(1:6),rbx,rby,rbz,rbvx,rbvy,rbvz,summass

!sorry, but open mp reduction just works for scalars...
rbx=0._double
rby=0._double
rbz=0._double
rbvx=0._double
rbvy=0._double
rbvz=0._double

!$omp parallel do default(shared)   &
!$omp private(j) &
!$omp reduction(+:rbx,rby,rbz,rbvx,rbvy,rbvz)
do j=1,body
   rbx=rbx+mass(j)*rv(j,1)
   rby=rby+mass(j)*rv(j,2)
   rbz=rbz+mass(j)*rv(j,3)
   rbvx=rbvx+mass(j)*rv(j,4)
   rbvy=rbvy+mass(j)*rv(j,5)
   rbvz=rbvz+mass(j)*rv(j,6)
end do
!$omp end parallel do 

summass=Sum(mass)

rb(1)=rbx/summass
rb(2)=rby/summass
rb(3)=rbz/summass
rb(4)=rbvx/summass
rb(5)=rbvy/summass
rb(6)=rbvz/summass

!$omp parallel do default(NONE)   &
!$omp shared(rv,body,rb) private(j)
do j=1,body
   rv(j,1:6)=rv(j,1:6)-rb(1:6)
end do
!$omp end parallel do

return
end subroutine bakoo
!*****************************************************************
SUBROUTINE binbakoo(rv,mass,body)
!-----------------------------------------------------
! Transformation ins Baryzentrische Koordinatensystem
!-----------------------------------------------------
implicit none
integer(kind=intdble)::body,j
real(kind=double),intent(inout)::rv(1:body,1:6)
real(kind=double),intent(in)::mass(1:body)
real(kind=double)::rb(1:6),rbx,rby,rbz,rbvx,rbvy,rbvz,summass,rbbin(1:6)


do j=1,body
write(*,*)'rv init',j,rv(j,:)
end do

rbbin(1:6)=0._double

do j=1,2
   rbbin=rbbin+mass(j)*rv(j,:)
end do

summass=Sum(mass(1:2))
rbbin=rbbin/summass
do j=1,2
   rv(j,:)=rv(j,:)-rbbin
end do
!back up old binary barycentric rv data
!rvbinold(1:2,1:6)=rv(1:2,1:6)

if(body.ge.3) then
!sorry, but open mp reduction just works for scalars...
rbx=rbbin(1)
rby=rbbin(2)
rbz=rbbin(3)
rbvx=rbbin(4)
rbvy=rbbin(5)
rbvz=rbbin(6)

!$omp parallel do default(shared)   &
!$omp private(j) &
!$omp reduction(+:rbx,rby,rbz,rbvx,rbvy,rbvz)
do j=3,body
   rbx=rbx+mass(j)*rv(j,1)
   rby=rby+mass(j)*rv(j,2)
   rbz=rbz+mass(j)*rv(j,3)
   rbvx=rbvx+mass(j)*rv(j,4)
   rbvy=rbvy+mass(j)*rv(j,5)
   rbvz=rbvz+mass(j)*rv(j,6)
end do
!$omp end parallel do 

summass=Sum(mass(1:body))

rb(1)=rbx/summass
rb(2)=rby/summass
rb(3)=rbz/summass
rb(4)=rbvx/summass
rb(5)=rbvy/summass
rb(6)=rbvz/summass

!$omp parallel do default(NONE)   &
!$omp shared(rv,body,rb) private(j)
do j=1,body
   rv(j,1:6)=rv(j,1:6)-rb(1:6)
end do
!$omp end parallel do
! 
! rbbin=rbbin-rb
! rv(1,:)=rbbin*sum(mass(1:2))/mass(1)-mass(2)/mass(1)*rvbinold(2,:)
! rv(2,:)=rbbin*sum(mass(1:2))/mass(2)-mass(1)/mass(2)*rvbinold(1,:)


call bakoo(rv,mass,body)

end if
! write(*,*)'rbbin',rbbin
! 
! write(*,*)'rb',rb
! 
! do j=1,body
! write(*,*)'rv',j,rv(j,:)
! end do
! 
! write(*,*)'rbbin after',(rv(1,:)*mass(1)+rv(2,:)*mass(2))/(mass(1)+mass(2))
! 
! write(*,*)'rb after',((rv(1,:)*mass(1)+rv(2,:)*mass(2))+rv(3,:)*mass(3))/sum(mass(1:3))


return
end subroutine binbakoo


!**************************************************************
SUBROUTINE jacoo(rv,mass,body)
!-----------------------------------------------------
! Jacobian coordinates
!-----------------------------------------------------
implicit none
integer(kind=intdble)::body,j
real(kind=double),intent(inout)::rv(1:body,1:6)
real(kind=double),intent(in)::mass(1:body)
real(kind=double)::msum(1:body),rbbin(1:body,1:6)

rbbin(:,:)=0._double

rbbin(1,1:6)=rv(2,:)-rv(1,:)
msum(1)=mass(1)

do j=2,body
 msum(j)=msum(j-1)+mass(j)
 rbbin(j,1:6)=(rbbin(j,1:6)+mass(j)*rv(j,1:6))/msum(j)
end do

!$omp parallel do default(NONE)   &
!$omp shared(rv,body,rb) private(j)
do j=1,body
   rv(j,1:6)=rbbin(j,1:6)
end do
!$omp end parallel do

return
end subroutine jacoo
!*****************************************************************
SUBROUTINE binhekoo(rv,mass,body)
!------------------------------------------------------
! Transformation into the coordinate system of the binary's barycenter
! secondary will be given with respect to primary
!-------------------------------------------------------
implicit none
integer(kind=intdble)::body,j
real(kind=double),intent(in)::mass(1:body)
real(kind=double),intent(inout)::rv(1:body,1:6)
real(kind=double)::rbbin(1:6),summass

rbbin(1:6)=0._double

do j=1,2
   rbbin=rbbin+mass(j)*rv(j,:)
end do

summass=Sum(mass(1:2))

rbbin=rbbin/summass

!coordinates of secondary will be with respect to primary
 rv(2,1:6)=rv(2,1:6)-rv(1,1:6)

!coordinates of the rest will be with respect to binary barycenter

do j=3,body
   rv(j,1:6)=rv(j,1:6)-rbbin(1:6)
end do

!primary is at the center of coordinates
   rv(1,1:6)=0._double

return
end subroutine binhekoo
!************************************************************

subroutine crossp3d(a,b,c)
!--------------------------------
!  3 dimensional crossproduct
!  a,b .... input vectors
!  c   .... output vector
!  c = a x b
!------------------------------- 
implicit none
real(kind=double),dimension(1:3),intent(in)::a,b
real(kind=double),dimension(1:3),intent(out)::c

c(1)=a(2)*b(3)-b(2)*a(3)
c(2)=a(3)*b(1)-b(3)*a(1)
c(3)=a(1)*b(2)-b(1)*a(2)

return
end subroutine crossp3d
!////////////////////////////////////////////////////////////////
end module transform_m
