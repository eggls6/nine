module output_m
  use global_m                      
  use transform_m

  implicit none

  public::OpenFile
  public::CloseFile
  public::Out

  character(len=32),private::filename(1:12)


!/////////////////////////////////////////////////////////////////////////////////////////////////
!                                           Ausgabe und Auswertung
!
!                                         by Siegfried Eggl  20080110 
!------------------------------------------------------------------------------
!**********************************************************************
contains
  
!**********************************************************************
  subroutine OpenFile(p,input,Dateiname)
  use global_m
                       implicit none
                        type(Problem)::p
                        integer(kind=intdble)::tdigit
                        real(kind=double)::dtd,nn
                        
                        character(len=5),intent(in)::Dateiname
                        character(len=2),intent(in)::input
                        character(len=4)::fileend1 = '.hco'
                        character(len=4)::fileend2 = '.hel'
                        character(len=4)::fileend3 = '.egy'
                        character(len=4)::fileend4='.bel'
                        character(len=4)::fileend5='.bco'
                        character(len=4)::fileend6='.def'
                       
                        character(len=4)::fileend7 = '.bhc'
                        character(len=4)::fileend8 = '.bhe'
                        character(len=4)::fileend9='.bbc'

                        character(len=4)::fileend10='.cut'

                        character(len=4)::fileend11='.mrg'
                        
                        character(len=4)::fileend12='.cen'
                 
                        character(len=20)::acc,form

!------------------------------------------------------------------------
!For Binary output: use gfortran or ifort(not supported???)

      acc = 'stream'
      form = 'unformatted'
!---------------------------------------------------------------------

                        write(unit=*,fmt=*)'output file(s) will be:'

                                               
                        filename(1) = Dateiname//fileend1
if(p%outhco) then
        open(unit=22, file=filename(1), status="replace", action="write")   
                        write(unit=22,fmt=1001)
                        write(unit=*,fmt=*)filename(1)
end if 
 
                         filename(2) = Dateiname//fileend2
if(p%outhele) then
        open(unit=23, file=filename(2), status="replace", action="write")
                if(input.eq.'ye'.or. input.eq.'Ye'.or.input.eq.'YE') then
                        write(unit=23,fmt=1007)
                else  
                        write(unit=23,fmt=1002)
                end if
                        write(unit=*,fmt=*)filename(2)
end if

                        filename(3) = Dateiname//fileend3
                     
if(p%outegy) then
        open(unit=24, file=filename(3), status="replace", action="write")
                        write(unit=24,fmt=1003)
                        write(unit=*,fmt=*)filename(3)
end if     

                         filename(4) = Dateiname//fileend4
if(p%outbele) then
        open(unit=25, file=filename(4), status="replace", action="write")
                        write(unit=25,fmt=1002)
                        write(unit=*,fmt=*)filename(4)
end if
                 
                         filename(5) = Dateiname//fileend5
if(p%outbco) then
        open(unit=26, file=filename(5), status="replace", action="write")
                        write(unit=26,fmt=1001)
                        write(unit=*,fmt=*)filename(5)
end if
               
                       filename(6) = Dateiname//fileend6
if(p%outdef) then
        open(unit=27, file=filename(6), status="replace", action="write")
                        write(unit=*,fmt=*)filename(6)
end if

!--------------------------Binary output-------------------------------------------

                        filename(7) = Dateiname//fileend7
if(p%outbinhco) then
        open(unit=28, file=filename(7), status='replace',form=form,Access=acc)   
                 
                        write(unit=*,fmt=*)filename(7)
end if 
 
                         filename(8) = Dateiname//fileend8
if(p%outbinhele) then
        open(unit=29, file=filename(8),  status='replace',form=form,Access=acc)
                 
                        write(unit=*,fmt=*)filename(8)
end if

                        filename(9) = Dateiname//fileend9
if(p%outbinbco) then
        open(unit=30, file=filename(9), status='replace',form=form,Access=acc)   
                     
                        write(unit=*,fmt=*)filename(9)
end if 
!-----------------------Cutoff output--------------------------------------------


                        filename(10) = Dateiname//fileend10
if(p%cutoff.gt.0.d0) then
        open(unit=31, file=filename(10), status='replace',action="write")   
                     
                        write(unit=*,fmt=*)filename(10)   
                        write(unit=31,fmt=1004)
end if 

!-----------------------Merge output-----------------------------------------------
                        filename(11) = Dateiname//fileend11
if(p%merge) then
        open(unit=32, file=filename(11), status='replace',action="write")   
                     
                        write(unit=*,fmt=*)filename(11)   
                        write(unit=32,fmt=1005)
end if 

!---------------------Close encouter output
                       filename(12) = Dateiname//fileend12
if(p%cen) then                      
        open(unit=33, file=filename(12), status='replace',action="write")   
                     
                        write(unit=*,fmt=*)filename(12)   
                        write(unit=33,fmt=1006)
end if 



if(showprog) then
write(unit=*,fmt=*)'progress:'
else
 write(unit=*,fmt=*)'working...'      
end if


!---------------------Fix time output format----------------

dtd=huge(dtd)
nn=1._double
tdigit=0
if(p%tAusg.lt.1._double) then
 do while(dtd.gt.1.d-13)
  tdigit=tdigit+1
  nn=nn*10._double
  dtd=nn*p%tAusg-real(int(nn*p%tAusg,kind=intdble),kind=double)
 end do 
else
 do while(dtd.gt.1.d-13)
   tdigit=tdigit+1
   nn=nn/10._double
  dtd=nn*p%tAusg-real(int(nn*p%tAusg,kind=intdble),kind=double)
 end do 
end if

 if (tdigit.gt.0.and.tdigit.lt.100) then
   write(tdigc,"(I2)")tdigit
 else 
   tdigc="00"
 end if
!define global output format vars tdigc, tendc
write(tdigc,"(I2)")int(log10(p%tEnde))+tdigit
write(tendc,"(I2)")int(log10(p%tEnde))+tdigit+7


!------------output format statements---------------------------------------

1001 format(10X,'time',15X,'body',20X,'x [AU]',20X,'y [AU]',20X,'z [AU]',20X, &
             'vx [AU/D]',20X,'vy [AU/D]',20X,'vz [AU/D]',20X,'mass [M_sun]',20X,'name')

1002 format(10X,'time',15X,'body',25X,'a',25X,'e',25X,'i',25X,'omega',25X,'OMEGA', &
              25X,'mean anomaly',20X,'mass [M_sun]',20X,'name')

1003 format( 18X, 'time',20X,'total energy',10X,'delta energy',10X,'kinetic energy',10X,'potential energy',& 
                             10X, 'sum of delta energy', &
                             10X,'total angular momentum', 10X,'delta t.a.m',10X,'sum of delta t.a.m.', &
                             20X,'t.a.m. x',20X,'t.a.m. y',20X,'t.a.m. z', & 
                             10X,'barycenter x',10X,'barycenter y',10X,'barycenter z',20X)
1004 format(10X,'time',15X,'body',25X,'x',25X,'y',25X,'z',25X,'vx',25X,'vy',25X,'vz',25X,'specific energy')

1005 format(10X,'time',15X,'body1',5X,'body2',10X,'x',15X,'y',15X,'z',15X,'v_relative/v_escape')

1006 format(10X,'time',15X,'body1',20X,'body2',15X,'distance [AU]',25X,'x1 [AU]',25X,'y1 [AU]',25X,'z1 [AU]',25X, &
                 'vx1 [AU/D]',20X,'vy1 [AU/D]', 20X,'vz1 [AU/D]',20X,'x2 [AU]',20X,'y2 [AU]',20X,'z2 [AU]',20X, &
                 'vx2 [AU/D]',25X,'vy2 [AU/D]',25X,'vz2 [AU/D]')
                            
1007 format(10X,'time',15X,'body',25X,'a',25X,'e',25X,'i',25X,'omega',25X,'OMEGA', &
              20X,'mean anomaly',5X,'Yark da/dt diurnal [AU/Myr]',3X, 'Yark < da/dt > diurnal [AU/Myr]', &
              5X,'mass [M_sun]',8X,'name')                            
 

end subroutine OpenFile

!*************************************************
  subroutine CloseFile
                       close(unit=21)
                       close(unit=22)
                       close(unit=23)
                       close(unit=24)
                       close(unit=25)
                       close(unit=26)
                       close(unit=27)
                       close(unit=28)
                       close(unit=29)
                       close(unit=30)
                       close(unit=31)
                       write(*,*)'all units closed'

  end subroutine CloseFile

!*************************************************
subroutine default(p,rv,ie,il)
!called, when no outputfiles have been specified
integer(kind=intdble)::i
type(Problem)::p
real(kind=double),dimension(1:p%NrKoerper,1:6)::rv
real(kind=double)::ie,il




  if(p%outdef) then
    write(unit=27,fmt='(A)')'Initial Conditions'
    write(unit=27,fmt=1020)
     do i=1,p%NrKoerper
          write(unit=27,fmt="(f15.2,3X,I8,3X)",advance="no") 0.d0,i
          write(unit=27,fmt="(7(f25.16,3X))",advance="yes") rv(i,1:3),rv(i,4:6)*k,bname(i)
    end do
    
    write(unit=27,fmt=*)'total energy = ',ie
    write(unit=27,fmt=*)'total angular momentum = ',il
end if

1020  format(18X,'time',5X,'body',35X,'x',35X,'y',35X,'z',35X,'vx',35X,'vy',35X,'vz')

end subroutine default
!***************************************************
subroutine Out(tt, rvin,massin,p)
use global_m

 !output routine for all integrators
real(kind=double),intent(in)::tt
!type(vec),dimension(:),intent(in)::r,v
real(kind=double),dimension(:),intent(in)::massin
type(Problem),intent(in)::p
real(kind=double),dimension(1:p%NrKoerper,1:6),intent(in)::rvin
real(kind=double),dimension(1:p%NrKoerper,1:6)::rv,ele,rvh,rvb
real(kind=double),dimension(1:p%NrKoerper)::m
real(kind=double)::engy,ekin,epot,lpsum,lpsum2(1:3),rb(1:6)
real(kind=double), dimension(:,:),allocatable::lp
!real(kind=double), dimension(1:p%NrKoerper,1:3)::lrl
real(kind=double)::lperror,engyerror,t
!real(kind=double)::omega1,omega2
!real(kind=double)::spezen,spezenel
integer(kind=intdble)::n


t=tt+p%tstart
if(tt.eq.toutprev) then
else

 if(p%outbco) then
          !restore initial body identities (resort bodies)
         if(intback) then
          rvb(bodypos(:),1:3)=rvin(:,1:3)
          rvb(bodypos(:),4:6)=-rvin(:,4:6)
         else
          rvb(bodypos(:),:)=rvin(:,:)
         end if
         m(bodypos(:))=massin(:)
          ! transformation into barycentric coordinates
      do n=1,p%NrKoerper
         !Barycentric
        write(unit=26,fmt="(ES"//tendc//"."//tdigc//",3X,I8,3X)",advance="no") t,n      
        write(unit=26,fmt="(7(f25.16,3X),A)",advance="yes") rvb(n,1:3),rvb(n,4:6)*k,m(n),bname(n)                     
      end do       
   end if
  !OUTPUT IN ORBITAL ELEMENTS WITH RESPECT TO THE BINARIES' BARYCENTER
   if(p%outbele) then  
           !restore initial body identities (resort bodies)
         !rvh(bodypos(:),:)=rvin(:,:)
         if(intback) then
          rvh(bodypos(:),1:3)=rvin(:,1:3)
          rvh(bodypos(:),4:6)=-rvin(:,4:6)
         else
          rvh(bodypos(:),:)=rvin(:,:)
         end if
         m(bodypos(:))=massin(:)
     ! transformation into binary elements

     call binhekoo(rvh,m,p%NrKoerper)
     call btrnsel( p%NrKoerper,rvh,m,ele)  
     do n=2,p%NrKoerper
         !output 
        write(unit=25,fmt="(ES"//tendc//"."//tdigc//",3X,I8,3X)",advance="no") t,n
        write(unit=25,fmt="(6(f25.16,3X),A)",advance="yes")ele (n,:),bname(n)                     
     end do
        if(cutoffstop(1)) then
     if (maxval(ele(2:totnbody-m0count,2)).ge.1.d0) then
                 call CloseFile
                      write(*,*)'eccentricity of massive body >= 1'
             STOP
             end if
      end if
      if(cutoffstop(2))then
            if (maxval(ele(totnbody-m0count+1:totnbody,2)).ge.1.d0) then
                 call CloseFile
                      write(*,*)'eccentricity of massless body >= 1'
             STOP
             end if
      end if  
   end if


   if(p%outhele) then
      !restore initial body identities (resort bodies)
         !rvh(bodypos(:),:)=rvin(:,:)
         if(intback) then
          rvh(bodypos(:),1:3)=rvin(:,1:3)
          rvh(bodypos(:),4:6)=-rvin(:,4:6)
         else
          rvh(bodypos(:),:)=rvin(:,:)
         end if
         m(bodypos(:))=massin(:)
      ! transformation into heliocentric elements
      call hekoo(rvh,p%NrKoerper)

      if(p%outtroj) then
        call trojtrnsel( p%NrKoerper,rvh,m,ele)
      else
        call htrnsel( p%NrKoerper,rvh,m,ele)
      end if
      do n=2,p%NrKoerper                                          
         !Ausgabe Bahnelemente
         write(unit=23,fmt="(ES"//tendc//"."//tdigc//",3X,I8,3X)",advance="no") t,n
         if(p%nrint.eq.8.or.p%nrint.eq.6) then
         write(unit=23,fmt="(9(f25.16,3X),A)",advance="yes")ele (n,:),dadt_yark(n,:),m(n),bname(n)
         else
          write(unit=23,fmt="(7(f25.16,3X),A)",advance="yes")ele (n,:),m(n),bname(n)       
         end if           
      end do
      
      if(cutoffstop(1)) then
             if (maxval(ele(2:totnbody-m0count,2)).ge.1.d0) then
                 call CloseFile
                 write(*,*)'eccentricity of massive body >= 1'
             STOP
             end if
      end if
      if(cutoffstop(2))then
            if (maxval(ele(totnbody-m0count+1:totnbody,2)).ge.1.d0) then
                 call CloseFile
                      write(*,*)'eccentricity of massless body >= 1'
             STOP
             end if
      end if       
   end if
   
  

   if(p%outhco) then
        !restore initial body identities (resort bodies)
         !rvh(bodypos(:),:)=rvin(:,:)
         if(intback) then
          rvh(bodypos(:),1:3)=rvin(:,1:3)
          rvh(bodypos(:),4:6)=-rvin(:,4:6)
         else
          rvh(bodypos(:),:)=rvin(:,:)
         end if
         m(bodypos(:))=massin(:)
      ! transformation into heliocentric coordinates 
         call hekoo(rvh,p%NrKoerper)
      do n=2,p%NrKoerper       
         !Augsgabe Koordinaten und Geschwindigkeiten (außer der der Sonne die ja heliozentisch = 0.d0 sind)
         write(unit=22,fmt="(ES"//tendc//"."//tdigc//",3X,I8,3X)",advance="no") t,n
         write(unit=22,fmt="(7(f25.16,3X),A)",advance="yes") rvh(n,1:3),rvh(n,4:6)*k,m(n),bname(n)       
      end do
   end if
       

!!!!!!!!!!BINARY OUTPUT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(p%outbinbco) then
          !restore initial body identities (resort bodies)
         !rvb(bodypos(:),:)=rvin(:,:)
         if(intback) then
          rvb(bodypos(:),1:3)=rvin(:,1:3)
          rvb(bodypos(:),4:6)=-rvin(:,4:6)
         else
          rvb(bodypos(:),:)=rvin(:,:)
         end if
         m(bodypos(:))=massin(:)
          ! transformation into barycentric coordinates
      do n=1,p%NrKoerper
         !Augsgabe Koordinaten und Geschwindigkeiten
         write(unit=30)t,n,rvb(n,1:3),rvb(n,4:6)*k                
      
      end do       
   end if
   
!!!!OUTPUT IN BARYCENTRIC ORBITAL ELEMENTS!!!!!!!   

   if(p%outbinhele) then
      !restore initial body identities (resort bodies)
         !rvh(bodypos(:),:)=rvin(:,:)
         if(intback) then
          rvh(bodypos(:),1:3)=rvin(:,1:3)
          rvh(bodypos(:),4:6)=-rvin(:,4:6)
         else
          rvh(bodypos(:),:)=rvin(:,:)
         end if
         m(bodypos(:))=massin(:)
      ! transformation into heliocentric elements
      call hekoo(rvh,p%NrKoerper)
      
       if(p%outtroj) then
        call trojtrnsel( p%NrKoerper,rvh,m,ele)
      else
        call htrnsel( p%NrKoerper,rvh,m,ele)
      end if

      do n=2,p%NrKoerper                                          
         !Ausgabe Bahnelemente
         write(unit=29)t,n,ele(n,:)                   
      end do
      if(cutoffstop(1)) then
             if (maxval(ele(2:totnbody-m0count,2)).ge.1.d0) then
                 call CloseFile
                 write(*,*)'eccentricity of massive body >= 1'
             STOP
             end if
      end if
      if(cutoffstop(2))then
            if (maxval(ele(totnbody-m0count+1:totnbody,2)).ge.1.d0) then
                 call CloseFile
                    write(*,*)'eccentricity of massless body >= 1'
             STOP
             end if
      end if  
   end if
     !restore initial body identities (resort bodies)
         !rvh(bodypos(:),:)=rvin(:,:)
         if(intback) then
          rvh(bodypos(:),1:3)=rvin(:,1:3)
          rvh(bodypos(:),4:6)=-rvin(:,4:6)
         else
          rvh(bodypos(:),:)=rvin(:,:)
         end if
         m(bodypos(:))=massin(:)
      ! transformation into heliocentric coordinates
   if(p%outbinhco) then
         call hekoo(rvh,p%NrKoerper)
      do n=2,p%NrKoerper       
         !Augsgabe Koordinaten und Geschwindigkeiten (außer der der Sonne die ja heliozentisch = 0.d0 sind)
         write(unit=28)t,n,rvh(n,1:3),rvh(n,4:6)*k       
      end do
   end if
   

if(p%outegy) then     

! Particles are NOT being reordered to their inital numbers in order to calculate energies and angular momentum!
!rv(:,4:6)=  rvin(:,4:6)*k !Geschwindigkeitsrücktransformation (wegen v=dr/(dt*k))
!rv(:,1:3)=  rvin(:,1:3)

 if(intback) then
          rv(:,1:3)=rvin(:,1:3)
          rv(:,4:6)=-rvin(:,4:6)*k
         else
          rv(:,1:3)=rvin(:,1:3)
          rv(:,4:6)=rvin(:,4:6)*k
 end if

!Berechnung des Drehimpulses (Variable initlp wurde im Modul global_m definiert und in main berechnet)
     allocate(lp(1:p%NrKoerper,1:3))

     call drehimpuls(rv,p%NrKoerper,massin,lp,lpsum,lpsum2)
     lperror=log10(abs((lpsum-initlp+lLsum)/initlp))
     if (initlp.eq.0.d0) then
        lperror=1d14
     end if
     dLtot=dLtot+abs((lpsum-initlp+lLsum)/initlp)
     

!Berechnung der Energien    (Variable initegy wurde im Modul global_m definiert und in main berechnet)  
     call energy(rv,massin,p%NrKoerper,engy,ekin,epot) 
     if((engy-initegy+lEtot).eq.0.d0) then
        engyerror=0.d0
     else
        engyerror=log10(abs((engy-initegy+lEtot)/initegy))
     end if
     if (initegy.eq.0.d0) then
        engyerror=1d14
     end if
     dEtot=dEtot+abs((engy-initegy+lEtot)/initegy)


!Calculate Vector of Barycenter
     call barycenter(rv,m,p%NrKoerper,rb)
     
!Output in .egy file
     write(unit=24,fmt="(ES"//tendc//"."//tdigc//",1X)",advance="no") t
     write(unit=24,fmt="(5(f25.16,1X))",advance="no")engy,engyerror,ekin, &
          epot,dEtot          
     write(unit=24,fmt="(12(f25.16,1X),A)",advance="yes")lpsum,lperror,dLtot,&
          lpsum2(:)+lLtot(:),rb(1:3),Sum(m)    
                 
!     do n=1,p%NrKoerper
         !Augsgabe Drehimpuls und Laplace Runge Lenz Vektor
!         write(unit=27,fmt="(f15."//tdigc//",3X,I4,1X)",advance="no") t,n
!         write(unit=27,fmt="(3(f20.14,1X))",advance="yes") lp(n,:)
!      write(unit=27,fmt="(3(f20.14,1X))",advance="yes")lrl(n,:)                            
!      end do

      deallocate(lp)
end if         

!Prozentanzeige im Terminal
if(showprog) then
       call  progress(tt,p%tEnde)
end if
!Ausgabe wenn keine spezifischenDateien gewählt wurden
       if (t.ge.p%tEnde.and.p%outdef) then
             write(unit=27,fmt='(A)')'Results'
             write(unit=27,fmt=1004)
              do n=1,p%NrKoerper      
                 write(unit=27,fmt="(ES"//tendc//"."//tdigc//",3X,I8,3X)",advance="no") t,n
                 write(unit=27,fmt="(6(f25.16,3X))",advance="yes") rv(n,:)
              end do

                call drehimpuls(rv,p%NrKoerper-m0count,m,lp,lpsum,lpsum2)
                call energy(rv,m,p%NrKoerper-m0count,engy,ekin,epot)  

                 
                write(unit=27,fmt=*)'total energy = ',engy
                write(unit=27,fmt=*)'total angular momentum = ',lpsum

                write(unit=27,fmt='(A)')' '
       end if
toutprev=tt       
end if !toutprev       

1004  format(18X,'time',5X,'body',35X,'x',35X,'y',35X,'z',35X,'vx',35X,'vy',35X,'vz')
 return
  end subroutine Out
!*********************************************************
subroutine energy(rv,mass,body,engy,ekin,epot)

implicit none
integer(kind=intdble)::j,l,body
real(kind=double),intent(in)::rv(:,:),mass(:)
real(kind=double),intent(out)::ekin,epot,engy
real(kind=double)::rjl(1:3),visq(1:body)

!$omp parallel do default(shared)   &
!$omp private(j)
     do j=1,body
        visq(j)=Dot_Product(rv(j,4:6),rv(j,4:6))
     end do
!$omp end parallel do

     epot=0.d0
   
   !$omp parallel do default(shared)   &
   !$omp private(j,l,rjl) &
   !$omp reduction(+:epot)
     do j=1,body-1    
          do l=j+1,body
            if(mass(j).eq.0.d0) then
             else
              rjl(1:3)=rv(l,1:3)-rv(j,1:3)
              if(all(rjl(:).eq.0.d0)) then
                else
                epot=epot+mass(j)*mass(l)/SQRT(Dot_Product(rjl(1:3),rjl(1:3)))
              end if
            end if
          end do
     end do
    !$omp end parallel do

     epot=epot*k*k

     ekin=sum(mass(1:body)*visq(1:body))/2.d0
     engy=ekin-epot


     return

     end subroutine energy
!*******************************************************
subroutine drehimpuls(rv,body,mass,lp,lpsum,lpsum2)
!***************************************************
! Berechnung des Drehimpulses
!                    lp=r x p
! des Gesamtdrehimpulses
!                    lpsum2=sum(lp)
! und dessen Betrag
!                    lpsum=|lpsum2|
!dependencies: crossproduct2
!**************************************************
!!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
implicit none
  integer(kind=intdble)::j,body
  real(kind=double),intent(in):: rv(:,:),mass(:)
  real(kind=double),intent(out)::lpsum2(1:3),lpsum,lp(:,:)

  lpsum=0.d0
  lpsum2(:)=0.d0


  call crossproduct2(rv,body,lp)
  

  do j=1,body
     if(mass(j).ne.0.d0) then 
          lpsum2(:)=lpsum2(:)+mass(j)*lp(j,:)
     end if
  end do
     lpsum=Sqrt(Dot_Product(lpsum2(1:3),lpsum2(1:3)))
  return
end subroutine drehimpuls
!*****************************************************************
 SUBROUTINE crossproduct2(a,body,c)

!---------------------------------------------------
!problemspezifisches Programm
! berechnet das Kreuzprodukt zweier 3 dimensionaler Vektoren 
! a(1:3) x a(4:6) für  N (body) - Körper
!
! dependencies: none
!--------------------------------------------------
        
        implicit none
        integer(kind=intdble)::body,i
        real(kind=double),intent(in)::a(:,:)
        real(kind=double)::b(1:body,1:3),c(:,:)

!$omp parallel do default(shared)   &
!$omp private(i)
        do i=1,body
         b(i,1:3)=a(i,4:6)
           c(i,1)=a(i,2)*b(i,3)-a(i,3)*b(i,2)
           c(i,2)=a(i,3)*b(i,1)-a(i,1)*b(i,3)
           c(i,3)=a(i,1)*b(i,2)-a(i,2)*b(i,1)
         end do
!$omp end parallel do
  
        return

      end subroutine crossproduct2

!*******************************************************
SUBROUTINE barycenter(rv,mass,body,rb)
!-----------------------------------------------------
! Transformation ins Baryzentrische Koordinatensystem
!-----------------------------------------------------
implicit none
integer(kind=intdble)::body,j
real(kind=double),intent(in)::rv(1:body,1:6)
real(kind=double),intent(in)::mass(1:body)
real(kind=double)::rb(1:6)

rb(:)=0.d0

do j=1,body
   rb(:)=rb(:)+mass(j)*rv(j,:)
end do

  rb(:)=rb(:)/Sum(mass)
  return
end SUBROUTINE barycenter


!*****************************************************
! ! subroutine f2(rv,n,mass,force)
! ! implicit none
! ! integer(kind=intdble)::n,i,j,l
! ! real(kind=double)::rv(1:n,1:6),force(1:n,1:3),mass(1:n),rij(1:n,1:n,1:3),d(1:n,1:n)
! ! 
! ! force(:,:)=0.d0
! ! 
! ! do i=1,n-1
! !    do j=i+1,n
! !      do l=1,3
! !          rij(i,j,l)=rv(j,l)-rv(i,l)
! !          rij(j,i,l)=rv(i,l)-rv(j,l)
! !        end do
! !     end do
! ! end do
! ! 
! ! do i=1,n
! !    rij(i,i,:)=0.d0
! ! end do
! ! 
! ! call dist3(rij,n,d)
! ! 
! ! 
! !    do i=1,n
! !       do j=1,n
! !        if(i.ne.j) then
! !           force(i,:)=force(i,:)+mass(j)*rij(i,j,:)/d(i,j)**3.d0
! !        end if
! !       end do
! !    end do
! !    
! ! return
! ! end subroutine f2
! ! !*********************************************************************************
! ! subroutine dist3(rij,n,d)
! ! !----------------------------------
! ! ! problemspezifisches Programm
! ! ! berechnet die Paarabstände zwischen den Körpern
! ! ! dependencies: none
! ! !----------------------------------
! ! 
! ! implicit none
! ! integer(kind=intdble)::i,j,n
! ! real(kind=double)::rij(1:n,1:n,1:3),d(1:n,1:n)
! ! 
! ! do i=1,n-1
! !   do j=i+1,n
! !    d(i,j) = SQRT(rij(i,j,1)*rij(i,j,1)+rij(i,j,2)*rij(i,j,2 )+rij(i,j,3)*rij(i,j,3))
! !    d(j,i)=d(i,j)
! !   end do
! ! end do
! ! 
! ! do i=1,n
! !    d(i,i)=0.d0
! ! end do
! ! 
! ! return
! ! 
! ! end subroutine dist3
!****************************************
subroutine progress(now,total)
!--------------------------------------
! view the progress of a calculation in %
! now[real] ... current number of iterations 
! total [real]... total number of iterations
!
! dependencies: none
!---------------------------------------
real(kind=double)::now,total
real(kind=double)::dum
character::back

!backspace character in ASCII Code
back=achar(08)
!percentage of progress
dum=now/total*100.d0
!repositioning the cursor
write(unit=*,fmt='(A)',advance='no')back//back//back//back//back//back//back
!write the result and suppress the new line ($)
write(unit=*,fmt='(F6.2,A)',advance='no')dum,'%'

!final carriage return
if (now.ge.total) then
   write(unit=*,fmt='(A)')' '
end if


end subroutine progress
!////////////////////////////////////////////////////////////////////
end module output_m


