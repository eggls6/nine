module cutoffmerge_m
  use global_m
  use output_m
  use transform_m

  implicit none

public::cutoffr
public::detectmerge
public::scattermerge
public::closeenc

contains
!*********************************************************************************************

subroutine cutoffr(m0flag,dim,time,totbody,body,m0body,cr,cutbody,m0cutbody,rv,mass,order)
!------------------------------------------------------------------------------------------
! subroutine cutoffr will check whether massive or massless paricles have moved beyond the 
! user supplied cutoffradius cr, reorder the transgressing particles to the end of the position 
! and velocitiy array rv, and return the reorderd rv-array together with adapted numbers of massive and massless paricles
!
! m0flag [logical].... does the 


! cutbody and m0cutbody have to be initialized with 0
!------------------------------------------------------------------------------------------
implicit none
integer(kind=intdble),intent(in)::totbody,dim
integer(kind=intdble),intent(inout)::body,m0body,order(1:totbody),cutbody,m0cutbody
integer(kind=intdble)::i,dumord,bmc0,dimx2,dimp1
logical::m0flag,outflag(1:totbody)
real(kind=double),intent(in)::time
real(kind=double),intent(inout),dimension(:,:)::rv
real(kind=double)::r2,v2,en(1:totbody),dumrv(1:2*dim),mass(1:totbody),cr,dummass

!real(kind=double)::Ldum(1:3)  for updating lost angular momentum 


dimx2=dim*2
dimp1=dim+1
!check massive particles for candiates that are beyond cutoff radius 

outflag=.false.

!$omp parallel do default(NONE)   &
!$omp shared(body,tstart,rv,dim,dimp1,mass,cr,dimx2,en,outflag,time,order,showprog,cutoffstop) &
!$omp private(i,r2,v2)
 massive:    do i=1,body  
    r2=Dot_Product(rv(i,1:dim),rv(i,1:dim))

    if (r2.gt.cr*cr) then
       !check the particles' specific orbital energy en = v^2/2 - mass*k^2/r       
        v2=Dot_Product(rv(i,dimp1:dimx2),rv(i,dimp1:dimx2))
          en(i)=v2/2.d0-(Sum(mass(:))-mass(i))/Sqrt(r2)
       
         !if specific energy en >0 get rid of them (system escaper)
         if(en(i).gt.0.d0) then
            outflag(i)=.true.
                if(showprog) then !defined in global_m
                  write(unit=*,fmt=*)'massive particle',order(i),' beyond cutoff radius'
                  write(unit=*,fmt=*)'specific energy of particle = ', en(i), '  => system escape at time',time 
               end if
                  if(intback) then
                write(31,2002)time+tstart,order(i),rv(i,1:3),-rv(i,4:6),en(i),en(i)*mass(i),sum(mass) 
                  else
                write(31,2002)time+tstart,order(i),rv(i,:),en(i),en(i)*mass(i),sum(mass) 
                  end if
               if(cutoffstop(1)) then
                    write(unit=*,fmt=*)'massive particle',order(i),' beyond cutoff radius'
                    write(unit=*,fmt=*)'integration stopped at',time
                    call CloseFile
                    STOP
               end if                                
         end if
    end if    
 end do massive
!$omp end parallel do         


!conserve lost energy and angular momentum?
! ! do i=1,body
! ! if(outflag(i)) then
! ! !compute energy and angular momentum loss
! !              lEtot=lEtot+en(i)*mass(i)
! !              call crossp3d(rv(i,1:3),mass(i)*rv(i,4:6),Ldum(1:3))
! !              lLtot(:)=lLtot(:)+Ldum(:)
! !              lLsum=lLsum+Sqrt(Dot_product(Ldum(:),Ldum(:)))
! ! end if
! ! end do


!relocate massive particles to reside at the end of the rv matirx and adapt particlenumbers in order to exclude the last ones from calculations. 

if(body.le.1) then
            call CloseFile
            write(unit=*,fmt=*)'all massive bodies except for one have left the system. Integration stopped...'
            STOP
else

 do i=body,1,-1

    if(outflag(i)) then
          
         !relocate location and veolcities
         dumrv(:)=rv(i,:)
         rv(i,:)=rv(body,:)
         rv(body,:)=dumrv(:)

          !relocate masses
         dummass=mass(i)
         mass(i)=mass(body)
         mass(body)=dummass

          !keep track of the original paritcle number
         dumord=order(i)
         order(i)=order(body)
         order(body)=dumord
         
          !adapt counters for the number of massive particles
         cutbody=cutbody+1
         body=body-1
        
       if(showprog) then !defined in global_m
         write(unit=*,fmt=*)'number of massive particles reduced to:',body
         write(unit=*,fmt=*)'Particle ID:',order(i), ' mass:',mass(i)
       !  write(*,*)'body,cutbody,totbody,order(body),m0cutbody',body,cutbody,totbody,order(body),m0cutbody
       end if
         if (intback) then
         write(31,fmt=*)time,order(i),rv(i,1:3),-rv(i,4:6),en(i),en(i)*mass(i),sum(mass)   
         else
         write(31,fmt=*)time,order(i),rv(i,:),en(i),en(i)*mass(i),sum(mass)  
         end if
      end if
end do 

end if

!If the integrator will seperate between massless and massive bodies...   
   if(m0flag) then 
  
!check massless particles for candiates that are beyond cutoff radius   
      outflag=.false.

      if(m0body.gt.0) then
 
         bmc0=totbody-m0cutbody

!$omp parallel do default(NONE)   &
!$omp shared(body,tstart,cutbody,bmc0,rv,dim,dimp1,cr,dimx2,en,outflag,time,order,mass,showprog,cutoffstop) &
!$omp private(i,r2,v2)
         do i=body+cutbody,bmc0   
            r2=Dot_Product(rv(i,1:dim),rv(i,1:dim))
            if (r2.gt.cr*cr) then

               !check the particles' specific orbital energy en = v^2/2 - mass*k^2/r
          
               v2=Dot_Product(rv(i,dimp1:dimx2),rv(i,dimp1:dimx2))
               en(i)=v2/2.d0-Sum(mass(:))/Sqrt(r2)
          
               !if specific energy en <0 get rid of them (system escaper)
               if(en(i).gt.0.d0) then
                  outflag(i)=.true.
                 if(showprog) then
                      write(*,*)'massless Particle',order(i), 'beyond cutoff radius'
                      write(*,*)'specific energy of particle = ', en(i), '  => system escape at time',time
                 end if
                      if(intback) then
                      write(31,2002)time+tstart,order(i),rv(i,1:3),-rv(i,4:6),en(i),en(i)*mass(i),sum(mass)
                      else
                       write(31,2002)time+tstart,order(i),rv(i,:),en(i),en(i)*mass(i),sum(mass)
                      end if
                      
                 if(cutoffstop(2)) then
                    write(unit=*,fmt=*)'massless particle',order(i),' beyond cutoff radius'
                    write(unit=*,fmt=*)'integration stopped at',time
                    call CloseFile
                    STOP
                end if     
               end if
            end if            
         end do
!$omp end parallel do 

!conserve lost energy and angular momentum?
! ! do i=body+cutbody,bmc0
! ! if(outflag(i)) then
! ! !compute energy and angular momentum loss
! !              lEtot=lEtot+en(i)*mass(i)
! !              call crossp3d(rv(i,1:3),mass(i)*rv(i,4:6),Ldum(1:3))
! !              lLtot(:)=lLtot(:)+Ldum(:)
! !              lLsum=lLsum+Sqrt(Dot_product(Ldum(:),Ldum(:)))
! ! end if
! ! end do
  
!relocate massless particles (not their masses of course) to reside at the end of the rv matirx in order to exclude them from calculations.
         massless:      do i=bmc0,body+cutbody,-1   !not from totbody, cause the last ones (totbody-m0body) left the system already
  
            if(outflag(i)) then
       
               !relocate location and veolcities
               dumrv(:)=rv(i,:)
               rv(i,:)=rv(bmc0,:)
               rv(bmc0,4:6)=dumrv(4:6)
               rv(bmc0,1:3)=0.d0
             
               !keep track of the original paritcle number
               dumord=order(i)
               order(i)=order(bmc0)
               order(bmc0)=dumord
               
               !adapt counters for the number of massles particles
               m0cutbody=m0cutbody+1
               m0body=m0body-1
               !for sequential relocation during one timestep
               bmc0=totbody-m0cutbody
            end if
         end do massless
         
      end if

 !update barycenter?
 !call bakoo(rv(1:bmc0,:),mass(1:bmc0),bmc0)
else
 !call bakoo(rv(1:body,:),mass(1:body),body)
end if !m0flag

2002 format (F15.4,5X,I10,5X,9(E15.7,5X))
   return

 end subroutine cutoffr

!***************************************************************************

subroutine hillrm(totbody,mbody,rho,rmin,rmax,mass,hillrmin)
implicit none

 integer(kind=intdble),intent(in) :: totbody,mbody
 real(kind=double),intent(in) :: rho(:,:),rmin,rmax
 real(kind=double),intent(in) :: mass(1:totbody)
 
 integer(kind=intdble)::locprim(1:1),locsecond(1:1),i,j
 real(kind=double)::onethird

 logical::himass(1:totbody)

 real(kind=double)::hillrmin(1:mbody)

himass=.true.

locprim=maxloc(mass(:))

himass(locprim)=.false.

locsecond=maxloc(mass(:),himass)

hillrmin(:)=rmin

onethird=1.d0/3.d0

!$omp parallel default(NONE)   &
!$omp shared(rho,onethird,mbody,hillrmin,locprim,locsecond,mass,rmin,rmax,totbody) &
!$omp private(i,j)
!$omp do
do i=1,mbody

!calculating hill's radius with respect to the heaviest body but preferring any smaller value with respect to another body. 
if(i.eq.locprim(1)) then !when dealing with heaviest body set smallest hillradius to rmax in order to avoid rmin being used  
   hillrmin(i)=max(min(rho(i,locsecond(1))*(onethird*mass(i)/(mass(locsecond(1))+mass(i)))**onethird, rmax),rmin)
else
   hillrmin(i)=max(min(rho(i,locprim(1))*(onethird*mass(i)/(mass(locprim(1))+mass(i)))**onethird, rmax),rmin)
end if

    do j=1,mbody
     !take the largest Hill's radius of the two interacting bodies, but compare it to the heaviest body in the system in order to guarantee merging
     if(mass(i).eq.0.d0.or.mass(j).eq.0.d0.or.i.eq.j) then
     else
        hillrmin(i)=min(rho(i,j)*(onethird*mass(i)/(mass(j)+mass(i)))**onethird,hillrmin(i))
     end if
   end do
end do
!$omp end do
!$omp end parallel
end subroutine

!*******************************************************************************
 subroutine detectmerge(totbody,mbody,mbodyp1,rmin,rmax,rho,mass,mflag,merger)
 
 implicit none

 integer(kind=intdble),intent(in) :: totbody,mbody,mbodyp1
 real(kind=double),intent(in) :: rho(:,:),rmin,rmax
 real(kind=double),intent(in) :: mass(1:totbody)
 
 integer(kind=intdble)::i,j
 real(kind=double)::hillrmin(1:mbody)
 logical,intent(inout) :: mflag,merger(1:totbody,1:mbody)


call hillrm(totbody,mbody,rho,rmin,rmax,mass,hillrmin)

mflag=.false.

!$omp parallel default(NONE)   &
!$omp shared(rho,mbody,hillrmin,mbodyp1,totbody,merger,mflag) &
!$omp private(i,j)
!$omp do
do i=1,mbody
   do j=i+1,mbody
     if (rho(j,i).le.(hillrmin(i))) then
        mflag=.true.
        merger(j,i)=.true.
     else
        merger(j,i)=.false.
     end if     
   end do

    do j=mbodyp1,totbody
     if (rho(j,i).le.(hillrmin(i))) then
        mflag=.true.
        merger(j,i)=.true.
     else
        merger(j,i)=.false.
     end if     
   end do
end do
!$omp end do
!$omp end parallel


return
end subroutine detectmerge

!**************************************************************************************
subroutine scattermerge(rmin,time,dim,totbody,body,m0body,cutbody,cutm0body,rv,rhom,mass,merger,order)

implicit none

 integer(kind=intdble),intent(in)::totbody,dim
 real(kind=double),intent(in)::rhom(:,:),time,rmin

 real(kind=double),intent(inout) ::rv(1:totbody,1:2*dim),mass(1:totbody)
 integer(kind=intdble),intent(inout)::body,m0body,cutbody,cutm0body,order(1:totbody)
 logical,intent(inout) :: merger(1:totbody,1:body)

 integer(kind=intdble)::mcount,i,j,dumord,bm0c
 real(kind=double)::vrelij(1:dim),vdiff2,vesc2

 real(kind=double)::rad(1:2)



mergemassive: do i=body-1,1,-1
 !merging of two massive bodies

  mcount=count(merger(i,:))

  if(mcount.ge.1) then
    
   do j=body,i+1,-1

      if(merger(i,j)) then
!elastic collision and scattering ? No way -> Harz4 Modell unjustified
 
         vrelij(:)=rv(j,4:6)-rv(i,4:6)
         vdiff2 = DOT_PRODUCT(vrelij(:),vrelij(:))
         vesc2 = MAX(2.d0*mass(j)/rhom(i,j),2.d0*mass(i)/rhom(i,j))    
        
        
         !calculate physical radii of bodies

        
         call radius(mass(i),rad(1))
         call radius(mass(j),rad(2))
           

         if (vdiff2<vesc2.or.rhom(i,j).le.rmin.or.rhom(i,j).le.sum(rad)) then
         !merging takes place

               rv(i,1:6)=(rv(i,1:6)*mass(i)+rv(j,1:6)*mass(j))/(mass(i)+mass(j))
               mass(i)=mass(i)+mass(j)
            
               write(*,*)' '
               write(*,*)'massive particles',order(i),' and ',order(j),'merged'
               write(32,2003)time+tstart,order(i),order(j),rv(j,1:3),Sqrt(vdiff2/vesc2)

         !relocate location and veolcities
          !merged bodies are being placed at system's barycenter
             !  dumrv(:)=rv(j,:)
               rv(j,:)=rv(body,:)
               ! rv(body,:)=dumrv(:)
               rv(body,:)=0.d0
               
               !relocate masses
            !merging accounted for the merged particle's mass by adding it to the other particle
            !   dummass=mass(j)
               mass(j)=mass(body)
            !   mass(body)=mass(j)
               mass(body)=0.d0
  

          !keep track of the original paritcle number
               dumord=order(j)
               order(j)=order(body)
               order(body)=dumord
         
          !adapt counters for the number of massive particles
               cutbody=cutbody+1
               body=body-1

             
                write(*,*)'number of massive particles reduced to:',body
                
                if(body.le.1.and.m0body.le.0) then
                   call CloseFile
                    write(*,*)'all bodies except for one have left the system. Integration stopped...'
                     STOP
                end if

             end if
         !maybe we dont have to run through all the j-loop


             mcount=mcount-1
             merger(i,j)=.false.
          end if
    
      if(mcount.le.0) then
           !all mergers are accounted for, exit loop
         exit 
      end if

   end do                   
end if
end do mergemassive

bm0c=totbody-cutm0body

!-----------------------------------------------------

mergemml: do i=bm0c,body+cutbody,-1
!merging of a massive with massless body

   mcount=count(merger(i,:))

   if(mcount.ge.1) then

      do j=body,1,-1

         if(merger(i,j)) then
!scattering ? No way -> Harz4 Modell unjustified
 
         !calculate physical radii of bodies        
         call radius(mass(i),rad(1))
         call radius(mass(j),rad(2))


            vrelij(:)=rv(j,4:6)-rv(i,4:6)
            vdiff2 = DOT_PRODUCT(vrelij(:),vrelij(:))
            vesc2 = 2.d0*mass(j)/rhom(i,j)

    
         if (vdiff2<vesc2.or.rhom(i,j).le.rmin.or.rhom(i,j).le.sum(rad)) then
               !merging takes place

  write(*,*)' '    
  write(*,*)'massive particle',order(j),' and massless particle ' ,order(i),'merged'
  write(32,2003)time+tstart,order(i),order(j),rv(j,1:3),Sqrt(vdiff2/vesc2)

         !impact changes positions and velocities of massive particle
               rv(j,1:6)=(rv(i,1:6)*mass(i)+rv(j,1:6)*mass(j))/(mass(i)+mass(j))
               mass(j)=mass(i)+mass(j)

        
         !relocate location and veolcities (exchange coordinates and masses with the 
         !last particle still encompassed by the calculation)
         !merged bodies are being placed at system's barycenter
               
               rv(i,:)=rv(bm0c,:)
               rv(bm0c,:)=0.d0
               
               !masses are being updated
                mass(i)=mass(bm0c)
                mass(bm0c)=0.d0

          !keep track of the original paritcle number
               dumord=order(i)
               order(i)=order(bm0c)
               order(bm0c)=dumord
         
          !adapt counters for the number of massive particles
               cutm0body=cutm0body+1
               m0body=m0body-1
               bm0c=totbody-cutm0body

               write(*,*)'number of massless particles reduced to:',m0body
            end if
 
 !           write(*,*)'particles ',i,'and ',j,' should have merged, but v_rel/v_esc:',Sqrt(vdiff2/vesc2) 
            
                 mcount=mcount-1
                 merger(i,j)=.false.
              end if
         if(mcount.le.0) then
            !all mergers are accounted for, exit loop
            exit 
         end if

      end do
   end if   
end do mergemml

!update barycenter?
!call bakoo(rv(1:bm0c,:),mass(1:bm0c),bm0c)


2003 format (F15.4,5X,I10,5X,I10,5X,4(F15.10,5X))
return
end subroutine scattermerge
!**************************************************************************
subroutine radius(mass,rad)
!calculates physical radius rad [AU] of bodies with mass [Msun]
         implicit none
         real(kind=double),intent(in)::mass
         real(kind=double),intent(inout)::rad
         real(kind=double),parameter::fourpibythree=4.d0*PI/3.d0
         real(kind=double)::density

      !    fourpibythree=4.d0*PI/3.d0

            if (mass.le.1.d-5) then
             !terrestrial planets' mean density ~ 5 [g/cm^3] = 9.25812d6[Msun/AU^3]
                density=9.25812d6
                rad=(mass/density/fourpibythree)**(1.d0/3.d0)
            elseif(mass.gt.1.d-5.and.mass.le.1.d-3) then
               ! Gasgiants mean density~1.4 [g/cm^3] =2.3d6[Msun/AU^3]
                density=2.35661d6
                 rad=(mass/density/fourpibythree)**(1.d0/3.d0)  
            elseif(mass.gt.1.d-3.and.mass.lt.1.d-2) then
              !brown dwarfs have almost all the same radius as jupiter
              ! constant is (r_jup [AU])   
                rad=4.77888d-7
            elseif(mass.ge.1.d-2) then
                !rsun[AU]=0.0465241
               rad=mass**0.78d0*0.0465241d0
            end if

        return
        end subroutine
!*****************************************************************************************
subroutine closeenc(time,totbody,body,bodyp1,cedist,cemlim,rho,rv,mass,order)
!identifies close encounters
    implicit none
 integer(kind=intdble),intent(in)::totbody,body,bodyp1
 real(kind=double),intent(in)::rho(:,:),time,cedist,cemlim
 real(kind=double),intent(in) ::rv(1:totbody,1:6),mass(1:totbody)
 integer(kind=intdble)::order(1:totbody)
 integer(kind=intdble)::i,j   
 real(kind=double)::hillr(1:body),kk
 !Gaussian constant
 kk=k
 
  !calculate Hill's radius for all massive bodies  
  call hillsr(totbody,body,rho,mass,hillr)
    
  
 do i=1,body-1
   if(mass(i).le.cemlim) then
   do j=i+1,body   
    if(rho(j,i).le.cedist) then    
     if(intback) then
     write(unit=33,fmt=*) time+tstart, order(i),order(j),rho(j,i),rv(i,1:3),-rv(i,4:6)*kk,rv(j,1:3),-rv(j,4:6)*kk
     else
     write(unit=33,fmt=*) time+tstart, order(i),order(j),rho(j,i),rv(i,1:3),rv(i,4:6)*kk,rv(j,1:3),rv(j,4:6)*kk
     end if
    end if
   end do
   else
   do j=i+1,body    
     if(rho(j,i).le.hillr(i)) then 
      if(intback) then
     write(unit=33,fmt=*) time+tstart, order(i),order(j),rho(j,i),rv(i,1:3),-rv(i,4:6)*kk,rv(j,1:3),-rv(j,4:6)*kk
     else
     write(unit=33,fmt=*) time+tstart, order(i),order(j),rho(j,i),rv(i,1:3),rv(i,4:6)*kk,rv(j,1:3),rv(j,4:6)*kk
     end if
    end if
   end do
   end if
 end do

 do i=1,body
  if(mass(i).le.cemlim) then
  do j=bodyp1,totbody
   if(rho(j,i).le.cedist) then
    if(intback) then
     write(unit=33,fmt=*)time+tstart, order(i),order(j),rho(j,i),rv(i,1:3),-rv(i,4:6)*kk,rv(j,1:3),-rv(j,4:6)*kk
    else
      write(unit=33,fmt=*)time+tstart, order(i),order(j),rho(j,i),rv(i,1:3),rv(i,4:6)*kk,rv(j,1:3),rv(j,4:6)*kk
    end if
   end if 
  end do
  else 
   do j=bodyp1,totbody
    if(rho(j,i).le.hillr(i)) then
     if(intback) then
     write(unit=33,fmt=*)time+tstart, order(i),order(j),rho(j,i),rv(i,1:3),-rv(i,4:6)*kk,rv(j,1:3),-rv(j,4:6)*kk
     else
      write(unit=33,fmt=*)time+tstart, order(i),order(j),rho(j,i),rv(i,1:3),rv(i,4:6)*kk,rv(j,1:3),rv(j,4:6)*kk
     end if
    end if 
   end do
  endif
 end do
 
end subroutine  
!********************************************************************
subroutine hillsr(totbody,mbody,rho,mass,hillr)
!calculates hills radius with respect to most massive body in the system 
implicit none

 integer(kind=intdble),intent(in) :: totbody,mbody
 real(kind=double),intent(in) :: rho(:,:)
 real(kind=double),intent(in) :: mass(1:totbody)
 
 integer(kind=intdble)::locprim(1:1),i
 real(kind=double)::onethird,rad

 real(kind=double)::hillr(1:mbody)


locprim=maxloc(mass(:))
hillr(:)=0.d0
onethird=1.d0/3.d0


do i=1,mbody
!calculating hill's radius with respect to the heaviest body but preferring any smaller value with respect to another body. 
if(i.eq.locprim(1)) then !when dealing with heaviest body use its physical radius  
   call radius(mass(i),rad)
   hillr(i)=rad
else
   hillr(i)=rho(i,locprim(1))*(onethird*mass(i)/(mass(locprim(1))+mass(i)))**onethird
end if  
end do

end subroutine
!////////////////////////////////////////////////////////////////////////////////////////
end module cutoffmerge_m


