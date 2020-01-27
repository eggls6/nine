module input_m
  use global_m

  implicit none

  public::ReadConfig
  public::ReadIniRV
  public::ReadIniEL


contains

!/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!                                                          Eingabe 
!
!                                         by Siegfried Eggl  20080110 
!------------------------------------------------------------------------------
!********************************************************************
  subroutine ReadConfig(p,input,add,plt)
      implicit none
                       integer(kind=intdble)::i
                       logical::infile
                       type(Problem)::p
                       type(planetoids)::plt
                       !Variablen zur Verwaltung von Input und Outputsteuerung
                       type(additional),intent(out)::add
                       !Input Steuerungsvariable
                        character(len=2),intent(out)::input
                       !Input of planetoid discs
                        character::ptoyn,mergeyn
                       !Output Steuerungsvariable
                       character(len=3),dimension(1:8)::output
                       character(len=1)::dum
                        
                       
        inquire(file="config.inn",exist=infile)

        if(infile) then

        open(unit=21, file="config.inn", status="old", action="read")

        !Nr. des Integrators
        read(unit=21, fmt=*) p%NrInt
                       !start time, end time, output time intervall
        read(unit=21, fmt=*)p%tstart, tend, p%tAusg
        !putting the start time value into a global variable
        tstart=p%tstart
        p%tEnde=abs(tend)-tstart
        if(abs(tend)<tstart) then
          write(*,*) 'integration start time > integration end time' 
       end if 
       
       if(tend.lt.0.d0) then
           intback=.true.
          else
          intback=.false.
        end if  
        
        
        read(unit=21, fmt=*) p%Ns
                       !Anzahl der Koerper
        read(unit=21, fmt=*) p%NrKoerper
                        !one step accuracy
                        read(unit=21,fmt=*)p%eps
                        !koorinate system
                        read(unit=21,fmt=*)p%center
                        !input format (keplerian elements, rv vectors)
                        read(unit=21,fmt=*)input
                        !Order of Lie - Series
                        !read(unit=21,fmt=*)add%lterms
                        add%lterms=12
                        !Minimale Schrittweite
                        read(unit=21,fmt=*)add%lstepmin

                      !Cutoffradius (if = 0, no cuttoff), stop calcualtion when massive (1)/massless (2) particle is beyond cutoff (logical)
                        read(unit=21,fmt=*)p%cutoff,cutoffstop(1:2)
!-------------------------------------------------------------------------------------------------------------------------
!                     !close encounter output?
                       read(unit=21,fmt=*)p%cen,p%cend,p%cemlim 
!--------------------------------------------------------------------------------------------------------
                      !merging on?
                        read(unit=21,fmt=*)mergeyn,p%mergermin,p%mergermax
                        if(mergeyn.eq.'y') then
                           p%merge=.true.
                        else
                           p%merge=.false.
                        endif
!-------------------------------------------------------------------------------------------------------
!                      limit mass for mutual planetesimal interaction (NICE Model)

                       read(unit=21,fmt=*)p%mlimit
!--------------------------------------------------------------------------------------------------------        
                       read(unit=21,fmt=*)showprog
!--------------------------------------------------------------------------------------------------------               
                      !Outputsteuerung
                       read(unit=21,fmt=*)output(:)
                       
                       p%outegy=.false.
                       p%outhele=.false.
                       p%outbele=.false.
                       p%outhco=.false.
                       p%outbco=.false.
                       p%outdef=.false.
                       p%outbinhele=.false.
                       p%outbinhco=.false.
                       p%outbinbco=.false.
                     

                       do i=1,8
                          select case (output(i))
                          case ('egy','en')
                            p% outegy=.true.
                         case('hel')
                           p% outhele=.true.
                         case('bel')
                           p% outbele=.true.   
                         case('hco')
                            p%outhco=.true.
                          case('bco')
                             p%outbco=.true.
                          case('bbc')
                             p%outbinbco=.true.
                          case('bhe')
                             p%outbinhele=.true. 
                           case('bhc')
                             p%outbinhco=.true.   
                         end select
                      end do


                      if (input.eq.'te') then      
                        p%outtroj=.true.
                      else
                        p%outtroj=.false.
                      end if

!no outputfiles given:
                      if (.not.p%outegy.and..not.p% outhele.and..not.p% outbele &
                           .and..not.p%outhco.and..not.p%outbco.and..not. &
                           p%outbinbco.and..not.p%outbinhco.and..not.p%outbinhele) then
                              p%outdef=.true.
                              write(unit=*,fmt=*)' '
                              write(unit=*,fmt=*) 'no specific outputfiles selected'
                       end if
                    
                       read(unit=21,fmt=*)ptoyn

                       select case(ptoyn)
                          case ('y','Y','j','J')
                            plt%pto=.true.
                          case default
                             plt%pto=.false.
                       end select

                        read(unit=21,fmt=*) dum
                        
                        
!generate standard config.inn file
                    else
write(unit=*,fmt='(A)')'No inputfile found (config.inn)!'
write(unit=*,fmt='(A)')'Creating template inputfile (config.inn)...'
                     
 open(unit=25, file="config.inn", status="unknown")
                       
write(unit=25,fmt='(A)', advance='yes')'2                    ....id number of integrator (see end of file)'
write(unit=25,fmt='(A)', advance='no')'4.d2  9.d3  10.d0          '      
write(unit=25,fmt='(A)')'...start time / end time / outputtimes [days]'
write(unit=25,fmt='(A)',advance='no')'1.d-3                '
write(unit=25,fmt='(A)')'...stepsize for Candy and ordercontrolled Lie Series [days]'
write(unit=25,fmt='(A)',advance='no')'4                    '
write(unit=25,fmt='(A)')'...number of bodies in this file (without counting bodies in planetoid rings)'
write(unit=25,fmt='(A)',advance='no')'1.d-14               '
write(unit=25,fmt='(A)')'...error mark (Lie Series: eps = local error)'
write(unit=25,fmt='(A)',advance='no')'b                    '
write(unit=25,fmt='(A)')'...coordinates during integration: barycentric only -> b)'
write(unit=25,fmt='(A)',advance='no')'he                   '
write(unit=25,fmt='(A)')'...Inputformat (rv-vectors: (rv), heliocentric Keplerian elements: (he), Yarkovsky elements: (ye))'
!write(unit=25,fmt='(A)',advance='no')'11                   '
!write(unit=25,fmt='(A)')'...Lie-Itegrator  SW  (=Anzahl der (minimalen) Lieterme) (Lie specific: number of sequence terms) > 2!!'
!write(unit=25,fmt='(A)',advance='no')'8  14              '
!write(unit=25,fmt='(A)')'...Lie-Integrator ST (order window for term control: minimum order / maximum order)'
write(unit=25,fmt='(A)',advance='no')'1d-14                '
write(unit=25,fmt='(A)')'...minimum stepsize [days * Gaussian Gravitational Constant]'
write(unit=25,fmt='(A)',advance='no')'100.d0  .false. .false.                   '
write(unit=25,fmt='(A)',advance='no')'...Cutoff-radius [AU] / cutoff <= 0.d0  means no cuttoff, stop calculation' 
write(unit=25,fmt='(A)')' when massive, massless body beyond cutoff? (.true./.false. .true./false.)'
write(unit=25,fmt='(A)',advance='no')'no   0.001d0              '
write(unit=25,fmt='(A)')'...Close Encounter output file? (yes/no) / miminum pair distance for output[AU]'
write(unit=25,fmt='(A)',advance='no')'no  1.d-10  1.d0     '
write(unit=25,fmt='(A)')'...Merging? (yes/no),   minimum merging radius [AU], maximum merging radius [AU]'
write(unit=25,fmt='(A)',advance='no')'1.d-8     '
write(unit=25,fmt='(A)')'...limit mass for mutual planetesimal interaction [Msun] (NICE Model)'
write(unit=25,fmt='(A)',advance='no')'.true.       '
write(unit=25,fmt='(A)')'...show progress in terminal'
write(unit=25,fmt='(A)',advance='no')'bco, hel,en          ...Outputfiles(en=energy conservation/'
write(unit=25,fmt='(A)',advance='no')'hel=heliocentric elements/hco=heliocentric coordinates/bco= barycentric coordinates/'
write(unit=25,fmt='(A)')'bbc=barycentric coordinates binary format/bhe=heliocentric elements binary output)'
write(unit=25,fmt='(A)',advance='no')'no                   '  
write(unit=25,fmt='(A)')'... Planetoid rings? (yes/no)'
write(unit=25,fmt='(A)')'0.d0	0.d0	0.d0	0.d0	0.d0	0.d0	1.d0'
write(unit=25,fmt='(A)')'5.20336301d0	0.04839266d0	1.30530d0	274.1977d0	100.55615d0	19.65053d0 1.8986d-3'
write(unit=25,fmt='(A)')'9.53707032d0	0.05415060d0	2.48446d0	338.7169d0	113.71504d0	317.51238d0 5.6846d-4'
write(unit=25,fmt='(A)')'1.d0           0.d0            0.d0            0.d0            0.d0            0.d0        0.d0'
write(unit=25,fmt='(A)')''
write(unit=25,fmt='(A)')''
write(unit=25,fmt='(A)')'Input format heliocentric Keplerian elements:'
write(unit=25,fmt='(A)',advance='no')'a   e   i  argument of perihelion (omega)  argument of'
write(unit=25,fmt='(A)')'the ascending node (Omega)  mean anomaly (M) mass (Solar masses)'
write(unit=25,fmt='(A)')''
write(unit=25,fmt='(A)')'Input format heliocentric Keplerian elements:'
write(unit=25,fmt='(A)')'rx   ry   rz    vx    vy    vz     mass (Solar masses)'
write(unit=25,fmt='(A)')''
write(unit=25,fmt='(A)',advance='no')'Input format Yarkovsky elements:'
write(unit=25,fmt='(A)',advance='no')'semi-major axis (a) [AU] eccentricity (e) [] inclination (i) [deg] '
write(unit=25,fmt='(A)',advance='no')' argument of perihelion (omega) [deg]'
write(unit=25,fmt='(A)',advance='no')'argument of the ascending node (Omega) [deg]  mean anomaly (M) [deg]  mass [Solar masses]  '
write(unit=25,fmt='(A)',advance='no')' prograde (0) / retrograde (180) rotation '
write(unit=25,fmt='(A)',advance='no')' thermal capacity [J/kg/K]'
write(unit=25,fmt='(A)',advance='no')'  k_0 and k_1 parameters of the surface thermalconductivity [K(T) = k_0 + k_1 T_av^3]   '
write(unit=25,fmt='(A)',advance='no')'density of the surface layer [kg/m^3]  radius of the body [km]  rotation frequency [Hz]  '
write(unit=25,fmt='(A)')'surface absorptivity=emissivity  bulk density [kg/m^3]'
write(unit=25,fmt='(A)')''
write(unit=25,fmt='(A)')'Integrators:'
write(unit=25,fmt='(A)')''
write(unit=25,fmt='(A)')'0......Candy (symplectic 4th order, fixed stepsize)'
write(unit=25,fmt='(A)')'1......Yoshida (symplectic & symmetric 8th order, fixed stepsize)'
write(unit=25,fmt='(A)')'2......standard Lie series Integrator with adaptive stepsize control (LieSW)'
write(unit=25,fmt='(A)')'3..... Bulirsch Stoer with adaptive stepsize(Mercury6)'
write(unit=25,fmt='(A)')'4..... Gauss Radau with adaptive stepsize'
write(unit=25,fmt='(A)')'5..... Gauss Radau with adaptive stepsize and GR (spheric heliocentric metric)'
!write(unit=25,fmt='(A)')'6..... Gauss Radau with adaptive stepsize, GR (spheric heliocentric metric) and Yarkovsky Thermal effect'
write(unit=25,fmt='(A)')'6..... Gauss Radau with adaptive stepsize, GR (EIH barycentric metric) '
!write(unit=25,fmt='(A)')'8..... Gauss Radau with adaptive stepsize, GR (EIH barycentric metric) and Yarkovsky Thermal effect'

   close(25)
   write(unit=*,fmt='(A)')'...done'
   write(unit=*,fmt='(A)')'Launch the program once again, please!'
   STOP
   
end if

                      return
  end subroutine ReadConfig
!*************************************************************  
  subroutine skipcomments (unit, line)
   integer, intent(in):: unit
   character*80, intent(inout):: line

10 continue
   read (unit, *) line
   if (line(1:1) .eq. '#') goto 10
   return
end subroutine skipcomments
  
!*************************************************************
 subroutine ReadPto(plt)
 implicit none
        logical::infile
         type(planetoids)::plt
         integer::i
         character::dumchar

 inquire(file="planetoid.inn",exist=infile)

 if(infile) then

    open(unit=22, file="planetoid.inn", status="old", action="read")
 
 else 
!create example planetoid input file

     open(unit=22, file="planetoid.inn", status="new")

write(unit=22,fmt='(A)')'2            ...total number of heavy bodies with planetoid rings'
write(unit=22,fmt='(A)', advance='no')'1 3          ... which bodies should the planetoid rings be around' 
write(unit=22,fmt='(A)')'(see input file: e.g. 1 3 6) (0: barycenter)'            
write(unit=22,fmt='(A)')'500 1000      ... number of planetoids of respective ring'
write(unit=22,fmt='(A)')'0.d0 0.d0  ... horizontal orientation angles of respective ring to referrence plane [deg]'
write(unit=22,fmt='(A)')'0.d0 0.d0  ... vertical orientation angles of respective ring to referrence plane [deg]'
write(unit=22,fmt='(A)')'3.d0 1.d0      ... minimum initial separation of planetoids [r_hill]  '
write(unit=22,fmt='(A)')'50 50           ... radial grid size (integer: 50 default)  '
!***************distribution parameters for elements and mass*******************
write(unit=22,fmt='(A)')'!********** distribution parameters for semi-major axis a ***************'
write(unit=22,fmt='(A)',advance='no')'0    1      ... distribution function for a of planetoids ' 
write(unit=22,fmt='(A)')                         '(0...c1*(a-a0)^c2; 1...exp(c1*(a-a0)^c2)'
write(unit=22,fmt='(A)')'1.d0  1.d0    ... distribution parameter c1'
write(unit=22,fmt='(A)')'-1.5d0 2.d0    ... distribution parameter c2'
write(unit=22,fmt='(A)')'0.5d0 0.5d0    ... distribution parameter a0'
write(unit=22,fmt='(A)')'0.5d0 0.3d0  ... lower limit of a of respective ring to its central body [AU]'
write(unit=22,fmt='(A)')'1.d0 0.7d0   ... upper limit of a of respective ring to its central body [AU]'

write(unit=22,fmt='(A)')'!********** distribution parameters for numeric eccentricity e ***************'
write(unit=22,fmt='(A)',advance='no')'0    1      ... distribution function for e of planetoids ' 
write(unit=22,fmt='(A)')                         '(0...c1*(e-e0)^c2; 1...exp(c1*(e-e0)^c2)'
write(unit=22,fmt='(A)')'1.d0  1.d0    ... distribution parameter c1'
write(unit=22,fmt='(A)')'0.d0 0.d0    ... distribution parameter c2'
write(unit=22,fmt='(A)')'0.d0 0.d0    ... distribution parameter e0'
write(unit=22,fmt='(A)')'0.1d0 0.3d0  ... lower limit of e of respective ring []'
write(unit=22,fmt='(A)')'0.5d0 0.4d0   ... upper limit of e of respective ring []'

write(unit=22,fmt='(A)')'!********** distribution parameters for inclination i ***************'
write(unit=22,fmt='(A)',advance='no')'0    0      ... distribution function for i of planetoids ' 
write(unit=22,fmt='(A)')                         '(0...c1*(i-i0)^c2; 1...exp(c1*(r-r0)^c2)'
write(unit=22,fmt='(A)')'1.d0  1.d0      ... distribution parameter c1'
write(unit=22,fmt='(A)')'0.d0 0.d0       ... distribution parameter c2'
write(unit=22,fmt='(A)')'0.d0 0.d0       ... distribution parameter i0'
write(unit=22,fmt='(A)')'0.d0 0.d0  ... lower limit of i of respective ring [deg]'
write(unit=22,fmt='(A)')'20.d0 0.5d0   ... upper limit of i of respective ring [deg]'

write(unit=22,fmt='(A)')'!********** distribution parameters for argument of pericenter omega ***************'
write(unit=22,fmt='(A)',advance='no')'0    0      ... distribution function for omega of planetoids ' 
write(unit=22,fmt='(A)')                         '(0...c1*(omega-omega0)^c2; 1...exp(c1*(omega-omega0)^c2)'
write(unit=22,fmt='(A)')'1.d0  1.d0      ... distribution parameter c1'
write(unit=22,fmt='(A)')'0.d0 0.d0       ... distribution parameter c2'
write(unit=22,fmt='(A)')'0.d0 0.d0       ... distribution parameter omega0'
write(unit=22,fmt='(A)')'0.d0 0.3d0  ... lower limit of omega of respective ring [deg]'
write(unit=22,fmt='(A)')'90.d0 0.4d0   ... upper limit of omega of respective ring [deg]'

write(unit=22,fmt='(A)')'!********** distribution parameters for argument of the ascending node Omega *******'
write(unit=22,fmt='(A)',advance='no')'0    0      ... distribution function for Omega of planetoids ' 
write(unit=22,fmt='(A)')                         '(0...c1*(Omega-Omega0)^c2; 1...exp(c1*(Omega-Omega0)^c2)'
write(unit=22,fmt='(A)')'1.d0  1.d0      ... distribution parameter c1'
write(unit=22,fmt='(A)')'0.d0 0.d0       ... distribution parameter c2'
write(unit=22,fmt='(A)')'0.d0 0.d0       ... distribution parameter Omega0'
write(unit=22,fmt='(A)')'0.d0 0.3d0  ... lower limit of Omega of respective ring [deg]'
write(unit=22,fmt='(A)')'360.d0 0.4d0   ... upper limit of Omega of respective ring [deg]'

write(unit=22,fmt='(A)')'!********** distribution parameters for mean anomaly M ***************'
write(unit=22,fmt='(A)',advance='no')'0    0      ... distribution function for M of planetoids ' 
write(unit=22,fmt='(A)')                         '(0...c1*(M-M0)^c2; 1...exp(c1*(M-M0)^c2)'
write(unit=22,fmt='(A)')'1.d0  1.d0      ... distribution parameter c1'
write(unit=22,fmt='(A)')'0.d0 0.d0       ... distribution parameter c2'
write(unit=22,fmt='(A)')'0.d0 0.d0       ... distribution parameter M0'
write(unit=22,fmt='(A)')'0.d0 0.3d0  ... lower limit of M of respective ring [deg]'
write(unit=22,fmt='(A)')'360.d0 0.4d0   ... upper limit of M of respective ring [deg]'

write(unit=22,fmt='(A)')'!********** distribution parameters for planetoid masses m ***************'
write(unit=22,fmt='(A)',advance='no')'0    0      ... distribution function for m of planetoids  ' 
write(unit=22,fmt='(A)')         ' (2...equidistribution of log10(m); 3... equidistribution of m)'
write(unit=22,fmt='(A)')'0.d0 1.d-8  ... lower limit of m of respective ring [M_sun]'
write(unit=22,fmt='(A)')'1.d-5 1.d-7   ... upper limit of m of respective ring [M_sun] '
write(unit=22,fmt='(A)',advance='no')'1.75d-3 3.5d-5    ...maximum mass of respecitve ring '
write(unit=22,fmt='(A)')'(when reached, the rest of particles will have mass=0) [M_sun]'
     close(22)

     write(unit=*,fmt='(A)')'samle planetoid file created'
     write(unit=*,fmt='(A)')'please adapt file "planetoid.inn" and restart program'
     STOP


 end if   
     

 read(unit=22,fmt=*)plt%n
 
 allocate(plt%bodyn(plt%n,1:2), & 
    plt%distinn(plt%n,1:7),plt%distout(plt%n,1:7),plt%incl(plt%n),plt%phi(plt%n), &
    plt%exp_or_r(plt%n,1:7),plt%gridsize(plt%n),plt%a(plt%n,1:7),plt%b(plt%n,1:7), &
    plt%xrhill(plt%n),plt%r0(plt%n,1:7),plt%massmax(plt%n))

 read(unit=22,fmt=*)plt%bodyn(:,1)
 read(unit=22,fmt=*)plt%bodyn(:,2)
 read(unit=22,fmt=*)plt%incl(:)
 read(unit=22,fmt=*)plt%phi(:)
 read(unit=22,fmt=*)plt%xrhill(:)
 read(unit=22,fmt=*)plt%gridsize(:)
do i=1,6
 read(unit=22,fmt=*)dumchar
 read(unit=22,fmt=*)plt%exp_or_r(:,i)
 read(unit=22,fmt=*)plt%a(:,i)
 read(unit=22,fmt=*)plt%b(:,i)
 read(unit=22,fmt=*)plt%r0(:,i)
 read(unit=22,fmt=*)plt%distinn(:,i)
 read(unit=22,fmt=*)plt%distout(:,i)
end do
 read(unit=22,fmt=*)dumchar
 read(unit=22,fmt=*)plt%exp_or_r(:,7)
 read(unit=22,fmt=*)plt%distinn(:,7)
 read(unit=22,fmt=*)plt%distout(:,7)
 read(unit=22,fmt=*)plt%massmax(:)

 close(22)

 end subroutine ReadPto
!*************************************************************
 subroutine ReadIniRV(rv, m, NrK)
        implicit none
        integer(kind=intdble),intent(in)::NrK
        integer(kind=intdble)::j
        real(kind=double),dimension(:,:),intent(out)::rv
        real(kind=double),dimension(:),intent(out)::m
        

        !Anfangsdaten
        do j=1,NrK
                read(unit=21, fmt=*) rv(j,1),rv(j,2),rv(j,3), rv(j,4),rv(j,5),rv(j,6),m(j),&
                bname(j)
        end do

        close(unit=21)

        return
end subroutine ReadIniRV
!*************************************************************
 subroutine ReadIniEL(ele,m,NrK)
           implicit none         
                        integer(kind=intdble),intent(in)::NrK
                        integer(kind=intdble)::j
                        type(kepele),dimension(1:NrK),intent(out)::ele
                        real(kind=double),dimension(1:NrK),intent(out)::m                       

                        do j=1,NrK
                           read(unit=21,fmt=*)ele(j)%a,ele(j)%e,ele(j)%i,ele(j)%kom,ele(j)%gom,ele(j)%man,m(j),&
                           bname(j)
                        end do                    
                        close(unit=21)
                        return
 end subroutine ReadIniEL

!*************************************************************
 subroutine ReadIniYEL(ele,m,NrK)
 use global_m
 use yark_m
   implicit none         
   type(kepele),dimension(:),intent(out)::ele
   real(kind=double),dimension(:),intent(out)::m
   integer(kind=intdble),intent(in)::NrK
   integer(kind=intdble)::j
  
   
   do j=1,NrK
    read(unit=21,fmt=*)ele(j)%a,ele(j)%e,ele(j)%i,ele(j)%kom,ele(j)%gom, &
                         ele(j)%man,m(j),yark_p(j,2:10), bname(j)
   ! call init_yark(yark_p(j,1:3))
    !given in km in the input file
    yark_p(j,7)=yark_p(j,7)*1000._double     
    yark_p(j,1)=0._double 
   end do         
                                           
   close(unit=21)
   return
 end subroutine ReadIniYEL
!*************************************************************
 subroutine ReadIniELT(ele,m,NrK)
 use global_m
   implicit none         
   type(kepele),dimension(:),intent(out)::ele
   real(kind=double),dimension(:),intent(out)::m
   integer(kind=intdble),intent(in)::NrK
   integer(kind=intdble)::j
 
   
   do j=1,NrK
    read(unit=21,fmt=*)ele(j)%a,ele(j)%e,ele(j)%i,ele(j)%kom,ele(j)%gom,&
                       ele(j)%man,m(j),trojan(j),bname(j)
   end do
                       
   close(unit=21)
   return
 end subroutine ReadIniELT
end module input_m
