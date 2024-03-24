!##############################################################################
Subroutine trunc_3d_vars (varn,nzp,nxp,nyp,a)

implicit none

character(len=32) :: varn
integer :: nzp,nxp,nyp
real :: a(nxp,nyp,nzp)
real :: t_mbud1,t_mbud2,t_mixr,t_numc,t_aeroc

!Common rounding thresholds for reducing precision and 
!enabling better data compression in LITE files
t_mbud1 = 1e8
t_mbud2 = 1e9
t_mixr  = 1e6
t_numc  = 1e1
t_aeroc = 1e1

!Application of rounding based on native units of variables within
!model runtime (e.g. m/s, kg/kg, #/kg, K, etc)

!Dynamic and general thermodynamic variables
if(    trim(varn)=='THETA')then
a = anint(a*1e3)/1e3
elseif(trim(varn)=='THP')then
a = anint(a*1e3)/1e3
elseif(trim(varn)=='UP')then
a = anint(a*1e2)/1e2
elseif(trim(varn)=='VP')then
a = anint(a*1e2)/1e2
elseif(trim(varn)=='WP')then
a = anint(a*1e3)/1e3
elseif(trim(varn)=='PP')then
a = anint(a*1e3)/1e3
elseif(trim(varn)=='PI')then
a = anint(a*1e3)/1e3
elseif(trim(varn)=='DN0')then
a = anint(a*1e3)/1e3
!Vapor and hydrometeor mixing ratio variables (kg/kg)
elseif(trim(varn)=='RTP')then
a = anint(a*t_mixr)/t_mixr
elseif(trim(varn)=='RV')then
a = anint(a*t_mixr)/t_mixr
elseif(trim(varn)=='RCP')then
a = anint(a*t_mixr)/t_mixr
elseif(trim(varn)=='RRP')then
a = anint(a*t_mixr)/t_mixr
elseif(trim(varn)=='RPP')then
a = anint(a*t_mixr)/t_mixr
elseif(trim(varn)=='RSP')then
a = anint(a*t_mixr)/t_mixr
elseif(trim(varn)=='RAP')then
a = anint(a*t_mixr)/t_mixr
elseif(trim(varn)=='RGP')then
a = anint(a*t_mixr)/t_mixr
elseif(trim(varn)=='RHP')then
a = anint(a*t_mixr)/t_mixr
elseif(trim(varn)=='RDP')then
a = anint(a*t_mixr)/t_mixr
!Hydrometeor number concentrations (#/kg)
elseif(trim(varn)=='CCP')then
a = anint(a*t_numc)/t_numc
elseif(trim(varn)=='CRP')then
a = anint(a*t_numc)/t_numc
elseif(trim(varn)=='CPP')then
a = anint(a*t_numc)/t_numc
elseif(trim(varn)=='CSP')then
a = anint(a*t_numc)/t_numc
elseif(trim(varn)=='CAP')then
a = anint(a*t_numc)/t_numc
elseif(trim(varn)=='CGP')then
a = anint(a*t_numc)/t_numc
elseif(trim(varn)=='CHP')then
a = anint(a*t_numc)/t_numc
elseif(trim(varn)=='CDP')then
a = anint(a*t_numc)/t_numc
!Radiation variables
elseif(trim(varn)=='FTHRD')then
a = anint(a*1e5)/1e5
elseif(trim(varn)=='SWUP')then
a = anint(a*1e2)/1e2
elseif(trim(varn)=='SWDN')then
a = anint(a*1e2)/1e2
elseif(trim(varn)=='LWUP')then
a = anint(a*1e2)/1e2
elseif(trim(varn)=='LWDN')then
a = anint(a*1e2)/1e2
elseif(trim(varn)=='BEXT')then
a = anint(a*1e4)/1e4
!Aerosol number concentrations
elseif(trim(varn)=='CN1NP')then
a = anint(a*t_aeroc)/t_aeroc
elseif(trim(varn)=='CN2NP')then
a = anint(a*t_aeroc)/t_aeroc
elseif(trim(varn)=='CN3NP')then
a = anint(a*t_aeroc)/t_aeroc
elseif(trim(varn)=='MD1NP')then
a = anint(a*t_aeroc)/t_aeroc
elseif(trim(varn)=='MD2NP')then
a = anint(a*t_aeroc)/t_aeroc
elseif(trim(varn)=='ABC1NP')then
a = anint(a*t_aeroc)/t_aeroc
elseif(trim(varn)=='ABC2NP')then
a = anint(a*t_aeroc)/t_aeroc
elseif(trim(varn)=='SALT_FILM_NP')then
a = anint(a*t_aeroc)/t_aeroc
elseif(trim(varn)=='SALT_JET_NP')then
a = anint(a*t_aeroc)/t_aeroc
elseif(trim(varn)=='SALT_SPUM_NP')then
a = anint(a*t_aeroc)/t_aeroc
elseif(trim(varn)=='REGEN_AERO1_NP')then
a = anint(a*t_aeroc)/t_aeroc
elseif(trim(varn)=='REGEN_AERO2_NP')then
a = anint(a*t_aeroc)/t_aeroc
!Level-1 mass process rates
elseif(trim(varn)=='NUCCLDRT')then
a = anint(a*t_mbud1)/t_mbud1
elseif(trim(varn)=='CLD2RAINT')then
a = anint(a*t_mbud1)/t_mbud1
elseif(trim(varn)=='ICE2RAINT')then
a = anint(a*t_mbud1)/t_mbud1
elseif(trim(varn)=='NUCICERT')then
a = anint(a*t_mbud1)/t_mbud1
elseif(trim(varn)=='VAPLIQT')then
a = anint(a*t_mbud1)/t_mbud1
elseif(trim(varn)=='VAPICET')then
a = anint(a*t_mbud1)/t_mbud1
elseif(trim(varn)=='EVAPLIQT')then
a = anint(a*t_mbud1)/t_mbud1
elseif(trim(varn)=='EVAPICET')then
a = anint(a*t_mbud1)/t_mbud1
elseif(trim(varn)=='FREEZINGT')then
a = anint(a*t_mbud1)/t_mbud1
elseif(trim(varn)=='MELTINGT')then
a = anint(a*t_mbud1)/t_mbud1
elseif(trim(varn)=='MELTICET')then
a = anint(a*t_mbud1)/t_mbud1
elseif(trim(varn)=='RIMECLDT')then
a = anint(a*t_mbud1)/t_mbud1
elseif(trim(varn)=='RAIN2ICET')then
a = anint(a*t_mbud1)/t_mbud1
elseif(trim(varn)=='AGGREGATET')then
a = anint(a*t_mbud1)/t_mbud1
!WP & Temperature process rates
elseif(trim(varn)=='LATHEATVAPT')then
a = anint(a*1e5)/1e5
elseif(trim(varn)=='LATHEATFRZT')then
a = anint(a*1e5)/1e5
elseif(trim(varn)=='LATHEATVAP')then
a = anint(a*1e5)/1e5
elseif(trim(varn)=='LATHEATFRZ')then
a = anint(a*1e5)/1e5
elseif(trim(varn)=='WP_BUOY_THETA')then
a = anint(a*1e4)/1e4
elseif(trim(varn)=='WP_BUOY_COND')then
a = anint(a*1e4)/1e4
elseif(trim(varn)=='WP_ADVDIF')then
a = anint(a*1e4)/1e4
!Level-2 mass process rates
elseif(trim(varn)=='INUCHOMRT')then
a = anint(a*t_mbud2)/t_mbud2
elseif(trim(varn)=='INUCCONTRT')then
a = anint(a*t_mbud2)/t_mbud2
elseif(trim(varn)=='INUCIFNRT')then
a = anint(a*t_mbud2)/t_mbud2
elseif(trim(varn)=='INUCHAZRT')then
a = anint(a*t_mbud2)/t_mbud2
elseif(trim(varn)=='VAPCLDT')then
a = anint(a*t_mbud2)/t_mbud2
elseif(trim(varn)=='VAPRAINT')then
a = anint(a*t_mbud2)/t_mbud2
elseif(trim(varn)=='VAPPRIST')then
a = anint(a*t_mbud2)/t_mbud2
elseif(trim(varn)=='VAPSNOWT')then
a = anint(a*t_mbud2)/t_mbud2
elseif(trim(varn)=='VAPAGGRT')then
a = anint(a*t_mbud2)/t_mbud2
elseif(trim(varn)=='VAPGRAUT')then
a = anint(a*t_mbud2)/t_mbud2
elseif(trim(varn)=='VAPHAILT')then
a = anint(a*t_mbud2)/t_mbud2
elseif(trim(varn)=='VAPDRIZT')then
a = anint(a*t_mbud2)/t_mbud2
elseif(trim(varn)=='EVAPCLDT')then
a = anint(a*t_mbud2)/t_mbud2
elseif(trim(varn)=='EVAPRAINT')then
a = anint(a*t_mbud2)/t_mbud2
elseif(trim(varn)=='EVAPPRIST')then
a = anint(a*t_mbud2)/t_mbud2
elseif(trim(varn)=='EVAPSNOWT')then
a = anint(a*t_mbud2)/t_mbud2
elseif(trim(varn)=='EVAPAGGRT')then
a = anint(a*t_mbud2)/t_mbud2
elseif(trim(varn)=='EVAPGRAUT')then
a = anint(a*t_mbud2)/t_mbud2
elseif(trim(varn)=='EVAPHAILT')then
a = anint(a*t_mbud2)/t_mbud2
elseif(trim(varn)=='EVAPDRIST')then
a = anint(a*t_mbud2)/t_mbud2
elseif(trim(varn)=='MELTPRIST')then
a = anint(a*t_mbud2)/t_mbud2
elseif(trim(varn)=='MELTSNOWT')then
a = anint(a*t_mbud2)/t_mbud2
elseif(trim(varn)=='MELTAGGRT')then
a = anint(a*t_mbud2)/t_mbud2
elseif(trim(varn)=='MELTGRAUT')then
a = anint(a*t_mbud2)/t_mbud2
elseif(trim(varn)=='MELTHAILT')then
a = anint(a*t_mbud2)/t_mbud2
elseif(trim(varn)=='RIMECLDSNOWT')then
a = anint(a*t_mbud2)/t_mbud2
elseif(trim(varn)=='RIMECLDAGGRT')then
a = anint(a*t_mbud2)/t_mbud2
elseif(trim(varn)=='RIMECLDGRAUT')then
a = anint(a*t_mbud2)/t_mbud2
elseif(trim(varn)=='RIMECLDHAILT')then
a = anint(a*t_mbud2)/t_mbud2
elseif(trim(varn)=='RAIN2PRT')then
a = anint(a*t_mbud2)/t_mbud2
elseif(trim(varn)=='RAIN2SNT')then
a = anint(a*t_mbud2)/t_mbud2
elseif(trim(varn)=='RAIN2AGT')then
a = anint(a*t_mbud2)/t_mbud2
elseif(trim(varn)=='RAIN2GRT')then
a = anint(a*t_mbud2)/t_mbud2
elseif(trim(varn)=='RAIN2HAT')then
a = anint(a*t_mbud2)/t_mbud2
elseif(trim(varn)=='AGGRSELFPRIST')then
a = anint(a*t_mbud2)/t_mbud2
elseif(trim(varn)=='AGGRSELFSNOWT')then
a = anint(a*t_mbud2)/t_mbud2
elseif(trim(varn)=='AGGRPRISSNOWT')then
a = anint(a*t_mbud2)/t_mbud2
endif

return
END SUBROUTINE trunc_3d_vars

!##############################################################################
Subroutine trunc_2d_vars (varn,nxp,nyp,a)

implicit none

character(len=32) :: varn
integer :: nxp,nyp
real :: a(nxp,nyp)
real :: t_accp,t_pcpr

!Common rounding thresholds for reducing precision and 
!enabling better data compression in LITE files
t_accp = 1e3
t_pcpr = 1e7
!Application of rounding based on native units of variables within
!model runtime (e.g. mm/sec, etc)

!Dynamic and general thermodynamic variables
if(    trim(varn)=='ACCPR'.or. &
       trim(varn)=='ACCPP'.or. &
       trim(varn)=='ACCPS'.or. &
       trim(varn)=='ACCPA'.or. &
       trim(varn)=='ACCPG'.or. &
       trim(varn)=='ACCPH'.or. &
       trim(varn)=='ACCPD')then
 a = anint(a*t_accp)/t_accp
elseif(trim(varn)=='PCPRR'.or. &
       trim(varn)=='PCPRP'.or. &
       trim(varn)=='PCPRS'.or. &
       trim(varn)=='PCPRA'.or. &
       trim(varn)=='PCPRG'.or. &
       trim(varn)=='PCPRH'.or. &
       trim(varn)=='PCPRD')then
 a = anint(a*t_pcpr)/t_pcpr
endif

return
END SUBROUTINE trunc_2d_vars

!##############################################################################
Subroutine rams_mm (indata,ni1,omin,omax)

implicit none

integer :: ni1
real :: indata(ni1),omin,omax
integer :: i

omax=indata(1)
omin=indata(1)
   do i=2,ni1
      omax=max(indata(i),omax)
      omin=min(indata(i),omin)
   enddo

return
END SUBROUTINE rams_mm

!##############################################################################
real Function walltime (wstart)

implicit none

real :: wstart
integer :: ii,ir

!Use lowercase "call" since this is a system call
call system_clock (count=ii,count_rate=ir)
walltime = float(ii)/float(ir) - wstart

return
END FUNCTION walltime

!##############################################################################
Subroutine rearrange (nzp,nxp,nyp,a,b)

implicit none

integer :: nzp,nxp,nyp
real :: a(nzp,nxp,nyp),b(nxp,nyp,nzp)
integer :: k,i,j

do i=1,nxp
   do j=1,nyp
      do k=1,nzp
         b(i,j,k)=a(k,i,j)
      enddo
   enddo
enddo

return
END SUBROUTINE rearrange

!##############################################################################
Subroutine unarrange (nzp,nxp,nyp,a,b)

implicit none

integer :: nzp,nxp,nyp
real :: a(nxp,nyp,nzp),b(nzp,nxp,nyp)
integer :: k,i,j

do i=1,nxp
   do j=1,nyp
      do k=1,nzp
         b(k,i,j)=a(i,j,k)
      enddo
   enddo
enddo

return
END SUBROUTINE unarrange

!##############################################################################
Subroutine rearrange_p (n2,n3,n4,n5,a,b)

implicit none

integer :: n2,n3,n4,n5
real :: a(n4,n2,n3,n5),b(n2,n3,n4,n5)
integer :: i,j,k,ip

do ip = 1,n5
   do k = 1,n4
      do j = 1,n3
         do i = 1,n2
            b(i,j,k,ip) = a(k,i,j,ip)
         enddo
      enddo
   enddo
enddo

return
END SUBROUTINE rearrange_p

!##############################################################################
Subroutine unarrange_p (n2,n3,n4,n5,a,b)

implicit none

integer :: n2,n3,n4,n5
real :: a(n2,n3,n4,n5),b(n4,n2,n3,n5)
integer :: i,j,k,ip

do ip = 1,n5
   do k = 1,n4
      do j = 1,n3
         do i = 1,n2
            b(k,i,j,ip) = a(i,j,k,ip)
         enddo
      enddo
   enddo
enddo

return
END SUBROUTINE unarrange_p

!##############################################################################
Subroutine makefnam (fname,prefix,tinc,iyr,imn,idy,itm,type,post,fmt)

! creates standard timestamped filename

implicit none

integer :: iyr,imn,idy,itm,oyr,omn,ody,otm,ib1,ib2
character(len=*) :: fname,prefix,post,fmt
character(len=1) :: type
real :: tinc
character(len=40) :: dstring

!   print*,iyr,imn,idy,itm,tinc
if(tinc == 0.) then
   oyr=iyr ; omn=imn ; ody=idy ; otm=itm
else
   CALL date_add_to (iyr,imn,idy,itm,tinc,'s',oyr,omn,ody,otm)
!   print*,oyr,omn,ody,otm
endif

write(dstring,100) '-',type,'-',oyr,'-',omn,'-',ody,'-',otm
100 format(3a1,i4.4,a1,i2.2,a1,i2.2,a1,i6.6)

ib1=len_trim(prefix)
fname=prefix(1:ib1)//dstring(1:20)
if (post(1:1) /= '$') then
   ib1=len_trim(fname)
   ib2=len_trim(post)
   fname=fname(1:ib1)//'-'//post(1:ib2)
endif
ib1=len_trim(fname)
fname=fname(1:ib1)//'.'//trim(fmt)

return
END SUBROUTINE makefnam

!##############################################################################
Subroutine rams_f_open (iunit,filenm,formt,stat,act,iclob)

! replaces old jclopen and jclget
! files are overwritten unless iclob (ICLOBBER) set to 1

use mem_grid, only:print_msg

implicit none

integer :: iunit,iclob
character(len=*) :: filenm,formt,stat,act
logical :: exans

!print*,'filenm,formt,stat1=',filenm,formt,stat

inquire(FILE=filenm,EXIST=exans)

if(exans.and.iclob.eq.0.and.  &
     (act(1:4).eq.'WRIT'.or.act(1:4).eq.'writ') .and. print_msg) then
   print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   print*,'!!!   trying to open file name :'
   print*,'!!!       ',filenm
   print*,'!!!   but it already exists. run is ended.'
   print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   stop 'rams_f_open - exists'
endif

!print*,'filenm,formt,stat2=',filenm(1:len_trim(filenm)),formt,stat
open(iunit,STATUS=stat,FILE=trim(filenm),FORM=formt)
if(print_msg) print*,'F_open - ',trim(filenm)

return
END SUBROUTINE rams_f_open
