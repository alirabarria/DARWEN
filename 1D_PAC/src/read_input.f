
       subroutine read_input

       include 'pac.common'
       integer i,j

       write(*,1000)inpfile
       i=nmaxcharfile
       do while (inpfile(i:i).ne.'.')
          i=i-1
       enddo
       ncharmodel=i-1
       namemodel(1:ncharmodel)=inpfile(1:ncharmodel)

c      ...read input file
       open(unit=2,err=100,file=inpfile,status='old')
       read(2,*)specfile
       read(2,*)zabfile
       read(2,*)nreacfile,i,warnreac
       if(warnreac.eq.0)write(*,1010)
       if(i.ne.0.and.i.ne.1)stop' E- iondip code must be 0/1'
       if(i.eq.1)iondip=.true.
       do i=1,nreacfile
          read(2,*)reacfile(i)
       enddo
       i=nmaxcharfile
       do while (reacfile(nreacfile)(i:i).ne.'.')
          i=i-1
       enddo
       if(reacfile(nreacfile)(i+1:i+5).eq.'therm')then
          therm=.true.
          thermfile(1:i+5)=reacfile(nreacfile)(1:i+5)
          nreacfile=nreacfile-1
          if(nreacfile.lt.1)stop ' E- No reaction files'
       endif
       read(2,*,err=110)eddyfile,i
       if(i.eq.1)longeddy=.true.
       read(2,*)starfile
       read(2,*)iphotfile
       read(2,*)pmax,pmin,nz
       pmax=pmax*1.0d6                           ! bar -> dyn cm-2
       pmin=pmin*1.0d6                           ! bar -> dyn cm-2
       write(*,1001)pmax,pmin,nz
       if(nz.gt.nmaxz)stop ' E- Too many requested z points'
       read(2,*)rstar
       rstar=rstar*rsun                          ! R(Sun)   -> cm
       read(2,*)rplanet,mplanet
       rplanet=rplanet*rearth                    ! R(Earth) -> cm
       mplanet=mplanet*mearth                    ! m(Earth) -> g
       read(2,*)dist,zenithangle,albedoplanet
       dist=dist*au                              ! AU       -> cm
       zenithangle=zenithangle*2.0d0*pi/360.0d0  ! deg      -> rad
       if(albedoplanet.lt.0.0d0.or.albedoplanet.gt.1.0d0)
     # stop ' E- Invalid planet surface albedo'
       read(2,*)i,j
       if((i.ne.0.and.i.ne.1).or.(j.ne.0.and.j.ne.1))
     # stop ' E- Diffusion, Photochemistry must be set at 0/1'
       if(i.eq.1)diffusion=.true.
       if(j.eq.1)photochem=.true.
       read(2,*)nummethod
       if(nummethod.eq.-1)then                   ! solid body rotation activated
          sbrotation=.true.
          read(2,*)lptfile                       ! longitude-pressure-temperature file
          read(2,*)nlong,nperiod                 ! No longitudes per period, No rotation periods
          nlong=nlong+1                          ! to account for the first substellar longitude
          if(nlong.gt.nmaxlong)stop' E- Many longitudes per period'
       endif
       if(nummethod.eq.2)write(*,1002)
       read(2, *)sensi,drate,stsensi,endsensi,max_nout
       read(2, *)rtol,maxerrsteady
       read_varch=0
       if(sensi.eq.3)read_varch=1
       close(2)

       return

100    write(*,*)' E- Error opening model input file: ',inpfile
       stop
110    write(*,*)' E- Error reading eddy file ', eddyfile,
     # ' and/or dimension (1d or 2d)'
       stop

1000   format(2x,'I- Running model ',a)
1001   format(2x,'I- Pressure range from '1pg10.2' to '1pg10.2
     # ' with 'i4' height layers')
1002   format(2x,'I- DLSODES selected as numerical method')
1010   format(2x,'W- Reaction warnings deactivated')

       end



c_______________________________________________________________________

       subroutine read_spec
c      ...read species included in the model

       use pac_data
       include 'pac.common'
       integer i,j,k,l,iwk,wknat(nmaxspec,nmaxelemxs)
       real*8 dumab,dumdhf
       character*130 line
       character*2 cwk(2*nmaxelemxs),wkelem(nmaxspec,nmaxelemxs)
       logical usedelem,idem

c      ...open species file
       open(unit=2,err=100,file=specfile,status='old')
       write(*,1000)specfile

c      ...read species
       nspec=0
       do while(nspec.le.nmaxspec+1)
          line(1:130)=''
          read(2,'(a130)',end=10)line
          if(line(1:1).eq.'!')cycle

c      ...read species name and properties
          nspec=nspec+1
          if(nspec.gt.nmaxspec) stop ' E- Too many species'
          read(line(1:85),*,err=110)spec(nspec),dumab,dumdhf,
     #    dipmom(nspec),polary(nspec),ionp(nspec),charge(nspec),
     #    nattot(nspec),nelemspec(nspec)

          if(nelemspec(nspec).gt.nmaxelemxs)goto 120
          i=nmaxcharspec
          do while (spec(nspec)(i:i).eq.' ')
             i=i-1
          enddo
          ilenspec(nspec)=i
          if(ilenspec(nspec).gt.nmaxcharspec)goto 130

c      ...read containing elements
          do i=1,2*nmaxelemxs
             cwk(i)(1:2)=''
          enddo
          read(line(90:130),*,err=110)(cwk(i),i=1,2*nelemspec(nspec))
          do i=1,nelemspec(nspec)
             j=2*(i-1)+1
             read(cwk(j)(1:2),*,err=110)wkelem(nspec,i)
             read(cwk(j+1)(1:2),*,err=110)wknat(nspec,i)
          enddo

c      ...verify consistency of number of atoms (total and of each element)
          iwk=0
          do i=1,nelemspec(nspec)
             iwk=iwk+wknat(nspec,i)
          enddo
          if(iwk.ne.nattot(nspec))goto 140
       enddo
10     close(2)

c      ...get elements included and evaluate species mass
       write(*,1020)
       nelem=0
       do i=1,nspec
          mass(i)=0.0d0
       enddo
       do i=1,nmaxelem
          usedelem=.false.
          do j=1,nspec
             do k=1,nelemspec(j)
                call same_string(2,element(i),wkelem(j,k),idem)
                if(idem)then
                   if(.not.usedelem)then
                      nelem=nelem+1
                      elem(nelem)(1:2)=element(i)(1:2)
                      usedelem=.true.
                      write(*,1030)element(i)
                   endif
                   nat(j,nelem)=wknat(j,k)
                   mass(j)=mass(j)+nat(j,nelem)*atomic_weight(i)
                endif
             enddo
          enddo
       enddo
       do i=1,nspec
          mass(i)=mass(i)*amu             ! amu -> g
       enddo

c      ...verify that each species appears just once
       do i=1,nspec-1
          k=ilenspec(i)
          do j=i+1,nspec
             l=ilenspec(j)
             if(k.ne.l)cycle
             if(spec(i)(1:k).eq.spec(j)(1:k))goto 150
          enddo
       enddo

c      ...evaluate if there are charged species
       do i=1,nspec
          if(charge(i).ne.0)ion=.true.
       enddo

c      ...identify key species (H2, H, e-, ...)
       do i=1,nspec
          if(ilenspec(i).eq.2.and.
     #    (spec(i)(1:2).eq.'h2'.or.spec(i)(1:2).eq.'H2'))  id_h2=i
          if(ilenspec(i).eq.1.and.
     #    (spec(i)(1:1).eq.'h'.or.spec(i)(1:1).eq.'H'))    id_h=i
          if(ilenspec(i).eq.2.and.
     #    (spec(i)(1:2).eq.'he'.or.spec(i)(1:2).eq.'He'.or.
     #    spec(i)(1:2).eq.'HE'))                           id_he=i
          if(ilenspec(i).eq.3.and.
     #    (spec(i)(1:3).eq.'h2o'.or.spec(i)(1:3).eq.'H2O'))id_h2o=i
          if(ilenspec(i).eq.3.and.
     #    (spec(i)(1:3).eq.'co2'.or.spec(i)(1:3).eq.'CO2'))id_co2=i
          if(ilenspec(i).eq.3.and.
     #    (spec(i)(1:3).eq.'ch4'.or.spec(i)(1:3).eq.'CH4'))id_ch4=i
          if(nattot(i).eq.0.and.charge(i).eq.-1)           id_electr=i
       enddo
       if(ion.and.id_electr.eq.0)stop ' E- electron is missing'

       write(*,1010)nspec

       return

100    write(*,*)' E- Error opening species file ',specfile
       stop
110    write(*,*)' E- Error reading species file at line:'
       write(*,*)' ',line
       stop
120    write(*,*)' E- Too many elements for ',spec(nspec)
       stop
130    write(*,*)' E- Too many characters for ',spec(nspec)
       stop
140    write(*,*)' E- Wrong number of atoms for ',spec(nspec)
       stop
150    write(*,*)' E- Species ',spec(i)(1:k),' included twice'
       stop

1000   format(1x,' I- Reading species from ',a45)
1010   format(1x,' I- Data read for ',i4,' species')
1020   format(1x,' I- Elements included:')
1030   format(1x,'                       ',a2)

       end




c_______________________________________________________________________

       subroutine read_zstruct_initab
c      ...read (z,p,T) structure and initial abundances

       include 'pac.common'
       integer nmaxz_read
       parameter(nmaxz_read=10000)

       integer i,j,k,l,xnz,nzspec,idab(nmaxspec),j1,j2,js,jj
       real*8 ab(nmaxspec),chtot,g,slope,zt,zb,zs
       real*8, allocatable, dimension(:)   :: xz,xpr,xtk
       real*8, allocatable, dimension(:,:) :: xab
       character*3 ext
       character*20 text(3)
       character*(nmaxcharspec) zspec(nmaxspec)
       character*100 line
       logical fnd_spec(nmaxspec),fnd_zspec(nmaxspec),idem

       i=nmaxcharfile
       do while (zabfile(i:i).ne.'.')
          i=i-1
       enddo
       ext(1:3)=zabfile(i+1:i+3)
       allocate(xz(nmaxz_read))
       allocate(xpr(nmaxz_read))
       allocate(xtk(nmaxz_read))
       allocate(xab(nspec,nmaxz_read))

c      ...initial abundances are provided
       if(ext(1:3).eq.'zab')then
          write(*,1000)zabfile
          open(unit=2,err=100,file=zabfile,status='old')
          do                                     ! read header of .zab file
             line(1:100)=''                      ! and species names
             read(2,'(a100)',err=110)line
             if(line(1:1).eq.'!')cycle
             exit
          enddo
          read(line(1:100),*,err=110)xnz
          if(xnz.gt.nmaxz_read) stop ' E- Too many input z points'
          read(2,*,err=110)nzspec
          if(nzspec.gt.nmaxspec) stop ' E- Too many .zab species'
          read(2,*,err=110)(text(j),j=1,3),(zspec(i),i=1,nzspec)
          do i=1,nspec                           ! compare species names in .zab file
             fnd_spec(i)=.false.                 ! with those read in .spec file
          enddo
          do j=1,nzspec
             l=nmaxcharspec
             do while (zspec(j)(l:l).eq.' ')
                l=l-1
             enddo
             fnd_zspec(j)=.false.
          enddo
          do i=1,nspec
             do j=1,nzspec
                if(fnd_zspec(j))cycle
                call same_string(nmaxcharspec,spec(i),zspec(j),idem)
                if(idem)then
                   idab(i)=j
                   fnd_zspec(j)=.true.
                   fnd_spec(i)=.true.
                   exit
                endif
             enddo
          enddo
          k=0
          do j=1,nzspec
             if(fnd_zspec(j))k=k+1
          enddo
          if(k.ne.nzspec)then
             do j=1,nzspec
                if(.not.fnd_zspec(j))write(*,1020)zspec(j)
             enddo
          endif

          do j=1,xnz                                ! read z,p,T and initial abundances
             read(2,*,err=110)xz(j),xpr(j),xtk(j),(ab(i),i=1,nzspec)
             xz(j)=xz(j)*1.0d5                      ! km  -> cm
             xpr(j)=xpr(j)*1.0d6                    ! bar -> dyn cm-2
             do i=1,nspec
                if(fnd_spec(i))then
                   xab(i,j)=10.0d0**ab(idab(i))     ! molar fraction
                else
                   xab(i,j)=0.0d0
                endif
             enddo
          enddo
          close(2)
          do i=1,nspec
             if(.not.fnd_spec(i))
     #       write(*,*)' W- Initial abundance set to zero for ',spec(i)
          enddo

c      ...initial abundances are not provided (everything atomic and same abundance for all atoms)
       elseif(ext(1:3).eq.'apt')then
          write(*,1010)zabfile
          open(unit=2,err=100,file=zabfile,status='old')
          j=0
          do
             line(1:100)=''
             read(2,'(a100)',err=110,end=10)line
             if(line(1:1).eq.'!')cycle
             j=j+1
             if(xnz.gt.nmaxz_read) stop ' E- Too many input z points'
             read(line(1:100),*,err=110)xz(j),xpr(j),xtk(j)
          enddo
10        close(2)
          xnz=j
          do j=1,xnz
             do i=1,nspec
                xab(i,j)=0.0d0
                if(nattot(i).eq.1.and.charge(i).eq.0)
     #          xab(i,j)=1.0d0/dble(nelem)
             enddo
          enddo

       else
          stop ' E- (z,p,T) file extension not recognized'
       endif

c      ...interpolate to (z,p,T) grid
       if(xz(1).gt.xz(xnz))then                  ! input (z,p,T) ordered from top to bottom
          j1=1
          j2=xnz
          js=1
       else                                      ! input (z,p,T) ordered from bottom to top
          j1=xnz
          j2=1
          js=-1
       endif
       if(xpr(j1).gt.pmin)         goto 120
       if(xpr(j2).lt.pmax)         goto 130
       do j=j1,j2-js,js
          if(xz(j).le.xz(j+js))    goto 140
          if(xpr(j).gt.xpr(j+js))  goto 150
       enddo
       j=j1                                      ! locate top z
       do while(xpr(j).le.pmin.and.(js*j).lt.(js*j2))
          j=j+js
       enddo
       slope=(xz(j)-xz(j-js))/
     # (dlog10(xpr(j))-dlog10(xpr(j-js)))
       zt=xz(j-js)+slope*(dlog10(pmin)-dlog10(xpr(j-js)))
       j=j2                                      ! locate bottom z
       do while(xpr(j).ge.pmax.and.(js*j).gt.(js*j1))
          j=j-js
       enddo
       slope=(xz(j)-xz(j+js))/
     # (dlog10(xpr(j))-dlog10(xpr(j+js)))
       zb=xz(j+js)+slope*(dlog10(pmax)-dlog10(xpr(j+js)))
       zs=(zt-zb)/dble(nz-1)
       do j=1,nz                                 ! create (z,p,T,n,abun) grid
          z(j)=zt-(j-1)*zs
          if(j.eq.nz)z(j)=zb
          jj=j1
          do while(xz(jj).ge.z(j).and.(js*jj).lt.(js*j2))
             jj=jj+js
          enddo
          slope=(dlog10(xpr(jj))-dlog10(xpr(jj-js)))/
     #    (xz(jj)-xz(jj-js))
          pressure(j)=10.0d0**
     #    (dlog10(xpr(jj-js))+slope*(z(j)-xz(jj-js)))
          slope=(xtk(jj)-xtk(jj-js))/(xz(jj)-xz(jj-js))
          tk(j)=xtk(jj-js)+slope*(z(j)-xz(jj-js))
          density(j)=pressure(j)/kboltz/tk(j)    ! cm-3
          do i=1,nspec
             slope=(dlog10(xab(i,jj))-dlog10(xab(i,jj-js)))/
     #       (xz(jj)-xz(jj-js))
             abun(i,j)=10.0d0**
     #       (dlog10(xab(i,jj-js))+slope*(z(j)-xz(jj-js)))
          enddo
       enddo
       deallocate(xz)
       deallocate(xpr)
       deallocate(xtk)
       deallocate(xab)

c      ...normalize abundances (Sum_i[abun_i]=1) for each layer
       call normalize_abun

c     ...check if charge is balanced
       do j=1,nz
          chtot=0.0d0
          do i=1,nspec
             chtot=chtot+charge(i)*abun(i,j)
          enddo
       enddo
       if(dabs(chtot).gt.1.0d-3)goto 160

c     ...compute size of ODEs = number of layers * number of species
       nzs=nz*nspec

c     ...evaluate scale height for each species and layer
       if(diffusion)then
          do j=1,nz
             g=ggravity*mplanet/(rplanet+z(j))**2.0d0
             do i=1,nspec
                hspe(i,j)=kboltz*tk(j)/(mass(i)*g)
             enddo
          enddo
       endif

       return

100    write(*,*)' E- Error opening (z,p,T) file: ',zabfile
       stop
110    write(*,*)' E- Error reading (z,p,T) file at line:'
       write(*,*)line
       stop
120    stop ' E- Requested pressure at top outside input grid'
130    stop ' E- Requested pressure at bottom outside input grid'
140    write(*,*)' E- Discontinuity in z at j=',j
       stop
150    write(*,*)' E- Discontinuity in pressure at j=',j
       stop
160    write(*,*)' E- Global charge is not zero but ',chtot
       stop

1000   format(1x,' I- Reading vertical structure',/,
     #        1x,'    and initial abundances from ',a45)
1010   format(1x,' I- Reading vertical structure from ',a45,/,
     #        1x,'    (initial abundances set as atomic)')
1020   format(1x,' W- Species: ',a12,' in .zab file not included')

       end




c_______________________________________________________________________

       subroutine read_reac
c      ...read reactions included in the model

       include 'pac.common'
       integer i,j,k,l,irfile,ir,ip,ilen,il,nrf,nrd,nrfnd,npfnd,
     # idk,idl,iwk1,iwk2,balr,balp,ipr,nrefpr,ind1
       real*8 rwk, refpr, a, b, c
       character*256 line, trash
       character*(nmaxcharspec) cwk
       character*(nmaxcharfile) usedreacfile
       logical pass,inclu,rfnd(nmaxreag),pfnd(nmaxprod),wkcoll,idem

       ilen=0                                    ! to avoid warning
       nrplog=0                                  ! number of plog dependent reactions
       nrcoll=0                                  ! number of collisional reactions

c     ...open file to write the reactions used
       usedreacfile(1:nmaxcharfile)=''
       usedreacfile(1:ncharmodel)=namemodel(1:ncharmodel)
       usedreacfile(ncharmodel+1:ncharmodel+14)='_reac_used.dat'
       open(unit=3,file=usedreacfile,status='unknown')

c      ...open reaction files
       i=0
       do irfile=1,nreacfile
          write(*,1000)reacfile(irfile)
          open(unit=2,file=reacfile(irfile),status='old')

          do while(i.le.nmaxreac)
c           write(*,*) 'i=', i

c      ...read reaction file lines; format: a50, 48 char for data, # error + source # ref
             line(1:256)=''
             read(2,'(a256)',end=10)line
             if(line(1:1).eq.'!'.or.line(1:256).eq.' ')cycle
             i=i+1
             read(line(1:50),'(a50)',err=100)reac(i)
             if(read_varch.eq.0)then
               read(line(51:83),*,err=100)alpha(i),beta(i),gamma(i),
     #       rtype(i)
             else
               read(line(51:83),*,err=100)alpha(i),beta(i),gamma(i),
     #       var_rate(i)
             read(line(87:108),1010,err=100)iwk1,iwk2,tkcode(i), ! @Arturo is this deprecated?
     #       lerror(i),ksource(i)
             read(line(111:256),'(a146)',err=100)comment(i)
             tkmin(i)=dble(iwk1)
             tkmax(i)=dble(iwk2)
             if(tkcode(i).ne.'<fix'.and.tkcode(i).ne.'>fix'.and. ! @Arturo is this deprecated?
     #          tkcode(i).ne.'<>fx'.and.tkcode(i).ne.'<zer'.and.
     #          tkcode(i).ne.'>zer'.and.tkcode(i).ne.'<>zr')
     #          tkcode(i)(1:4)='    '
             if(lerror(i).ne.'A'.and.lerror(i).ne.'B'.and.
     #          lerror(i).ne.'C'.and.lerror(i).ne.'D'.and.
     #          lerror(i).ne.'E')
     #          lerror(i)(1:1)=' '
             if(ksource(i).ne.'M'.and.ksource(i).ne.'E'.and.
     #          ksource(i).ne.'L'.and.ksource(i).ne.'C')
     #          ksource(i)(1:1)=' '
             if(lerror(i).eq.'A')kerror(i)=1.25d0
             if(lerror(i).eq.'B')kerror(i)=1.50d0
             if(lerror(i).eq.'C')kerror(i)=2.00d0
             if(lerror(i).eq.'D')kerror(i)=10.0d0
             if(lerror(i).eq.'E')kerror(i)=10.0d0 ! highly uncertain => factor 10
             if(lerror(i).eq.' ')kerror(i)=2.00d0 ! if unknown error => factor 2

c      ...decode reaction into reagents and products
             ir=1
             ip=1
             pass=.false.
             do j=1,50
                if(reac(i)(j:j).eq.' ')then
                   cycle
                elseif(reac(i)(j-1:j+1).eq.' + ')then
                   if(.not.pass)ir=ir+1
                   if(pass)     ip=ip+1
                elseif(reac(i)(j-1:j+1).eq.' = ')then
                   pass=.true.
                else
                   if(j.eq.1.or.reac(i)(j-1:j-1).eq.' ')then
                      ilen=1
                      if(.not.pass)then
                         reag(i,ir)(1:nmaxcharspec)=''
                         reag(i,ir)(ilen:ilen)=reac(i)(j:j)
                         ilenreag(i,ir)=ilen
                      else
                         prod(i,ip)(1:nmaxcharspec)=''
                         prod(i,ip)(ilen:ilen)=reac(i)(j:j)
                         ilenprod(i,ip)=ilen
                      endif
                   else
                      ilen=ilen+1
                      if(.not.pass)then
                         reag(i,ir)(ilen:ilen)=reac(i)(j:j)
                         ilenreag(i,ir)=ilen
                      else
                         prod(i,ip)(ilen:ilen)=reac(i)(j:j)
                         ilenprod(i,ip)=ilen
                      endif
                   endif
                endif
             enddo
             nreag(i)=ir
             nprod(i)=ip

c      ...do some verifications
             if(reac(i)(1:50).eq.' '.or.ir.le.0.or.ip.le.0)goto 110
             if(ir.gt.nmaxreag.or.ip.gt.nmaxprod)          goto 120
             do j=1,ir
                if(ilenreag(i,j).gt.nmaxcharspec)          goto 130
             enddo
             do j=1,ip
                if(ilenprod(i,j).gt.nmaxcharspec)          goto 130
             enddo
             if( iabs(rtype(i)).ne.  3 .and.     ! neutral-neutral reaction
     #          (iabs(rtype(i)).lt. 30 .or.      ! neutral-neutral reaction with plog dependence
     #           iabs(rtype(i)).gt. 45).and.     ! ! @Arturo changed max number of pressures allowed
     #           iabs(rtype(i)).ne.300 .and.     ! various-body reaction
     #           iabs(rtype(i)).ne.  4 .and.     ! ion-neutral reaction
     #          (iabs(rtype(i)).lt.400 .or.      ! ion-neutral reaction with special treatment
     #           iabs(rtype(i)).gt.420).and.     !
     #           iabs(rtype(i)).ne.  5 .and.     ! radiative association reaction
     #           iabs(rtype(i)).ne.  6 .and.     ! dissociative recombination reaction
     #           iabs(rtype(i)).ne.  7 .and.     ! mutual neutralization reaction
     #           iabs(rtype(i)).ne.  8 .and.     ! three body  reaction (explicitly written)
     #           iabs(rtype(i)).ne.800 .and.     ! three body  reaction (general Troe form: k_0/k_infinite + M enhancements)
     #           iabs(rtype(i)).ne.801 .and.     ! three body  reaction (k_0 + M enhancements)
     #           iabs(rtype(i)).ne.802 .and.     ! three body  reaction (k_infinite)
     #           iabs(rtype(i)).ne.810 .and.     ! three body  reaction (general SRI  form: k_0/k_infinite + M enhancements)
     #           iabs(rtype(i)).ne.811 .and.     ! three body  reaction (general Troe form: k_0/k_infinite, M=H2O only)
     #           iabs(rtype(i)).ne.812 .and.     ! three body  reaction (general Troe form: k_0/k_infinite, M=CO2 only)
     #           iabs(rtype(i)).ne.  9 .and.     ! thermal dissociation/isomerization (explicitly written)
     #           iabs(rtype(i)).ne.900 .and.     ! thermal dissociation/isomerization (general Troe form: k_0/k_infinite + M enhancements)
     #           iabs(rtype(i)).ne.901 .and.     ! thermal dissociation/isomerization (k_0 + M enhancements)
     #           iabs(rtype(i)).ne.902 .and.     ! thermal dissociation/isomerization (k_infinite)
     #          (iabs(rtype(i)).lt. 902 .or.     ! thermal dissociation/isomerization (k_infinite) @Thomas: thermal diss. with plog
     #           iabs(rtype(i)).gt. 909).and.    !
     #           iabs(rtype(i)).ne.910 .and.     ! thermal dissociation/isomerization (general SRI  form: k_0/k_infinite + M enhancements)
     #           iabs(rtype(i)).ne.911 .and.     ! thermal dissociation/isomerization (general Troe form: k_0/k_infinite, M=CH4 only)
     #           iabs(rtype(i)).ne. 20)          ! radiative de-excitation
     #                                                     goto 140

c     ...neutral-neutral reaction with plog dependence
             idrplog(i)=0
             if((iabs(rtype(i)).gt.30.and.
     #          iabs(rtype(i)).lt.45))then
                nrplog=nrplog+1
                idrplog(i)=nrplog
                nrefpr=iabs(rtype(i))-30         ! number of reference pressures to be read
                
                read(comment(i),*)refpr          ! for first line, convert comment into ref. pressure (string-->double)
                refpr=refpr*1.0d6                ! bar --> dyne cm-2
                prefplog(nrplog,1)=refpr         ! save plog-specific info (reference pressure and rate coefficients)
                alphaplog(nrplog,1)=alpha(i)
                betaplog(nrplog,1)=beta(i)
                gammaplog(nrplog,1)=gamma(i)
                do ipr=2,nrefpr                  ! for remaining reference pressures, read data and save
                   read(2,*)a,b,c,trash,refpr
                   refpr=refpr*1.0e6             ! bar --> dyne cm-2
                   if(refpr.lt.prefplog(nrplog,ipr-1))goto 200
                   prefplog(nrplog,ipr)=refpr
                   alphaplog(nrplog,ipr)=a
                   betaplog(nrplog,ipr)=b
                   gammaplog(nrplog,ipr)=c
                enddo
                
             endif


c     ...@Thomas: thermal diss. with plog dependance
c             idrplog(i)=0 this makes that the previous indexes get lost @Arturo
             if((iabs(rtype(i)).gt.902.and.
     #          iabs(rtype(i)).lt.910))then
                nrplog=nrplog+1
                idrplog(i)=nrplog
                nrefpr=iabs(rtype(i))-902         ! number of reference pressures to be read
                read(comment(i),*)refpr          ! for first line, convert comment into ref. pressure (string-->double)
                refpr=refpr*1.0d6                ! bar --> dyne cm-2
                prefplog(nrplog,1)=refpr         ! save plog-specific info (reference pressure and rate coefficients)
                alphaplog(nrplog,1)=alpha(i)
                betaplog(nrplog,1)=beta(i)
                gammaplog(nrplog,1)=gamma(i)
                do ipr=2,nrefpr                  ! for remaining reference pressures, read data and save
                   read(2,*)a,b,c,trash,refpr
                   refpr=refpr*1.0e6             ! bar --> dyne cm-2
                   if(refpr.lt.prefplog(nrplog,ipr-1))goto 200
                   prefplog(nrplog,ipr)=refpr
                   alphaplog(nrplog,ipr)=a
                   betaplog(nrplog,ipr)=b
                   gammaplog(nrplog,ipr)=c
                enddo
             endif

c     ...collisional reaction with additional parameters
             wkcoll=.false.
             idrcoll(i)=0
             if(iabs(rtype(i)).eq.800.or.
     #          iabs(rtype(i)).eq.810.or.
     #          iabs(rtype(i)).eq.801.or.
     #          iabs(rtype(i)).eq.811.or.
     #          iabs(rtype(i)).eq.812.or.
     #          iabs(rtype(i)).eq.900.or.
     #          iabs(rtype(i)).eq.910.or.
     #          iabs(rtype(i)).eq.901.or.
     #          iabs(rtype(i)).eq.911)then
                wkcoll=.true.
                nrcoll=nrcoll+1
                idrcoll(i)=nrcoll
                if(iabs(rtype(i)).eq.800.or.
     #             iabs(rtype(i)).eq.810.or.
     #             iabs(rtype(i)).eq.811.or.
     #             iabs(rtype(i)).eq.812.or.
     #             iabs(rtype(i)).eq.900.or.
     #             iabs(rtype(i)).eq.910.or.
     #             iabs(rtype(i)).eq.911)then
                   read(2,*)alphainf(nrcoll),betainf(nrcoll),
     #             gammainf(nrcoll)
                endif
                if(iabs(rtype(i)).eq.800.or.
     #             iabs(rtype(i)).eq.811.or.
     #             iabs(rtype(i)).eq.812.or.
     #             iabs(rtype(i)).eq.900.or.
     #             iabs(rtype(i)).eq.911)then
                   read(2,*)afoff(nrcoll),bfoff(nrcoll),cfoff(nrcoll),
     #             dfoff(nrcoll),iwk1
                elseif(iabs(rtype(i)).eq.810.or.
     #                 iabs(rtype(i)).eq.910)then
                   read(2,*)afoff(nrcoll),bfoff(nrcoll),cfoff(nrcoll),
     #             dfoff(nrcoll),efoff(nrcoll),iwk1
                elseif(iabs(rtype(i)).eq.801.or.
     #                 iabs(rtype(i)).eq.901)then
                   read(2,*)iwk1
                endif
                nmrcoll(nrcoll)=0
                if(iwk1.gt.0)then
                   do j=1,iwk1
                      read(2,*)cwk,rwk
                      do k=1,nspec
                         call same_string
     #                        (nmaxcharspec,cwk,spec(k),idem)
                         if(idem)then
                            nmrcoll(nrcoll)=nmrcoll(nrcoll)+1
                            idmrcoll(nrcoll,nmrcoll(nrcoll))=k
                            fmrcoll(nrcoll,nmrcoll(nrcoll))=rwk
                         endif
                      enddo
                   enddo
                endif
             endif

c      ...check whether the species involved are actually considered in the model
             inclu=.true.
             loop40: do j=1,ir
                il=ilenreag(i,j)
                do k=1,nspec
                   call same_string
     #                  (nmaxcharspec,reag(i,j),spec(k),idem)
                   if(idem)then
                      idreag(i,j)=k
                      cycle loop40
                   endif
                enddo
                inclu=.false.
                goto 20
             enddo loop40
             loop50: do j=1,ip
                il=ilenprod(i,j)
                do k=1,nspec
                   call same_string
     #                  (nmaxcharspec,prod(i,j),spec(k),idem)
                   if(idem)then
                      idprod(i,j)=k
                      cycle loop50
                   endif
                enddo
                inclu=.false.
                goto 20
             enddo loop50
             if(ir.eq.1.and.warnreac.eq.1)write(*,1020)reac(i) ! @Arturo hid warnings

c      ...do some more verifications
             k=0
             do j=1,ir
                k=k+charge(idreag(i,j))
             enddo
             if(iabs(rtype(i)).eq.3)then
                if(ir.gt.2)                                goto 140 ! changed ne for gt @Arturo
                if(k.ne.0)                                 goto 140
             elseif(iabs(rtype(i)).eq.300)then
                if(ir.le.2)                                goto 140
             elseif(iabs(rtype(i)).eq.4.or.
     #       (iabs(rtype(i)).ge.400.and.iabs(rtype(i)).le.420))then
                if(ir.ne.2)                                goto 140
                if(k.ne.1.and.k.ne.-1)                     goto 140
             elseif(iabs(rtype(i)).eq.5)then
                if(ir.ne.2.or.ip.ne.1)                     goto 140
             elseif(iabs(rtype(i)).eq.6)then
                if(ir.ne.2)                                goto 140
                if(charge(idreag(i,1)).ne.1)               goto 140
                if(idreag(i,2).ne.id_electr)               goto 140
             elseif(iabs(rtype(i)).eq.7)then
                if(ir.ne.2)                                goto 140
                if(charge(idreag(i,1)).eq.0)               goto 140
                if(charge(idreag(i,2)).eq.0)               goto 140
                if(k.ne.0)                                 goto 140
             elseif(iabs(rtype(i)).eq.8)then
                if(ir.ne.3.or.ip.ne.2)                     goto 140
             elseif(iabs(rtype(i)).ge.800.and.
     #       iabs(rtype(i)).le.812)then
                if((ir.ne.2.or.ip.ne.1).and.
     #           iabs(rtype(i)).ne.802)                    goto 140 ! @Arturo added 802 as an special case, for Veillet CHV and OHV
             elseif(iabs(rtype(i)).eq.9)then
                if(ir.ne.2.or.ip.ne.3)                     goto 140
             elseif(iabs(rtype(i)).ge.900.and.
     #       iabs(rtype(i)).le.911)then
                if(ir.ne.1.and.
     #           iabs(rtype(i)).ne.901)                    goto 140
             elseif(iabs(rtype(i)).eq.20)then
                if(ir.ne.1.or.ip.ne.1)                     goto 140
             endif

c      ...check the reaction mass and charge balance
             do k=0,nelem
                balr=0
                do j=1,ir
                   if(k.eq.0)balr=balr+charge(idreag(i,j))
                   if(k.gt.0)balr=balr+nat(idreag(i,j),k)
                enddo
                balp=0
                do j=1,ip
                   if(k.eq.0)balp=balp+charge(idprod(i,j))
                   if(k.gt.0)balp=balp+nat(idprod(i,j),k)
                enddo
                if(balr.ne.balp)goto 150
             enddo

c      ...check that the reaction did not appear before
             if(i.eq.1)goto 20
30           do j=1,i-1
                if(nreag(i).eq.nreag(j).and.nprod(i).eq.nprod(j))then
                   do l=1,nreag(j)
                      rfnd(l)=.false.
                   enddo
                   nrfnd=0
                   do k=1,nreag(i)
                      idk=idreag(i,k)
                      do l=1,nreag(j)
                         idl=idreag(j,l)
                         if(idk.eq.idl.and..not.rfnd(l))then
                            rfnd(l)=.true.
                            nrfnd=nrfnd+1
                            exit
                         endif
                      enddo
                   enddo
                   do l=1,nprod(j)
                      pfnd(l)=.false.
                   enddo
                   npfnd=0
                   do k=1,nprod(i)
                      idk=idprod(i,k)
                      do l=1,nprod(j)
                         idl=idprod(j,l)
                         if(idk.eq.idl.and..not.pfnd(l))then
                            pfnd(l)=.true.
                            npfnd=npfnd+1
                            exit
                         endif
                      enddo
                   enddo
                   if(nrfnd.eq.nreag(i).and.npfnd.eq.nprod(i))then
                      if((tkcode(j).eq.'>zer'.or.tkcode(j).eq.'<>zr')
     #                .and.(tkcode(i).eq.'<zer'.or.tkcode(i).eq.'<>zr')
     #                .and.tkmax(j).eq.tkmin(i))then
c      ...reactions duplicated for various temperature intervals must be in order of increasing temperature
                         write(*,1030)reac(i)
                      else
                         if(warnreac.eq.1)then ! @Arturo hid warnings
                          write(*,*)' E- These reactions are repeated:'
                          write(*,*)reac(j)
                          write(*,*)reac(i)
                         endif
                      endif
                   endif
                endif
             enddo

c      ...depending whether reaction is included or not ...
20           if(inclu)then
                write(3,220)i,line
             else
                i=i-1
                if(wkcoll)nrcoll=nrcoll-1
             endif

c     ...include reverse reaction if needed
             if(inclu.and.rtype(i).lt.0)then
                i=i+1
                idrvreac(i)=i-1
                nreag(i)=nprod(i-1)
                if(nreag(i).gt.nmaxreag)goto 120
                reac(i)=''
                ilen=0
                reac(i)(1:1)=' '
                ilen=1
                do j=1,nreag(i)
                   reag(i,j)=prod(i-1,j)
                   ilenreag(i,j)=ilenprod(i-1,j)
                   idreag(i,j)=idprod(i-1,j)
                   if(idthspec(idreag(i,j)).le.0)goto 160
                   reac(i)(ilen+1:ilen+ilenreag(i,j))=
     #             reag(i,j)(1:ilenreag(i,j))
                   ilen=ilen+ilenreag(i,j)
                   reac(i)(ilen+1:ilen+3)=' + '
                   ilen=ilen+3
                enddo
                ilen=ilen-3
                reac(i)(ilen+1:ilen+3)=' = '
                ilen=ilen+3
                nprod(i)=nreag(i-1)
                if(nprod(i).gt.nmaxprod)goto 120
                do j=1,nprod(i)
                   prod(i,j)=reag(i-1,j)
                   ilenprod(i,j)=ilenreag(i-1,j)
                   idprod(i,j)=idreag(i-1,j)
                   if(idthspec(idprod(i,j)).le.0)goto 160
                   reac(i)(ilen+1:ilen+ilenprod(i,j))=
     #             prod(i,j)(1:ilenprod(i,j))
                   ilen=ilen+ilenprod(i,j)
                   reac(i)(ilen+1:ilen+3)=' + '
                   ilen=ilen+3
                enddo
                ilen=ilen-3
                reac(i)(ilen+1:ilen+3)='   '
                rtype(i)=0
                tkmin(i)=tkmin(i-1)
                tkmax(i)=tkmax(i-1)
                tkcode(i)=tkcode(i-1)
                lerror(i)=lerror(i-1)
                ksource(i)=ksource(i-1)
                comment(i)=' detailed balance'
                line=''
                write(line(2:50),'(a49)')reac(i)(1:49)
                write(line(52:83),'(1pg9.2,1x,f8.2,1x,i8,1x,i4)')
     #          alpha(i),beta(i),idint(gamma(i)),rtype(i)
                write(line(85:85),'(a1)')':'
                write(line(86:108),'(1x,i6,1x,i6,1x,a4,1x,a1,1x,a1)')
     #          idint(tkmin(i)),idint(tkmax(i)),tkcode(i),lerror(i),
     #          ksource(i)
                write(line(110:110),'(a1)')':'
                write(line(111:256),'(a146)')comment(i)
                goto 30
             endif
          enddo
      
10        close(2)
       enddo
       nreac=i
       if(nreac.ge.nmaxreac)   stop ' E- Too many reactions'
       if(nrcoll.ge.nmaxrcoll) stop ' E- Too many collision reactions'

       close(3) ! close ...reactions_used.dat file

c      ...check that all species have at least one formation/destruction reaction
       do i=1,nspec
          nrd=0
          nrf=0
          if(spec(i).eq.'He')then
           nrd=1
           nrf=1
          endif
          do j=1,nreac
             do k=1,nreag(j)
                if(idreag(j,k).eq.i)nrd=nrd+1
             enddo
             do k=1,nprod(j)
                if(idprod(j,k).eq.i)nrf=nrf+1
             enddo
          enddo
          if(nrd.eq.0)goto 170
          if(nrf.eq.0)goto 180
          if(nrd.lt.2.and.warnreac.eq.1)then ! @Arturo hid warnings
             write(*,1040)spec(i),nrd
          endif
          if(nrf.lt.2.and.warnreac.eq.1)then ! @Arturo hid warnings
             write(*,1050)spec(i),nrf
          endif
       enddo

c      ...build rate equations, i.e. dn_i/dt=...
       do i=1,nspec
          nterm(i)=0
          do j=1,nreac
             do k=1,nreag(j)             ! add sink terms
                if(idreag(j,k).eq.i)then
                   nterm(i)=nterm(i)+1
                   idterm(i,nterm(i))=j
                   sterm(i,nterm(i))=-1
                   fterm(i,nterm(i))=1
                   if(k.eq.nreag(j))goto 80
                   do l=k+1,nreag(j)
                      if(idreag(j,l).eq.i)
     #                fterm(i,nterm(i))=fterm(i,nterm(i))+1
                   enddo
                   goto 80
                endif
             enddo

80           do k=1,nprod(j)             ! add source terms
                if(idprod(j,k).eq.i)then
                   nterm(i)=nterm(i)+1
                   idterm(i,nterm(i))=j
                   sterm(i,nterm(i))=1
                   fterm(i,nterm(i))=1
                   if(k.eq.nprod(j))exit
                   do l=k+1,nprod(j)
                      if(idprod(j,l).eq.i)
     #                fterm(i,nterm(i))=fterm(i,nterm(i))+1
                   enddo
                   exit
                endif
             enddo
          enddo
          if(nterm(i).gt.nmaxrs)goto 190
       enddo

       write(*,2000)nreac

       return



100    write(*,*)' E- Error reading reactions file at line:'
       write(*,*)'  ',line
       stop
110    write(*,*)' E- ',reac(i)
       stop      '    wrong syntax in above reaction'
120    write(*,*)' E- ',reac(i)
       stop      '    too many reagents/products in above reaction'
130    write(*,*)' E- ',reac(i)
       stop      '    too many characters species in above reaction'
140    write(*,*)' E- ',reac(i)
       stop      '    wrong rtype for above reaction'
150    write(*,*)' E- ',reac(i)
       stop      '    Balance not fulfilled for above reaction'
160    write(*,*)' E- ',reac(i-1)
       stop      '    Lack of thermochemical data for above reaction'
170    write(*,*)' E- No destruction reaction for species:',spec(i)
       stop
180    write(*,*)' E- No formation reaction for species:',spec(i)
       stop
190    write(*,*)' E- Species ',spec(i)
       stop      '    involved in too many reactions'
200    write(*,*)' E- ',reac(i)
       stop      '    PLOG pressures should be monotonously increasing'
220    format(i4,3x,(a256))

1000   format(1x,' I- Reading reactions from ',a45)
1010   format(i6,1x,i6,1x,a4,1x,a1,1x,a1)
1020   format(1x,' W- Unimolecular reaction: ',a50)
1030   format(1x,' I- >1 temp intervals for: ',a50)
1040   format(1x,' W- Species ',a12,
     # ' has just ',i2,' destruction reaction')
1050   format(1x,' W- Species ',a12,
     # ' has just ',i2,' formation reaction')
2000   format(1x,' I- ',i4,'       reactions included in the model')


       end




c_______________________________________________________________________

       subroutine read_therm
c      ...read therm data

       include 'pac.common'
       integer i,j,k,l,ilen,read_nattot,read_nelem,read_charge,
     # read_thtype,read_nat(nmaxelemxs),wknat,read_nthtemp
       real*8 read_a(9,nmaxthtemp),read_thtemp(nmaxthtemp+1)
       character*256 line
       character*(nmaxcharspec) read_spec
       character*2 read_txt(2*nmaxelemxs),read_elem(nmaxelemxs)
       logical fnd_spec(nmaxspec),idem

       if(.not.therm)return

       write(*,1000)thermfile
       open(unit=2,file=thermfile,status='old')
       nthspec=0
       fnd_spec(1:nmaxspec)=.false.

c Big loop over therm species
       do
          line(1:256)=''
          read(2,'(a256)',end=10)line
          if(line(1:1).eq.'!')cycle

          read(line,*,err=100,end=10)            ! read therm parameters
     #    read_spec,read_nattot,read_nelem,read_charge,read_thtype
          read(2,*,err=110)(read_txt(j),j=1,2*read_nelem)
          do i=1,2*read_nelem
             if(((-1)**i).lt.0)read_elem((i+1)/2)(1:2)=read_txt(i)(1:2)
             if(((-1)**i).gt.0)read(read_txt(i)(1:2),*)read_nat(i/2)
          enddo
          if(read_thtype.eq.1)then
             read_nthtemp=2
             read(2,*,err=120)(read_thtemp(j),j=1,3)
             read(2,*,err=120)(read_a(j,1),j=1,7)
             read(2,*,err=120)(read_a(j,2),j=1,7)
          elseif(read_thtype.eq.2)then
             read(2,*,err=120)read_nthtemp,
     #       (read_thtemp(j),j=1,read_nthtemp+1)
             do i=1,read_nthtemp
                read(2,*,err=120)(read_a(j,i),j=1,9)
             enddo
          else
             write(*,*)' E- unknown therm data type for: ',read_spec
             stop
          endif

          i=nmaxcharspec                         ! verify whether species is included or not
          do while (read_spec(i:i).eq.' ')
             i=i-1
          enddo
          ilen=i
          do i=1,nspec
             if(fnd_spec(i))cycle
             if(ilen.ne.ilenspec(i))cycle
             call same_string(nmaxcharspec,read_spec,spec(i),idem)

             if(idem)then                        ! species is included
                if(read_nattot.ne.nattot(i))   goto 130
                if(read_nelem.ne.nelemspec(i)) goto 130
                if(read_charge.ne.charge(i))   goto 130
                do j=1,nelem
                   l=0
                   do k=1,read_nelem
                      call same_string(2,read_elem(k),elem(j),idem)
                      if(idem)l=k
                   enddo
                   wknat=0
                   if(l.ne.0)wknat=read_nat(l)
                   if(wknat.ne.nat(i,j))       goto 130
                enddo
                nthspec=nthspec+1                ! store therm data
                idthspec(i)=nthspec
                thtype(nthspec)=read_thtype

                nthtemp(nthspec)=read_nthtemp
                do k=1,nthtemp(nthspec)+1
                   thtemp(nthspec,k)=read_thtemp(k)
                   if(k.eq.nthtemp(nthspec)+1)exit
                   do j=1,9
                      tha(nthspec,j,k)=read_a(j,k)
                   enddo
                enddo
                fnd_spec(i)=.true.
                exit
             endif
          enddo

c End of big loop over therm species
       enddo
10     close(2)

       return

100    write(*,*)' E- Error reading therm file at line:'
       write(*,*)'    ',line
       stop
110    write(*,*)' E- Error reading therm file elements for: ',
     # read_spec
       stop
120    write(*,*)' E- Error reading thermo data for: ',read_spec
       stop
130    write(*,*)' E- Wrong data for species: ',read_spec
       stop      '    either in species file or in therm file'

1000   format(1x,' I- Reading therm data from ',a45)

       end




c_______________________________________________________________________

       subroutine read_eddy
c      ...read and create grid for eddy diffusion coefficient profile

       include 'pac.common'
       integer nmaxeddy_read
       parameter(nmaxeddy_read=10000)
       integer i,j,xnz,j1,j2,js,jj
       real*8 slope
       real*8, allocatable, dimension(:) :: xpr,xeddy
       character*100 line

       if(.not.diffusion)return

       allocate(xpr(nmaxeddy_read))
       allocate(xeddy(nmaxeddy_read))

c      ...read eddy diffusion coefficient profile
       open(unit=2,file=eddyfile,status='old')
       i=0
       do
          line(1:100)=''
          read(2,'(a100)',err=100,end=10)line
          if(line(1:1).eq.'!')cycle

          i=i+1
          if(i.gt.nmaxeddy_read)stop ' E- Eddy grid too fine'
          read(line(1:100),*,err=100)xpr(i),xeddy(i)       ! pressure [bar], K_eddy [cm2 s-1]
          xpr(i)=xpr(i)*1.0d6                   ! bar -> dyn cm-2
       enddo
10     close(2)
       xnz=i

c      ...interpolate to (z,p,T) grid
       if(xpr(1).lt.xpr(xnz))then                ! input (p,eddy) ordered from top to bottom
          j1=1
          j2=xnz
          js=1
          if(xpr(1).gt.pmin)       goto 110
          if(xpr(xnz).lt.pmax)     goto 120
          do j=1,xnz-1
             if(xpr(j).ge.xpr(j+1))goto 130
          enddo
       else                                      ! input (p,eddy) ordered from bottom to top
          j1=xnz
          j2=1
          js=-1
          if(xpr(xnz).gt.pmin)     goto 110
          if(xpr(1).lt.pmax)       goto 120
          do j=1,xnz-1
             if(xpr(j).le.xpr(j+1))goto 130
          enddo
       endif
       do j=1,nz                                 ! create K_eddy grid
          do jj=j1,j2,js
             if(xpr(jj).gt.pressure(j))exit
          enddo
          slope=(dlog10(xeddy(jj))-dlog10(xeddy(jj-js)))/
     #    (dlog10(xpr(jj))-dlog10(xpr(jj-js)))
          kdif(j)=xeddy(jj-js)*(pressure(j)/xpr(jj-js))**slope
       enddo

       deallocate(xpr)
       deallocate(xeddy)

       return

100    write(*,*)' E- Error reading eddy diffusion file: ',eddyfile
       stop
110    stop ' E- Requested pressure at top outside input eddy grid'
120    stop ' E- Requested pressure at bottom outside input eddy grid'
130    write(*,*)' E- Discontinuity in eddy grid at j=',j
       stop

       end




c_______________________________________________________________________

       subroutine read_photon
c      ...create wavelength grid and read stellar spectrum and photo cross sections

       include 'pac.common'
       integer i,j,k,l,ii,s,nxy,ir,ip,ilen,il,wkilenreag(2),
     # balr,balp,npfnd,idk,idl,nphotochannel
       real*8 dwlth
       real*8, allocatable, dimension(:)   :: x,y
       real*8, allocatable, dimension(:,:) :: qy
       character*256 line
       character*(nmaxcharfile) phfile
       character*(nmaxcharspec) wkreag(nmaxphotochannel)
       character*50 wkphotoreac(nmaxphotochannel)
       logical pass,pfnd(nmaxprod),idem,inclu

       if(.not.photochem)return

       ilen=0                                    ! to avoid warning
       nwlth=1000                                ! number of wavelengths
       wlthmin=1.0d0                             ! low wavelegth cutoff [nm]
       wlthmax=750.0d0                           ! up  wavelegth cutoff [nm]
       if(nwlth.gt.nmaxwlth)stop 'E- Wavelength grid is too fine'

c      ...create wavelength grid
       dwlth=(wlthmax-wlthmin)/dble(nwlth-1)
       do i=1,nwlth
          wlth(i)=wlthmin+(i-1)*dwlth
       enddo


c Stellar spectrum.........................................

c      ...read stellar spectrum from starfile
       allocate(x(nmaxxy),y(nmaxxy))
       write(*,1000)starfile
       open(unit=2,file=starfile,status='old')
       i=0
       do
          line(1:256)=''
          read(2,'(a256)',err=100,end=10)line
          if(line(1:1).eq.'!')cycle

          i=i+1
          if(i.gt.nmaxxy)stop ' E- Stellar spectrum too fine'
          read(line(1:256),*,err=100)x(i),y(i)             ! wavelength [nm], intensity [erg s-1 cm-2 Hz-1 sr-1]
          y(i)=y(i)*1.0d7*clight/x(i)/x(i)                 ! I [=] erg s-1 cm-2 nm-1 sr-1
          y(i)=y(i)*pi*(rstar/dist)**2.0d0                 ! I [=] erg s-1 cm-2 nm-1
          y(i)=y(i)*x(i)*1.0d-7/hplanck/clight             ! I [=] photon s-1 cm-2 nm-1
       enddo
10     close(2)
       nxy=i

c      ...resample stellar flux to the wavelength grid
       call resample_to_wlth(x,y,nxy,fphot)


c Photo cross sections.....................................

c      ...read photo cross section info file
       allocate(qy(nmaxphotochannel,nmaxxy))
       write(*,1010)iphotfile
       open(unit=2,file=iphotfile,status='old')
       s=0                                       ! number of photospecies
       k=0                                       ! number of photoreactions
       do
          line(1:256)=''
          read(2,'(a256)',err=110,end=20)line
          if(line(1:1).eq.'!')cycle

c      ...open and read photo cross section file
          phfile(1:nmaxcharfile)=''
          read(line(1:256),*,err=110)phfile
          open(unit=3,file='photoreac/'//phfile,status='old')
          j=-1
          do
             line(1:256)=''
             read(3,'(a256)',err=120,end=30)line
             if(line(1:1).eq.'!')cycle

c     ...read photo reaction(s)
             if(j.eq.-1)then
                read(line,*,err=120)nphotochannel
                if(nphotochannel.gt.nmaxphotochannel)      goto 130
                do ii=1,nphotochannel
                   wkphotoreac(ii)(1:50)=''
                   read(3,'(a50)',err=120)wkphotoreac(ii)
                enddo
                j=0
                cycle
             endif

c      ...read photo cross section: wavelength[nm], cross section[cm2], quantum yield(s)
             j=j+1
             if(j.gt.nmaxxy)                               goto 140
             read(line(1:256),*,err=120)x(j),y(j),
     #       (qy(ii,j),ii=1,nphotochannel)
          enddo
30        close(3)
          nxy=j

c >>>  Absorption cross sections of species

c      ...verify that reagent is included and store absorption cross section
          do ii=1,nphotochannel
             wkreag(ii)(1:nmaxcharspec)=''
             read(wkphotoreac(ii),*,err=120)wkreag(ii)
             i=nmaxcharspec
             do while(wkreag(ii)(i:i).eq.' ')
                i=i-1
             enddo
             if(ii.gt.1)then
                if(i.ne.j)                                 goto 150
                if(wkreag(ii)(1:i).ne.wkreag(ii-1)(1:i))   goto 150
             endif
             j=i
          enddo
          do i=1,nspec
             call same_string(nmaxcharspec,wkreag(1),spec(i),idem)
             if(idem)then
                s=s+1
                if(s.gt.nmaxphotoreac)stop ' E- Many photo species'
                idphotospec(s)=i
                call resample_to_wlth(x,y,nxy,csphotabs(s,:))
                exit
             endif
          enddo
          if(.not.idem)then
             write(*,1020)wkreag(1)
             cycle
          endif

c >>>  Cross sections of photo reactions (photochemical channels)

          do ii=1,nphotochannel
             k=k+1
             if(k.gt.nmaxphotoreac)stop ' E- Many photo processes'
             photoreac(k)(1:50)=''
             photoreac(k)(1:50)=wkphotoreac(ii)(1:50)

c      ...decode photo reaction into reagents and products
             ir=1
             ip=1
             pass=.false.
             wkreag(1)(1:nmaxcharspec)=''
             wkreag(2)(1:nmaxcharspec)=''
             photoreag(k)(1:nmaxcharspec)=''
             do i=1,nmaxprod
                photoprod(k,i)(1:nmaxcharspec)=''
             enddo
             do j=1,50
                if(photoreac(k)(j:j).eq.' ')then
                   cycle
                elseif(photoreac(k)(j-1:j+1).eq.' + ')then
                   if(.not.pass)ir=ir+1
                   if(pass)     ip=ip+1
                elseif(photoreac(k)(j-1:j+1).eq.' = ')then
                   pass=.true.
                else
                   ilen=ilen+1
                   if(j.eq.1.or.photoreac(k)(j-1:j-1).eq.' ')ilen=1
                   if(.not.pass)then
                      wkreag(ir)(ilen:ilen)=photoreac(k)(j:j)
                      wkilenreag(ir)=ilen
                   else
                      photoprod(k,ip)(ilen:ilen)=photoreac(k)(j:j)
                      ilenphotoprod(k,ip)=ilen
                   endif
                endif
             enddo
             ilenphotoreag(k)=wkilenreag(1)
             photoreag(k)(1:nmaxcharspec)=wkreag(1)(1:nmaxcharspec)
             nphotoprod(k)=ip

c      ...check for possible errors in the photo reaction syntax
             if(ir.ne.2.or.ip.le.0)                        goto 160
             if(wkilenreag(2).ne.6)                        goto 160
             if(wkreag(2)(1:6).ne.'photon')                goto 160
             if(ilenphotoreag(k).gt.nmaxcharspec)          goto 160
             if(ip.gt.nmaxprod)                            goto 160
             do j=1,ip
                if(ilenphotoprod(k,j).gt.nmaxcharspec)     goto 160
             enddo

c      ...verify that species involved are included in the model
             inclu=.true.
             il=ilenphotoreag(k)
             do i=1,nspec
                call same_string
     #               (nmaxcharspec,photoreag(k),spec(i),idem)
                if(idem)then
                   idphotoreag(k)=i
                   exit
                endif
             enddo
             if(idphotoreag(k).eq.0)inclu=.false.
             do j=1,ip
                il=ilenphotoprod(k,j)
                do i=1,nspec
                   call same_string
     #                  (nmaxcharspec,photoprod(k,j),spec(i),idem)
                   if(idem)then
                      idphotoprod(k,j)=i
                      exit
                   endif
                enddo
                if(idphotoprod(k,j).eq.0)inclu=.false.
             enddo
             if(.not.inclu)then
                k=k-1
                cycle
             endif

c      ...check the photo reaction mass and charge balance
             do i=0,nelem
                balr=0
                if(i.eq.0)balr=charge(idphotoreag(k))
                if(i.gt.0)balr=nat(idphotoreag(k),i)
                balp=0
                do j=1,ip
                   if(i.eq.0)balp=balp+charge(idphotoprod(k,j))
                   if(i.gt.0)balp=balp+nat(idphotoprod(k,j),i)
                enddo
                if(balr.ne.balp)                           goto 170
             enddo

c      ...check that the photo reaction did not appear before
             if(k.gt.1)then
                do j=1,k-1
                   if(idphotoreag(k).ne.idphotoreag(j))cycle
                   if(nphotoprod(k).ne.nphotoprod(j))cycle
                   pfnd(1:nmaxprod)=.false.
                   npfnd=0
                   do i=1,nphotoprod(k)
                      idk=idphotoprod(k,i)
                      do l=1,nphotoprod(j)
                         idl=idphotoprod(j,l)
                         if(idk.eq.idl.and..not.pfnd(l))then
                            pfnd(l)=.true.
                            npfnd=npfnd+1
                            exit
                         endif
                      enddo
                   enddo
                   if(npfnd.eq.nphotoprod(k))              goto 180
                enddo
             endif

c      ...convert quantum yield to cross section and resample to the wavelength grid
             do j=1,nxy
                qy(ii,j)=qy(ii,j)*y(j)
             enddo
             call resample_to_wlth(x,qy(ii,:),nxy,csphot(k,:))

          enddo
       enddo
20     close(2)
       nphotospec=s
       nphotoreac=k
       write(*,1030)nphotoreac
       deallocate(x,y,qy)

c      ...build rate equations, i.e. dn_i/dt=...
       do i=1,nspec
          nphototerm(i)=0
          do j=1,nphotoreac
             if(idphotoreag(j).eq.i)then ! add sink terms
                nphototerm(i)=nphototerm(i)+1
                idphototerm(i,nphototerm(i))=j
                sphototerm(i,nphototerm(i))=-1
                fphototerm(i,nphototerm(i))=1
             endif

             do k=1,nphotoprod(j)        ! add source terms
                if(idphotoprod(j,k).eq.i)then
                   nphototerm(i)=nphototerm(i)+1
                   idphototerm(i,nphototerm(i))=j
                   sphototerm(i,nphototerm(i))=1
                   fphototerm(i,nphototerm(i))=1
                   if(k.eq.nphotoprod(j))exit
                   do l=k+1,nphotoprod(j)
                      if(idphotoprod(j,l).eq.i)
     #                fphototerm(i,nphototerm(i))=
     #                fphototerm(i,nphototerm(i))+1
                   enddo
                   exit
                endif
             enddo
          enddo
          if(nphototerm(i).gt.nmaxps)then
             write(*,*)' E- Species ',spec(i)
             stop      '    involved in too many photo reactions'
          endif
       enddo

       return

100    write(*,*)' E- Error reading stellar spectrum file: ',starfile
       stop
110    write(*,*)' E- Error reading photo cross section info file: ',
     # iphotfile
120    write(*,*)' E- Error reading photo cross section file: ',phfile
       stop
130    write(*,*)' E- Many points in photo cross section file: ',phfile
       stop
140    write(*,*)' E- Many photochemical channels in: ',phfile
       stop
150    write(*,*)' E- Different reagents in : ',phfile
       stop
160    write(*,*)' E- In photo reaction: ',photoreac(k)(1:50)
       stop      '    There is an error in the syntax'
170    write(*,*)' E- In photo reaction: ',photoreac(k)(1:50)
       stop      '    Balance not fullfilled'
180    write(*,*)' E- : ',photoreac(j)
       write(*,*)'    : ',photoreac(k)
       stop      '    The above photo reactions are repeated'

1000   format(1x,' I- Reading stellar spectrum from ',a45)
1010   format(1x,' I- Reading photo cross section list from ',a35)
1020   format(1x,' W- Photo-species ',a10,' not included in the model')
1030   format(1x,' I- ',i4,' photo reactions included in the model')

       end




c_______________________________________________________________________

       subroutine same_string(nmaxchar,string1,string2,idem)
c      ...evaluate if two character strings are similar (not-case-sensitive)
c            idem=.false. --> different strings
c            idem=.true.  --> similar strings

       implicit none
       integer nmaxchar,l1,l2,k,l
       character*(nmaxchar) string1,string2,cwk1,cwk2
       logical idem

       idem=.false.

       l1=nmaxchar
       do while(string1(l1:l1).eq.' ')
          l1=l1-1
       enddo
       l2=nmaxchar
       do while(string2(l2:l2).eq.' ')
          l2=l2-1
       enddo
       if(l1.eq.l2)then
          cwk1(1:nmaxchar)=''
          cwk2(1:nmaxchar)=''
          do k=1,l1
             l=ichar(string1(k:k))
             if(l.ge.65.and.l.le.90)l=l+32
             cwk1(k:k)=char(l)
             l=ichar(string2(k:k))
             if(l.ge.65.and.l.le.90)l=l+32
             cwk2(k:k)=char(l)
          enddo
          if(cwk1(1:l1).eq.cwk2(1:l1))idem=.true.
       endif

       return
       end




c_______________________________________________________________________

       subroutine resample_to_wlth(x,y,nxy,yout)
c      ...resample y(x) to wavelenth grid

       include 'pac.common'
       integer nxy,i,j,jl,ju
       real*8 x(*),y(*),yout(nmaxwlth),wl,wu,slope,yl,yu,x1,x2,y1,y2

       do i=1,nwlth
          if(i.eq.1)then
             wl=wlth(i)
             wu=wlth(i)+0.5d0*(wlth(i+1)-wlth(i))
          elseif(i.eq.nwlth)then
             wl=wlth(i)-0.5d0*(wlth(i)-wlth(i-1))
             wu=wlth(i)
          else
             wl=wlth(i)-0.5d0*(wlth(i)-wlth(i-1))
             wu=wlth(i)+0.5d0*(wlth(i+1)-wlth(i))
          endif

          j=1
          do while(x(j).lt.wl.and.j.lt.nxy)
             j=j+1
          enddo
          j=max(j-1,1)
          jl=j
          do while(x(j).lt.wu.and.j.lt.nxy)
             j=j+1
          enddo
          ju=min(max(j-1,1),nxy-1)
          slope=(y(jl+1)-y(jl))/(x(jl+1)-x(jl))
          yl=y(jl)+slope*(wl-x(jl))
          slope=(y(ju+1)-y(ju))/(x(ju+1)-x(ju))
          yu=y(ju)+slope*(wu-x(ju))

          yout(i)=0.0d0
          if(x(jl).ge.wu.or.x(ju+1).le.wl)cycle
          do j=jl,ju
             x1=x(j)
             y1=y(j)
             x2=x(j+1)
             y2=y(j+1)
             if(j.eq.jl.and.x(jl).lt.wl)then
                x1=wl
                y1=yl
             endif
             if(j.eq.ju.and.x(ju+1).gt.wu)then
                x2=wu
                y2=yu
             endif
             yout(i)=yout(i)+0.5d0*(y1+y2)*(x2-x1)
          enddo
          yout(i)=yout(i)/(wu-wl)
       enddo

       return
       end
