
       subroutine write_results

       include 'pac.common'
       integer i,j,k,l,ilen,noutfile,i1,i2,n,ihou,imin,isec,aux_rv
       real*8 sum
       character*(nmaxcharfile) outfile,outfile1,outfile2,photofile
       character*(nmaxcharfile*nmaxreacfile) rfiletxt
       character*3 difftxt,phottxt,stetxt
       character*15 numtxt

c      ...prepare character strings to write in header
       rfiletxt(1:nmaxcharfile*nmaxreacfile)=''
       ilen=1
       do i=1,nreacfile
          j=nmaxcharfile
          do while(reacfile(i)(j:j).eq.' ')
             j=j-1
          enddo
          rfiletxt(ilen:ilen+j)=reacfile(i)(1:j)
          ilen=ilen+j+4
       enddo
       difftxt='No '
       if(diffusion)difftxt='Yes'
       phottxt='No '
       if(photochem)phottxt='Yes'
       numtxt(1:15)=''
       if(nummethod.eq.1)then
          numtxt='Newton Raphson '
       elseif(nummethod.eq.2)then
          numtxt='DLSODES        '
          stetxt='No '
          if(steady)stetxt='Yes'
       elseif(nummethod.eq.3)then
          numtxt='Crank-Nicolson '
          stetxt='No '
          if(steady)stetxt='Yes'
       endif
c      ...write cpu time in header
       call get_cputime(ihou,imin,isec)
c      ...write out abundances
       noutfile=nspec/nspecxfile
       if(noutfile*nspecxfile.lt.nspec)noutfile=noutfile+1
       if(noutfile.gt.9) stop 'E- Too many output files'
       n=ncharmodel
       do k=1,noutfile
          outfile(1:nmaxcharfile)=''
          outfile(1:n)=namemodel(1:n)
          outfile1=outfile
          outfile2=outfile
          outfile(n+1:n+13)='_out_00_0.dat'
          outfile1(n+1:n+13)='_che_00_0.dat'
          outfile2(n+1:n+13)='_dif_00_0.dat'

          if(sensi.gt.0)then
            outfile(n+10:n+11)='_i'
            if(irate.lt.10)then
              write(outfile(n+12:n+13),'(i1)')irate
              outfile(n+13:n+17)='.dat'
            endif
            if(irate.ge.10.and.irate.lt.100)then
              write(outfile(n+12:n+14),'(i2)')irate
              outfile(n+14:n+18)='.dat'
            endif
            if(irate.ge.100.and.irate.lt.1000)then
              write(outfile(n+12:n+15),'(i3)')irate
              outfile(n+15:n+19)='.dat'
            endif
            if(irate.ge.1000.and.irate.lt.10000)then
              write(outfile(n+12:n+16),'(i4)')irate
              outfile(n+16:n+20)='.dat'
            endif

          endif
          if(itime.lt.10)then
            write(outfile(n+7:n+7),'(i1)')itime
            write(outfile1(n+7:n+7),'(i1)')itime
            write(outfile2(n+7:n+7),'(i1)')itime
          endif

          if(itime.ge.10)write(outfile(n+6:n+7),'(i2)')itime
          write(outfile(n+9:n+9),'(i1)')k
          
          open(unit=4,file=outfile,status='unknown')
          if(writekrate)then
           open(unit=16,file=outfile1,status='unknown')
           open(unit=17,file=outfile2,status='unknown')
          endif
c         Output for Newton-Raphson
          if(nummethod.eq.1)then
           write(4,1001)specfile,zabfile,nreacfile,rfiletxt,starfile,
     #     iphotfile,dist/au,rstar/rsun,rplanet/rearth,mplanet/mearth,
     #     zenithangle*360.0d0/2.0d0/pi,difftxt,phottxt,numtxt,
     #     errsteady
           write(4,1000)
          endif
c         Output for DLSODES and Crank-Nicholson
          if(nummethod.gt.1)then 
           write(4,1002)specfile,zabfile,nreacfile,rfiletxt,starfile,
     #    iphotfile,dist/au,rstar/rsun,rplanet/rearth,mplanet/mearth,
     #    zenithangle*360.0d0/2.0d0/pi,difftxt,phottxt,numtxt,time,
     #    errsteady,stetxt,ihou,imin,isec
           if(nummethod.eq.2)then
            write(4,1003)rtol,errflag ! Only for DLSODES
           endif
           if(irate.gt.0.and.sensi.gt.0)then ! Adding output relevant for Sensitivity Analysis @Arturo
             write(4,1022)irate,old_alpha(irate),alpha(irate),
     #            drate,rtype(irate)
             if(rvirate.ne.0)write(4,1023)rvirate
           endif
           write(4,1000)
          endif


          i1=nspecxfile*(k-1)+1
          i2=min(nspecxfile*k,nspec)
          
          if(sensi.gt.0.and.irate.gt.0)then
           write(4,1025)(i,i=1,i2-i1+1+5)
           write(4,1035)(spec(i),i=i1,i2)
          else
           write(4,1020)(i,i=1,i2-i1+1+3)
           write(4,1030)(spec(i),i=i1,i2)
           
           if(writekrate)then ! To check the right side of equations @Arturo
            write(16,1020)(i,i=1,i2-i1+1+3)
            write(16,1030)(spec(i),i=i1,i2)
            write(17,1020)(i,i=1,i2-i1+1+3)
            write(17,1030)(spec(i),i=i1,i2)
           endif
          endif

          do j=nz,1,-1
             sum=0.0d0                           ! normalize abundances
             do i=1,nspec
                sum=sum+abun(i,j)
             enddo
             if(sensi.gt.0.and.irate.gt.0)then
              aux_rv=irate
              if(rvirate.gt.0)aux_rv=rvirate
              write(4,1045)z(j)/1.0d5,pressure(j)/1.0d6,tk(j),
     #       dlog10(krate(irate,j)),dlog10(krate(aux_rv,j)),
     #       (dlog10(max(abunmin,abun(i,j)/sum)),i=i1,i2)
             else
              write(4,1040)z(j)/1.0d5,pressure(j)/1.0d6,tk(j),
     #       (dlog10(max(abunmin,abun(i,j)/sum)),i=i1,i2)
              if(writekrate)then
               write(16,1041)z(j)/1.0d5,pressure(j)/1.0d6,tk(j),
     #       (chem_ctbn(i,j),i=i1,i2)
               write(17,1041)z(j)/1.0d5,pressure(j)/1.0d6,tk(j),
     #       (diff_ctbn(i,j),i=i1,i2)
              endif

             endif
          enddo
          close(4)
          if(writekrate)then
           close(16)
           close(17)
          endif
       enddo

c      ...write location where opacity is 1 as a function of wavelength
       if(photochem)then
          photofile=outfile
          photofile(1:nmaxcharfile)=''
          photofile(1:n)=namemodel(1:n)
          photofile(n+1:n+9)='_tau1.dat'
          open(unit=7,file=photofile,status='unknown')
          write(7,1050)time
          do l=1,nwlth
             write(7,1060)wlth(l),p_tau1(l)/1.0d6,z_tau1(l)/1.0d5
          enddo
          close(7)
       endif

c      ...write coolin/heating rates
       call cool_heat_rates

       return

1000   format( 
     # '! .......................................................',/,
     # '! Mixing ratios : ')

1001   format(
     # '!',/,
     # '! 1D atmosphere chemical model with parameters:',/,
     # '! species file : ',a,/,
     # '! (z,p,T) structure and initial abundances file : ',a,/,
     # '! reaction rate constant files = ',i2,' : ',a,/,
     # '! stellar spectrum file : ',a,/,
     # '! photo cross section processes list file : ',a,/,
     # '! Star-planet distance                =',1pg12.4,' AU',/,
     # '! Star radius                         =',1pg12.4,' Rsun',/,
     # '! Planet radius                       =',1pg12.4,' REarth',/,
     # '! Planet mass                         =',1pg12.4,' MEarth',/,
     # '! Zenith angle                        =',1pg12.4,' deg',/,
     # '! Diffusion included ?                : ',a,/,
     # '! Photochemistry included ?           : ',a,/,
     # '! Numerical method used               : ',a,/,
     # '! Relative convergence error          =',1pg12.4,/)

1002   format(
     # '!',/,
     # '! 1D atmosphere chemical model with parameters:',/,
     # '! species file : ',a,/,
     # '! (z,p,T) structure and initial abundances file : ',a,/,
     # '! reaction rate constant files = ',i2,' : ',a,/,
     # '! stellar spectrum file : ',a,/,
     # '! photo cross section processes list file : ',a,/,
     # '! Star-planet distance                 =',1pg12.4,' AU',/,
     # '! Star radius                          =',1pg12.4,' Rsun',/,
     # '! Planet radius                        =',1pg12.4,' REarth',/,
     # '! Planet mass                          =',1pg12.4,' MEarth',/,
     # '! Zenith angle                         =',1pg12.4,' deg',/,
     # '! Diffusion included ?                 : ',a,/,
     # '! Photochemistry included ?            : ',a,/,
     # '! Numerical method used                : ',a,/,
     # '! Continuity equation integrated up to =',1pg12.4,' s ',/,
     # '!   with molar fraction relative error =',1pg12.4,/,
     # '! Steady state reached ?               : ',a,/,
     # '! CPU time: ',i3,'h 'i2,'m ',i2,'s')

1003   format('! RTOL for DLSODES:',1pg12.4,/,
     # '! Conservation error?: ',i4)

1022   format(
     # '! SENSITIVITY ANALYSIS ACTIVATED in reaction: ',i4,/,
     # '! Changed alpha from :',1pg12.4,' to ',1pg12.4,/,
     # '! Variation of reaction: ',1pg12.4,/,
     # '! Reaction type ID: ', i4)

1023   format(
     # '! Also changed in reaction:',i4)


1020   format('!',3(7x,i2,7x),4x,999(3x,i3,9x))
1025   format('!',5(7x,i2,7x),4x,999(3x,i3,9x))

1030   format('!','   height[km]   ',' pressure[bar]  ',
     # ' temperature[K] ',7x,999(a15))

1035   format('!','   height[km]   ',' pressure[bar]  ',
     # ' temperature[K] ', '     kf_varied     ', '  kr_varied  ', 
     # 7x,999(a15))
1040   format(3x,1pg12.4,4x,1pg12.4,4x,0pf12.3,4x,999(3x,0pf12.8))
1041   format(3x,1pg12.4,4x,1pg12.4,4x,0pf12.3,4x,999(3x,1pg12.5))
1045   format(3x,1pg12.4,4x,1pg12.4,4x,0pf12.3,4x,1pg12.4,4x,
     # 1pg12.4,4x,999(3x,0pf12.8))
1050   format('! tau=1 layer at time:',1pg12.4,' s',/,
     #        '!','  wlth[nm]   ',
     #         1x,'pressure[bar]',
     #         1x,' height[km]  ')
1060   format(2x,1pg12.4,2x,1pg12.4,2x,1pg12.4)




       end




c_______________________________________________________________________

       subroutine cool_heat_rates

       include 'pac.common'
       integer i,j,l,id,n
       real*8 tkin,sum,dh,h,s,cp,dn,chrate(nmaxz),cpmean
       character*(nmaxcharfile) outfile

c Loop over height
       do j=1,nz
          sum=0.0d0                              ! normalize abundances
          do i=1,nspec
             sum=sum+abun(i,j)
          enddo
          tkin=tk(j)
          chrate(j)=0.0d0

c ---  Loop over chemical reactions: add chemical contribution to cooling/heating rate
          do i=1,nreac
             dh=0.0d0 ! DeltaH/RT                !...evaluate reaction enthalpy
             do l=1,nreag(i)
                id=idthspec(idreag(i,l))
                call evaluate_therm(id,tkin,h,s,cp)
                dh=dh-h
             enddo
             do l=1,nprod(i)
                id=idthspec(idprod(i,l))
                call evaluate_therm(id,tkin,h,s,cp)
                dh=dh+h
             enddo

             dn=krate(i,j)                       !...evaluate reaction rate
             do l=1,nreag(i)
                dn=dn*max(abunmin,abun(idreag(i,l),j)/sum)*density(j)
             enddo

             chrate(j)=chrate(j)-                ! add contribution: cooling(-), heating(+)
     #       dn*dh*kboltz*tkin                   ! [=] erg cm-3 s-1
          enddo

c ---  Loop over photochemical reactions:  add photochemical contribution to cooling/heating rate
c      A break of a molecule by photodissociation does not cool the atmosphere
c      since the heat taken for breaking the bond(s) comes from UV photons
c          do i=1,nphotoreac
c             dh=0.0d0 ! DeltaH/RT                !...evaluate photoreaction enthalpy
c             id=idthspec(idphotoreag(i))
c             call evaluate_therm(id,tkin,h,s,cp)
c             dh=dh-h
c             do l=1,nphotoprod(i)
c                id=idthspec(idphotoprod(i,l))
c                call evaluate_therm(id,tkin,h,s,cp)
c                dh=dh+h
c             enddo
c
c             dn=jphot(i,j)*                      !...evaluate photoreaction rate
c     #       max(abunmin,abun(idreag(i,l),j)/sum)
c
c             chrate(j)=chrate(j)-                ! add contribution: cooling(-), heating(+)
c     #       dn*dh*kboltz*tkin                   ! [=] erg cm-3 s-1
c          enddo

c      ...convert cooling/heating rate to units of [K s-1]
          cpmean=0.0d0 ! Cp/R                    ! evaluate mean heat capacity
          do i=1,nspec
             id=idthspec(i)
             call evaluate_therm(id,tkin,h,s,cp)
             cpmean=cpmean+cp*max(abunmin,abun(i,j)/sum)
          enddo
          chrate(j)=chrate(j)/
     #    (cpmean*kboltz*density(j))             ! [=] K s-1
       enddo

c      ...write out cooling/heating rate
       n=ncharmodel
       outfile(1:nmaxcharfile)=''
       outfile(1:n)=namemodel(1:n)
       outfile(n+1:n+11)='_chrate.dat'
       open(unit=9,file=outfile,status='unknown')
       write(9,1000)time
       do j=nz,1,-1
          write(9,1010)z(j)/1.0d5,pressure(j)/1.0d6,tk(j),chrate(j)
       enddo
       close(9)

1000   format('! cooling/heating rate at time:',1pg12.4,' s',/,
     #        '! cooling (-) and heating (+)',/,
     #        '!','   height[km]   ',
     #            ' pressure[bar]  ',
     #            ' temperature[K] ',
     #            ' cool/heat rate[K s-1]')
1010   format(3x,1pg12.4,4x,1pg12.4,4x,0pf12.3,4x,1pg12.4)

       end
