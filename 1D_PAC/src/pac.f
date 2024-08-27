
       include 'pac.common'
       integer i,j
       write(*,1000)
       call get_command_argument(1, inpfile)
       open(unit=1,err=100,file=inpfile,status='old')
       
       do i=1,nmaxmodel        
          call reset_variables
          write(*,1020)iday0,imon0,iyea0
          write(*,1001)
          read(1,*,err=110,end=10)inpfile
          if(inpfile(1:3).eq.'END'.or.inpfile(1:3).eq.'end')exit
          call read_input
          call read_spec
          call read_zstruct_initab
          call read_therm
          call read_reac
          call read_photon

          if(.not.sbrotation)then
             call read_eddy
             call solve_continuity_equation
             call write_results
          else
             if(.not.longeddy)then
                call read_eddy
             endif
             call compute_solidbodyrotation
          endif
       enddo
10     close(1)

       write(*,1010)
       stop

100    write(*,*)' E- Error opening pac.inp'
       stop
110    write(*,*)' E- Error reading pac.inp'
       stop


1000   format(/,
     # '---------------------------------------------------------',/,
     # '                       PAC program                       ',/,
     # '           1D Planetary Atmosphere Chemistry             ',/,
     # '                    version March 2023                   ',/)

1001   format(/,
     # '---------------------------------------------------------',/)
     
1010   format(/,
     # '---------------------------------------------------------',/,
     # '                     Calculations done                   ',/,
     # '                                                         ',/,
     # '                        Exiting PAC                      ',/,
     # '---------------------------------------------------------',/)

1020   format(1x,'           Calculations date: ',i2,'/'i2,'/',i4 )


1111   format(1pg12.4)
       end




c_______________________________________________________________________

       subroutine reset_variables

       include 'pac.common'
       integer time_array(8), ind0

c      ...reset CPU time
       call date_and_time(values=time_array)
       iyea0=time_array(1)
       imon0=time_array(2)
       iday0=time_array(3)
       ihou0=time_array(5)
       imin0=time_array(6)
       isec0=time_array(7)

c      ...input files
       inpfile(1:nmaxcharfile)=''
       namemodel(1:nmaxcharfile)=''
       ncharmodel=0
       specfile(1:nmaxcharfile)=''
       zabfile(1:nmaxcharfile)=''
       nreacfile=0
       reacfile(:)(1:nmaxcharfile)=''
       thermfile(1:nmaxcharfile)=''
       eddyfile(1:nmaxcharfile)=''

c      ...model parameters
       ion=.false.
       therm=.false.
       diffusion=.false.
       longeddy=.false.
       photochem=.false.
       iondip=.false.

c      ...physical parameters
       z(:)=0.0d0
       pressure(:)=0.0d0
       tk(:)=0.0d0
       density(:)=0.0d0

c      ...element data
       nelem=0
       elem(:)(1:2)=''
       cdelem(:)=0.0d0
       cdelem0(:)=0.0d0
       elab(:,:)=0.0d0
       elab0(:,:)=0.0d0

c      ...species data
       nspec=0
       spec(:)(1:nmaxcharspec)=''
       dipmom(:)=0.0d0
       polary(:)=0.0d0
       ionp(:)=0.0d0
       charge(:)=0
       nattot(:)=0
       nelemspec(:)=0
       mass(:)=0.0d0
       nat(:,:)=0
       abun(:,:)=0.0d0
       idthspec(:)=0
       nthtemp(:)=0
       thtemp(:,:)=0.0d0
       tha(:,:,:)=0.0d0
       id_h2=0
       id_h=0
       id_he=0
       id_h2o=0
       id_co2=0
       id_ch4=0
       id_electr=0

c      ...reaction data
       nreac=0
       reac(:)(1:50)=''
       alpha(:)=0.0d0
       beta(:)=0.0d0
       gamma(:)=0.0d0
       rtype(:)=0
       tkmin(:)=0.0d0
       tkmax(:)=0.0d0
       tkcode(:)=''
       lerror(:)=''
       ksource(:)=''
       comment(:)=''
       kwarning(:)=.false.
       ilenreag(:,:)=0
       reag(:,:)(1:nmaxcharspec)=''
       idreag(:,:)=0
       ilenprod(:,:)=0
       prod(:,:)(1:nmaxcharspec)=''
       idprod(:,:)=0
       idrvreac(:)=0
       krate(:,:)=0.0d0
       

c      ...collisional reaction data
       nrcoll=0
       idrcoll(:)=0
       nmrcoll(:)=0
       alphainf(:)=0.0d0
       betainf(:)=0.0d0
       gammainf(:)=0.0d0
       afoff(:)=0.0d0
       bfoff(:)=0.0d0
       cfoff(:)=0.0d0
       dfoff(:)=0.0d0
       efoff(:)=0.0d0
       idmrcoll(:,:)=0
       fmrcoll(:,:)=0.0d0

c      ...diffusion data
       hatm(:)=0.0d0
       kdif(:)=0.0d0
       hspe(:,:)=0.0d0
       ddif(:,:)=0.0d0

c      ...photo reaction data
       nphotoreac=0
       photoreac(:)(1:50)=''
       ilenphotoreag(:)=0
       photoreag(:)(1:nmaxcharspec)=''
       idphotoreag(:)=0
       idphotospec(:)=0
       ilenphotoprod(:,:)=0
       photoprod(:,:)(1:nmaxcharspec)=''
       idphotoprod(:,:)=0
       csphot(:,:)=0.0d0
       csphotabs(:,:)=0.0d0
       jphot(:,:)=0.0d0

c      ...stellar flux data
       nwlth=0
       fphot(:)=0.0d0

c      ...continuity equation data
       time=0.0d0
       steady=.false.
       errsteady=0.0d0
       itime=0
       f(:,:)=0.0d0
       df_jm_dif(:,:)=0.0d0
       df_jp_dif(:,:)=0.0d0
       df(:,:,:)=0.0d0

c      ...solid body rotation data
       sbrotation=.false.
       lptfile(1:nmaxcharfile)=''
       windspeed=0.0d0
       nperiod=0
       nlong=0
       long(:)=0.0d0
       tk_grid(:,:)=0.0d0
       kdif_grid(:,:)=0.0d0
       longitude=0.0d0
       tin_sbr=0.0d0
       ostat_sbr=0

       return
       end




c_______________________________________________________________________

       subroutine get_cputime(ihou,imin,isec)

       include 'pac.common'
       integer time_array(8),i,j,cputime,ihou,imin,isec

       call date_and_time(values=time_array)
       j=0
       do i=imon0,time_array(2)-1
          if(i.eq. 1)j=j+31                      ! january
          if(i.eq. 2)j=j+28                      ! february (non leap year)
          if(i.eq. 3)j=j+31                      ! march
          if(i.eq. 4)j=j+30                      ! april
          if(i.eq. 5)j=j+31                      ! may
          if(i.eq. 6)j=j+30                      ! june
          if(i.eq. 7)j=j+31                      ! july
          if(i.eq. 8)j=j+31                      ! august
          if(i.eq. 9)j=j+30                      ! september
          if(i.eq.10)j=j+31                      ! october
          if(i.eq.11)j=j+30                      ! november
          if(i.eq.12)j=j+31                      ! december
       enddo
       cputime=(time_array(1)-iyea0)*31536000+   ! years (non leap year)
     #         (time_array(2)-imon0)*j+          ! months
     #         (time_array(3)-iday0)*86400+      ! days
     #         (time_array(5)-ihou0)*3600+       ! hours
     #         (time_array(6)-imin0)*60+         ! minutes
     #          time_array(7)-isec0              ! seconds
       ihou=cputime/3600
       imin=cputime/60-ihou*60
       isec=mod(cputime,60)

       return
       end
