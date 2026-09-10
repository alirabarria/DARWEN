
       subroutine get_rate_constants
c      Get rate constants of chemical reactions
       USE IEEE_ARITHMETIC
       include 'pac.common'
       integer i,j,k,l,dv,iwk,id,thisid,ref,nrefpr,seedy,iseedy
       real*8 tkin,dh,ds,keq,mneu,dneu,pneu,klangevin,twk,kfr,x,
     # feff(nmaxspec),k0,kinf,m,fcent,c,n,d,fk,kmax,h,s,cp,klow,
     # khigh,logk,o_krate_i,aux_alpha,norm,mu,sigma,
     # random_normal
       character*(nmaxcharfile) kwarningfile
       character*(nmaxcharfile) kratefile !to print krate matrix @Arturo

c         Vary reaction constant if Sensitivity analysis is needed @Arturo
       if(sensi.gt.0.and.irate.gt.0)then
         write(*,*)'Checking sensitivity at rate',irate
         if(first_vary)then ! Only saving alpha before first variation @Arturo
          old_alpha=alpha
          old_alphaplog=alphaplog
          old_alphainf=alphainf
         endif
         if(change_rate)then
          alpha=old_alpha ! Restoring alpha to the way it was before @Arturo
          alphaplog=old_alphaplog
          alphainf=old_alphainf

          aux_alpha=old_alpha(irate) 
          alpha(irate)=(1+drate)*aux_alpha
         endif
         write(*,*)'Old alpha',old_alpha(irate),
     #   'New alpha',alpha(irate)
       endif


c Loop over height
       do j=1,nz

c Loop over reactions to evaluate rate constants
          do i=1,nreac

c         compute rate constant of reverse with therm data
             if(rtype(i).eq.0)then
                tkin=tk(j)
                dh=0.0d0             ! DeltaH/RT
                ds=0.0d0             ! DeltaS/R
                do l=1,nreag(i)
                   id=idthspec(idreag(i,l))
                   call evaluate_therm(id,tkin,h,s,cp)
                   dh=dh-h
                   ds=ds-s
                enddo
                do l=1,nprod(i)
                   id=idthspec(idprod(i,l))
                   call evaluate_therm(id,tkin,h,s,cp)
                   dh=dh+h
                   ds=ds+s
                enddo
                dv=nprod(i)-nreag(i)
                keq=(1.0d6/kboltz/tkin)**dv*dexp(ds-dh)
                krate(i,j)=keq*krate(idrvreac(i),j)
c                if(idrvreac(i).eq.1)then
c                   write(*,*)reac(i),z(j)
c                   write(*,*)ds,dh
c                   write(*,*)krate(idrvreac(i),j),krate(i,j)
c                   write(*,*)keq
c                endif

                if(idrvreac(i).eq.irate.and.change_rate)then
                   rvirate=0 !!!Changed rvirate
                   if(rtype(idrvreac(i)).lt.0)rvirate=i
                   if(j.eq.1)write(*,*)'Also changed rate:',i
                   krate(i,j)=(1+drate)*krate(i,j)
                endif
                cycle
             endif

c         compute rate constant of chemical reaction with kinetic data
             if(((tkcode(i).eq.'<zer'.or.tkcode(i).eq.'<>zr').and.
     #       tk(j).le.tkmin(i)).or.
     #       ((tkcode(i).eq.'>zer'.or.tkcode(i).eq.'<>zr').and.
     #       tk(j).gt.tkmax(i)))then
                krate(i,j)=0.0d0
                cycle
             endif
             if((tkcode(i).eq.'<fix'.or.tkcode(i).eq.'<>fx').and.
     #       tk(j).lt.tkmin(i))then
                tkin=tkmin(i)
             elseif((tkcode(i).eq.'>fix'.or.tkcode(i).eq.'<>fx').and.
     #       tk(j).gt.tkmax(i))then
                tkin=tkmax(i)
             else
                tkin=tk(j)
             endif
             krate(i,j)=alpha(i)*(tkin/300.0d0)**beta(i)*
     #       dexp(-gamma(i)/tkin)

c         update rate constant if interpolation of PLOG reaction is needed
             if(iabs(rtype(i)).gt.30.and.iabs(rtype(i)).lt.45)then ! @Arturo changed the maximum pressures limits
                thisid=idrplog(i)                                   ! plog-reaction identifier
                nrefpr=iabs(rtype(i))-30                            ! number of reference pressures
                if(pressure(j).le.prefplog(thisid,1))then            ! minimal pressure case

                 if(change_rate.and.i.eq.irate)then ! Varying PLOG reactions @Arturo
                  alphaplog=old_alphaplog
                  aux_alpha=alphaplog(thisid,1)
                  alphaplog(thisid,1)=(1+drate)*aux_alpha
                 endif

                   krate(i,j)=alphaplog(thisid,1)*
     #             (tkin/300.0d0)**betaplog(thisid,1)*
     #             dexp(-gammaplog(thisid,1)/tkin)

                if(IEEE_IS_NAN(krate(i,j)))then
                  write(*,*)'NaN found',krate(i,j),i,j
                  write(*,*)'alphaplog',alphaplog(thisid,1)
                  write(*,*)'betaplog',betaplog(thisid,1)
                  write(*,*)'gammaplog',gammaplog(thisid,1)
                  write(*,*)'tkin',tkin
                  write(*,*)(tkin/300.0d0)**betaplog(thisid,1)
                  write(*,*)'dexp',dexp(-gammaplog(thisid,1)/tkin)
                  errflag=2  ! @Arturo, I want this to stop if NaNs are found?
                  write(*,*)
                endif

                elseif(pressure(j).ge.prefplog(thisid,nrefpr))then  ! maximal pressure case

                  if(change_rate.and.i.eq.irate)then ! Varying PLOG reactions @Arturo
                   alphaplog=old_alphaplog
                   aux_alpha=alphaplog(thisid,nrefpr)
                   alphaplog(thisid,nrefpr)=(1+drate)*aux_alpha
                   write(*,*)aux_alpha,alphaplog(thisid,nrefpr)
                  endif

                   krate(i,j)=alphaplog(thisid,nrefpr)*
     #             (tkin/300.0d0)**betaplog(thisid,nrefpr)*
     #             dexp(-gammaplog(thisid,nrefpr)/tkin)
                  
                   


c                   if(thisid.eq.1)then
c                      write(*,*) '--> above ',prefplog(thisid,nrefpr)
c                   endif
                else
                   do ref=1,nrefpr-1
                      if(pressure(j).le.prefplog(thisid,ref+1))then ! pressure in between cases
                        if(change_rate.and.i.eq.irate)then !! Varying PLOG reactions @Arturo
                         alphaplog=old_alphaplog
                         aux_alpha=alphaplog(thisid,ref)
                         alphaplog(thisid,ref)=(1+drate)*aux_alpha

                         aux_alpha=alphaplog(thisid,ref+1)
                         alphaplog(thisid,ref+1)=(1+drate)*aux_alpha

                        endif

                         klow=alphaplog(thisid,ref)*                ! --> lower boundary
     #                   (tkin/300.0d0)**betaplog(thisid,ref)*
     #                   dexp(-gammaplog(thisid,ref)/tkin)
                         khigh=alphaplog(thisid,ref+1)*             ! --> upper boundary
     #                   (tkin/300.0d0)**betaplog(thisid,ref+1)*
     #                   dexp(-gammaplog(thisid,ref+1)/tkin)
                         logk=dlog(klow)+(dlog(khigh)-dlog(klow))*  ! logarithmic interpolation
     #                   (dlog(pressure(j))-dlog(prefplog(thisid,ref)))/
     #                   (dlog(prefplog(thisid,ref+1))-
     #                   dlog(prefplog(thisid,ref)))
                         krate(i,j)=dexp(logk)                      ! reaction rate
c                         if(thisid.eq.1)then
c                            write(*,*)'--> between ',
c     #                      prefplog(thisid,ref),' and ',
c     #                      prefplog(thisid,ref+1)
c                         endif
                         exit
                      endif
                   enddo
                endif
c                do ref=1,nrefpr
c                   write(*,*)prefplog(thisid,ref),alphaplog(thisid,ref),
c     #             betaplog(thisid,ref),gammaplog(thisid,ref)
c                enddo
             endif

c         update rate constant if needed if special treatment for ion-neutral reactions
             if(iondip.and.
     #       iabs(rtype(i)).ge.400.and.iabs(rtype(i)).lt.420)then
                if(charge(idreag(i,1)).ne.0.and.
     #          charge(idreag(i,2)).eq.0)then
                   iwk=2
                elseif(charge(idreag(i,1)).eq.0.and.
     #          charge(idreag(i,2)).ne.0)then
                   iwk=1
                else
                   write(*,*)' E- ',reac(i)
                   stop      '    wrong rtype 400-420 in above'
                endif
                mneu=mass(idreag(i,iwk))
                dneu=dipmom(idreag(i,iwk))
                pneu=polary(idreag(i,iwk))
                if(dneu.lt.-990.0d0)cycle
                if(pneu.le.0.0d0)cycle
                klangevin=2.3420d-9*dsqrt(pneu/mneu)
                if(iabs(rtype(i)).eq.400)then    ! scale with experimental data
                   twk=0.5d0*(tkmin(i)+tkmax(i))
                   x=3.4744d0*dneu*dsqrt(300.0d0/twk)/dsqrt(pneu)
                   if(x.ge.2.0d0)then
                      kfr=klangevin*(0.4767d0*x+0.6200d0)
                   else
                      kfr=klangevin*
     #                (((x+0.5090d0)**2.0d0)/10.526d0+0.9754d0)
                   endif
                   kfr=alpha(i)/kfr
                   if(kfr.gt.1.0d0.and..not.kwarning(i))then
                      kwarning(i)=.true.
                      if(warnreac.eq.1)then ! @Arturo hid warnings
                       write(7,1000)reac(i),idint(twk),idint(twk),
     #                 alpha(i),alpha(i)/kfr,reac(i)
                      endif
                   endif
                else                             ! experimental data not available
                   kfr=1.0d0/dble(iabs(rtype(i))-400)
                endif
                if(dneu.lt.0.1d0)then            !_ion-non polar
                   krate(i,j)=kfr*klangevin
                else                             !_ion-dipole
                   x=3.4744d0*dneu*dsqrt(300.0d0/tkin)/dsqrt(pneu)
                   if(x.ge.2.0d0)then
                      krate(i,j)=kfr*klangevin*(0.4767d0*x+0.6200d0)
                   else
                      krate(i,j)=kfr*klangevin*
     #                (((x+0.5090d0)**2.0d0)/10.526d0+0.9754d0)
                   endif
                endif
             endif

c         update rate constant if collisional reaction
             if(iabs(rtype(i)).eq.800.or.        ! three-body or thermal dissociation
     #          iabs(rtype(i)).eq.810.or.        !  - general form Troe (k_0/k_infinite)
     #          iabs(rtype(i)).eq.811.or.        !                      (k_0/k_infinite + M=specific)
     #          iabs(rtype(i)).eq.812.or.        !  - general form SRI  (k_0/k_infinite)
     #          iabs(rtype(i)).eq.900.or.
     #          iabs(rtype(i)).eq.910.or.
     #          iabs(rtype(i)).eq.911)then

                feff(:)=1.0d0                    ! ...evaluate k_0, k_infinite, [M]
                id=idrcoll(i)
                do k=1,nmrcoll(id)
                   feff(idmrcoll(id,k))=fmrcoll(id,k)
                enddo
                k0=alpha(i)*(tkin/300.)**beta(i)*dexp(-gamma(i)/tkin)

                if(change_rate.and.i.eq.irate)then ! Varying alphainf reactions @Arturo
                  alphainf=old_alphainf
                  aux_alpha=alphainf(id)
                  alphainf(id)=(1+drate)*aux_alpha
                endif

                kinf=alphainf(id)*(tkin/300.)**betainf(id)*
     #          dexp(-gammainf(id)/tkin)
                m=0.0d0
                do k=1,nspec
                   m=m+feff(k)*abun(k,j)*density(j)
                enddo
                if(iabs(rtype(i)).eq.811)        ! ...o2 + h (+h2o) = ooh (+h2o)
     #             m=abun(id_h2o,j)*density(j)
                if(iabs(rtype(i)).eq.812)        ! ...m = co2
     #             m=abun(id_co2,j)*density(j)
                if(iabs(rtype(i)).eq.911)        ! ...ch4 (+ch4) = ch3 + h (+ch4)
     #             m=abun(id_ch4,j)*density(j)

                if(iabs(rtype(i)).eq.800.or.     ! ...evaluate F under Troe formulation
     #             iabs(rtype(i)).eq.811.or.
     #             iabs(rtype(i)).eq.812.or.
     #             iabs(rtype(i)).eq.900.or.
     #             iabs(rtype(i)).eq.911)then

                   fcent=(1.0d0-afoff(id))*dexp(-tkin/bfoff(id))+
     #             afoff(id)*dexp(-tkin/cfoff(id))+
     #             dexp(-dfoff(id)/tkin)
                   fcent=max(1.0d-300,fcent)
                   c=-0.4d0-0.67d0*dlog10(fcent)
                   n=0.75d0-1.27d0*dlog10(fcent)
                   d=0.14d0
                   if(m.gt.0.0d0)then
                      fk=dlog10(fcent)/(1.0d0+((dlog10(k0*m/kinf)+c)/
     #                (n-d*dlog10(k0*m/kinf)+c))**2.0d0)
                      fk=10.0d0**fk
                   else
                      fk=0.0d0                   ! ...avoid NaN if [M]=0
                   endif
                elseif(iabs(rtype(i)).eq.810.or. ! ...evaluate F under SRI formulation
     #                 iabs(rtype(i)).eq.910)then
                   x=1.0d0/(1.0d0+(dlog10(k0*m/kinf))**2.0d0)
                   fk=dfoff(id)*tkin**efoff(id)*(afoff(id)*
     #             dexp(-bfoff(id)/tkin)+dexp(-tkin/cfoff(id)))**x
                else
                   write(*,*)' E- ',reac(i)
                   stop      '    Need of factor F for above reaction'
                endif

                krate(i,j)=k0*m/(1.0d0+k0*m/kinf)*fk  ! ...evaluate rate constant


             elseif(iabs(rtype(i)).eq.801.or.    ! three-body or thermal dissociation
     #              iabs(rtype(i)).eq.901)then   !  - k_0 + M enhancements
                feff(:)=1.0d0
                id=idrcoll(i)
                do k=1,nmrcoll(id)
                   feff(idmrcoll(id,k))=fmrcoll(id,k)
                enddo
                k0=alpha(i)*(tkin/300.0)**beta(i)*dexp(-gamma(i)/tkin)
                m=0.0d0
                do k=1,nspec
                   m=m+feff(k)*abun(k,j)*density(j)
                enddo
                krate(i,j)=k0*m
c                                                !  - k_infinite (computed above as a normal reaction)
             
             endif


             if(IEEE_IS_NAN(krate(i,j)))then
              write(*,*)'NaN found',krate(i,j),i,j
              write(*,*)'reac ',i,reac(i),rtype(i)
              errflag=2  ! @Arturo, I want this to stop if NaNs are found?
              write(*,*)
             endif



c End of loop over reactions to evaluate rate constants
          enddo

c End of loop over height
       enddo

c      ...if rate constant is too high, impose upper limits and write to warning file
       kwarningfile(1:nmaxcharfile)=''
       kwarningfile(1:ncharmodel)=namemodel(1:ncharmodel)
       kwarningfile(ncharmodel+1:ncharmodel+18)='_krate_warning.dat'
       open(unit=2,file=kwarningfile,status='unknown')
       do j=1,nz
          do i=1,nreac
             kmax=huge
             if(nreag(i).eq.1)kmax=kratemax1
             if(nreag(i).eq.2)kmax=kratemax2
             if(nreag(i).eq.3)kmax=kratemax3
             if(krate(i,j).gt.kmax)then
                kwarning(i)=.true.
                if(warnreac.eq.1)then ! @Arturo hid warnings
                 write(2,1010)reac(i),idint(tk(j)),krate(i,j)
                endif
                if(rtype(i).eq.0)
     #          krate(idrvreac(i),j)=krate(idrvreac(i),j)*
     #          kmax/krate(i,j)
                krate(i,j)=kmax
             endif

             if(sensi.gt.0.and.i.eq.irate)then ! write reaction to warning file @Arturo
               write(2,1020)reac(i),idint(tk(j)),krate(i,j)
             endif
          enddo
       enddo
       close(2)
c      write krate to file for sensitivity analysis @Arturo
       if(writekrate)then
         kratefile(1:nmaxcharfile)=''
         kratefile(1:ncharmodel)=namemodel(1:ncharmodel)
         kratefile(ncharmodel+1:ncharmodel+18)='_krate.out'
         open(unit=10,file=kratefile,status='unknown')
         write(10,1030)(reac(i), i=1,nreac)
         do j=1,nz
            write(10,1040) (krate(i,j),i=1,nreac)
         enddo
         close(10)
       endif


       return

1000   format(1x,'W- ',a50,' k_measured(',i4',K) > k_collision(',i4,
     # 'K)','     : ',1pg9.2,' > ',1pg9.2)
1010   format(1x,'W- ',a50,' k(',i4',K) = ',1pg9.2)
1020   format(1x,'S-',a50,' k(',i4',K) = ',1pg9.2)
1030   format('!  ',2000(a50))
1040   format(2000(2x,1pg12.6))
1050   format(1x,'Reaction ',i4,',',i4,' changed from ',1pg9.2,
     # ' to ',1pg9.2)
1060   format(1pg10.10)

       end




c_______________________________________________________________________

       subroutine get_diffusion_terms
c      Get diffusion coefficients, and other diffusion stuff

       use pac_data
       include 'pac.common'
       real*8 p_a,p_b,p_c,p_d,p_e,p_f,p_g,p_h
       parameter(p_a=1.06036d0,p_b=0.15610d0,p_c=0.19300d0,
     # p_d=0.47635d0,p_e=1.03587d0,p_f=1.52996d0,p_g=1.76474d0,
     # p_h=3.89411d0)

       integer i,j,k,l,h
       real*8 norm,s,lj_lena,lj_lenb,lj_enea,lj_eneb,omega,musqrt,
     # lj_len2,ta,atd,r,dz,rm,rp,dzm,dzp,nm,np,tkm,tkp,hspem,hspep,
     # ddifm,ddifp,kdifm,kdifp,dtkm,dtkp,fdm,fdp,hatmm,hatmp

       logical lj_fnda,lj_fndb,idem

       if(.not.diffusion)return

c      ...to avoid warning
       lj_enea=0.0d0
       lj_eneb=0.0d0

c      ...evaluate atmospheric scale height
       call get_scale_height

c      ...evaluate coefficients of molecular diffusion
       do j=1,nz
          do i=1,nspec                            !......main species a
             lj_fnda=.false.
             do l=1,nmaxlj
                call same_string
     #          (nmaxcharspec,spec(i),lj_spec(l),idem)
                if(idem)then
                   lj_fnda=.true.
                   lj_lena=lj_length(l)*1.0d-8    ! A -> cm
                   lj_enea=lj_energy(l)
                   exit
                endif
             enddo
             if(.not.lj_fnda)lj_lena=4.0d-8       ! assumed characteristic length (4 A)
             s=0.0d0
             norm=0.0d0
             do k=1,nspec                         !......second collider b
                lj_fndb=.false.
                do l=1,nmaxlj
                   call same_string
     #             (nmaxcharspec,spec(k),lj_spec(l),idem)
                   if(idem)then
                      lj_fndb=.true.
                      lj_lenb=lj_length(l)*1.0d-8 ! A -> cm
                      lj_eneb=lj_energy(l)
                      exit
                   endif
                enddo
                if(.not.lj_fndb)lj_lenb=4.0d-8    ! assumed characteristic length (4 A)
                if(lj_fnda.and.lj_fndb.and.
     #          lj_enea.gt.0.0.and.lj_eneb.gt.0.0)then
                   ta=tk(j)/dsqrt(lj_enea*lj_eneb)
                   omega=p_a/(ta**p_b)+p_c/(p_d*ta)+
     #             p_e/(dexp(p_f*ta))+p_g/(dexp(p_h*ta))
                else
                   omega=1.0d0
                endif
                musqrt=dsqrt(mass(i)*mass(k)/(mass(i)+mass(k)))
                lj_len2=(0.5d0*(lj_lena+lj_lenb))**2.0d0
                s=s+abun(k,j)*musqrt*lj_len2*omega
                norm=norm+abun(k,j)
             enddo
             s=s*density(j)/norm
             ddif(i,j)=3.0d0/16.0d0*dsqrt(2.0d0*kboltz*tk(j)/pi)/s
          enddo
       enddo
       if(nummethod.eq.1)ddif(:,:)=0.0d0         ! neglect molecular diffusion if steady state Newton-Raphson solution

c      ...compute diffusion terms to evaluate function f in continuity equation
       do j=1,nz
          r=rplanet+z(j)
          if(j.eq.1)then
             dz=z(j)-z(j+1)
          elseif(j.eq.nz)then
             dz=z(j-1)-z(j)
          else
             dz=0.5d0*(z(j-1)-z(j+1))
          endif

          h=max(2,j)                             ! evaluate at j-1/2 (if j=1, then values are not used)
          rm=rplanet+0.5d0*(z(h-1)+z(h))
          dzm=z(h-1)-z(h)
          nm=0.5d0*(density(h-1)+density(h))
          tkm=0.5d0*(tk(h-1)+tk(h))
          dtkm=(tk(h-1)-tk(h))/(z(h-1)-z(h))
c          if(kdif(h).ne.0.0)write(*,*)kdif(h)
          kdifm=0.5d0*(kdif(h-1)+kdif(h))
          hatmm=0.5d0*(hatm(h-1)+hatm(h))
          h=min(j,nz-1)                          ! evaluate at j+1/2 (if j=nz, then values are not used)
          rp=rplanet+0.5d0*(z(h)+z(h+1))
          dzp=z(h)-z(h+1)
          np=0.5d0*(density(h)+density(h+1))
          tkp=0.5d0*(tk(h)+tk(h+1))
          dtkp=(tk(h)-tk(h+1))/(z(h)-z(h+1))
          kdifp=0.5d0*(kdif(h)+kdif(h+1))
          hatmp=0.5d0*(hatm(h)+hatm(h+1))

          do i=1,nspec
             if(i.eq.id_h2.or.i.eq.id_h.or.i.eq.id_he)then
                atd=-0.25d0                      ! thermal diffusion factor set to -0.25 for the light species H,H2,He (Bauer 1973, Physics of Planetary Ionospheres, p 10)
             else
                atd=0.0d0                        ! thermal diffusion factor set to zero for any other species
             endif

             h=max(2,j)                          ! evaluate at j-1/2 (if j=1, then values are not used)
             hspem=0.5d0*(hspe(i,h-1)+hspe(i,h))
             ddifm=0.5d0*(ddif(i,h-1)+ddif(i,h))
             fdm=1.0d0/hspem-1.0d0/hatmm+atd*dtkm/tkm
             h=min(j,nz-1)                       ! evaluate at j+1/2 (if j=nz, then values are not used)
             hspep=0.5d0*(hspe(i,h)+hspe(i,h+1))
             ddifp=0.5d0*(ddif(i,h)+ddif(i,h+1))
             fdp=1.0d0/hspep-1.0d0/hatmp+atd*dtkp/tkp

             df_dif(i,j)=0.0d0                   ! evaluate d[f_i^j]/d(y_i^j)
             if(j.gt.1)df_dif(i,j)=df_dif(i,j)+
     #       (rm/r)**2.0d0*nm/dz*
     #       (ddifm*(0.5d0*fdm-1.0d0/dzm)-kdifm/dzm)
             if(j.lt.nz)df_dif(i,j)=df_dif(i,j)-
     #       (rp/r)**2.0d0*np/dz*
     #       (ddifp*(0.5d0*fdp+1.0d0/dzp)+kdifp/dzp)

             df_jm_dif(i,j)=0.0d0                ! evaluate d[f_i^j]/d(y_i^j-1)
             if(j.gt.1)df_jm_dif(i,j)=(rm/r)**2.0d0*nm/dz*
     #       (ddifm*(0.5d0*fdm+1.0d0/dzm)+kdifm/dzm)

             df_jp_dif(i,j)=0.0d0                ! evaluate d[f_i^j]/d(y_i^j+1)
             if(j.lt.nz)df_jp_dif(i,j)=-(rp/r)**2.0d0*np/dz*
     #       (ddifp*(0.5d0*fdp-1.0d0/dzp)-kdifp/dzp)
          enddo
       enddo

       return
       end




c_______________________________________________________________________

       subroutine get_photorates
c      Get photorates doing the ultraviolet radiative transfer
c      with absorption due to photodissociations and Rayleigh scattering by all species
c      We use a two-stream iterative algorithm assuming that 1/2 of scattered photons go up and 1/2 down
c      Algorithm from Isaksen et al 1977 (see F. Selsis thesis) slightly modified

       include 'pac.common'
       integer nmaxit
       parameter(nmaxit=1000)
       real*8 varmin
       parameter(varmin=1.0d-2)                  ! 1 %

       integer i,j,k,l,id,it
       real*8 csrayleigh(nmaxspec),ta,ts,n,dz,sigma,ftop(nmaxwlth),
     # down,up,fx,o_fx,var,w1,w2,s1,s2,f1,f2,tau,tu,tl,pu,pl,slope,
     # evening,morning,r,path,plusangle
       real*8, allocatable, dimension(:,:) :: tau_abs,tau_sca,
     #                                        fdown,o_fdown,fup,o_fup

       if(.not.photochem)return
       if(sbrotation)then                        ! evaluate zenith angle
          plusangle=datan((rstar-rplanet)/dist)+0.5d0*pi/180.0d0
          evening=0.5d0*pi+plusangle             ! at limbs include the effect of the star apparent size
          morning=1.5d0*pi-plusangle             ! ...plus half degree due to atmospheric refraction (as in the Earth)
          if(longitude.lt.evening.or.longitude.gt.morning)then
             if(longitude.lt.evening)zenithangle=longitude
             if(longitude.gt.morning)zenithangle=2.0d0*pi-longitude
             zenithangle=min(zenithangle,0.5d0*pi-plusangle)
          else
             jphot(1:nphotoreac,1:nz)=0.0d0
             return
          endif
       endif

c      ...evaluate Rayleigh scattering cross section (part independent of wavelength)
       do i=1,nspec
          csrayleigh(i)=128.0d0/3.0d0*pi**5.0d0*
     #    (polary(i)*1.0d-24)**2.0d0
       enddo

c      ...evaluate flux above top layer assuming isothermal atmosphere and constant molar fractions
       do l=1,nwlth
          ta=0.0d0                               ! absorption
          do k=1,nphotospec
             id=idphotospec(k)
             n=abun(id,1)*density(1)
             dz=hspe(id,1)
             r=rplanet+z(1)+0.5d0*dz
c             path=dz/dcos(zenithangle)                    !...plane parallel atmosphere
             path=r*(dsqrt(dcos(zenithangle)**2.0d0+       !...spherical atmosphere; see http://en.wikipedia.org/wiki/Air_mass_(astronomy)
     #       (dz/r)**2.0d0+2.0d0*dz/r)-dcos(zenithangle))
             ta=ta+csphotabs(k,l)*n*path
          enddo
          ts=0.0d0                               ! scattering
          do i=1,nspec
             sigma=csrayleigh(i)/(wlth(l)*1.0d-7)**4.0d0
             n=abun(i,1)*density(1)
             dz=hspe(i,1)
             r=rplanet+z(1)+0.5d0*dz
c             path=dz/dcos(zenithangle)                    !...plane parallel atmosphere
             path=r*(dsqrt(dcos(zenithangle)**2.0d0+       !...spherical atmosphere
     #       (dz/r)**2.0d0+2.0d0*dz/r)-dcos(zenithangle))
             ts=ts+sigma*n*path
          enddo
          ftop(l)=fphot(l)*(dexp(-ta-ts)+0.5d0*(1.0d0-dexp(-ts)))
       enddo

c      ...compute opacity due to absorption and scattering for each layer
c                                      - tau(1)    is opacity between top of top layer and z(1)
c                                      - tau(j)    is opacity between z(j-1) and z(j)
c                                      - tau(nz+1) is opacity between z(nz) and bottom of bottom layer
       allocate(tau_abs(nz+1,nwlth),tau_sca(nz+1,nwlth))
       do j=1,nz+1
          r=rplanet+z(min(nz,j))
          if(j.eq.1)then
             dz=0.5d0*(z(j)-z(j+1))
          elseif(j.eq.nz+1)then
             dz=0.5d0*(z(j-2)-z(j-1))
             r=rplanet+z(j-1)-0.5d0*dz
          else
             dz=z(j-1)-z(j)
          endif
c          path=dz/dcos(zenithangle)                       !...plane parallel atmosphere
          path=r*(dsqrt(dcos(zenithangle)**2.0d0+          !...spherical atmosphere
     #    (dz/r)**2.0d0+2.0d0*dz/r)-dcos(zenithangle))
          do l=1,nwlth
             tau_abs(j,l)=0.0d0                  ! absorption
             do k=1,nphotospec
                id=idphotospec(k)
                if(j.eq.1)then
                   n=abun(id,j)*density(j)
                elseif(j.eq.nz+1)then
                   n=abun(id,j-1)*density(j-1)
                else
                   n=0.5d0*
     #             (abun(id,j-1)*density(j-1)+abun(id,j)*density(j))
                endif
                tau_abs(j,l)=tau_abs(j,l)+csphotabs(k,l)*n*path
             enddo
             tau_sca(j,l)=0.0d0                  ! scattering
             do i=1,nspec
                if(j.eq.1)then
                   n=abun(i,j)*density(j)
                elseif(j.eq.nz+1)then
                   n=abun(i,j-1)*density(j-1)
                else
                   n=0.5d0*
     #             (abun(i,j-1)*density(j-1)+abun(i,j)*density(j))
                endif
                sigma=csrayleigh(i)/(wlth(l)*1.0d-7)**4.0d0
                tau_sca(j,l)=tau_sca(j,l)+sigma*n*path
             enddo
          enddo
       enddo

c      ...compute flux at each layer with a two-stream iterative algorithm
       allocate(fdown(0:nz+1,nwlth),o_fdown(0:nz+1,nwlth))
       allocate(fup(0:nz+1,nwlth),o_fup(0:nz+1,nwlth))
       fdown(:,:)=0.0d0
       o_fdown(:,:)=0.0d0
       fup(:,:)=0.0d0
       o_fup(:,:)=0.0d0

c Loop over order of Isaksen algorithm
       do it=1,nmaxit

          var=0.0d0
          do l=1,nwlth
             do j=0,nz+1                         ! evaluate downward flux
                if(j.eq.0)then                   ! ...downward flux at top of top layer
                   fdown(j,l)=ftop(l)
                   cycle
                endif
                ta=tau_abs(j,l)
                ts=tau_sca(j,l)
                if(it.eq.1)then                  ! ...first iteration (zero order estimation)
                   fdown(j,l)=fdown(j-1,l)*dexp(-ta-ts)
                   cycle
                endif
                if(j.eq.1)then
                   down=ftop(l)
                else
                   down=o_fdown(j-1,l)
                endif
                up=o_fup(j,l)
                fdown(j,l)=(0.5d0*(1.0d0-dexp(-ts))*(down+up)+
     #          fdown(j-1,l))*dexp(-ta-ts)
             enddo
             do j=nz+1,0,-1                      ! evaluate upward flux
                if(j.eq.nz+1)then                ! ...upward flux at bottom of bottom layer (planet surface)
                   if(it.eq.1)fup(j,l)=fdown(j,l)*albedoplanet
                   if(it.gt.1)fup(j,l)=o_fdown(j,l)*albedoplanet
                   cycle
                endif
                ta=tau_abs(j+1,l)
                ts=tau_sca(j+1,l)
                if(it.eq.1)then                  ! ...first iteration (zero order estimation)
                   fup(j,l)=fup(j+1,l)*dexp(-ta-ts)
                   cycle
                endif
                if(j.eq.0)then
                   down=ftop(l)
                else
                   down=o_fdown(j,l)
                endif
                up=o_fup(j+1,l)
                fup(j,l)=(0.5d0*(1.0d0-dexp(-ts))*(down+up)+
     #          fup(j+1,l))*dexp(-ta-ts)
             enddo
             do j=0,nz+1                         ! evaluate convergence criterion
                if(it.eq.1)exit
                fx=fdown(j,l)+fup(j,l)
                o_fx=o_fdown(j,l)+o_fup(j,l)
                if(o_fx.eq.0)cycle
                if(dabs(fx-o_fx)/o_fx.gt.var)var=dabs(fx-o_fx)/o_fx
             enddo
             do j=0,nz+1                         ! store old flux values
                o_fdown(j,l)=fdown(j,l)
                o_fup(j,l)=fup(j,l)
             enddo
          enddo
          if(it.gt.1.and.var.lt.varmin)exit

c End of loop order of Isaksen algorithm
       enddo

       if(it.eq.nmaxit)write(*,1000)var

c      ...compute rate of photo reactions at each height
       do j=1,nz
          do k=1,nphotoreac
             jphot(k,j)=0.0d0
             do l=1,nwlth-1
                w1=wlth(l)
                w2=wlth(l+1)
                s1=csphot(k,l)
                s2=csphot(k,l+1)
                f1=fdown(j,l)+fup(j,l)
                f2=fdown(j,l+1)+fup(j,l+1)
                jphot(k,j)=jphot(k,j)+0.5d0*(s1*f1+s2*f2)*(w2-w1)
             enddo
          enddo
       enddo

c      ...locate tau=1 height as a function of wavelength
c         opacity is defined with respect to the non-attenuated incoming stellar flux
c         and with respect to the downward flux at each layer
       do l=1,nwlth
          do j=1,nz
             tau=dlog(fphot(l)/(fdown(j,l)))
             if(tau.gt.1.0d0)exit
          enddo
          if(j.eq.1.or.j.eq.nz)then
             p_tau1(l)=pressure(j)
             z_tau1(l)=z(j)
             cycle
          endif
          tu=dlog(fphot(l)/(fdown(j-1,l)))
          tl=dlog(fphot(l)/(fdown(j,l)))
          pu=pressure(j-1)
          pl=pressure(j)
          slope=(pu-pl)/(tu-tl)
          p_tau1(l)=pl+slope*(1.0d0-tl)
          slope=(z(j-1)-z(j))/(dlog10(pu)-dlog10(pl))
          z_tau1(l)=z(j)+slope*(dlog10(p_tau1(l))-dlog10(pl))
       enddo

       deallocate(tau_abs,tau_sca,fdown,o_fdown,fup,o_fup)

       return

1000   format(1x,
     # ' W- Rayleigh scattering two-stream algorithm not converged',/,
     # '    Reached accuracy = ',1pg10.3)

       end




c_______________________________________________________________________

       subroutine evaluate_therm(id,tkin,h,s,cp)
c      Evaluate H/RT, S/R, and Cp/R for requested species and temperature

       include 'pac.common'
       integer id,k
       real*8 tkin,a(9),h,s,cp

       if(thtype(id).eq.1)then                   ! 7-term NASA polynomials
          k=1
          if(tkin.lt.thtemp(id,3))k=2
          a(1:7)=tha(id,1:7,k)
          h=a(1)+a(2)*tkin/2.0d0+a(3)*tkin**2.0d0/3.0d0+
     #    a(4)*tkin**3.0d0/4.0d0+a(5)*tkin**4.0d0/5.0d0+a(6)/tkin
          s=a(1)*dlog(tkin)+a(2)*tkin+a(3)*tkin**2.0d0/2.0d0+
     #    a(4)*tkin**3.0d0/3.0d0+a(5)*tkin**4.0d0/4.0d0+a(7)
          cp=a(1)+a(2)*tkin+a(3)*tkin**2.0d0+a(4)*tkin**3.0d0+
     #    a(5)*tkin**4.0d0
       elseif(thtype(id).eq.2)then               ! 9-term NASA polynomials
          k=1
          do while(tkin.ge.thtemp(id,k+1).and.k.lt.nthtemp(id))
             k=k+1
          enddo
          a(1:9)=tha(id,1:9,k)
          h=-a(1)*tkin**(-2.0d0)+a(2)*dlog(tkin)/tkin+a(3)+
     #    a(4)*tkin/2.0d0+a(5)*tkin**2.0d0/3.0d0+
     #    a(6)*tkin**3.0d0/4.0d0+a(7)*tkin**4.0d0/5.0d0+a(8)/tkin
          s=-a(1)*tkin**(-2.0d0)/2.0d0-a(2)/tkin+a(3)*dlog(tkin)+
     #    a(4)*tkin+a(5)*tkin**2.0d0/2.0d0+a(6)*tkin**3.0d0/3.0d0+
     #    a(7)*tkin**4.0d0/4.0d0+a(9)
          cp=a(1)*tkin**(-2.0d0)+a(2)/tkin+a(3)+a(4)*tkin+
     #    a(5)*tkin**2.0d0+a(6)*tkin**3.0d0+a(7)*tkin**4.0d0
       endif

       return
       end




c_______________________________________________________________________

       subroutine get_scale_height
c      Evaluate atmospheric scale height at each layer

       include 'pac.common'
       integer i,j
       real*8 g,m,norm

       do j=1,nz
          g=ggravity*mplanet/(rplanet+z(j))**2.0d0
          m=0.0d0
          norm=0.0d0
          do i=1,nspec
             m=m+mass(i)*abun(i,j)
             norm=norm+abun(i,j)
          enddo
          m=m/norm
          hatm(j)=kboltz*tk(j)/(m*g)
       enddo

       return
       end





