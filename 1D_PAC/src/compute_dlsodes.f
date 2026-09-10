
       subroutine compute_dlsodes

       include 'pac.common'
       integer lrw_dlsodes,liw_dlsodes
       parameter(lrw_dlsodes=20+9*nmaxzs+
     #           2*nmaxnzjac+2*nmaxzs+(nmaxnzjac+9*nmaxzs)/2)
       parameter(liw_dlsodes=31+nmaxzs+nmaxnzjac)

       integer i,j,ii,neg,nnz,ihou,imin,isec,itask,istate,
     # iwork_dlsodes(liw_dlsodes),ia(nmaxzs+1),ja(nmaxnzjac),itsteady,
     # ntsteady,k,errspec,erralti,h
       real*8 rwork_dlsodes(lrw_dlsodes),t_total,o_t,tout,o_tout,tmax,
     # y(nmaxzs),atol(nmaxzs),tstep,tmin,wkerr,
     # tin,o_abun(nmaxspec,nmaxz),
     # stdy_abun(nmaxspec,nmaxz),stdy_o_abun(nmaxspec,nmaxz),
     # tloop,ini_tout,stdy_tout,stdy_time,abunminvary,abun_min
       external f_dlsodes,jac_dlsodes

c      ...Initialize variables
      
       t_total=0.0d0                               ! initial time of simulation @Arturo
       ini_tout=1.0d-10                          !! first output time, should it be that small? 1.0d-10
       tout=ini_tout                           
       tmin=day                                  ! minimum time to integrate
       tmax=13.75d9*year                         ! maximum time: age of the Universe
       tstep=10.0d0**0.2d0                       ! time step factor after tmin
       tiny_value=1e-80                         ! threshold to consider values = 0
c       rtol=1.0d-4                              ! relative tolerance for DLSODES (previously -8 RT=1.0d-3, -5 for @Arturo)
       abunminsteady=1.0d-20                  !! lower abundance cutoff to evaluate steady state (RT=1.0d-20) @Arturo
       abunminvary=1.0d-30
       writekrate=.false.
c       maxerrsteady=3.0*rtol                    ! relative tolerance to consider steady state @Arturo
       ntsteady=2                                ! minimum number of consecutive steps with dabun<maxerrsteady to consider steady state
       itsteady=0                                ! number of consecutive steps with dabun<maxerrsteady
       itime=0                                   ! counter of integration times to write out results
       if(sensi.gt.0)irate=0                     ! initializing irate to be 0 so it does noy vary at first @Arturo
       rvirate=0
       if(endsensi.eq.0)endsensi=nreac           ! setting end in nreac as maximum of indexes variation for sensitivity @Arturo
       itask=1                                   ! itask for DLSODES
       first_vary=.false.
       errflag=0
       nout=0
       abun_min=1.0d0
       call conservation(0)
       write(*, 1001)rtol
       write(*,*)
      


c Big loop of time (use DLSODES to integrate the system of ODES)
       do while(t_total.lt.tmax.and..not.steady)

c      ...store times and abundances as "old" values, assign abundances to array y
          o_t=t_total
          o_tout=tout
          do j=1,nz
             do i=1,nspec
                ii=(j-1)*nspec+i
                y(ii)=abun(i,j)
                o_abun(i,j)=abun(i,j)
             enddo
          enddo
          do i=1,nzs
             atol(i)=max(1.0d-78,1.0d-15*y(i))
c             atol(i)=1.0d-99
          enddo

c      ...get rate constants, photorates, and diffusion terms, why should they update them? @Arturo 
          call get_rate_constants
          call get_diffusion_terms
          call get_photorates
          change_rate=.false.
          first_vary=.false.
          

c      ...reset current initial time to zero, correcting bad naming @Arturo
c          tin=t
c          t=0.0d0
c          tout=tout-tin
          tloop=0.0d0
          tout=tout-t_total
          

c      ...small loop of time (normally just one step if DLSODES integrates up to tout)
          do while (tloop.lt.tout) ! We need this condition in case DLSODES do not manage to propagate until tout @Arturo
             istate=1                            ! call DLSODES
             iwork_dlsodes(:)=0
             rwork_dlsodes(:)=0.0d0
             iwork_dlsodes(6)=500               ! MXTSEP to DLSODES @Arturo
             call compute_iaja(nmaxzs+1,nmaxnzjac,ia,ja,nnz)
             iwork_dlsodes(30+1:30+nzs+1)=ia(1:nzs+1)
             iwork_dlsodes(31+nzs+1:31+nzs+nnz)=ja(1:nnz)
             call DLSODES(f_dlsodes,nzs,y,tloop,tout,2,rtol,atol,itask,
     #       istate,1,rwork_dlsodes,lrw_dlsodes,iwork_dlsodes,
     #       liw_dlsodes,jac_dlsodes,21) ! DLSODES does not know the history @Arturo
c             write(*,*)'t_total=',t_total,'tloop=',tloop,'tout=',tout
             neg=0
             do j=1,nz                           ! update abundances
                do i=1,nspec
                   ii=(j-1)*nspec+i
                   abun(i,j)=y(ii)
                   if(abun(i,j).lt.-atol(i))then !! That is why ATOL is not zero?
                      write(*,1030)j,z(j)/1.0d5,spec(i),abun(i,j)
                      neg=neg+1
                      abun(i,j)=abunmin
                   endif
                   if(abun(i,j).lt.abun_min)abun_min=abun(i,j)
                   if (abun(i,j).lt.tiny_value)then
                     abun(i, j)=tiny_value
                     write(*,*)'W-Zero abundance detected:',i,j
                   endif
                enddo
             enddo
c             if(t_total.gt.tmin)write(*,*)'abun last err',
c     #       abun(errspec,erralti)
c             if(neg.gt.0)stop 'E- Negative abundance found'
c      ...end of small loop of time
          enddo
         write(*,*)'abun_min = ',abun_min
c      ...return to real times
c          t=tin+t
c          tout=tin+tout
c          write(*,*)'tloop=',tloop,'tout=',tout
          tout=t_total+tout
          t_total=t_total+tloop

c      ...verify conservation of elements
          call conservation(2)

c      ...verify whether steady state has been reached
c          open(unit=70,file='outhack.out',action='write',
c     #    position='append')
c          write(70,1060)namemodel
          if(t_total.ge.tmin)then
             errsteady=tiny
             do j=1,nz
                do i=1,nspec
                   if(abun(i,j).lt.abunminsteady)cycle ! Discard abundances smaller than abunminsteady
                   wkerr=dabs((abun(i,j)-o_abun(i,j))/o_abun(i,j))
                   if(wkerr.gt.errsteady)then
                     errsteady=wkerr
                     errspec=i
                     erralti=j
                   endif
                enddo
             enddo
c             write(70,1050)errspec,spec(errspec),erralti,
c     #       abun(errspec,erralti)
c             write(*,*)t_total,spec(errspec),erralti,
c     #       abun(errspec,erralti)
             if(errsteady.gt.maxerrsteady)itsteady=0
             if(errsteady.le.maxerrsteady)itsteady=itsteady+1
             if(itsteady.ge.ntsteady)steady=.true.
             itime=itime+1
             time=t_total
             call write_results

             if(sensi.gt.0.and.steady)then
              if(irate.eq.0)then
               irate=stsensi
               stdy_abun=abun
               stdy_o_abun=o_abun
               stdy_time=t_total
               stdy_tout=tout
               first_vary=.true.
               abunminsteady=abunminvary
               if(sensi.eq.1)irate=stsensi-1
              endif
            

              if(irate.le.endsensi)then
               call write_results
               irate=irate+1
               if(.not.idrvreac(irate).eq.0)irate=irate+1
               abun=stdy_abun
               o_abun=stdy_o_abun
               t_total=0.0d0 !! Resetting time @Arturo
               tout=ini_tout 
               steady=.false.
               ntsteady=3 !! Resetting numbers of steps to call it convergence @Arturo
               change_rate=.true.
               itsteady=0
               itime=0
              endif
             endif
          else
             errsteady=0.0d0
          endif
          close(70)

c      ...write on screen current time
          call get_cputime(ihou,imin,isec)
          write(*,1000)t_total,ihou,imin,isec,
     #    errsteady

c      ...get next tout "tout=tout*tstep"
          tout=tout*tstep                      ! if t > tmin
          if(t_total.gt.tmin)nout=nout+1                          ! Cut the loop if too many files @Arturo
          if (max_nout.gt.0.and.nout.gt.max_nout)then
             steady=.true.
          endif
          if (istate.eq.-1)then ! Cut the loop if it starts looping
             steady=.true.
          endif
          if(t_total.lt.tmin)tout=o_tout*ftstep        ! if t < tmin
          if(t_total.lt.tmin.and.tout.gt.tmin)tout=tmin
          if(tout.gt.tmax)tout=tmax

c End of big loop of time
       enddo

       if(steady)then
         if(nout.gt.max_nout.and.max_nout.gt.0)then
           write(*,1313)t_total
         else
           write(*,1010)t_total
         endif
       else
          write(*,1020)t_total
       endif

       return

1000   format(2x,'time[s]=',1pg10.3,3x,
     #           'CPU=',i3,'h 'i2,'m ',i2,'s',3x,
     #           'err=',1pg10.3)
1001   format(2x,'I- RTOL for DLSODES is set to' 1pg10.2)
1010   format(1x,' I- Steady state reached at t=',1pg8.1,'s')
1313   format(1x,' I- Max number of outputs reached at t=',1pg8.1,'s')
1020   format(1x,' I- Integration done up to t=',1pg8.1,'s')
1030   format(1x,' E- In layer ',i3,' z=',1pg10.2,' km ',/,
     #        1x,'    molar fraction of ',a12,' is ',1pg10.2)
1040   format(1x, 1pg10.2, 1pg10.2, 1pg10.2)

1050   format(1x,'i,j:',i4,'(',a,')',i4,' abun=',1pg10.2)
1060   format(1x,a)

       end




c_______________________________________________________________________

       subroutine f_dlsodes(neq,t,y,ydot)
c      user subroutine for DLSODES to be used in evolution calculation
c      it evaluates the derivatives of the abundances

       include 'pac.common'
       integer i,j,ii,neq
       double precision t,y(*),ydot(*)

       i=int(t)                        ! to avoid warning
       if(neq.ne.nzs) stop ' E- Dimension error in f_dlsodes'
       do j=1,nz
          do i=1,nspec
             ii=(j-1)*nspec+i
             abun(i,j)=y(ii)
          enddo
       enddo

       call evaluate_derivatives

       do j=1,nz
          do i=1,nspec
             ii=(j-1)*nspec+i
             ydot(ii)=f(i,j)/density(j)
          enddo
       enddo

       return
       end




c_______________________________________________________________________

       subroutine jac_dlsodes(neq,t,y,jth,ian,jan,pdj)
c      user subroutine for DLSODES to be used in evolution calculation
c      it evaluates the jacobian of the system of ODEs

       include 'pac.common'
       integer i,j,l,g,gwk,id,ii,jth,neq,iwk,jwk
       double precision t,y(*),ian(*),jan(*),pdj(*)
       real*8 rwk,iswk

       i=int(t)                        ! to avoid warning
       i=int(ian(1))                   ! to avoid warning
       i=int(jan(1))                   ! to avoid warning
       if(neq.ne.nzs) stop ' E- Dimension error in jac_dlsodes'
       do j=1,nz
          do i=1,nspec
             ii=(j-1)*nspec+i
             abun(i,j)=y(ii)
          enddo
       enddo

       jwk=(jth-1)/nspec+1
       iwk=jth-(jwk-1)*nspec
       do j=1,nz
          do i=1,nspec
             ii=(j-1)*nspec+i
             if(jwk.eq.j)then

c      ...evaluate jacobian of formation/destruction: d[dn_i^j/dt]/dy_iwk^j = d[dn_i^j/dt]/dn_iwk^j * n^j
                df(i,j,iwk)=0.0d0
                do l=1,nterm(i)
                   id=idterm(i,l)
                   iswk=0.0d0
                   gwk=0
                   do g=1,nreag(id)
                      if(idreag(id,g).eq.iwk)then
                         iswk=iswk+1.0d0
                         gwk=g
                      endif
                   enddo
                   if(iswk.gt.0.0d0)then
                      rwk=1.0d0
                      do g=1,nreag(id)
                         if(g.ne.gwk)
     #                   rwk=rwk*abun(idreag(id,g),j)*
     #                   density(j)
                      enddo
                      df(i,j,iwk)=df(i,j,iwk)+density(j)*
     #                sterm(i,l)*fterm(i,l)*krate(id,j)*rwk*iswk
                   endif
                enddo
c     ...add diffusion contribution to jacobian: d[dn_i^j/dt]/dy_iwk^j
                if(diffusion)then
                   if(i.eq.iwk)then
                      df(i,j,iwk)=df(i,j,iwk)+df_dif(i,j)
                   endif
                endif
c     ...add photochemical contribution to jacobian: d[dn_i^j/dt]/dy_iwk^j
                if(photochem)then
                   do l=1,nphototerm(i)
                      id=idphototerm(i,l)
                      if(idphotoreag(id).eq.iwk)then
                         df(i,j,iwk)=df(i,j,iwk)+density(j)*
     #                   sphototerm(i,l)*fphototerm(i,l)*jphot(id,j)
                      endif
                   enddo
                endif
c     ...assign to pdj
                pdj(ii)=df(i,j,iwk)/density(j)
             elseif(jwk.eq.(j-1).and.iwk.eq.i)then
                pdj(ii)=df_jm_dif(i,j)/density(j)
             elseif(jwk.eq.(j+1).and.iwk.eq.i)then
                pdj(ii)=df_jp_dif(i,j)/density(j)
             else
                pdj(ii)=0.0d0
             endif
          enddo
       enddo

       return
       end




c_______________________________________________________________________

       subroutine compute_iaja(nia,nja,ia,ja,nnz)
c      ...evaluate ia,ja for DLSODES

       include 'pac.common'
       integer nia,nja,ia(nia),ja(nja),nnz,i,j,icol,jcol,irow,jrow

       call evaluate_jacobian
       nnz=1
       do jcol=1,nz
          do icol=1,nspec
             j=(jcol-1)*nspec+icol
             ia(j)=nnz
             do jrow=1,nz
                do irow=1,nspec
                   i=(jrow-1)*nspec+irow
                   if(jcol.eq.jrow)then
                      if(df(irow,jrow,icol).ne.0.0d0)then
                         ja(nnz)=i
                         nnz=nnz+1
                      endif
                   endif
                   if(.not.diffusion)cycle
                   if(jcol.eq.(jrow-1).and.icol.eq.irow)then
                      if(df_jm_dif(irow,jrow).ne.0.0d0)then
                         ja(nnz)=i
                         nnz=nnz+1
                      endif
                   elseif(jcol.eq.(jrow+1).and.icol.eq.irow)then
                      if(df_jp_dif(irow,jrow).ne.0.0d0)then
                         ja(nnz)=i
                         nnz=nnz+1
                      endif
                   endif
                enddo
             enddo
          enddo
       enddo
       ia(nzs+1)=nnz

       return
       end
