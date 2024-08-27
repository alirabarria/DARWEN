
       subroutine compute_dlsodes

       include 'pac.common'
       integer lrw_dlsodes,liw_dlsodes
       parameter(lrw_dlsodes=20+9*nmaxzs+
     #           2*nmaxnzjac+2*nmaxzs+(nmaxnzjac+9*nmaxzs)/2)
       parameter(liw_dlsodes=31+nmaxzs+nmaxnzjac)

       integer i,j,ii,neg,nnz,ihou,imin,isec,itask,istate,
     # iwork_dlsodes(liw_dlsodes),ia(nmaxzs+1),ja(nmaxnzjac),itsteady,
     # ntsteady,k,errspec,erralti,iszero
       real*8 rwork_dlsodes(lrw_dlsodes),t,o_t,tout,o_tout,tmax,
     # y(nmaxzs),rtol,atol(nmaxzs),tstep,tmin,wkerr,errsensi,
     # maxerrsteady,abunminsteady,tin,o_abun(nmaxspec,nmaxz),
     # stdy_abun(nmaxspec,nmaxz), stdy_o_abun(nmaxspec,nmaxz),
     # stdy_time
       external f_dlsodes,jac_dlsodes

c      ...Initialize variables
       t=0.0d0                                   ! initial time
       tout=1.0d-10                              ! first output time
       tmin=day                                  ! minimum time to integrate
       tmax=13.75d9*year                         ! maximum time: age of the Universe
       tstep=10.0d0**0.2d0                       ! time step factor after tmin
       abunminsteady=1.0d-20                     ! lower abundance cutoff to evaluate steady state
       maxerrsteady=1.0d-4                       ! relative tolerance to consider steady state
       ntsteady=2                                ! minimum number of consecutive steps with dabun<maxerrsteady to consider steady state
       itsteady=0                                ! number of consecutive steps with dabun<maxerrsteady
       itime=0                                   ! counter of integration times to write out results
       errsensi=1.0e-4                           ! error when to start changing the rate constant to do sensitivity analysis @Arturo
       if(sensi.eq.2)irate=0                     ! initializing irate to be 0 so it does noy vary at first @Arturo
       rvirate=0
       if(endsensi.eq.0)endsensi=nreac           ! setting end in nreac as maximum of indexes variation for sensitivity @Arturo
c       rtol=1.0d-8                               ! relative tolerance for DLSODES
       rtol=1.0d-3                               ! relative tolerance for DLSODES
       itask=1                                   ! itask for DLSODES
       call conservation(0)
       write(*,*)

c Big loop of time (use DLSODES to integrate the system of ODES)
       do while(t.lt.tmax.and..not.steady)

c      ...store times and abundances as "old" values, assign abundances to array y
          o_t=t
          o_tout=tout
          do j=1,nz
             do i=1,nspec
                ii=(j-1)*nspec+i
                y(ii)=abun(i,j)
                o_abun(i,j)=abun(i,j)
             enddo
          enddo
          do i=1,nzs
c             atol(i)=max(1.0d-20,1.0d-15*y(i))
             atol(i)=1.0d-99
          enddo

c      ...get rate constants, photorates, and diffusion terms
          call get_rate_constants
          call get_diffusion_terms
          call get_photorates

c      ...reset current initial time to zero
          tin=t
          t=0.0d0
          tout=tout-tin

c      ...small loop of time (normally just one step if DLSODES integrates up to tout)
          do while (t.lt.tout) ! Why do we need this condition, is this not handled at DLSODES itself? @Arturo
             istate=1                            ! call DLSODES
             iwork_dlsodes(:)=0
             rwork_dlsodes(:)=0.0d0
             iwork_dlsodes(6)=1000
c      ...apparently, DLSODES does not know which epoch is integrating @Arturo, only jac & tout changes?        
             call compute_iaja(nmaxzs+1,nmaxnzjac,ia,ja,nnz)
             iwork_dlsodes(30+1:30+nzs+1)=ia(1:nzs+1)
             iwork_dlsodes(31+nzs+1:31+nzs+nnz)=ja(1:nnz)
             call DLSODES(f_dlsodes,nzs,y,t,tout,2,rtol,atol,itask,
     #       istate,1,rwork_dlsodes,lrw_dlsodes,iwork_dlsodes,
     #       liw_dlsodes,jac_dlsodes,21)
             write(*,*)'inside small loop',t,tout,y(32)
             neg=0
             do j=1,nz                           ! update abundances
                do i=1,nspec
                   ii=(j-1)*nspec+i
                   abun(i,j)=y(ii)
                   if(abun(i,j).lt.-atol(i))then
                      write(*,1030)j,z(j)/1.0d5,spec(i),abun(i,j)
                      neg=neg+1
                      abun(i,j)=abunmin
                   endif
                enddo
             enddo
c             if(neg.gt.0)stop 'E- Negative abundance found'
c      ...end of small loop of time
          enddo

c      ...return to real times
          t=tin+t
          tout=tin+tout

c      ...verify conservation of elements
          call conservation(2)

c      ...verify whether steady state has been reached
          open(unit=70,file='outhack.out',action='write',
     #    position='append')
          write(70,1060)namemodel
          if(t.ge.tmin)then
             errsteady=tiny
             do j=1,nz
                do i=1,nspec
                   if(abun(i,j).lt.abunminsteady)cycle
                   wkerr=dabs((abun(i,j)-o_abun(i,j))/o_abun(i,j))
                   if(wkerr.gt.errsteady)then
                     errsteady=wkerr
                     errspec=i
                     erralti=j
                   endif
                enddo
             enddo
             write(70,1050)errspec,spec(errspec),erralti,
     #       abun(errspec,erralti)
             write(*,*)t,spec(errspec),erralti,abun(errspec,erralti)
             if(errsteady.gt.maxerrsteady)itsteady=0
             if(errsteady.le.maxerrsteady)itsteady=itsteady+1
             if(itsteady.ge.ntsteady)steady=.true.
             itime=itime+1
             time=t
             call write_results

             if(sensi.eq.2.and.errsteady.le.errsensi)then
              if(irate.eq.0)then
               irate=stsensi
               stdy_abun=abun
               stdy_o_abun=o_abun
              endif
              if(irate.le.endsensi)then
               call write_results
               irate=irate+1
               if(.not.idrvreac(irate).eq.0)irate=irate+1
               abun=stdy_abun
               o_abun=stdy_o_abun
               t=stdy_time
               steady=.false.
              endif
             endif
          else
             errsteady=0.0d0
          endif
          close(70)

c      ...write on screen current time
          call get_cputime(ihou,imin,isec)
          iszero=0
          if(errsteady.eq.0)iszero=1
          write(*,1000)t,ihou,imin,isec,errsteady,
     #    iszero

c      ...get next tout "tout=tout*tstep"
          tout=tout*tstep                      ! if t > tmin
          if(t.lt.tmin)tout=o_tout*ftstep        ! if t < tmin
          if(t.lt.tmin.and.tout.gt.tmin)tout=tmin
          if(tout.gt.tmax)tout=tmax

c End of big loop of time
       enddo

       if(steady)then
          write(*,1010)t
       else
          write(*,1020)t
       endif

       return

1000   format(2x,'time[s]=',1pg10.3,3x,
     #           'CPU=',i3,'h 'i2,'m ',i2,'s',3x,
     #           'err=',1pg10.3, 'err=0? ',i1) ! @Arturo asked if err.eq.0
1010   format(1x,' I- Steady state reached at t=',1pg8.1,'s')
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
