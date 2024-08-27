
       subroutine compute_newtonraphson

       include 'pac.common'
       integer i,j,k,ktop,niter,nmaxiter,nminiter,ihou,imin,isec,
     # ibad,jbad,kbad
       real*8 abwk,minabconv,p(nmaxspec,nmaxz),err,o_err,psum,
     # slow,o_abun(nmaxspec,nmaxz),abunmax(nmaxspec),wkerr,errbad
       logical converge,specfnd(nmaxspec)

c      ...Initialize variables
       nmaxiter=2000
       nminiter=10
       niter=1
       converge=.false.
c       minabconv=1.0d-20
       minabconv=1.0d-15
       slow=0.7d0                                ! slow convergence factor
       rtol=1.0d-2
       call conservation_layer(0)
       write(*,*)

c      ...verify that there is not initial elemental abundance vertical gradient
       call conservation(1)                      ! locate most abundant element in atmosphere
       ktop=0
       abwk=0.0d0
       do k=1,nelem
          if(cdelem(k).gt.abwk)then
             abwk=cdelem(k)
             ktop=k
          endif
       enddo
       if(ktop.eq.0)stop ' E- Most abundant element not found'
       do k=1,nelem                              ! compare among layers
          do j=2,nz                              ! ...element abundance relative to most abundant element
             err=dabs((elab0(k,j)/  elab0(ktop,j)-
     #                 elab0(k,j-1)/elab0(ktop,j-1))/
     #                (elab0(k,j-1)/elab0(ktop,j-1)))
             if(err.gt.5.0d-3)goto 100
          enddo
       enddo

c      ...estimate maximum abundance reachable by each species
       do i=1,nspec
          abunmax(i)=0.0d0
          do k=1,nelem
             if(nat(i,k).eq.0)cycle
             do j=1,nz
                if(elab0(k,j)/dble(nat(i,k)).gt.abunmax(i))
     #          abunmax(i)=elab0(k,j)/dble(nat(i,k))
             enddo
          enddo
       enddo
       if(ion.and.id_electr.gt.0)abunmax(id_electr)=1.0d0


c Big loop to converge towards steady state
       do while (niter.le.nmaxiter.and..not.converge)

c          if(niter.eq.1.or.mod(niter,5).eq.0)then ! at first iteration and each 5

c      ...compute mean values to scale dy/dt and conservation equations
             call get_rate_constants
             call get_diffusion_terms
             call get_photorates
             call evaluate_derivatives
             dydtscale=0.0d0
             do j=1,nz
                do i=1,nspec
                   dydtscale=dydtscale+f(i,j)/density(j)
                enddo
             enddo
             dydtscale=dydtscale/dble(nzs)
             call conservation_layer(1)
             elabscale=0.0d0
             do k=1,nelem
                do j=1,nz
                   elabscale=elabscale+elab(k,j)
                enddo
             enddo
             elabscale=elabscale/dble(nelem*nz)

c      ...find most abundant species locking each element in each layer for conservation equations
             do j=1,nz
                specfnd(:)=.false.
                do k=1,nelem
                   abwk=0.0d0
                   do i=1,nspec
                      if(specfnd(i))cycle
                      if(nat(i,k).eq.0)cycle
                      if(nat(i,k)*abun(i,j).gt.abwk)then
                         abwk=nat(i,k)*abun(i,j)
                         idcons(k,j)=i
                      endif
                   enddo
                   specfnd(idcons(k,j))=.true.
                enddo
             enddo
c          endif

c      ...save abundances
          do j=1,nz
             do i=1,nspec
                o_abun(i,j)=abun(i,j)
             enddo
          enddo

c      ...call Newton-Raphson subroutine, 1 iteration and back
          call newtonraphson(p)
c          call newtonraphson2(p)

c      ...recover abundance corrections
          o_err=err
          err=0.0d0
          psum=0.0d0
          errbad=0.0d0
          ibad=0
          jbad=0
          do j=1,nz
             do i=1,nspec
                psum=psum+dabs(p(i,j))
                if((abun(i,j)+p(i,j)).lt.minabconv)cycle   ! only molar fractions > threshold
                wkerr=dabs(p(i,j)/(o_abun(i,j)+tiny))      ! ...count for convergence
                err=err+wkerr
                if(wkerr.gt.errbad)then
                   errbad=wkerr
                   ibad=i
                   jbad=j
                endif
             enddo
          enddo

c      ...have we converged after a minimum number of iterations?
          if(err.lt.rtol.and.o_err.lt.rtol.and.niter.gt.nminiter)then
             do j=1,nz
                do i=1,nspec
                   abun(i,j)=abun(i,j)+p(i,j)
                enddo
             enddo
             converge=.true.
             errsteady=err
c      ...if not, continue iterating after correcting abundances (robusteness corrections)
          else
             do j=1,nz
                do i=1,nspec
                   if(psum.gt.1.0d0)p(i,j)=p(i,j)/psum     ! scale residuals if too big
                   abun(i,j)=o_abun(i,j)+slow*p(i,j)       ! get new molar fractions
                   if(abun(i,j).gt.abunmax(i))             ! correct if too large molar fraction
     #                abun(i,j)=0.5d0*(o_abun(i,j)+abunmax(i))
                   if(abun(i,j).lt.0.0d0)                  ! correct if negative molar fraction
     #                abun(i,j)=0.1d0*o_abun(i,j)
                   if(mod(niter,2).eq.0)                   ! take the mean each 2 iterations
     #                abun(i,j)=0.5d0*(o_abun(i,j)+abun(i,j))
                enddo
             enddo
          endif

c      ...write on screen current time
          call get_cputime(ihou,imin,isec)
          write(*,1000)niter,nmaxiter,ihou,imin,isec,err,
     #    spec(ibad),jbad,pressure(jbad)/1.0d6

          niter=niter+1
       enddo
c End of big loop to converge towards steady state


       if(converge)then
          write(*,1010)niter-1,err
       else
          write(*,1020)niter-1,err
          stop
       endif

c      ...write on screen maximum element conservation error
       err=0.0d0
       kbad=0
       jbad=0
       do j=1,nz
          do k=1,nelem
             abwk=0.0d0
             do i=1,nspec
                abwk=abwk+nat(i,k)*abun(i,j)
             enddo
             abwk=(abwk-elab0(k,j))/elab0(k,j)
             if(dabs(abwk).gt.dabs(err))then
                err=abwk
                kbad=k
                jbad=j
             endif
          enddo
       enddo
       write(*,1030)err,elem(kbad),jbad,pressure(jbad)/1.0d6

       return

100    write(*,*)' E- Initial vertical gradient for: ',elem(k)
       write(*,*)'    molar fraction sum in layer ',j,elab(k,j)
       write(*,*)'    molar fraction sum in layer ',j-1,elab(k,j-1)
       stop

1000   format(2x,'it=',i4,'/',i4,3x,
     #           'CPU=',i3,'h 'i2,'m ',i2,'s',3x,
     #           'err=',1pg10.3,3x,
     #           'bad=',a6,1x,i3,1x,1pg8.1,'bar')
1010   format(/,2x,'I- Converged after ',i5,' iterations',/,
     #          2x,'   relative abundance error = ',1pg10.3)
1020   format(/,2x,'E- Not converged after ',i4,' iterations',/,
     #          2x,'   relative abundance error = ',1pg10.3)
1030   format(/,2x,'I- Maximum element conservation relative error =',
     #        1pg12.4,/,2x,' element =',a2,
     #        2x,'layer ',i3,1x,1pg8.1,'bar')

       end




c_______________________________________________________________________

       subroutine newtonraphson(p)
c      subroutine equivalent to "mnewt" from Numerical Recipes (just one iteration)

       include 'pac.common'
       integer i,j,k
       real*8 abwk,p(nmaxspec,nmaxz),ajac(nmaxspec,nmaxspec,nmaxz),
     # bjac(nmaxspec,nmaxspec,nmaxz-1),
     # cjac(nmaxspec,nmaxspec,nmaxz-1),yfwd(nmaxspec,nmaxz)

c Get derivatives and jacobian

c      ...get dy/dt terms, and jacobian
       call evaluate_derivatives
       call evaluate_jacobian

c     ...scale and assign values to P (right-hand terms) and J (jacobian)
c        J * DY = P ;  J * DY = - f_i^j / n^j
       do j=1,nz                                 ! assign values of righ-hand terms P
          do i=1,nspec
             p(i,j)=-f(i,j)/density(j)/dydtscale
          enddo
       enddo
       ajac(:,:,:)=0.0d0                         ! assign set of matrix A in J main diagonal
       do j=1,nz
          do i=1,nspec
             do k=1,nspec
                ajac(i,k,j)=df(i,j,k)/density(j)/dydtscale
             enddo
          enddo
       enddo
       bjac(:,:,:)=0.0d0                         ! assign set of matrix B in J subdiagonal
       do j=2,nz
          do i=1,nspec
             bjac(i,i,j-1)=df_jm_dif(i,j)/density(j)/dydtscale
          enddo
       enddo
       cjac(:,:,:)=0.0d0                         ! assign set of matrix C in J superdiagonal
       do j=1,nz-1
          do i=1,nspec
             cjac(i,i,j)=df_jp_dif(i,j)/density(j)/dydtscale
          enddo
       enddo

c      ...substitute by the corresponding conservation equations
       do j=1,nz
          do k=1,nelem
             abwk=0.0d0
             do i=1,nspec
                abwk=abwk+nat(i,k)*abun(i,j)
                ajac(idcons(k,j),i,j)=nat(i,k)/elabscale
             enddo
             p(idcons(k,j),j)=-(abwk-elab0(k,j))/elabscale
          enddo
          if(ion)then
             abwk=0.0d0
             do i=1,nspec
                abwk=abwk+charge(i)*abun(i,j)
                ajac(id_electr,i,j)=charge(i)/elabscale
             enddo
             p(id_electr,j)=-abwk/elabscale
          endif
       enddo

c      ...LU decompose full (tridiagonal by blocks) matrix
       call lubtridag(nspec,nmaxspec,nz,nmaxz,ajac,bjac,cjac)

c      ...solve system of equations (get Newton step p_i^j = dy_i^j)
       do j=1,nz                                 ! forward substitution
          do i=1,nspec
             yfwd(i,j)=p(i,j)
             if(j.gt.1)then
                do k=1,nspec
                   yfwd(i,j)=yfwd(i,j)-bjac(i,k,j-1)*yfwd(k,j-1)
                enddo
             endif
             do k=1,i-1
                yfwd(i,j)=yfwd(i,j)-ajac(i,k,j)*yfwd(k,j)
             enddo
          enddo
       enddo
       do j=nz,1,-1                              ! back substitution
          do i=nspec,1,-1
             p(i,j)=yfwd(i,j)
             do k=i+1,nspec
                p(i,j)=p(i,j)-ajac(i,k,j)*p(k,j)
             enddo
             if(j.lt.nz)then
                do k=1,nspec
                   p(i,j)=p(i,j)-cjac(i,k,j)*p(k,j+1)
                enddo
             endif
             p(i,j)=p(i,j)/ajac(i,i,j)
          enddo
       enddo

       return
       end




c_______________________________________________________________________

       subroutine newtonraphson2(p)
c      subroutine equivalent to "mnewt" from Numerical Recipes (just one iteration)

       include 'pac.common'
       integer i,j,k,ii,jj,irow,icol,id,imax,ll
       real*8 p(nmaxzs),abwk,aamax,sum,dum,sdp
       integer, allocatable, dimension(:)   :: indx
       real*8,  allocatable, dimension(:)   :: fvec,fvec_save,vv,r
       real*8,  allocatable, dimension(:,:) :: fjac,fjac_save

c      ...assign dimension of variables
       if(nzs.gt.10000)stop'E- Too large size in newtonraphson'
       allocate(fvec(nzs))
       allocate(fvec_save(nzs))
       allocate(fjac(nzs,nzs))
       allocate(fjac_save(nzs,nzs))
       allocate(vv(nzs))
       allocate(r(nzs))
       allocate(indx(nzs))

c Get derivatives and jacobian

c      ...get dy/dt terms, and jacobian
       call evaluate_derivatives
       call evaluate_jacobian

c     ...scale and assign values to dy/dt (fvec) and jacobian (fjac)
       do j=1,nz
          do i=1,nspec
             irow=(j-1)*nspec+i
             fvec(irow)=f(i,j)/density(j)/dydtscale
             do jj=1,nz
                do ii=1,nspec
                   icol=(jj-1)*nspec+ii
                   if(jj.eq.j)then
                      fjac(irow,icol)=df(i,j,ii)/density(j)/dydtscale
                   elseif(jj.eq.j-1.and.ii.eq.i)then
                      fjac(irow,icol)=
     #                df_jm_dif(i,j)/density(j)/dydtscale
                   elseif(jj.eq.j+1.and.ii.eq.i)then
                      fjac(irow,icol)=
     #                df_jp_dif(i,j)/density(j)/dydtscale
                   else
                      fjac(irow,icol)=0.0d0
                   endif
                enddo
             enddo
          enddo
       enddo

c      ...substitute by the corresponding conservation equations
       do j=1,nz
          do k=1,nelem
             id=(j-1)*nspec+idcons(k,j)
             abwk=0.0d0
             do i=1,nspec
                ii=(j-1)*nspec+i
                abwk=abwk+nat(i,k)*abun(i,j)
                fjac(id,ii)=nat(i,k)/elabscale
             enddo
             fvec(id)=(abwk-elab0(k,j))/elabscale
          enddo
c no ions
       enddo

c Save fvec and fjac, and assign values to p (right-hand side of linear equations)
        do i=1,nzs
          fvec_save(i)=-fvec(i)
          do j=1,nzs
             fjac_save(i,j)=fjac(i,j)
          enddo
          p(i)=-fvec(i)
       enddo

c ludcmp *** LU decompose jacobian matrix

       do i=1,nzs
          aamax=0.0d0
          do j=1,nzs
             if(dabs(fjac(i,j)).gt.aamax)aamax=dabs(fjac(i,j))
          enddo
          if(aamax.eq.0.0d0)stop 'singular matrix in newtonraphson'
          vv(i)=1.0d0/aamax
       enddo
       do j=1,nzs
          do i=1,j-1
             sum=fjac(i,j)
             do k=1,i-1
                sum=sum-fjac(i,k)*fjac(k,j)
             enddo
             fjac(i,j)=sum
          enddo
          aamax=0.0d0
          imax=0
          do i=j,nzs
             sum=fjac(i,j)
             do k=1,j-1
                sum=sum-fjac(i,k)*fjac(k,j)
             enddo
             fjac(i,j)=sum
             dum=vv(i)*dabs(sum)
             if(dum.ge.aamax)then
                imax=i
                aamax=dum
             endif
          enddo
          if(j.ne.imax)then
             do k=1,nzs
                dum=fjac(imax,k)
                fjac(imax,k)=fjac(j,k)
                fjac(j,k)=dum
             enddo
             vv(imax)=vv(j)
          endif
          indx(j)=imax
          if(fjac(j,j).eq.0.0d0)fjac(j,j)=tiny
          if(j.ne.nzs)then
             dum=1.0d0/fjac(j,j)
             do i=j+1,nzs
                fjac(i,j)=fjac(i,j)*dum
             enddo
          endif
       enddo

c lubksb *** forward and back substitution
       ii=0
       do i=1,nzs
          ll=indx(i)
          sum=p(ll)
          p(ll)=p(i)
          if(ii.ne.0)then
             do j=ii,i-1
                sum=sum-fjac(i,j)*p(j)
             enddo
          elseif(sum.ne.0.0d0)then
             ii=i
          endif
          p(i)=sum
       enddo
       do i=nzs,1,-1
          sum=p(i)
          do j=i+1,nzs
             sum=sum-fjac(i,j)*p(j)
          enddo
          p(i)=sum/fjac(i,i)
       enddo

c mprove *** improve solution
       do i=1,nzs
          sdp=-fvec_save(i)
          do j=1,nzs
             sdp=sdp+fjac_save(i,j)*p(j)
          enddo
          r(i)=sdp
       enddo
       ii=0                                      ! call lubksb to solve for the error term
       do i=1,nzs
          ll=indx(i)
          sum=r(ll)
          r(ll)=r(i)
          if(ii.ne.0)then
             do j=ii,i-1
                sum=sum-fjac(i,j)*r(j)
             enddo
          elseif(sum.ne.0.0d0)then
             ii=i
          endif
          r(i)=sum
       enddo
       do i=nzs,1,-1
          sum=r(i)
          do j=i+1,nzs
             sum=sum-fjac(i,j)*r(j)
          enddo
          r(i)=sum/fjac(i,i)
       enddo
       do i=1,nzs                                ! subtract error term from old solution
          p(i)=p(i)-r(i)
       enddo

c      ...deallocate variables
       deallocate(fvec)
       deallocate(fvec_save)
       deallocate(fjac)
       deallocate(fjac_save)
       deallocate(vv)
       deallocate(r)
       deallocate(indx)

       return
       end
