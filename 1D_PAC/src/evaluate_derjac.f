
       subroutine evaluate_derivatives
c      Evaluate function f in continuity equation: f_i^j = d[n_i^j]/dt  
       USE IEEE_ARITHMETIC
       include 'pac.common'       
       integer i,j,k,l,id
       real*8 rwk
c Loop over height and species
       do j=1,nz
          do i=1,nspec
c            if (j.eq.1)then
c             write(*,*)'abun',abun(idreag(id,1),j),j,spec(i)
c             write(*,*)
c            endif
c      ...evaluate chemical rate of formation/destruction: f_i^j = dn_i^j/dt 
             f(i,j)=0.0d0
             do k=1,nterm(i)
                id=idterm(i,k)
                rwk=1.0d0
                do l=1,nreag(id)
                   rwk=rwk*abun(idreag(id,l),j)*density(j)
                enddo   
                f(i,j)=f(i,j)+sterm(i,k)*fterm(i,k)*krate(id,j)*rwk

                if(IEEE_IS_NAN(f(i,j)))then     
                 write(*,*)'NaN found',f(i,j)
                 write(*,*)'spec',spec(i),i
                 write(*,*)'krate ',krate(id,j)
                 write(*,*)'rwk',rwk
                 write(*,*)'density',density(j),j
                 write(*,*)'abun',abun(idreag(id,1),j)
                 errflag=2  ! @Arturo, I want this to stop if NaNs are found?
                 write(*,*)
               endif

             enddo
             if (writekrate)chem_ctbn(i,j)=f(i,j)


c      ...add diffusion contribution to f_i^j 
             if(diffusion)then
                           f(i,j)=f(i,j)+df_dif(i,j)   *abun(i,j)
                if(j.gt.1) f(i,j)=f(i,j)+df_jm_dif(i,j)*abun(i,j-1)
                if(j.lt.nz)f(i,j)=f(i,j)+df_jp_dif(i,j)*abun(i,j+1)
                if (writekrate)diff_ctbn(i,j)=f(i,j)-chem_ctbn(i,j)
             endif

c      ...add photochemical contribution to f_i^j 
             if(photochem)then
                do k=1,nphototerm(i)
                   id=idphototerm(i,k)
                   f(i,j)=f(i,j)+sphototerm(i,k)*fphototerm(i,k)*
     #             jphot(id,j)*abun(idphotoreag(id),j)*density(j)
                enddo
             endif

c End of loop over height and species
          enddo
       enddo

       return
       end




c_______________________________________________________________________

       subroutine evaluate_jacobian
c      Evaluate derivatives of function f with respect to molar fraction of species within same layer 
c      df(i,j,k) = d[f_i^j]/dy_k^j = d[dn_i^j/dt]/dy_k^j
       
       include 'pac.common'       
       integer i,j,k,l,g,gwk,id
       real*8 rwk,iswk

c Loop over height and species i,k
       do j=1,nz
          do i=1,nspec
             do k=1,nspec

c      ...evaluate jacobian of formation/destruction: d[dn_i^j/dt]/dy_k^j = d[dn_i^j/dt]/dn_k^j * n^j 
                df(i,j,k)=0.0d0
                do l=1,nterm(i)
                   id=idterm(i,l)
                   iswk=0.0d0
                   gwk=0
                   do g=1,nreag(id)
                      if(idreag(id,g).eq.k)then
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
                      df(i,j,k)=df(i,j,k)+density(j)*
     #                sterm(i,l)*fterm(i,l)*krate(id,j)*rwk*iswk
                   endif
                enddo

c     ...add diffusion contribution to jacobian: d[dn_i^j/dt]/dy_k^j 
               if(diffusion)then
                   if(i.eq.k)then
                      df(i,j,k)=df(i,j,k)+df_dif(i,j)
                   endif
                endif

c     ...add photochemical contribution to jacobian: d[dn_i^j/dt]/dy_k^j 
                if(photochem)then
                   do l=1,nphototerm(i)
                      id=idphototerm(i,l)
                      if(idphotoreag(id).eq.k)then
                         df(i,j,k)=df(i,j,k)+density(j)*
     #                   sphototerm(i,l)*fphototerm(i,l)*jphot(id,j)
                      endif
                    enddo
                 endif

c End of loop over height and species i,k
             enddo
          enddo
       enddo

       return
       end
