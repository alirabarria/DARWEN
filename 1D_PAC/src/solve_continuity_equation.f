
       subroutine solve_continuity_equation

       include 'pac.common'

       if(diffusion)                    write(*,1000)
       if(diffusion.and.nummethod.eq.1) write(*,1010)
       if(.not.diffusion)               write(*,1020)
       if(photochem)                    write(*,1030)
       if(.not.photochem)               write(*,1040)
       if(sensi.gt.0)                   write(*,1050)
       if(sensi.eq.0)                   write(*,1060)

       if(nummethod.eq.1)then
c          call compute_newton
          call compute_newtonraphson
c          call compute_newtonraphson_oldstyle
       elseif(nummethod.eq.2)then
          call compute_dlsodes
       elseif(nummethod.eq.3)then
          call compute_cranknicolson
       else
          stop' E- Wrong requested numerical method'
       endif

       return

1000   format(1x,' I- Diffusion included')
1010   format(1x,' I- Molecular diffusion neglected (nummethod=1)')
1020   format(1x,' I- Diffusion not included')
1030   format(1x,' I- Photochemistry included')
1040   format(1x,' I- Photochemistry not included')
1050   format(1x,' I- Sensitivity analysis included')
1060   format(1x,' I- Sensitivity analysis not included')

       end




c_______________________________________________________________________

       subroutine normalize_abun
c      ...normalize abundances so that (Sum_i[abun_i]=1 for each layer)

       include 'pac.common'
       integer i,j
       real*8 rwk

       do j=1,nz
          rwk=0.0d0
          do i=1,nspec
             rwk=rwk+abun(i,j)
          enddo
          do i=1,nspec
             abun(i,j)=abun(i,j)/rwk
          enddo
       enddo

       return
       end




c_______________________________________________________________________

       subroutine conservation(ido)
c      ...compute column densities of atoms as
c         the total number of atoms in the atmosphere divided by planet surface area

       include 'pac.common'
       real*8 constol                            ! tolerance for conservation
       parameter(constol=0.05d0)                 ! 5 %

       integer i,j,k,ido
       real*8 dz,rm,rp,err

       do k=1,nelem
          cdelem(k)=0.0d0
          do j=1,nz
             if(j.eq.1)then
                dz=z(j)-z(j+1)
             elseif(j.eq.nz)then
                dz=z(j-1)-z(j)
             else
                dz=dabs(0.5d0*(z(j-1)+z(j+1)))
             endif
             rm=rplanet+z(j)-0.5d0*dz
             rp=rplanet+z(j)+0.5d0*dz
             do i=1,nspec
                if(nat(i,k).eq.0)cycle
                cdelem(k)=cdelem(k)+nat(i,k)*abun(i,j)*density(j)*
     #          ((rp/rplanet)**3.0d0-(rm/rplanet)**3.0d0)
             enddo
          enddo
          cdelem(k)=cdelem(k)*rplanet/3.0d0
       enddo

       if(ido.eq.0)then                          ! compute initial column densities
          do k=1,nelem
             cdelem0(k)=cdelem(k)
          enddo
       elseif(ido.eq.1)then                      ! simply compute current column densities and return
          return
       elseif(ido.eq.2)then                      ! compare with initial column densities
          i=0                                    ! ...and if different stop
          do k=1,nelem
             err=dabs((cdelem(k)-cdelem0(k))/cdelem0(k))
             if(err.gt.constol)then
                write(*,*)' E- Error in conservation for ',elem(k)
                write(*,*)'    Current-to-initial column density =',
     #          cdelem(k)/cdelem0(k)
                errflag=1
                i=i+1
             endif
          enddo
c          if(i.gt.0)stop
       else
          stop ' E- Wrong input integer in conservation'
       endif

       return
       end




c_______________________________________________________________________

       subroutine conservation_layer(ido)
c      ...compute sum of molar fraction of each atom in any form, for each layer

       include 'pac.common'
       integer i,j,k,ido
       real*8 err

       do k=1,nelem
          do j=1,nz
             elab(k,j)=0.0d0
             do i=1,nspec
                if(nat(i,k).eq.0)cycle
                elab(k,j)=elab(k,j)+nat(i,k)*abun(i,j)
             enddo
          enddo
       enddo

       if(ido.eq.0)then                          ! compute initial values
          do k=1,nelem
             do j=1,nz
                elab0(k,j)=elab(k,j)
             enddo
          enddo
       elseif(ido.eq.1)then                      ! simply compute current column densities and return
          return
       elseif(ido.eq.2)then                      ! compare with initial column densities
          i=0                                    ! ...and if different stop
          do j=1,nz
             do k=1,nelem
                err=dabs((elab(k,j)-elab0(k,j))/elab(k,j))
                if(err.gt.5.0d-3)then
                   write(*,*)' E- Element: ',elem(k)
                   write(*,*)'    layer:   ',j,z(j)/1.0d5,' km '
                   write(*,*)'    current-to-initial abundances sum =',
     #             elab(k,j)/elab0(k,j)
                   i=i+1
                endif
             enddo
          enddo
          if(i.gt.0)stop
       else
          stop ' E- Wrong input integer in conservation_layer'
       endif

       return
       end
