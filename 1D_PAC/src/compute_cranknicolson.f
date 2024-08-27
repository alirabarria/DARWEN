
       subroutine compute_cranknicolson
       
       include 'pac.common'       
       integer i,j,nt,nmaxt,ntsteady,itsteady,ihou,imin,isec
       real*8 t,tmin,tmax,maxdy,wkerr

c      ...Set initial conditions
       t=0.0d0                                   ! initial time
       cn_dt=1.0d-12                             ! initial time step
       tmin=day                                  ! minimum time to integrate
       tmax=13.75d9*year                         ! maximum time: age of the Universe
       nt=1                                      ! number of time integration steps
       nmaxt=10000                               ! maximum number of time integration steps
       abunminsteady=1.0d-20                     ! lower abundance cutoff to evaluate steady state
       maxerrsteady=1.0d-2                       ! relative tolerance to consider steady state
       ntsteady=3                                ! minimum number of consecutive steps with dabun<maxerrsteady to consider steady state
       itsteady=0                                ! number of consecutive steps with dabun<maxerrsteady
       itime=0                                   ! counter of integration times to write out results
       call conservation(0)
       write(*,*)

c      ...get rate constants, photorates, and diffusion terms
       call get_rate_constants
       call get_diffusion_terms
       call get_photorates

c Start time integration of equation of continuity
       do while(t.lt.tmax.and.nt.le.nmaxt.and..not.steady)

c      ...evaluate derivatives and jacobian
          call evaluate_derivatives
          call evaluate_jacobian

c      ...integrate continuity equation from t to t+dt
          call solve_matrix_system

c      ...evaluate maximum relative molar fraction variation
          maxdy=0.0d0
          do j=1,nz
             do i=1,nspec
                if(cn_dy(i,j)/abun(i,j).gt.maxdy)
     #          maxdy=cn_dy(i,j)/abun(i,j)
             enddo
          enddo

c      ...variation > 10%, decrease time step and go back
          if(maxdy.gt.0.1d0)then
             cn_dt=0.75d0*cn_dt
             nt=nt+1
             cycle
          endif

c      ...update time step and abundances
          if(maxdy.gt.0.05d0)then             ! ...update time step dt
             cn_dt=cn_dt
          elseif(maxdy.le.0.05d0.and.maxdy.gt.0.01d0)then
             cn_dt=1.2d0*cn_dt
          elseif(maxdy.le.0.01d0.and.maxdy.gt.0.001d0)then
             cn_dt=1.5d0*cn_dt
          else
             cn_dt=2.0d0*cn_dt
          endif
          t=t+cn_dt                           ! ...update time
          do j=1,nz                           ! ...update abundances at t+dt
             do i=1,nspec
                abun(i,j)=abun(i,j)+cn_dy(i,j)
             enddo
          enddo

c      ...verify conservation of elements
          call conservation(2)

c      ...verify whether steady state has been reached
          if(t.gt.tmin)then
             errsteady=tiny
             do j=1,nz
                do i=1,nspec
                   if(abun(i,j).lt.abunminsteady)cycle
                   wkerr=dabs(cn_dy(i,j)/(abun(i,j)-cn_dy(i,j)))
                   if(wkerr.gt.errsteady)errsteady=wkerr
                enddo
             enddo
             if(errsteady.ge.maxerrsteady)itsteady=0
             if(errsteady.lt.maxerrsteady)itsteady=itsteady+1
             if(itsteady.ge.ntsteady)steady=.true.
             itime=itime+1
             time=t
             call write_results
          else
             errsteady=0.0d0
          endif

c      ...write on screen current time
          call get_cputime(ihou,imin,isec)
          write(*,1000)t,ihou,imin,isec,errsteady

          nt=nt+1                                ! a further time integration step
c Finish time integration of equation of continuity
       enddo

       if(steady)then
          write(*,1010)t
       else
          write(*,1020)t
       endif

       return

1000   format(2x,'time[s]=',1pg10.3,3x,
     #           'CPU=',i3,'h 'i2,'m ',i2,'s',3x,
     #           'err=',1pg10.3) 
1010   format(1x,' I- Steady state reached at t=',1pg8.1,'s')
1020   format(1x,' I- Integration done up to t=',1pg8.1,'s')

       end




c_______________________________________________________________________

       subroutine solve_matrix_system
c      Solve system of linear algebraic equations to get abundances at time (t+dt)
       
       include 'pac.common'       
       integer i,j,k
       real*8 a(nmaxspec,nmaxspec,nmaxz),b(nmaxspec,nmaxspec,nmaxz-1),
     # c(nmaxspec,nmaxspec,nmaxz-1),yfwd(nmaxspec,nmaxz)

c      ...create jacobian matrix (tridiagonal by blocks)
       a(:,:,:)=0.0d0                            ! assign set of matrix A in main diagonal
       do j=1,nz
          do i=1,nspec
             do k=1,nspec
                if(i.eq.k)then
                   a(i,k,j)=density(j)/cn_dt-cn_theta*df(i,j,k)
                else
                   a(i,k,j)=-cn_theta*df(i,j,k)
                endif
             enddo
          enddo
       enddo
       b(:,:,:)=0.0d0                            ! assign set of matrix B in subdiagonal
       do j=2,nz
          do i=1,nspec
             b(i,i,j-1)=-cn_theta*df_jm_dif(i,j)
          enddo
       enddo
       c(:,:,:)=0.0d0                            ! assign set of matrix C in superdiagonal
       do j=1,nz-1
          do i=1,nspec
             c(i,i,j)=-cn_theta*df_jp_dif(i,j)
          enddo
       enddo

c      ...LU decompose full (tridiagonal by blocks) matrix
       call lubtridag(nspec,nmaxspec,nz,nmaxz,a,b,c)

c      ...solve system of equations (get dy_i^j) by forward and back substitution
       do j=1,nz                                 ! forward substitution
          do i=1,nspec
             yfwd(i,j)=f(i,j)
             if(j.gt.1)then
                do k=1,nspec
                   yfwd(i,j)=yfwd(i,j)-b(i,k,j-1)*yfwd(k,j-1)
                enddo
             endif
             do k=1,i-1
                yfwd(i,j)=yfwd(i,j)-a(i,k,j)*yfwd(k,j)
             enddo
          enddo
       enddo
       do j=nz,1,-1                              ! back substitution
          do i=nspec,1,-1
             cn_dy(i,j)=yfwd(i,j)
             do k=i+1,nspec
                cn_dy(i,j)=cn_dy(i,j)-a(i,k,j)*cn_dy(k,j)
             enddo
             if(j.lt.nz)then
                do k=1,nspec
                   cn_dy(i,j)=cn_dy(i,j)-c(i,k,j)*cn_dy(k,j+1)
                enddo
             endif
             cn_dy(i,j)=cn_dy(i,j)/a(i,i,j)
          enddo
       enddo

       return
       end




c_______________________________________________________________________

       subroutine lubtridag(n,nmax,m,mmax,a,b,c)
c      Block LU decomposition of a block tridiagonal matrix

c      m: dimension of full tridiagonal matrix by blocks, that we call F 
c      n: dimension of block matrix, that we call A,B,C
c      a(i,j,k): block matrix A(i,j) in main diagonal of F, in position F(k,k)
c      b(i,j,k): block matrix B(i,j) in   subdiagonal of F, in position F(k+1,k)
c      c(i,j,k): block matrix C(i,j) in superdiagonal of F, in position F(k,k+1)
c
c      The resulting L and U matrix are overwritten in A,B,C
c      The full L matrix is contained in B and in the lower triangular values of A (diagonal is 1 and is not stored)
c      The full U matrix is contained in C and in the diagonal and upper triangular values of A

       implicit none
       integer n,nmax,m,mmax,i,j,k,l
       real*8 a(nmax,nmax,mmax),b(nmax,nmax,mmax-1),
     # c(nmax,nmax,mmax-1),aux,minv(nmax,nmax),maux(nmax,nmax)

       do i=1,m

c      ...LU decompose block matrix A in position F(i,i)
          if(i.gt.1)then
             do j=1,n
                do k=1,n
                  aux=0.0d0
                  do l=1,n
                     aux=aux+b(j,l,i-1)*c(l,k,i-1)
                  enddo
                  a(j,k,i)=a(j,k,i)-aux
                enddo
             enddo
          endif
          call ludecomp(a(:,:,i),n,nmax)
          if(i.eq.m)exit

c      ...get L(i+1,i) and store it in B in position F(i+1,i)
          call invtri('U',n,nmax,a(:,:,i),minv)
          do j=1,n
             do k=1,n
                maux(j,k)=0.0d0
                do l=1,n
                   maux(j,k)=maux(j,k)+b(j,l,i)*minv(l,k)
                enddo
             enddo
          enddo
          b(:,:,i)=maux(:,:)

c      ...get U(i,i+1) and store it in C in position F(i,i+1)
          call invtri('L',n,nmax,a(:,:,i),minv)
          do j=1,n
             do k=1,n
                maux(j,k)=0.0d0
                do l=1,n
                   maux(j,k)=maux(j,k)+minv(j,l)*c(l,k,i)
                enddo
             enddo
          enddo
          c(:,:,i)=maux(:,:)
       enddo

       return
       end




c_______________________________________________________________________

       subroutine invtri(lu,n,nmax,a,inva)
c      Calculate the inverse of a Lower/Upper triangular matrix
c      Lower: only lower triangular elements of A (a_ij, j<i) are considered
c             diagonal values are not considered but assumed to be 1
c      Upper: only upper triangular elements of A (a_ij, j>=i) are considered

       implicit none
       integer n,nmax,i,j,k
       real*8 a(nmax,nmax),inva(nmax,nmax)
       character*1 lu

       if(lu.ne.'L'.and.lu.ne.'U') stop ' LU error in invtri'
       inva(:,:)=0.0d0

c      ...calculate the inverse of a lower triangular matrix by forward substitution
       if(lu.eq.'L')then
          do j=n,1,-1
             inva(j,j)=1.0d0
             do i=j+1,n
                do k=1,i-1
                   inva(i,j)=inva(i,j)-a(i,k)*inva(k,j)
                enddo
             enddo
          enddo

c      ...calculate the inverse of an upper triangular matrix by back substitution
       elseif(lu.eq.'U')then
          do j=1,n
             inva(j,j)=1.0d0/a(j,j)
             do i=j-1,1,-1
                do k=i+1,n
                   inva(i,j)=inva(i,j)-a(i,k)*inva(k,j)
                enddo
                inva(i,j)=inva(i,j)/a(i,i)
             enddo
          enddo
       endif

       return
       end




c_______________________________________________________________________

       subroutine ludecomp(a,n,np)
c      From Fortran Numerical Recipes (ludcmp) but without permutation

       implicit none
       integer n,np,i,j,k
       real*8 a(np,np),aamax,dum,sum,vv(np)

       do i=1,n
          aamax=0.0d0
          do j=1,n
             if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
          enddo
          if (aamax.eq.0.0d0) stop ' E- singular matrix in ludecomp'
          vv(i)=1.0d0/aamax
       enddo
       do j=1,n
          do i=1,j-1
             sum=a(i,j)
             do k=1,i-1
                sum=sum-a(i,k)*a(k,j)
             enddo
             a(i,j)=sum
          enddo
          do i=j,n
             sum=a(i,j)
             do k=1,j-1
                sum=sum-a(i,k)*a(k,j)
             enddo
             a(i,j)=sum
          enddo
          if(a(j,j).eq.0.0d0) stop ' E- singularity in ludcmp'
          if(j.ne.n)then
             dum=1.0d0/a(j,j)
             do i=j+1,n
                a(i,j)=a(i,j)*dum
             enddo
          endif
       enddo

       return
       end
