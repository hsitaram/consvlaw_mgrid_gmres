program convdiffmain

      use solvergmres
      use convdiff
      implicit none

      integer :: i,nv,j
      integer :: np
      integer :: m,n,it
      real*8, allocatable :: x0(:),x(:),b(:),xold(:)
      real*8, allocatable :: dcoeff(:),vel(:),reac(:)
      real*8, allocatable :: source(:),ax(:)
      real*8 :: timederivfactor
      logical :: dirc_bc_flags(2)
      logical :: flux_bc_flags(2)
      real*8 :: dircvals(2)
      real*8 :: fluxvals(2)
      real*8 :: dt,xcoord
      integer :: ntsteps
      real*8,allocatable :: res(:)
      real*8 :: resnorm

      print *,"**********************************************"
      print *,"Solving convection diffusion reaction source equations &
		      using upwind discretization and semi-implicit time"
      print *,"**********************************************"
      
      open(unit=5,file='soln.dat')

      !mesh and time step
      np = 257
      h = 1.d0/(np-1)
      n = np
      m = 5
      it = 40
      ntsteps=1
      dt = 10000.d0

      allocate(x0(n))
      allocate(x(n))
      allocate(b(n))
      allocate(xold(n))
      allocate(ax(n))
      
      allocate(dcoeff(n))
      allocate(vel(n))
      allocate(reac(n))
      allocate(source(n))
      allocate(res(n))

!=============================================================      
!set values for variables      
!=============================================================     
      !dcoeff=1.d0
      !vel=0.d0
      timederivfactor=0.d0
      do i=1,n
         xcoord    = (i-1)*h
         dcoeff(i) = xcoord
	 vel(i)    = 1.d0
         reac(i)   = 0.d0
      	 source(i) = -2.d0*xcoord
      enddo
  
      dirc_bc_flags(1) = .false.
      dirc_bc_flags(2) = .true.
      flux_bc_flags(1) = .true.
      flux_bc_flags(2) = .false.

      dircvals(1) = 0.d0
      dircvals(2) = 0.d0
      fluxvals(1) = 0.d0
      fluxvals(2) = 0.d0
     
      ax   = 0.d0 
      xold = 0.d0
      x0 = 0.d0
      x  = 0.d0
      b  = 0.d0
      do i=1,n
        xcoord=(i-1)*h
      	!xold(i)=xcoord*xcoord-xcoord
	!xold(i)=-0.25
	!x = xcoord
     enddo
!===============================================================

      do i=1,ntsteps

      	call findrhs(b,xold,timederivfactor,source,dirc_bc_flags,&
			 flux_bc_flags,dircvals,fluxvals,h,dt,n)
	x=xold

	!call gauss_seidel_smoothing(res,b,x,timederivfactor,vel,dcoeff,reac,&
	!	dirc_bc_flags,flux_bc_flags,dircvals,&
	!			fluxvals,h,dt,n,500)

	!do nv=1,1

	!	call dovcycle(x,b,timederivfactor,vel,dcoeff,reac,&
	!	dirc_bc_flags,flux_bc_flags,dircvals,&
	!			fluxvals,h,dt,n)

!		call findAX(ax,x,timederivfactor,vel,dcoeff,reac,dirc_bc_flags,&
!			flux_bc_flags,dircvals,fluxvals,h,dt,n)
!		res=b-ax
!		resnorm=0.d0
!		do j=1,n
!			resnorm=resnorm+res(j)*res(j)
!		enddo
!		resnorm=sqrt(resnorm)
!		print *,"resnorm:",resnorm
!	enddo
 	call performgmres(b,xold,x,timederivfactor,vel,dcoeff,reac,dirc_bc_flags,&
			flux_bc_flags,dircvals,fluxvals,h,dt,&
			m,n,it,findAX,mgridprecond)
	
	!residual norm calculation
	resnorm=0.d0
	call findAX(ax,x,timederivfactor,vel,dcoeff,reac,dirc_bc_flags,&
		flux_bc_flags,dircvals,fluxvals,h,dt,n)
	res = b - ax
	do j=1,n
		resnorm=resnorm+res(j)*res(j)
	enddo
	resnorm=sqrt(resnorm)
	print *,"resnorm:",resnorm
	
	xold=x
	print *,"************Finished timestep:",i

      enddo

      do i=1,n
      	write(5,'(F10.5 F10.5)'),(i-1)*h,x(i)
      enddo
      
      close(5)	

end program convdiffmain 
