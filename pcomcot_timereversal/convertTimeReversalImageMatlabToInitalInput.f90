!!    Input:     _xcoordinates01_interp.dat, _ycoordinates01_interp.dat, _InitialElevation_interp.dat
!!               _xcoordinates01.dat, _ycoordinates01.dat
!!    Output:    InitialElevation.xyz

      implicit NONE
      integer*4               ::  ios, NX, NY
      real*8, allocatable     ::  H(:,:), x(:), y(:)
      integer*4               ::  i, j, ii
      real*8                  ::  h0,tmp

      integer*4               ::  NXF, NYF
      real*8, allocatable     ::  HF(:,:), xF(:), yF(:)


      NX = 0; NY = 0
      open(1,file='_xcoordinate01_interp.dat',form='formatted',status='old')
      do
          read(1,*,iostat=ios) tmp
          if(ios.ne.0) exit
          NX = NX+1
      enddo
      close(1)
      open(1,file='_ycoordinate01_interp.dat',form='formatted',status='old')
      do
          read(1,*,iostat=ios) tmp
          if(ios.ne.0) exit
          NY = NY+1
      enddo
      close(1)

      ALLOCATE(H(NX,NY));
      ALLOCATE(x(NX));    ALLOCATE(y(NY))

      open(1,file='_xcoordinate01_interp.dat',form='formatted',status='old')
      do i = 1,NX
          read(1,*) x(i)
      enddo
      close(1)

      open(1,file='_ycoordinate01_interp.dat',form='formatted',status='old')
      do j = 1,NY
          read(1,*) y(j)
      enddo
      close(1)

      open(1,file='_InitialElevation_interp.dat',form='formatted')
      do j = 1,NY
          read(1,*) (H(i,j), i=1,NX)
      enddo
      close(1)



      NXF = 0; NYF = 0
      open(1,file='_xcoordinate01.dat',form='formatted',status='old')
      do
          read(1,'(f15.5)',advance='no',iostat=ios) tmp
          if(ios.ne.0) exit
          NXF = NXF+1
      enddo
      close(1)
      open(1,file='_ycoordinate01.dat',form='formatted',status='old')
      do
          read(1,'(f15.5)',advance='no',iostat=ios) tmp
          if(ios.ne.0) exit
          NYF = NYF+1
      enddo
      close(1)

      ALLOCATE(HF(NXF,NYF)); ALLOCATE(xF(NXF)); ALLOCATE(yF(NYF))
      
      open(1,file='_xcoordinate01.dat',form='formatted',status='old')
      read(1,*) (xF(i), i=1,NXF)
      close(1)

      open(1,file='_ycoordinate01.dat',form='formatted',status='old')
      read(1,*) (yF(j), j=1,NYF)
      close(1)

      do j = 1,NYF
          do i = 1,NXF
              call interpH(x,y,H,NX,NY,xF(i),yF(j),h0)
              HF(i,j) = h0
          enddo
      enddo


      write(*,*) 'writing to file...'
      open(1,file='InitialElevation.xyz',form='formatted')
      do j = 1,NYF
          do i = 1,NXF
              write(1,*) xF(i), yF(j), HF(i,j)
          enddo
      enddo
      close(1)

      end



      subroutine interpH(x,y,H,NX,NY,x0,y0,z0)

      implicit NONE
      real*8     ::  x(NX), y(NY), H(NX,NY), x0, y0, z0
      integer*4  ::  NX, NY, i, j, i0, j0
      real*8     ::  x1, x2, y1, y2, z1, z2, z3, z4

      if(x0.lt.x(1).or.x0.gt.x(NX).or.y0.lt.y(1).or.y0.gt.y(NY)) then
          z0 = 0.0d0
          return
      endif

      i0 = -1
      do i = 1,NX-1
          if(x0.ge.x(i).and.x0.le.x(i+1)) then
              i0 = i; exit
          endif
      enddo

      j0 = -1
      do j = 1,NY-1
          if(y0.ge.y(j).and.y0.le.y(j+1)) then
              j0 = j; exit
          endif
      enddo

      if(i0.eq.-1.or.j0.eq.-1) then
          z0 = 0.0d0
          return
      endif

      x1 = x(i0); x2 = x(i0+1)
      y1 = y(j0); y2 = y(j0+1)
      z1 = H(i0,j0); z2 = H(i0+1,j0);
      z3 = H(i0,j0+1); z4 = H(i0+1,j0+1)
      z0 = z1*(x2-x0)*(y2-y0)+z2*(x0-x1)*(y2-y0)+z3*(x2-x0)*(y0-y1)+z4*(x0-x1)*(y0-y1)
      z0 = z0/(x2-x1)/(y2-y1)

      end subroutine interpH

