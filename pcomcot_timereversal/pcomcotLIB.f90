
      subroutine readConfig(GP)
      
      use VariableDefination
      implicit NONE
      include "mpif.h"
      type(GlobalParameters)  ::  GP
      integer*4               ::  ios, i, ierror
      character(999)          :: s

      open(666,file=TRIM(ADJUSTL(GP%COMCOTParametersFileName)),status='old',form='formatted',iostat=ios)
      if(ios.ne.0) then
          write(*,*)
          write(*,*) 'ERROR: cann''t find file comcot.ctl'
          write(*,*)
          call MPI_ABORT(MPI_COMM_WORLD, ierror)
      endif

      write(*,*)
      write(*,*) '##################################################'
      write(*,*) '#                                                #'
      write(*,'(a,f4.1,a)') ' #         PCOMCOT Simulation Version',GP%Version,'         #'
      write(*,*) '#                                                #'
      write(*,*) '##################################################'
      do i = 1,5
          read(666, '(a)') s
      enddo
      read(666, '(2/)')
      read(666, '(49x,i30)')   GP%CoordinatesType
      read(666, '(49x,i30)')   GP%GoverningEquationType
      read(666, '(49x,f30.6)') GP%WaterDepthLimit
      read(666, '(49x,i30)')   GP%InitialConditionType
      read(666, '(49x,i30)')   GP%BoundaryConditionType
      read(666, '(49x,f30.6)') GP%TotalTime
      read(666, '(49x,f30.6)') GP%dt
      read(666, '(49x,f30.6)') GP%DTSaveSTART
      read(666, '(49x,f30.6)') GP%DTSaveData
      read(666, '(49x,i30)')   GP%SaveFlux
      read(666, '(49x,i30)')   GP%ComputeGreen
      read(666, '(49x,i30)')   GP%MinGridsPerNode
      read(666, '(49x,i30)')   GP%FeedbackToParentLayer
      if(GP%InitialConditionType.eq.0.and.GP%ComputeGreen.eq.1) then
          write(*,*) "ERROR: to compute Green's fucntions, initial condition type must be 'fault'."
          call MPI_ABORT(MPI_COMM_WORLD, ierror)
      endif

      close(666)

      GP%nCalculations = 1
      GP%TotalTimeSteps = MAX(1,NINT(GP%TotalTime/GP%dt))
      GP%NDTSaveData = MAX(1,NINT(GP%DTSaveData/GP%dt))
      if(GP%BoundaryConditionType.ne.1) then
          write(*,*) "ERROR: only wall boundary is supported."
          call MPI_ABORT(MPI_COMM_WORLD, ierror)
      endif
      if(GP%GoverningEquationType.ne.0) then
          write(*,*) "ERROR: only linear equation is supported."
          call MPI_ABORT(MPI_COMM_WORLD, ierror)
      endif

      end subroutine ReadConfig



      subroutine checkFiles(GP)
      
      use VariableDefination
      implicit NONE
      include "mpif.h"
      type(GlobalParameters)  ::  GP
      integer*4               ::  iLayer, ios, ierror
      character(999)          ::  fn1, fn2, s
      logical                 ::  fe, fe2

      write(*,*)
      write(*,*) 'checking files that are needed ...'
      write(*,*)

      write(*,'(a)',advance='no') '         Bathymetry Data File:  '
      GP%NumLayers = 0
      do iLayer=1, 99
          write(fn1,'(a,i2.2,a)') TRIM(ADJUSTL(GP%BathymetryFilePrefix)),iLayer,'.nf'
          inquire(file=fn1, exist=fe)
          write(fn2,'(a,i2.2,a)') TRIM(ADJUSTL(GP%BathymetryFilePrefix)),iLayer,'.xyz'
          inquire(file=fn2, exist=fe2)
          if(fe.or.fe2) then
          GP%NumLayers = GP%NumLayers+1
          if(GP%NumLayers.ne.1.and.MOD(GP%NumLayers-1,5).eq.0) then
              write(*,*)
              write(*,'(a)',advance='no'),'                                '
          elseif(GP%NumLayers.ne.1) then
              write(*,'(a)',advance='no') ', '
          endif
          if(fe) then
              write(*,'(a11)',advance='no') ADJUSTL(fn1)
          else
              write(*,'(a11)',advance='no') ADJUSTL(fn2)
          endif
          endif
      enddo
      write(*,*)
      if(GP%NumLayers.eq.0) then
          write(*,*) 'ERROR: cann''t find any bathymetry data files:  layerXX.nf or layerXX.xyz'
          call MPI_ABORT(MPI_COMM_WORLD, ierror)
      endif

      if(GP%InitialConditionType.eq.0) then
          fn1 = TRIM(ADJUSTL(GP%InitialConditionFilePrefix))//'.nf'
          inquire(file=fn1, exist=fe)
          fn2 = TRIM(ADJUSTL(GP%InitialConditionFilePrefix))//'.xyz'
          inquire(file=fn2, exist=fe2)
          if((.not.fe).and.(.not.fe2)) then
              write(*,*) 'ERROR: cann''t find initial water elevation file:  ',TRIM(ADJUSTL(fn1)),' or ',TRIM(ADJUSTL(fn2))
              call MPI_ABORT(MPI_COMM_WORLD, ierror)
          endif
          if(fe) then
              GP%InitialConditionFileName = fn1
              GP%InitialConditionFileFormat = 1
          else if(fe2) then
              GP%InitialConditionFileName = fn2
              GP%InitialConditionFileFormat = 2
          endif
          write(*,*) '        Initial Water Elevation Data File:  ',TRIM(ADJUSTL(GP%InitialConditionFileName))
      endif
      
      if(GP%InitialConditionType.eq.1) then
          inquire(file=TRIM(ADJUSTL(GP%FaultParametersFileName)),exist=fe)
          if(.not.fe) then
              write(*,*) 'ERROR: cann''t find muli-fault information file:  ', TRIM(ADJUSTL(GP%FaultParametersFileName))
              call MPI_ABORT(MPI_COMM_WORLD, ierror)
          endif
          write(*,*) '        Muli-Fault Information File:  ',TRIM(ADJUSTL(GP%FaultParametersFileName))
      endif

      inquire(file=TRIM(ADJUSTL(GP%StationFileName)),exist=fe)
      if(.not.fe) then
          write(*,*) 'ERROR: cann''t find stations information file:  ',TRIM(ADJUSTL(GP%StationFileName))
          call MPI_ABORT(MPI_COMM_WORLD, ierror)
      endif
      write(*,*) '        Stations coordinates file:  ',TRIM(ADJUSTL(GP%StationFileName))

      end subroutine checkFiles



      subroutine getBathymetryDataSize(GP, LP)
      
      use VariableDefination
      implicit NONE
      include "mpif.h"
      type(GlobalParameters)  ::  GP
      type(LayerParameters)   ::  LP(100)
      integer*4               ::  iLayer, ierror, iLayerCheck
      character(999)          ::  fn1, fn2
      logical                 ::  fe, fe2
      real*8     ::  x1,x2,x,y1,y
      integer*4  ::  i, ios, endoffile


      write(*,*)
      write(*,*) 'getting bathymetry data size ...'
      write(*,*)

      iLayer = 0
      do iLayerCheck=1, 99
          write(fn1,'(a,i2.2,a)') TRIM(ADJUSTL(GP%BathymetryFilePrefix)),iLayerCheck,'.nf'
          inquire(file=fn1, exist=fe)
          write(fn2,'(a,i2.2,a)') TRIM(ADJUSTL(GP%BathymetryFilePrefix)),iLayerCheck,'.xyz'
          inquire(file=fn2, exist=fe2)
          if(fe.or.fe2) then
          iLayer = iLayer+1
          if(fe) then
              LP(iLayer)%BathymetryFileName = fn1
              LP(iLayer)%BathymetryFileFormat = 1
          elseif(fe2) then
              LP(iLayer)%BathymetryFileName = fn2
              LP(iLayer)%BathymetryFileFormat = 2
          endif
          if(LP(iLayer)%BathymetryFileFormat.eq.1) then
              call getBathymetryDataSizeNetCDF(GP, LP, iLayer)
          else
              open(23,file=TRIM(ADJUSTL(LP(iLayer)%BathymetryFileName)),status='old',form='formatted',iostat=ios)
              read(23,*) x1,y1
              if(iLayer.eq.1) then
                  GP%StartEastWest = 0
                  if((GP%CoordinatesType.eq.0).and.(x1.ge.0).and.(x1.le.180))then
                      GP%StartEastWest = 1
                  else if((GP%CoordinatesType.eq.0).and.(x1.gt.180.or.x1.lt.0.0))then
                      GP%StartEastWest = 2
                  endif
              endif
              if((GP%StartEastWest.eq.1).and.(x1.lt.0)) x1 = x1+360
              if((GP%StartEastWest.eq.2).and.(x1.gt.180)) x1 = x1-360
              LP(iLayer)%xmin = x1; LP(iLayer)%ymin = y1
              read(23,*) x2
              if((GP%StartEastWest.eq.1).and.(x2.lt.0)) x2 = x2+360
              if((GP%StartEastWest.eq.2).and.(x2.gt.180)) x2 = x2-360
              LP(iLayer)%dx = x2-x1
              LP(iLayer)%NX = 2; LP(iLayer)%NY = 1
              do
                  read(23,*,iostat=ios) x,y
                  if((GP%StartEastWest.eq.1).and.(x.lt.0)) x = x+360
                  if((GP%StartEastWest.eq.2).and.(x.gt.180)) x = x-360
                  if(ABS(x-x1).gt.0.1*LP(iLayer)%dx) then
                      LP(iLayer)%NX = LP(iLayer)%NX+1;  LP(iLayer)%xmax = x
                  else
                      LP(iLayer)%NY = LP(iLayer)%NY+1;  LP(iLayer)%dy = y-y1
                      exit
                  endif
              enddo
              do i=1, LP(iLayer)%NX-1
                  read(23,*) x, LP(iLayer)%ymax
              enddo
              do
                  endoffile = 0
                  do i=1,LP(iLayer)%NX
                      read(23,*,iostat=ios) x, LP(iLayer)%ymax
                      if(ios.ne.0) then
                          endoffile = 1
                          exit
                      endif
                  enddo
                  if(endoffile.eq.1) exit
                  LP(iLayer)%NY = LP(iLayer)%NY + 1
              enddo
              close(23)
          endif

          ! adjust dx and dy
          LP(iLayer)%dx = (LP(iLayer)%xmax - LP(iLayer)%xmin) / LP(iLayer)%NX
          LP(iLayer)%dy = (LP(iLayer)%ymax - LP(iLayer)%ymin) / LP(iLayer)%NY

          write(*,'(a,a11)',advance='no') '        ',ADJUSTL(LP(iLayer)%BathymetryFileName)
          write(*,*)                        &
              ':  xmin: ',LP(iLayer)%xmin,'    xmax: ',LP(iLayer)%xmax, &
              '    dx:',LP(iLayer)%dx,'    nx: ',LP(iLayer)%NX
          write(*,*)                      &
              '                   ',                                    &
              '   ymin: ',LP(iLayer)%ymin,'    ymax: ',LP(iLayer)%ymax, &
              '    dy:',LP(iLayer)%dy,'    ny: ',LP(iLayer)%NY
          write(*,*)
          endif
      enddo

      end subroutine getBathymetryDataSize



      subroutine determineLayerDependency(GP, LP)

      use VariableDefination
      implicit NONE
      include "mpif.h"
      type(GlobalParameters)  ::  GP
      type(LayerParameters)   ::  LP(100)
      integer*4               ::  ierror
      integer*4               ::  i, j, k, NumMotherLayers, NumTopLayers, IsOverlap
      real*8                  ::  p1(4,2), p2(4,2)

      GP%NumLayerLevels = 0
      do i=1, GP%NumLayers
          NumMotherLayers = 0
          do j=1, GP%NumLayers
              if(j.ne.i) then
                  if((LP(i)%xmin.ge.LP(j)%xmin).and.(LP(i)%xmax.le.LP(j)%xmax).and. &
                      (LP(i)%ymin.ge.LP(j)%ymin).and.(LP(i)%ymax.le.LP(j)%ymax)) &
                      NumMotherLayers = NumMotherLayers+1
              endif
          enddo
          LP(i)%Level = NumMotherLayers+1
          GP%NumLayerLevels = MAX(GP%NumLayerLevels, LP(i)%Level)
      enddo
      
      GP%TopLayer = 0; NumTopLayers = 0
      do i=1, GP%NumLayers
          if(LP(i)%Level.eq.1) then
              GP%TopLayer = i
              NumTopLayers = NumTopLayers+1
          endif
      enddo
      if(NumTopLayers.gt.1) then
          write(*,'(a)',advance='no') 'ERROR: more than 1 top layer found: '
          do i=1, GP%NumLayers
              if(LP(i)%Level.eq.1) write(*,'(a,a)',advance='no') TRIM(ADJUSTL(LP(i)%BathymetryFileName)),'  '
          enddo
          write(*,*)
          call MPI_ABORT(MPI_COMM_WORLD, ierror)
      endif

      do i=1, GP%NumLayers
          do j=1, GP%NumLayers
              if(i.ne.j.and.LP(i)%Level.eq.LP(j)%Level) then
                  IsOverlap = 0
                  p1(1,1)=LP(i)%xmin; p1(1,2)=LP(i)%ymin; p1(2,1)=LP(i)%xmax; p1(2,2)=LP(i)%ymin
                  p1(3,1)=LP(i)%xmax; p1(3,2)=LP(i)%ymax; p1(4,1)=LP(i)%xmin; p1(4,2)=LP(i)%ymax
                  p2(1,1)=LP(j)%xmin; p2(1,2)=LP(j)%ymin; p2(2,1)=LP(j)%xmax; p2(2,2)=LP(j)%ymin
                  p2(3,1)=LP(j)%xmax; p2(3,2)=LP(j)%ymax; p2(4,1)=LP(j)%xmin; p2(4,2)=LP(j)%ymax
                  do k=1, 4
                      if(p1(k,1).ge.p2(1,1).and.p1(k,1).le.p2(3,1).and. &
                          p1(k,2).ge.p2(1,2).and.p1(k,2).le.p2(3,2)) then
                          IsOverlap = 1
                          exit
                      endif
                  enddo
                  do k=1, 4
                      if(p2(k,1).ge.p1(1,1).and.p2(k,1).le.p1(3,1).and. &
                          p2(k,2).ge.p1(1,2).and.p2(k,2).le.p1(3,2)) then
                          IsOverlap = 1
                          exit
                      endif
                  enddo
                  if(IsOverlap.eq.1) then
                      write(*,*) "ERROR: layers overlaped: ",               &
                          TRIM(ADJUSTL(LP(i)%BathymetryFileName)), " and ", &
                          TRIM(ADJUSTL(LP(j)%BathymetryFileName))
                      call MPI_ABORT(MPI_COMM_WORLD, ierror)
                  endif
              endif
          enddo
      enddo

      do i=1, GP%NumLayers
          LP(i)%Parent = 0
          do j=1, GP%NumLayers
              if(j.ne.i) then
                  if((LP(i)%xmin.ge.LP(j)%xmin).and.(LP(i)%xmax.le.LP(j)%xmax).and. &
                      (LP(i)%ymin.ge.LP(j)%ymin).and.(LP(i)%ymax.le.LP(j)%ymax)) &
                      LP(i)%Parent = MAX(LP(i)%Parent, j)
              endif
          enddo
      enddo

      do i=1, GP%NumLayers
          write(*,'(a,a11)',advance='no') '        ',ADJUSTL(LP(i)%BathymetryFileName)
          write(*,'(a,i3,a)',advance='no') ':  Level: ',LP(i)%Level,'    Parent: '
          if(LP(i)%Parent.ne.0) then
              write(*,'(a11)') ADJUSTL(LP(LP(i)%Parent)%BathymetryFileName)
          else
              write(*,*)
          endif
      enddo
      write(*,*)
      do i=1, GP%NumLayers
      if(LP(i)%Level.gt.1) then
          if(LP(i)%dx.gt.LP(LP(i)%Parent)%dx) &
              write(*,'(a,a,a,f7.4,a,a,a,f7.4)')                               &
                  'WARNING: ',TRIM(ADJUSTL(LP(i)%BathymetryFileName)),         &
                  ' dx: ',LP(i)%dx,' larger than Parent ',                     &
                  TRIM(ADJUSTL(LP(LP(i)%Parent)%BathymetryFileName)),' dx: ', LP(LP(i)%Parent)%dx
          if(LP(i)%dy.gt.LP(LP(i)%Parent)%dy) &
              write(*,'(a,a,a,f7.4,a,a,a,f7.4)')                               &
                  'WARNING: ',TRIM(ADJUSTL(LP(i)%BathymetryFileName)),         &
                  ' dy: ',LP(i)%dy,' larger than Parent ', &
                  TRIM(ADJUSTL(LP(LP(i)%Parent)%BathymetryFileName)),' dy: ', LP(LP(i)%Parent)%dy
      endif
      enddo


      end subroutine determineLayerDependency



      subroutine readBathymetry(GP, LP)
      
      use VariableDefination
      implicit NONE
      include "mpif.h"
      type(GlobalParameters)  ::  GP
      type(LayerParameters)   ::  LP(100)
      integer*4  ::  iLayer, i, j, ios
      real*8     ::  x, y
      character(999)  ::  s

      write(*,*)
      write(*,*) 'reading bathymetry data ...'
      write(*,*)
      write(*,'(a)',advance='no') '        '

      do iLayer=1, GP%NumLayers
          if(iLayer.ne.1.and.MOD(iLayer-1,5).eq.0) then
              write(*,*)
              write(*,'(a)',advance='no'),'        '
          elseif(iLayer.ne.1) then
              write(*,'(a)',advance='no') ', '
          endif
          write(*,'(a11)',advance='no') ADJUSTL(LP(iLayer)%BathymetryFileName)
          LP(iLayer)%zmax = -20000; LP(iLayer)%zmin = 20000
          if(LP(iLayer)%BathymetryFileFormat.eq.1) then
              call readBathymetryNetCDF(GP, LP, iLayer)
          else
              open(23,file=TRIM(ADJUSTL(LP(iLayer)%BathymetryFileName)), &
                  status='old',form='formatted',iostat=ios)
              do j=1, LP(iLayer)%NY
                  do i=1, LP(iLayer)%NX
                      read(23,*) x,y,LP(iLayer)%Z(i,j)
                      LP(iLayer)%zmax = MAX(LP(iLayer)%zmax, LP(iLayer)%Z(i,j))
                      LP(iLayer)%zmin = MIN(LP(iLayer)%zmin, LP(iLayer)%Z(i,j))
                      if((GP%StartEastWest.eq.1).and.(x.lt.0)) x = x+360
                      if((GP%StartEastWest.eq.2).and.(x.gt.180)) x = x-360
                      if(j.eq.1) LP(iLayer)%X(i) = x
                      if(i.eq.1) LP(iLayer)%Y(j) = y
                  enddo
              enddo
              close(23)
          endif
if(1.eq.0) then
          do j=1, LP(iLayer)%NY
              if(LP(iLayer)%Z(1,j).gt.GP%WaterDepthLimit.and. &
                  LP(iLayer)%Z(2,j).le.GP%WaterDepthLimit) then
                  LP(iLayer)%Z(1,j) = 0.0d0
              endif
              if(LP(iLayer)%Z(LP(iLayer)%NX,j).gt.GP%WaterDepthLimit.and. &
                  LP(iLayer)%Z(LP(iLayer)%NX-1,j).le.GP%WaterDepthLimit) then
                  LP(iLayer)%Z(LP(iLayer)%NX,j) = 0.0d0
              endif
          enddo
          do i=1, LP(iLayer)%NX
              if(LP(iLayer)%Z(i,1).gt.GP%WaterDepthLimit.and. &
                  LP(iLayer)%Z(i,2).le.GP%WaterDepthLimit) then
                  LP(iLayer)%Z(i,1) = 0.0d0
              endif
              if(LP(iLayer)%Z(i,LP(iLayer)%NY).gt.GP%WaterDepthLimit.and. &
                  LP(iLayer)%Z(i,LP(iLayer)%NY-1).le.GP%WaterDepthLimit) then
                  LP(iLayer)%Z(i,LP(iLayer)%NY) = 0.0d0
              endif
          enddo
endif

          do j=1, LP(iLayer)%NY
              do i=1, LP(iLayer)%NX
                  if(LP(iLayer)%Z(i,j).le.GP%WaterDepthLimit) then
                      call smoothBathymetry(GP,LP,iLayer,i,j)
                  endif
if(1.eq.0) then
                  if(LP(iLayer)%Z(i,j).le.GP%WaterDepthLimit.and.    &
                      LP(iLayer)%Z(i+1,j).gt.GP%WaterDepthLimit.and. &
                      LP(iLayer)%Z(i-1,j).gt.GP%WaterDepthLimit.and. &
                      LP(iLayer)%Z(i,j+1).gt.GP%WaterDepthLimit.and. &
                      LP(iLayer)%Z(i,j-1).gt.GP%WaterDepthLimit) then
                      
                      LP(iLayer)%Z(i,j)=0.25*(LP(iLayer)%Z(i+1,j)+LP(iLayer)%Z(i-1,j)+ &
                          LP(iLayer)%Z(i,j+1)+LP(iLayer)%Z(i,j-1))
                  endif
endif
              enddo
          enddo

          write(s,'(a,i2.2,a)') '_xcoordinate',iLayer,'.dat'
          open(23,file=TRIM(ADJUSTL(s)),form='formatted')
          do i=1,LP(iLayer)%NX
              write(23,'(f15.5)',advance='no') LP(iLayer)%X(i)
          enddo
          close(23)
          write(s,'(a,i2.2,a)') '_ycoordinate',iLayer,'.dat'
          open(23,file=TRIM(ADJUSTL(s)),form='formatted')
          do i=1,LP(iLayer)%NY
              write(23,'(f15.5)',advance='no') LP(iLayer)%Y(i)
          enddo
          close(23)
          write(s,'(a,i2.2,a)') '_bathymetry',iLayer,'.dat'
          open(23,file=TRIM(ADJUSTL(s)),form='unformatted')
          write(23) LP(iLayer)%NX, LP(iLayer)%NY
          do j=1,LP(iLayer)%NY
              write(23) (LP(iLayer)%Z(i,j),i=1,LP(iLayer)%NX)
          enddo
          close(23)
      enddo
      write(*,*)

      end subroutine readBathymetry



      recursive subroutine smoothBathymetry(GP, LP, iLayer, m, n)
      
      use VariableDefination
      implicit NONE
      type(GlobalParameters)  ::  GP
      type(LayerParameters)   ::  LP(100)
      integer*4  ::  iLayer, m, n, i, j

      if(m.le.GP%nRowBathymetryBoundary) then
          do i=m-1, 1, -1
              if(LP(iLayer)%Z(i,n).gt.GP%WaterDepthLimit) then
                  LP(iLayer)%Z(i,n) = 0.0d0
                  call smoothBathymetry(GP,LP,iLayer,i,n)
              endif
          enddo
      endif
      if(m.ge.LP(iLayer)%NX-GP%nRowBathymetryBoundary+1) then
          do i=m+1, LP(iLayer)%NX
              if(LP(iLayer)%Z(i,n).gt.GP%WaterDepthLimit) then
                  LP(iLayer)%Z(i,n) = 0.0d0
                  call smoothBathymetry(GP,LP,iLayer,i,n)
              endif
          enddo
      endif
      if(n.le.GP%nRowBathymetryBoundary) then
          do i=n-1, 1, -1
              if(LP(iLayer)%Z(m,i).gt.GP%WaterDepthLimit) then
                  LP(iLayer)%Z(m,i) = 0.0d0
                  call smoothBathymetry(GP,LP,iLayer,m,i)
              endif
          enddo
      endif
      if(n.ge.LP(iLayer)%NY-GP%nRowBathymetryBoundary+1) then
          do i=n+1, LP(iLayer)%NY
              if(LP(iLayer)%Z(m,i).gt.GP%WaterDepthLimit) then
                  LP(iLayer)%Z(m,i) = 0.0d0
                  call smoothBathymetry(GP,LP,iLayer,m,i)
              endif
          enddo
      endif
      if(m.ge.GP%nRowBathymetryBoundary+1) then
          do i=m-1, m-GP%nRowBathymetry, -1
              if(LP(iLayer)%Z(i,n).le.GP%WaterDepthLimit) then
                  do j=m-1, i+1, -1
                      if(LP(iLayer)%Z(j,n).gt.GP%WaterDepthLimit) then
                          LP(iLayer)%Z(j,n) = 0.0d0
                          call smoothBathymetry(GP,LP,iLayer,j,n)
                      endif
                  enddo
              endif
          enddo
      endif
      if(m.le.LP(iLayer)%NX-GP%nRowBathymetryBoundary) then
          do i=m+1, m+GP%nRowBathymetry
              if(LP(iLayer)%Z(i,n).le.GP%WaterDepthLimit) then
                  do j=m+1, i-1
                      if(LP(iLayer)%Z(j,n).gt.GP%WaterDepthLimit) then
                          LP(iLayer)%Z(j,n) = 0.0d0
                          call smoothBathymetry(GP,LP,iLayer,j,n)
                      endif
                  enddo
              endif
          enddo
      endif
      if(n.ge.GP%nRowBathymetryBoundary+1) then
          do i=n-1, n-GP%nRowBathymetry, -1
              if(LP(iLayer)%Z(m,i).le.GP%WaterDepthLimit) then
                  do j=n-1, i+1, -1
                      if(LP(iLayer)%Z(m,j).gt.GP%WaterDepthLimit) then
                          LP(iLayer)%Z(m,j) = 0.0d0
                          call smoothBathymetry(GP,LP,iLayer,m,j)
                      endif
                  enddo
              endif
          enddo
      endif
      if(n.le.LP(iLayer)%NY-GP%nRowBathymetryBoundary) then
          do i=n+1, n+GP%nRowBathymetry
              if(LP(iLayer)%Z(m,i).le.GP%WaterDepthLimit) then
                  do j=n+1, i-1
                      if(LP(iLayer)%Z(m,j).gt.GP%WaterDepthLimit) then
                          LP(iLayer)%Z(m,j) = 0.0d0
                          call smoothBathymetry(GP,LP,iLayer,m,j)
                      endif
                  enddo
              endif
          enddo
      endif

if(1.eq.0) then
      if(m.ge.3.and.LP(iLayer)%Z(m-1,n).gt.GP%WaterDepthLimit.and. &
          LP(iLayer)%Z(m-2,n).le.GP%WaterDepthLimit) then
          LP(iLayer)%Z(m-1,n)=0.5*(LP(iLayer)%Z(m,n)+LP(iLayer)%Z(m-2,n))
          call smoothBathymetry(GP,LP,iLayer,m-1,n)
      endif
      if(m.le.LP(iLayer)%NX-2.and.LP(iLayer)%Z(m+1,n).gt.GP%WaterDepthLimit.and. &
          LP(iLayer)%Z(m+2,n).le.GP%WaterDepthLimit) then
          LP(iLayer)%Z(m+1,n)=0.5*(LP(iLayer)%Z(m,n)+LP(iLayer)%Z(m+2,n))
          call smoothBathymetry(GP,LP,iLayer,m+1,n)
      endif
      if(n.ge.3.and.LP(iLayer)%Z(m,n-1).gt.GP%WaterDepthLimit.and. &
          LP(iLayer)%Z(m,n-2).le.GP%WaterDepthLimit) then
          LP(iLayer)%Z(m,n-1)=0.5*(LP(iLayer)%Z(m,n)+LP(iLayer)%Z(m,n-2))
          call smoothBathymetry(GP,LP,iLayer,m,n-1)
      endif
      if(n.le.LP(iLayer)%NY-2.and.LP(iLayer)%Z(m,n+1).gt.GP%WaterDepthLimit.and. &
          LP(iLayer)%Z(m,n+2).le.GP%WaterDepthLimit) then
          LP(iLayer)%Z(m,n+1)=0.5*(LP(iLayer)%Z(m,n)+LP(iLayer)%Z(m,n+2))
          call smoothBathymetry(GP,LP,iLayer,m,n+1)
      endif
endif

      end subroutine smoothBathymetry



      subroutine cflCheck(GP, LP)
      
      use VariableDefination
      implicit NONE
      type(GlobalParameters)  ::  GP
      type(LayerParameters)   ::  LP(100)
      integer*4               ::  iLayer, i, j
      real*8                  ::  zmax, dmin, cr

      write(*,*)
      write(*,*) 'check CFL and determine dt for each layer ...'
      zmax = LP(GP%TopLayer)%zmax
      if(GP%CoordinatesType.eq.0) then !spehrical coordinate
          dmin = MIN(COS(LP(GP%TopLayer)%ymin*GP%PI/180.0),COS(LP(GP%TopLayer)%ymax*GP%PI/180.0))*GP%R_EARTH
          dmin = MIN(GP%R_EARTH*LP(GP%TopLayer)%DX, dmin*LP(GP%TopLayer)%DY)*GP%PI/180.0
      else
          dmin = MIN(LP(GP%TopLayer)%DX, LP(GP%TopLayer)%DY)
      endif
      cr = GP%dt*SQRT(GP%GRAV*zmax)/dmin
      write(*,*)
      write(*,'(a,f6.3)',advance='no') '        CFL number:',cr
      if(GP%CoordinatesType.eq.0) then
          write(*,*) '  for linear non-dispersive equations ...'
      else if(GP%CoordinatesType.eq.1) then
          write(*,*) '  for nonlinear non-dispersive equations ...'
      endif
      if(cr.gt.0.5) then
          write(*,*) 'WARNING: CFL number is greater than 0.5 ...'
      endif
      write(*,*) '       determine dt for each layer according to cfl number:'
      write(*,*)
      do iLayer=1, GP%NumLayers
          if(iLayer.eq.GP%TopLayer) then
              LP(iLayer)%dt = GP%dt
              LP(iLayer)%nStepsPerTimeStep = 1
              LP(iLayer)%dtratio = 1.0d0
          else
              zmax = LP(iLayer)%zmax
              if(GP%CoordinatesType.eq.0) then !spehrical coordinate
                  dmin = MIN(COS(LP(iLayer)%ymin*GP%PI/180.0),COS(LP(iLayer)%ymax*GP%PI/180.0))*GP%R_EARTH
                  dmin = MIN(GP%R_EARTH*LP(iLayer)%DX, dmin*LP(iLayer)%DY)*GP%PI/180.0
              else
                  dmin = MIN(LP(iLayer)%DX, LP(iLayer)%DY)
              endif
              LP(iLayer)%dt = cr*dmin/SQRT(GP%GRAV*zmax)
              i = 1
              do
                  if(LP(iLayer)%dt*i.ge.GP%dt.or.ABS(LP(iLayer)%dt*i-GP%dt).lt.1e-6) then
                      LP(iLayer)%nStepsPerTimeStep = i
                      if(ABS(LP(iLayer)%dt*i-GP%dt).lt.1e-6) then
                          LP(iLayer)%dtratio = 1.0d0
                      else
                          LP(iLayer)%dtratio = (GP%dt-(i-1)*LP(iLayer)%dt)/LP(iLayer)%dt
                      endif
                      exit
                  endif
                  i = i + 1
              enddo
          endif
          write(*,'(a,a11,a,f7.4)') '        ',ADJUSTL(LP(iLayer)%BathymetryFileName),':  dt = ',LP(iLayer)%dt
      enddo
      write(*,*)

      end subroutine cflCheck



      subroutine readStations(GP, LP, SP)
      
      use VariableDefination
      implicit NONE
      type(GlobalParameters)    ::  GP
      type(LayerParameters)     ::  LP(100)
      type(StationParameters)   ::  SP(999)
      character(999)            ::  s
      integer*4  ::  i, i0, j, ios
      real*8     ::  x, y

      write(*,*) 'reading location of stations from file ', &
          TRIM(ADJUSTL(GP%StationFileName)),' ...'
      open(23,file=TRIM(ADJUSTL(GP%StationFileName)), &
          status='old',form='formatted',iostat=ios)
      i = 0; i0 = 0
      do
          read(23,'(a)',iostat=ios) s
          if(ios.ne.0) exit
          if(LEN(TRIM(ADJUSTL(s))).eq.0) exit
          i = i+1
          read(s,*) x,y
          if((GP%StartEastWest.eq.1).and.(x.lt.0.0)) x=x+360.0
          if((GP%StartEastWest.eq.2).and.(x.gt.180.0)) x=x-360.0
          if((x.ge.LP(GP%TopLayer)%xmin+LP(GP%TopLayer)%dx).and. &
              (x.le.LP(GP%TopLayer)%xmax-LP(GP%TopLayer)%dx).and. &
              (y.ge.LP(GP%TopLayer)%ymin+LP(GP%TopLayer)%dy).and. &
              (y.le.LP(GP%TopLayer)%ymax-LP(GP%TopLayer)%dy)) then
              i0 = i0 + 1
              SP(i0)%X = x; SP(i0)%Y = y
          else
              if(i0.eq.i-1)  write(*,'(a)',advance='no') &
                  '        Station(s) not in computing domain: #'
              write(*,'(i4)',advance='no') i
          endif
      enddo
      close(23)
      if(i0.ne.GP%NumStations) then
          write(*,*)
          GP%NumStations = i0
      endif
      do i=1, GP%NumStations
          SP(i)%nLayer = 0
          do j=1, GP%NumLayers
              if(SP(i)%X.ge.LP(j)%xmin.and.SP(i)%X.le.LP(j)%xmax.and.  &
                  SP(i)%Y.ge.LP(j)%ymin.and.SP(i)%Y.le.LP(j)%ymax) then
                  if(SP(i)%nLayer.eq.0) then
                      SP(i)%nLayer = j
                  elseif(LP(j)%Level.gt.LP(SP(i)%nLayer)%Level) then
                      SP(i)%nLayer = j
                  endif
              endif
          enddo
          write(*,'(a,i3,a,f9.4,a,f9.4,a,i3)')  &
              '        station ',i,':    x:',SP(i)%x,',    y:',SP(i)%y,',    layer:',SP(i)%nLayer
      enddo
      write(*,*)

      end subroutine readStations



      subroutine readFaultParameters(GP, FP)
      
      use VariableDefination
      implicit NONE
      type(GlobalParameters)    ::  GP
      type(FaultParameters)     ::  FP(999)
      integer*4                 ::  i, ios
      character(999)            ::  s

      write(*,*) 'reading fault parameters from file: ', &
          TRIM(ADJUSTL(GP%FaultParametersFileName))
      open(23,file=TRIM(ADJUSTL(GP%FaultParametersFileName)), &
          status='old',form='formatted',iostat=ios)
      read(23, '(3/)')
      i=1
      do
          read(23, '(a)', iostat=ios) s
          if(ios.ne.0) exit
          read(23, '(a)', iostat=ios) s
          if(ios.ne.0) exit
          read(23, '(a)', iostat=ios) s
          if(ios.ne.0) exit
          if(INDEX(s,'Parameters for Fault Segment').eq.0) exit
          read(23, '(a)', iostat=ios) s
          read(23, '(49x,f30.6)') FP(i)%T0
          if((FP(i)%T0).lt.0.0.or.FP(i)%T0.gt.GP%TotalTime) then
              FP(i)%NT = -1
          else
              FP(i)%NT = NINT(FP(i)%T0/GP%dt)
          endif
          read(23, '(49x,f30.6)') FP(i)%Depth
          read(23, '(49x,f30.6)') FP(i)%Length
          read(23, '(49x,f30.6)') FP(i)%Width
          read(23, '(49x,f30.6)') FP(i)%Slip
          read(23, '(49x,f30.6)') FP(i)%Rake
          read(23, '(49x,f30.6)') FP(i)%Strike
          read(23, '(49x,f30.6)') FP(i)%Dip
          read(23, '(49x,f30.6)') FP(i)%Y0
          read(23, '(49x,f30.6)') FP(i)%X0
          if(GP%StartEastWest.eq.1.and.FP(i)%X0.lt.0.0) FP(i)%X0 = FP(i)%X0+360.0
          if(GP%StartEastWest.eq.2.and.FP(i)%X0.gt.180.0) FP(i)%X0 = FP(i)%X0-360.0
          i=i+1
      enddo
      GP%NumFaults = i-1
      close(23)
      write(*,*)
      write(*,'(a,i4,a)') '        There are ',GP%NumFaults,'  fault segments.'
      write(*,*)
      if(GP%ComputeGreen.eq.1) GP%nCalculations = GP%NumFaults

      end subroutine readFaultParameters



      subroutine partitionDomain(GP, LP, SP)
      
      use VariableDefination
      implicit NONE
      type(GlobalParameters)  ::  GP
      type(LayerParameters)   ::  LP(100)
      type(StationParameters) ::  SP(999)
      integer*4  ::  nsizeTotal, nsize
      integer*4  ::  iLayerLevel, iLayer, pLayer, iSta, i, j, ii, jj
      integer*4  ::  npartx, nparty, ndiff, ntmp, nlength
      integer*4  ::  nstartx, nendx, nstarty, nendy
      integer*4  ::  istart, iend, jstart, jend
      integer*4  ::  ixy, nBlocks, iBlock, AllBlocks(99,2), IsNewBlock, IsMergeBlock
      integer*4  ::  iMarginx, iMarginy, MaxGrids, MaxGridsOld, MinGrids, AvgGrids
      integer*4  ::  MaxNX, MinNX, MaxNY, MinNY
      integer*4  ::  nCnt, iNode, iCnt, iBoundary, iHMN
      integer*4  ::  iNodeFrom, iNodeTo, HasBoundary, InDomain, IsExist, IsNewSendRecv
      real*8     ::  x, y

      write(*,*) 'dividing each computational domain into subdomains ...'
      nsizeTotal = GP%nsizeTotal
      do iLayerLevel=1, GP%NumLayerLevels
      do iLayer=1, GP%NumLayers
      if(LP(iLayer)%Level.eq.iLayerLevel) then

          if(GP%MinGridsPerNode.le.0) then
              nsize = nsizeTotal
          else
              nsize = MIN(nsizeTotal,LP(iLayer)%NX*LP(iLayer)%NY/GP%MinGridsPerNode)
              nsize = MAX(1,nsize)
          endif
          npartx = 1; nparty = nsize
          ndiff = ABS(LP(iLayer)%NX/npartx-LP(iLayer)%NY/nparty)
          do i=1, nsize
              if(MOD(nsize,i).eq.0) then
                  ntmp = ABS(LP(iLayer)%NX/i-LP(iLayer)%NY/(nsize/i))
                  if(ntmp.lt.ndiff) then
                      npartx = i; nparty = nsize/i; ndiff = ntmp;
                  endif
              endif
          enddo
          LP(iLayer)%nsize = nsize
          LP(iLayer)%npartx = npartx; LP(iLayer)%nparty = nparty

          LP(iLayer)%PartitionInfo = -1
          do i=0,nsize-1
              nlength = LP(iLayer)%NX/npartx
              if(MOD(i,npartx)+1.le.(MOD(LP(iLayer)%NX,npartx))) then
                  LP(iLayer)%PartitionInfo(i+1,1) = MOD(i,npartx)*(nlength+1)+1
                  LP(iLayer)%PartitionInfo(i+1,2) = LP(iLayer)%PartitionInfo(i+1,1) + nlength
              else
                  LP(iLayer)%PartitionInfo(i+1,1) = MOD(LP(iLayer)%NX,npartx)*(nlength+1)+ &
                      (MOD(i,npartx)-MOD(LP(iLayer)%NX,npartx))*nlength + 1
                  LP(iLayer)%PartitionInfo(i+1,2) = LP(iLayer)%PartitionInfo(i+1,1)+nlength-1
              endif
              nlength = LP(iLayer)%NY/nparty
              if(i/npartx+1.le.MOD(LP(iLayer)%NY,nparty)) then
                  LP(iLayer)%PartitionInfo(i+1,3) = (i/npartx)*(nlength+1) + 1
                  LP(iLayer)%PartitionInfo(i+1,4) = LP(iLayer)%PartitionInfo(i+1,3) + nlength
              else
                  LP(iLayer)%PartitionInfo(i+1,3) = MOD(LP(iLayer)%NY,nparty)*(nlength+1)+ &
                      (i/npartx-MOD(LP(iLayer)%NY,nparty))*nlength+1
                  LP(iLayer)%PartitionInfo(i+1,4) = LP(iLayer)%PartitionInfo(i+1,3)+nlength-1
              endif
          enddo

          ! To save communation time, child layers are contained in one parent layer compute node !
          ! Partition of parent layers have to be optimized !
          if(GP%ComputeDivisionOpt.eq.2) then
          MaxGridsOld = 0
          if(LP(iLayer)%Level.ne.GP%NumLayerLevels) then
          do
          iMarginx = 0; iMarginy = 0;
          MaxNX = 0; MinNX = LP(iLayer)%NX; MaxNY = 0; MinNY = LP(iLayer)%NY
          do ixy=1, 2
              nBlocks = 0; AllBlocks = -1;
              do i=1, GP%NumLayers
              if(LP(i)%Level.eq.iLayerLevel+1) then
                  IsNewBlock = 1
                  if(ixy.eq.1) then
                      istart = NINT((LP(i)%xmin-LP(iLayer)%xmin)/LP(iLayer)%dx)+1
                      iend   = NINT((LP(i)%xmax-LP(iLayer)%xmin)/LP(iLayer)%dx)+1
                  else
                      istart = NINT((LP(i)%ymin-LP(iLayer)%ymin)/LP(iLayer)%dy)+1
                      iend   = NINT((LP(i)%ymax-LP(iLayer)%ymin)/LP(iLayer)%dy)+1
                  endif
                  do iBlock=1, nBlocks
                      if(istart.ge.AllBlocks(iBlock,1).and.istart.le.AllBlocks(iBlock,2)) then
                          IsNewBlock = 0; AllBlocks(iBlock,2) = MAX(AllBlocks(iBlock,2),iend); exit
                      elseif(iend.ge.AllBlocks(iBlock,1).and.iend.le.AllBlocks(iBlock,2)) then
                          IsNewBlock = 0; AllBlocks(iBlock,1) = MIN(AllBlocks(iBlock,1),istart); exit
                      elseif(istart.lt.AllBlocks(iBlock,1).and.iend.gt.AllBlocks(iBlock,2)) then
                          IsNewBlock = 0; AllBlocks(iBlock,1) = istart
                          AllBlocks(iBlock,2) = iend; exit
                      endif
                  enddo
                  if(IsNewBlock.eq.1) then
                      nBlocks = nBlocks+1; AllBlocks(nBlocks,1) = istart; ALlBlocks(nBlocks,2) = iend
                  endif
              endif
              enddo
              do i=1, nBlocks-1
                  do j=i+1, nBlocks
                      if(AllBlocks(i,1).gt.AllBlocks(j,1)) then
                          ntmp = AllBlocks(i,1); AllBlocks(i,1) = AllBlocks(j,1); AllBlocks(j,1) = ntmp
                          ntmp = AllBlocks(i,2); AllBlocks(i,2) = AllBlocks(j,2); AllBlocks(j,2) = ntmp
                      endif
                  enddo
              enddo
              IsMergeBlock = 0
              do
                  do iBlock=1, nBlocks-1
                      if(AllBlocks(iBlock,2).ge.AllBlocks(iBlock+1,1)) then
                          AllBlocks(iBlock,2) = MAX(AllBlocks(iBlock,2),AllBlocks(iBlock+1,2))
                          do i=iBlock+1, nBlocks-1
                              AllBlocks(i,1) = AllBlocks(i+1,1); AllBlocks(i,2) = AllBlocks(i+1,2)
                          enddo
                          IsMergeBlock = 1; exit
                      endif
                  enddo
                  if(IsMergeBlock.eq.0) exit
              enddo
              if(ixy.eq.1) then
                  nlength = LP(iLayer)%NX/npartx
              else
                  nlength = LP(iLayer)%NY/nparty
              endif
              if((ixy.eq.1.and.npartx.gt.nBlocks).or.(ixy.eq.2.and.nparty.gt.nBlocks)) then
                  nlength = 0
                  do iBlock=1, nBlocks
                      nlength = nlength+AllBlocks(iBlock,2)-AllBlocks(iBlock,1)+1
                  enddo
                  if(ixy.eq.1) then
                      nlength = (LP(iLayer)%NX-nlength)/(npartx-nBlocks)
                  else
                      nlength = (LP(iLayer)%NY-nlength)/(nparty-nBlocks)
                  endif
              endif
              do i=0, nsize-1
                  istart = 0; iend = 0
                  if(ixy.eq.1) then
                      if(MOD(i,npartx).eq.0) then
                          istart = 1
                      else
                          istart = LP(iLayer)%PartitionInfo(i,2)+1
                      endif
                  else
                      if(i/npartx.eq.0) then
                          istart = 1
                      else
                          istart = LP(iLayer)%PartitionInfo(i-npartx+1,4)+1
                      endif
                  endif
                  iend = istart + nlength
                  do iBlock=1, nBlocks
                      if(iend.ge.AllBlocks(iBlock,1).and.iend.le.AllBlocks(iBlock,2)) then
                          if(iend-AllBlocks(iBlock,1).lt.AllBlocks(iBlock,2)-iend.and. &
                              AllBlocks(iBlock,1)-3.ge.istart+3) then
                              iend = AllBlocks(iBlock,1)-3;
                          else
                              iend = AllBlocks(iBlock,2)+3;
                              if(ixy.eq.1) iMarginx = MAX(iMarginx, AllBlocks(iBlock,1)-istart)
                              if(ixy.eq.2) iMarginy = MAX(iMarginy, AllBlocks(iBlock,1)-istart)
                          endif
                      endif
                  enddo
                  if(ixy.eq.1.and.MOD(i,npartx).eq.npartx-1) then
                      iend = LP(iLayer)%NX
                      if(istart.le.AllBlocks(nBlocks,2)) iMarginx = MAX(iMarginx, iend-AllBlocks(nBlocks,2));
                  elseif(ixy.eq.2.and.i/npartx.eq.nparty-1) then
                      iend = LP(iLayer)%NY
                      if(istart.le.AllBlocks(nBlocks,2)) iMarginy = MAX(iMarginy, iend-AllBlocks(nBlocks,2));
                  endif
                  if(ixy.eq.1) then
                      LP(iLayer)%PartitionInfo(i+1,1) = istart; LP(iLayer)%PartitionInfo(i+1,2) = iend
                  else
                      LP(iLayer)%PartitionInfo(i+1,3) = istart; LP(iLayer)%PartitionInfo(i+1,4) = iend
                  endif
              enddo
          enddo !do loop: ixy
          do i=0, nsize-1
              MaxNX = MAX(MaxNX, (LP(iLayer)%PartitionInfo(i+1,2)-LP(iLayer)%PartitionInfo(i+1,1)+1)+1)
              MinNX = MIN(MinNX, (LP(iLayer)%PartitionInfo(i+1,2)-LP(iLayer)%PartitionInfo(i+1,1)+1)+1)
              MaxNY = MAX(MaxNY, (LP(iLayer)%PartitionInfo(i+1,4)-LP(iLayer)%PartitionInfo(i+1,3)+1)+1)
              MinNY = MIN(MinNY, (LP(iLayer)%PartitionInfo(i+1,4)-LP(iLayer)%PartitionInfo(i+1,3)+1)+1)
          enddo
          MaxGrids = MAX(MaxGrids,MaxNX*MaxNY)
          if((MaxGridsOld.eq.0.or.MaxGridsOld-MaxGrids.gt.0.2*MaxGridsOld).and. &
              MaxGrids.gt.3*GP%MinGridsPerNode.and.iMarginx.lt.iMarginy.and.nsize+npartx.lt.nsizeTotal) then
              nparty = nparty+1; nsize = npartx*nparty; LP(iLayer)%nsize = nsize
              LP(iLayer)%npartx = npartx; LP(iLayer)%nparty = nparty
          elseif((MaxGridsOld.eq.0.or.MaxGridsOld-MaxGrids.gt.0.2*MaxGridsOld).and. &
              MaxGrids.gt.3*GP%MinGridsPerNode.and.iMarginx.gt.iMarginy.and.nsize+nparty.lt.nsizeTotal) then
              npartx = npartx+1; nsize = npartx*nparty; LP(iLayer)%nsize = nsize
              LP(iLayer)%npartx = npartx; LP(iLayer)%nparty = nparty
          elseif(MinNX.lt.0.1*MaxNX) then
              npartx = npartx-1; nsize = npartx*nparty; LP(iLayer)%nsize = nsize
              LP(iLayer)%npartx = npartx; LP(iLayer)%nparty = nparty
          elseif(MinNY.lt.0.1*MaxNY) then
              nparty = nparty-1; nsize = npartx*nparty; LP(iLayer)%nsize = nsize
              LP(iLayer)%npartx = npartx; LP(iLayer)%nparty = nparty
          else
              exit
          endif
          MaxGridsOld = MaxGrids
          enddo !do loop: end looping for optimizing partition
          endif !if: layer is a parent layer
          endif !if: if child layer on ONE node

      endif
      enddo
      enddo

      !////// calculate MaxNX, MinNX, MaxNY, MinNY //////!
      do iLayer=1, GP%NumLayers
          LP(iLayer)%MaxNX = 0; LP(iLayer)%MinNX = LP(iLayer)%NX;
          LP(iLayer)%MaxNY = 0; LP(iLayer)%MinNY = LP(iLayer)%NY;
          do i=1, LP(iLayer)%nsize
              LP(iLayer)%MaxNX = MAX(LP(iLayer)%MaxNX,LP(iLayer)%PartitionInfo(i,2)-LP(iLayer)%PartitionInfo(i,1)+1)
              LP(iLayer)%MinNX = MIN(LP(iLayer)%MinNX,LP(iLayer)%PartitionInfo(i,2)-LP(iLayer)%PartitionInfo(i,1)+1)
              LP(iLayer)%MaxNY = MAX(LP(iLayer)%MaxNY,LP(iLayer)%PartitionInfo(i,4)-LP(iLayer)%PartitionInfo(i,3)+1)
              LP(iLayer)%MinNY = MIN(LP(iLayer)%MinNY,LP(iLayer)%PartitionInfo(i,4)-LP(iLayer)%PartitionInfo(i,3)+1)
          enddo
      enddo
      GP%MaxNX = LP(1)%MaxNX; GP%MinNX = LP(1)%MinNX;
      GP%MaxNY = LP(1)%MaxNY; GP%MinNY = LP(1)%MinNY;
      do iLayer=2, GP%NumLayers
          GP%MaxNX = MAX(GP%MaxNX, LP(iLayer)%MaxNX); GP%MinNX = MIN(GP%MinNX, LP(iLayer)%MinNX);
          GP%MaxNY = MAX(GP%MaxNY, LP(iLayer)%MaxNY); GP%MinNY = MIN(GP%MinNY, LP(iLayer)%MinNY);
      enddo

      !****** Calculate Boundary Send Receive Table for Each Layer ******!
      do iLayer=1, GP%NumLayers
          nsize = LP(iLayer)%nsize; npartx = LP(iLayer)%npartx; nparty = LP(iLayer)%nparty
          LP(iLayer)%BoundarySendRecvCount = 0; LP(iLayer)%BoundarySendRecv = -1; nCnt = 0
          do i=0, nsize-1
              nstartx = LP(iLayer)%PartitionInfo(i+1,1); nendx = LP(iLayer)%PartitionInfo(i+1,2)
              nstarty = LP(iLayer)%PartitionInfo(i+1,3); nendy = LP(iLayer)%PartitionInfo(i+1,4)
              !////// send H/M/N on left boundary of this compute domain //////!
              if(MOD(i,npartx).ne.0) then
                  nCnt = nCnt+1
                  LP(iLayer)%BoundarySendRecv(nCnt, 1) = i
                  LP(iLayer)%BoundarySendRecv(nCnt, 2) = i-1
                  LP(iLayer)%BoundarySendRecv(nCnt, 3) = nstartx
                  LP(iLayer)%BoundarySendRecv(nCnt, 4) = nstartx+GP%nRowBoundary-1
                  LP(iLayer)%BoundarySendRecv(nCnt, 5) = nstarty
                  LP(iLayer)%BoundarySendRecv(nCnt, 6) = nendy
                  LP(iLayer)%BoundarySendRecv(nCnt, 7) = nstartx
                  LP(iLayer)%BoundarySendRecv(nCnt, 8) = nstartx+GP%nRowBoundaryFlux-1
                  LP(iLayer)%BoundarySendRecv(nCnt, 9) = nstarty
                  LP(iLayer)%BoundarySendRecv(nCnt,10) = nendy
                  LP(iLayer)%BoundarySendRecv(nCnt,11) = nstartx
                  LP(iLayer)%BoundarySendRecv(nCnt,12) = nstartx+GP%nRowBoundaryFlux-1
                  LP(iLayer)%BoundarySendRecv(nCnt,13) = nstarty
                  LP(iLayer)%BoundarySendRecv(nCnt,14) = MIN(nendy,LP(iLayer)%NY-1)
              endif
              !////// send H/M/N on right boundary of this compute domain //////!
              if(MOD(i,npartx).ne.npartx-1) then
                  nCnt = nCnt+1
                  LP(iLayer)%BoundarySendRecv(nCnt, 1) = i
                  LP(iLayer)%BoundarySendRecv(nCnt, 2) = i+1
                  LP(iLayer)%BoundarySendRecv(nCnt, 3) = nendx-GP%nRowBoundary+1
                  LP(iLayer)%BoundarySendRecv(nCnt, 4) = nendx
                  LP(iLayer)%BoundarySendRecv(nCnt, 5) = nstarty
                  LP(iLayer)%BoundarySendRecv(nCnt, 6) = nendy
                  LP(iLayer)%BoundarySendRecv(nCnt, 7) = nendx-GP%nRowBoundaryFlux+1
                  LP(iLayer)%BoundarySendRecv(nCnt, 8) = nendx
                  LP(iLayer)%BoundarySendRecv(nCnt, 9) = nstarty
                  LP(iLayer)%BoundarySendRecv(nCnt,10) = nendy
                  LP(iLayer)%BoundarySendRecv(nCnt,11) = nendx-GP%nRowBoundaryFlux+1
                  LP(iLayer)%BoundarySendRecv(nCnt,12) = nendx
                  LP(iLayer)%BoundarySendRecv(nCnt,13) = nstarty
                  LP(iLayer)%BoundarySendRecv(nCnt,14) = MIN(nendy,LP(iLayer)%NY-1)
              endif
              !////// send H/M/N on bottom boundary of this compute domain //////!
              if(i/npartx.ne.0) then
                  nCnt = nCnt+1
                  LP(iLayer)%BoundarySendRecv(nCnt, 1) = i
                  LP(iLayer)%BoundarySendRecv(nCnt, 2) = i-npartx
                  LP(iLayer)%BoundarySendRecv(nCnt, 3) = nstartx
                  LP(iLayer)%BoundarySendRecv(nCnt, 4) = nendx
                  LP(iLayer)%BoundarySendRecv(nCnt, 5) = nstarty
                  LP(iLayer)%BoundarySendRecv(nCnt, 6) = nstarty+GP%nRowBoundary-1
                  LP(iLayer)%BoundarySendRecv(nCnt, 7) = nstartx
                  LP(iLayer)%BoundarySendRecv(nCnt, 8) = MIN(nendx,LP(iLayer)%NX-1)
                  LP(iLayer)%BoundarySendRecv(nCnt, 9) = nstarty
                  LP(iLayer)%BoundarySendRecv(nCnt,10) = nstarty+GP%nRowBoundaryFlux-1
                  LP(iLayer)%BoundarySendRecv(nCnt,11) = nstartx
                  LP(iLayer)%BoundarySendRecv(nCnt,12) = nendx
                  LP(iLayer)%BoundarySendRecv(nCnt,13) = nstarty
                  LP(iLayer)%BoundarySendRecv(nCnt,14) = nstarty+GP%nRowBoundaryFlux-1
              endif
              !////// send H/M/N on top boundary of this compute domain //////!
              if(i/npartx.ne.nparty-1) then
                  nCnt = nCnt+1
                  LP(iLayer)%BoundarySendRecv(nCnt, 1) = i
                  LP(iLayer)%BoundarySendRecv(nCnt, 2) = i+npartx
                  LP(iLayer)%BoundarySendRecv(nCnt, 3) = nstartx
                  LP(iLayer)%BoundarySendRecv(nCnt, 4) = nendx
                  LP(iLayer)%BoundarySendRecv(nCnt, 5) = nendy-GP%nRowBoundary+1
                  LP(iLayer)%BoundarySendRecv(nCnt, 6) = nendy
                  LP(iLayer)%BoundarySendRecv(nCnt, 7) = nstartx
                  LP(iLayer)%BoundarySendRecv(nCnt, 8) = MIN(nendx,LP(iLayer)%NX-1)
                  LP(iLayer)%BoundarySendRecv(nCnt, 9) = nendy-GP%nRowBoundaryFlux+1
                  LP(iLayer)%BoundarySendRecv(nCnt,10) = nendy
                  LP(iLayer)%BoundarySendRecv(nCnt,11) = nstartx
                  LP(iLayer)%BoundarySendRecv(nCnt,12) = nendx
                  LP(iLayer)%BoundarySendRecv(nCnt,13) = nendy-GP%nRowBoundaryFlux+1
                  LP(iLayer)%BoundarySendRecv(nCnt,14) = nendy
              endif
          enddo
          LP(iLayer)%BoundarySendRecvCount = nCnt
      enddo

      !****** Calculate Parent to Child Send Receive Table: Boundaries on each node ******!
      do iLayer=1, GP%NumLayers
      if(LP(iLayer)%Level.gt.1) then
          pLayer = LP(iLayer)%Parent
          LP(iLayer)%ParentToChildSendRecvCount = 0; LP(iLayer)%ParentToChildSendRecv = -1; nCnt = 0
          do iNodeTo=0, LP(iLayer)%nsize-1
              HasBoundary = 0
              nstartx = LP(iLayer)%PartitionInfo(iNodeTo+1,1)
              nendx   = LP(iLayer)%PartitionInfo(iNodeTo+1,2)
              nstarty = LP(iLayer)%PartitionInfo(iNodeTo+1,3)
              nendy   = LP(iLayer)%PartitionInfo(iNodeTo+1,4)
              if(nstartx.eq.1.or.nendx.eq.LP(iLayer)%NX.or.nstarty.eq.1.or.nendy.eq.LP(iLayer)%NY) then
                  HasBoundary = 1
              endif
              if(HasBoundary.eq.1) then
              do iBoundary=1, 4
                  istart = 0; iend = -1; jstart = 0; jend = -1
                  if(nstartx.eq.1.and.iBoundary.eq.1) then
                      istart = 1; iend = istart;
                      jstart = MAX(nstarty-GP%nRowBoundary,1); jend = MIN(nendy+GP%nRowBoundary,LP(iLayer)%NY)
                  elseif(nendx.eq.LP(iLayer)%NX.and.iBoundary.eq.2) then
                      istart = LP(iLayer)%NX; iend = istart;
                      jstart = MAX(nstarty-GP%nRowBoundary,1); jend = MIN(nendy+GP%nRowBoundary,LP(iLayer)%NY)
                  elseif(nstarty.eq.1.and.iBoundary.eq.3) then
                      istart = MAX(nstartx-GP%nRowBoundary,1); iend = MIN(nendx+GP%nRowBoundary,LP(iLayer)%NX)
                      jstart = 1; jend = jstart;
                  elseif(nendy.eq.LP(iLayer)%NY.and.iBoundary.eq.4) then
                      istart = MAX(nstartx-GP%nRowBoundary,1); iend = MIN(nendx+GP%nRowBoundary,LP(iLayer)%NX)
                      jstart = LP(iLayer)%NY; jend = jstart;
                  endif
                  do j=jstart, jend
                  do i=istart, iend
                      iNodeFrom = -1; x = LP(iLayer)%X(i); y = LP(iLayer)%Y(j)
                      ii = FLOOR((x-LP(pLayer)%X(1))/LP(pLayer)%dx)+1
                      jj = FLOOR((y-LP(pLayer)%Y(1))/LP(pLayer)%dy)+1
                      do iNode=0, LP(pLayer)%nsize-1
                          if(ii.ge.LP(pLayer)%PartitionInfo(iNode+1,1).and. &
                              ii.le.LP(pLayer)%PartitionInfo(iNode+1,2).and. &
                              jj.ge.LP(pLayer)%PartitionInfo(iNode+1,3).and. &
                              jj.le.LP(pLayer)%PartitionInfo(iNode+1,4)) then
                              iNodeFrom = iNode; exit
                          endif
                      enddo
                      IsNewSendRecv = 0
                      if(nCnt.eq.0) then
                          IsNewSendRecv = 1
                      elseif(iNodeFrom.ne.LP(iLayer)%ParentToChildSendRecv(nCnt,1).or. &
                          iNodeTo.ne.LP(iLayer)%ParentToChildSendRecv(nCnt,2).or. &
                          iBoundary.ne.LP(iLayer)%ParentToChildSendRecv(nCnt,3)) then
                          IsNewSendRecv = 1
                      endif
                      if(IsNewSendRecv.eq.1) then
                          nCnt = nCnt + 1
                          LP(iLayer)%ParentToChildSendRecv(nCnt,1) = iNodeFrom
                          LP(iLayer)%ParentToChildSendRecv(nCnt,2) = iNodeTo
                          LP(iLayer)%ParentToChildSendRecv(nCnt,3) = iBoundary
                          if(iBoundary.eq.1.or.iBoundary.eq.2) then
                              LP(iLayer)%ParentToChildSendRecv(nCnt,4) = j
                              LP(iLayer)%ParentToChildSendRecv(nCnt,5) = j
                          else
                              LP(iLayer)%ParentToChildSendRecv(nCnt,4) = i
                              LP(iLayer)%ParentToChildSendRecv(nCnt,5) = i
                          endif
                      else
                          if(iBoundary.eq.1.or.iBoundary.eq.2) then
                              LP(iLayer)%ParentToChildSendRecv(nCnt,5) = j
                          else
                              LP(iLayer)%ParentToChildSendRecv(nCnt,5) = i
                          endif
                      endif
                  enddo
                  enddo
              enddo !iBoundary=1,4
              endif !HasBoundary
          enddo !iNodeTo, all nodes with boundaries
          LP(iLayer)%ParentToChildSendRecvCount = nCnt
      endif
      enddo

!iLayer=9
!write(*,*) 'iLayer,ParentToChildCount====',iLayer,LP(iLayer)%ParentToChildSendRecvCount
!do iBoundary=1,4
!do iCnt=1, LP(iLayer)%ParentToChildSendRecvCount
!if(LP(iLayer)%ParentToChildSendRecv(iCnt,3).eq.iBoundary) then
!write(*,'(5i7)') (LP(iLayer)%ParentToChildSendRecv(iCnt,i),i=1,5)
!endif
!enddo
!enddo

      !****** Calculate Child to Parent Send Receive Table: blocks on child layer ******!
      do iLayer=1, GP%NumLayers
      if(LP(iLayer)%Level.gt.1) then
          pLayer = LP(iLayer)%Parent
          LP(iLayer)%ChildToParentSendRecvCount = 0; LP(iLayer)%ChildToParentSendRecv = -1; nCnt = 0
          do iNodeTo=0, LP(pLayer)%nsize-1
              nstartx = LP(pLayer)%PartitionInfo(iNodeTo+1,1)
              nendx   = LP(pLayer)%PartitionInfo(iNodeTo+1,2)
              nstarty = LP(pLayer)%PartitionInfo(iNodeTo+1,3)
              nendy   = LP(pLayer)%PartitionInfo(iNodeTo+1,4)
              do iHMN=1, 3 !H/M/N
                  if(iHMN.eq.1) then
                      istart = MAX(nstartx-GP%nRowBoundary,1); iend = MIN(nendx+GP%nRowBoundary,LP(pLayer)%NX)
                      jstart = MAX(nstarty-GP%nRowBoundary,1); jend = MIN(nendy+GP%nRowBoundary,LP(pLayer)%NY)
                  else
                      istart = MAX(nstartx-GP%nRowBoundaryFlux,1); iend = MIN(nendx+GP%nRowBoundaryFlux,LP(pLayer)%NX)
                      jstart = MAX(nstarty-GP%nRowBoundaryFlux,1); jend = MIN(nendy+GP%nRowBoundaryFlux,LP(pLayer)%NY)
                  endif
                  do j=jstart, jend
                  do i=istart, iend
                      iNodeFrom = -1; InDomain = 0
                      x = LP(pLayer)%X(i); y = LP(pLayer)%Y(j)
                      if(x.ge.LP(iLayer)%xmin.and.x.le.LP(iLayer)%xmax.and. &
                          y.ge.LP(iLayer)%ymin.and.y.le.LP(iLayer)%ymax.and. &
                          x+0.5*LP(pLayer)%dx.ge.LP(iLayer)%xmin+0.5*LP(iLayer)%dx.and.&
                          x+0.5*LP(pLayer)%dx.le.LP(iLayer)%xmax-0.5*LP(iLayer)%dx.and.&
                          y+0.5*LP(pLayer)%dy.ge.LP(iLayer)%ymin+0.5*LP(iLayer)%dy.and.&
                          y+0.5*LP(pLayer)%dy.le.LP(iLayer)%ymax-0.5*LP(iLayer)%dy) then
                          InDomain = 1
                      endif
                      if(iHMN.eq.2) x = x+0.5*LP(pLayer)%dx
                      if(iHMN.eq.3) y = y+0.5*LP(pLayer)%dy
if(1.eq.0) then
                      if(iHMN.eq.1.and.x.ge.LP(iLayer)%xmin.and.x.le.LP(iLayer)%xmax.and. &
                          y.ge.LP(iLayer)%ymin.and.y.le.LP(iLayer)%ymax) then
                          InDomain = 1
                      elseif(iHMN.eq.2.and. &
                          x.ge.LP(iLayer)%xmin+0.5*LP(iLayer)%dx.and.x.le.LP(iLayer)%xmax-0.5*LP(iLayer)%dx.and. &
                          y.ge.LP(iLayer)%ymin.and.y.le.LP(iLayer)%ymax) then
                          InDomain = 1
                      elseif(iHMN.eq.3.and. &
                          x.ge.LP(iLayer)%xmin.and.x.le.LP(iLayer)%xmax.and. &
                          y.ge.LP(iLayer)%ymin+0.5*LP(iLayer)%dy.and.y.le.LP(iLayer)%ymax-0.5*LP(iLayer)%dy) then
                          InDomain = 1
                      endif
endif
                      if(InDomain.eq.1) then
                      ii = FLOOR((x-LP(iLayer)%X(1))/LP(iLayer)%dx)+1
                      jj = FLOOR((y-LP(iLayer)%Y(1))/LP(iLayer)%dy)+1
                      do iNode=0, LP(iLayer)%nsize-1
                          if(ii.ge.LP(iLayer)%PartitionInfo(iNode+1,1).and.ii.le.LP(iLayer)%PartitionInfo(iNode+1,2).and. &
                              jj.ge.LP(iLayer)%PartitionInfo(iNode+1,3).and.jj.le.LP(iLayer)%PartitionInfo(iNode+1,4)) then
                              iNodeFrom = iNode; exit
                          endif
                      enddo
                      IsExist = 0
                      do iCnt=1, nCnt
                          if(iNodeFrom.eq.LP(iLayer)%ChildToParentSendRecv(iCnt,1).and. &
                              iNodeTo.eq.LP(iLayer)%ChildToParentSendRecv(iCnt,2).and. &
                              iHMN.eq.LP(iLayer)%ChildToParentSendRecv(iCnt,3)) then
                              IsExist = iCnt; exit
                          endif
                      enddo
                      if(IsExist.eq.0) then
                          nCnt = nCnt + 1
                          LP(iLayer)%ChildToParentSendRecv(nCnt,1) = iNodeFrom
                          LP(iLayer)%ChildToParentSendRecv(nCnt,2) = iNodeTo
                          LP(iLayer)%ChildToParentSendRecv(nCnt,3) = iHMN
                          LP(iLayer)%ChildToParentSendRecv(nCnt,4) = i
                          LP(iLayer)%ChildToParentSendRecv(nCnt,5) = i
                          LP(iLayer)%ChildToParentSendRecv(nCnt,6) = j
                          LP(iLayer)%ChildToParentSendRecv(nCnt,7) = j
                      else
                          LP(iLayer)%ChildToParentSendRecv(IsExist,5) = MAX(i,LP(iLayer)%ChildToParentSendRecv(IsExist,5))
                          LP(iLayer)%ChildToParentSendRecv(IsExist,7) = MAX(j,LP(iLayer)%ChildToParentSendRecv(IsExist,7))
                      endif
                      endif !If HasBoundary=1: parent point(x,y) in a child domain
                  enddo
                  enddo
              enddo !iBoudanry for H/M/N
          enddo !iNodeTo for all parent nodes
          LP(iLayer)%ChildToParentSendRecvCount = nCnt
      endif
      enddo

!iLayer=9
!write(*,*) 'iLayer,ChildToParentCount====',iLayer,LP(iLayer)%ChildToParentSendRecvCount
!do iCnt=1, LP(iLayer)%ChildToParentSendRecvCount
!write(*,'(7i7)') (LP(iLayer)%ChildToParentSendRecv(iCnt,i),i=1,7)
!enddo

      !////// Calculate on which node is a station located //////!
      do iSta=1, GP%NumStations
          iLayer  = SP(iSta)%nLayer; SP(iSta)%nNode = -1
          ii = FLOOR((SP(iSta)%X-LP(iLayer)%X(1))/LP(iLayer)%dx)+1
          jj = FLOOR((SP(iSta)%Y-LP(iLayer)%Y(1))/LP(iLayer)%dy)+1
          do i=0, LP(iLayer)%nsize-1
              nstartx = LP(iLayer)%PartitionInfo(i+1,1)
              nendx   = LP(iLayer)%PartitionInfo(i+1,2)
              nstarty = LP(iLayer)%PartitionInfo(i+1,3)
              nendy   = LP(iLayer)%PartitionInfo(i+1,4)
              if(ii.ge.nstartx.and.ii.le.nendx.and.jj.ge.nstarty.and.jj.le.nendy) then
                  SP(iSta)%nNode = i; exit
              endif
          enddo
      enddo

      !////// Display partition infomation //////!
      do iLayer=1, GP%NumLayers
          write(*,*)
          write(*,'(a,a11,a)',advance='no') '        ',ADJUSTL(LP(iLayer)%BathymetryFileName),': '
          write(*,'(a,i5,a,i5,a,i10,a,i4,a)') 'NX, NY, NX*NY = ', LP(iLayer)%NX, ',', &
              LP(iLayer)%NY, ',', LP(iLayer)%NX*LP(iLayer)%NY, ',  use', LP(iLayer)%nsize,' nodes...'
          MaxGrids = LP(iLayer)%MaxNX*LP(iLayer)%MaxNY
          MinGrids = LP(iLayer)%MinNX*LP(iLayer)%MinNY
          AvgGrids = NINT((LP(iLayer)%NX*LP(iLayer)%NY*1.0)/(LP(iLayer)%npartx*LP(iLayer)%nparty))
          write(*,'(a,i7,a,i7,a,i7,a,f6.2,a,f6.2,a,f6.2)')           &
              '        Max, Min, Average grids on compute node = ',  &
              MaxGrids, ',', MinGrids, ',', AvgGrids, ',  ',         &
              MaxGrids*1.0/AvgGrids, ',', MinGrids*1.0/AvgGrids, ',', 1.0d0
          !write(*,'(a,a)') '        ComputingNode    nstartx     nendx    ',&
          !    'nstarty     nendy    TotalGridPoints'
          !do i=0, LP(iLayer)%nsize-1
          !    write(*,'(i16,i15,i11,i11,i10,i19)') i,&
          !        LP(iLayer)%PartitionInfo(i+1,1),LP(iLayer)%PartitionInfo(i+1,2), &
          !        LP(iLayer)%PartitionInfo(i+1,3),LP(iLayer)%PartitionInfo(i+1,4), &
          !        (LP(iLayer)%PartitionInfo(i+1,2)-LP(iLayer)%PartitionInfo(i+1,1)+1)* &
          !        (LP(iLayer)%PartitionInfo(i+1,4)-LP(iLayer)%PartitionInfo(i+1,3)+1)
          !enddo
          if(LP(iLayer)%nsize.ne.nsizeTotal) then
              write(*,'(a,a,a,i4,a,i4,a)') &
                  'WARNING: ',TRIM(ADJUSTL(LP(iLayer)%BathymetryFileName)), &
                  ' uses ',LP(iLayer)%nsize,'  out of ',nsizeTotal,'  computing nodes.'
          endif
      enddo
      write(*,*)
      open(99,file='PartitionInfo.dat',form='formatted')
      do iLayer=1, GP%NumLayers
          do i=0, LP(iLayer)%nsize-1
              write(99,'(8i7)') LP(iLayer)%nsize,LP(iLayer)%npartx,LP(iLayer)%nparty, &
                  i,(LP(iLayer)%PartitionInfo(i+1,j),j=1,4)
          enddo
      enddo
      close(99)

      end subroutine partitionDomain



      subroutine bcastCommonInfo(GP, LP, SP, FP)
      
      use VariableDefination
      implicit NONE
      include "mpif.h"
      type(GlobalParameters)   ::  GP
      type(LayerParameters)    ::  LP(100)
      type(StationParameters)  ::  SP(999)
      type(FaultParameters)    ::  FP(999)
      integer*4  ::  irank, nsize, master
      integer*4  ::  ierror, istatus(MPI_STATUS_SIZE)
      integer*4  ::  iLayer, iSta, iFault, LocalData(999*4*14), i, j

      irank = GP%irank; nsize = GP%nsizeTotal; master = GP%master
      if(irank.eq.master) then
          write(*,*) 'broadcasting common info from master node ...'
          write(*,*)
      endif

      call MPI_BCAST(GP%Version,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(GP%CoordinatesType,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(GP%GoverningEquationType,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(GP%WaterDepthLimit,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(GP%TotalTime,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(GP%dt,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(GP%InitialConditionType,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(GP%BoundaryConditionType,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(GP%DTSaveData,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(GP%DTSaveSTART,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(GP%SaveFlux,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(GP%ComputeGreen,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(GP%MinGridsPerNode,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(GP%FeedbackToParentLayer,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)

      call MPI_BCAST(GP%nCalculations,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(GP%TotalTimeSteps,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(GP%NDTSaveData,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(GP%NumLayers,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(GP%NumLayerLevels,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(GP%TopLayer,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
      !**** InitialConditionFileName   not broadcasted
      !**** InitialConditionFileFormat not broadcasted
      call MPI_BCAST(GP%StartEastWest,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(GP%NumStations,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(GP%NumFaults,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)

      call MPI_BCAST(GP%MaxNX,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(GP%MinNX,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(GP%MaxNY,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(GP%MinNY,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)

      do iLayer = 1,GP%NumLayers
          call MPI_BCAST(LP(iLayer)%Level,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
          call MPI_BCAST(LP(iLayer)%Parent,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
          !**** BathymetryFileName    not broadcasted
          !**** BathymetryFileFormat  not broadcasted
          call MPI_BCAST(LP(iLayer)%xmin,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
          call MPI_BCAST(LP(iLayer)%dx,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
          call MPI_BCAST(LP(iLayer)%xmax,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
          call MPI_BCAST(LP(iLayer)%ymin,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
          call MPI_BCAST(LP(iLayer)%dy,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
          call MPI_BCAST(LP(iLayer)%ymax,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
          call MPI_BCAST(LP(iLayer)%NX,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
          call MPI_BCAST(LP(iLayer)%NY,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
          !**** X, Y, Z broadcasted after memory allocated
          call MPI_BCAST(LP(iLayer)%zmin,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
          call MPI_BCAST(LP(iLayer)%zmax,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
          call MPI_BCAST(LP(iLayer)%dt,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
          call MPI_BCAST(LP(iLayer)%dtratio,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
          call MPI_BCAST(LP(iLayer)%nStepsPerTimeStep,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
          call MPI_BCAST(LP(iLayer)%nsize,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
          call MPI_BCAST(LP(iLayer)%npartx,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
          call MPI_BCAST(LP(iLayer)%nparty,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
          call MPI_BCAST(LP(iLayer)%MaxNX,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
          call MPI_BCAST(LP(iLayer)%MinNX,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
          call MPI_BCAST(LP(iLayer)%MaxNY,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
          call MPI_BCAST(LP(iLayer)%MinNY,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
          if(irank.eq.master) then
              do i=1, LP(iLayer)%nsize
                  do j=1, 4
                      LocalData((i-1)*4+j) = LP(iLayer)%PartitionInfo(i,j)
                  enddo
              enddo
              do i=0, nsize-1 ! common info is sent to all nodes including unused
                  if(i.ne.master)  then
                      call MPI_SEND(LocalData,LP(iLayer)%nsize*4,MPI_INTEGER,i,2015,MPI_COMM_WORLD,ierror)
                  endif
              enddo
          else
              call MPI_RECV(LocalData,LP(iLayer)%nsize*4,MPI_INTEGER,master,2015,MPI_COMM_WORLD,istatus,ierror)
              do i=1, LP(iLayer)%nsize
                  do j=1, 4
                      LP(iLayer)%PartitionInfo(i,j) = LocalData((i-1)*4+j)
                  enddo
              enddo
          endif
          call MPI_BCAST(LP(iLayer)%BoundarySendRecvCount,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
          if(irank.eq.master) then
              do i=1, LP(iLayer)%BoundarySendRecvCount
                  do j=1, 14
                      LocalData((i-1)*14+j) = LP(iLayer)%BoundarySendRecv(i,j)
                  enddo
              enddo
              do i=0, nsize-1
                  if(i.ne.master) call MPI_SEND(LocalData, LP(iLayer)%BoundarySendRecvCount*14, &
                      MPI_INTEGER, i, 2015, MPI_COMM_WORLD,ierror)
              enddo
          else
              call MPI_RECV(LocalData, LP(iLayer)%BoundarySendRecvCount*14, &
                  MPI_INTEGER, master, 2015, MPI_COMM_WORLD, istatus, ierror)
              do i=1, LP(iLayer)%BoundarySendRecvCount
                  do j=1, 14
                      LP(iLayer)%BoundarySendRecv(i,j) = LocalData((i-1)*14+j)
                  enddo
              enddo
          endif
          call MPI_BCAST(LP(iLayer)%ParentToChildSendRecvCount,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
          if(irank.eq.master) then
              do i=1, LP(iLayer)%ParentToChildSendRecvCount
                  do j=1, 5
                      LocalData((i-1)*5+j) = LP(iLayer)%ParentToChildSendRecv(i,j)
                  enddo
              enddo
              do i=0, nsize-1
                  if(i.ne.master) call MPI_SEND(LocalData, LP(iLayer)%ParentToChildSendRecvCount*5, &
                      MPI_INTEGER, i, 2015, MPI_COMM_WORLD,ierror)
              enddo
          else
              call MPI_RECV(LocalData, LP(iLayer)%ParentToChildSendRecvCount*5, &
                  MPI_INTEGER, master, 2015, MPI_COMM_WORLD, istatus, ierror)
              do i=1, LP(iLayer)%ParentToChildSendRecvCount
                  do j=1, 5
                      LP(iLayer)%ParentToChildSendRecv(i,j) = LocalData((i-1)*5+j)
                  enddo
              enddo
          endif
          call MPI_BCAST(LP(iLayer)%ChildToParentSendRecvCount,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
          if(irank.eq.master) then
              do i=1, LP(iLayer)%ChildToParentSendRecvCount
                  do j=1, 7
                      LocalData((i-1)*7+j) = LP(iLayer)%ChildToParentSendRecv(i,j)
                  enddo
              enddo
              do i=0, nsize-1
                  if(i.ne.master) call MPI_SEND(LocalData, LP(iLayer)%ChildToParentSendRecvCount*7, &
                      MPI_INTEGER, i, 2015, MPI_COMM_WORLD,ierror)
              enddo
          else
              call MPI_RECV(LocalData, LP(iLayer)%ChildToParentSendRecvCount*7, &
                  MPI_INTEGER, master, 2015, MPI_COMM_WORLD, istatus, ierror)
              do i=1, LP(iLayer)%ChildToParentSendRecvCount
                  do j=1, 7
                      LP(iLayer)%ChildToParentSendRecv(i,j) = LocalData((i-1)*7+j)
                  enddo
              enddo
          endif
      enddo

      do iSta=1, GP%NumStations
          call MPI_BCAST(SP(iSta)%X,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
          call MPI_BCAST(SP(iSta)%Y,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
          call MPI_BCAST(SP(iSta)%nLayer,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
          call MPI_BCAST(SP(iSta)%nNode,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
      enddo

      if(GP%InitialConditionType.eq.1) then
          do iFault=1, GP%NumFaults
              call MPI_BCAST(FP(iFault)%T0,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
              call MPI_BCAST(FP(iFault)%NT,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
              call MPI_BCAST(FP(iFault)%Depth,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
              call MPI_BCAST(FP(iFault)%Length,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
              call MPI_BCAST(FP(iFault)%Width,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
              call MPI_BCAST(FP(iFault)%Slip,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
              call MPI_BCAST(FP(iFault)%Rake,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
              call MPI_BCAST(FP(iFault)%HSlip,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
              call MPI_BCAST(FP(iFault)%PSlip,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
              call MPI_BCAST(FP(iFault)%Strike,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
              call MPI_BCAST(FP(iFault)%Dip,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
              call MPI_BCAST(FP(iFault)%X0,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
              call MPI_BCAST(FP(iFault)%Y0,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
          enddo
      endif

      end subroutine bcastCommonInfo



      subroutine bcastBathymetry(GP, LP, SP, FP, LocalData, LocalDataLength)

      use VariableDefination
      implicit NONE
      include "mpif.h"
      type(GlobalParameters)   ::  GP
      type(LayerParameters)    ::  LP(100)
      type(StationParameters)  ::  SP(999)
      type(FaultParameters)    ::  FP(999)
      integer*4  ::  LocalDataLength
      real*8     ::  LocalData(LocalDataLength)
      integer*4  ::  irank, nsize, master
      integer*4  ::  ierror, istatus(MPI_STATUS_SIZE)
      integer*4  ::  iLayer, iNode, i, j
      integer*4  ::  istartx, iendx, istarty, iendy, nRowBoundary

      irank = GP%irank; nsize = GP%nsizeTotal; master = GP%master
      if(irank.eq.master) then
          write(*,*) 'broadcasting bathymetry data from master node ...'
          write(*,*)
      endif

      do iLayer=1, GP%NumLayers
          call MPI_BCAST(LP(iLayer)%X,LP(iLayer)%NX,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
          call MPI_BCAST(LP(iLayer)%Y,LP(iLayer)%NY,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
      enddo

      nRowBoundary = GP%nRowBoundary
      if(irank.eq.master) then
          do iLayer=1, GP%NumLayers
              do iNode=0, LP(iLayer)%nsize-1
                  if(iNode.ne.master) then
                      istartx = MAX(1,LP(iLayer)%PartitionInfo(iNode+1,1)-nRowBoundary)
                      iendx   = MIN(LP(iLayer)%NX,LP(iLayer)%PartitionInfo(iNode+1,2)+nRowBoundary)
                      istarty = MAX(1,LP(iLayer)%PartitionInfo(iNode+1,3)-nRowBoundary)
                      iendy   = MIN(LP(iLayer)%NY,LP(iLayer)%PartitionInfo(iNode+1,4)+nRowBoundary)
                      do i=istartx, iendx
                          do j=istarty, iendy
                              LocalData(i-istartx+1+(j-istarty)*(iendx-istartx+1)) = LP(iLayer)%Z(i,j)
                          enddo
                      enddo
                      call MPI_SEND(LocalData,(iendx-istartx+1)*(iendy-istarty+1), &
                          MPI_DOUBLE_PRECISION,iNode,2014,MPI_COMM_WORLD,ierror)
                  endif
              enddo
          enddo
      else
          do iLayer=1, GP%NumLayers
              if(irank.lt.LP(iLayer)%nsize) then
                  istartx = MAX(1,LP(iLayer)%PartitionInfo(irank+1,1)-nRowBoundary)
                  iendx   = MIN(LP(iLayer)%NX,LP(iLayer)%PartitionInfo(irank+1,2)+nRowBoundary)
                  istarty = MAX(1,LP(iLayer)%PartitionInfo(irank+1,3)-nRowBoundary)
                  iendy   = MIN(LP(iLayer)%NY,LP(iLayer)%PartitionInfo(irank+1,4)+nRowBoundary)
                  call MPI_RECV(LocalData,(iendx-istartx+1)*(iendy-istarty+1),  &
                      MPI_DOUBLE_PRECISION,master,2014,MPI_COMM_WORLD,istatus,ierror)
                  do i=istartx, iendx
                      do j=istarty, iendy
                          LP(iLayer)%Z(i,j) = LocalData(i-istartx+1+(j-istarty)*(iendx-istartx+1))
                      enddo
                  enddo
              endif
          enddo
      endif

      end subroutine bcastBathymetry



      subroutine computeParameters(GP, LP, SP, FP)

      use VariableDefination
      implicit NONE
      include "mpif.h"
      type(GlobalParameters)   ::  GP
      type(LayerParameters)    ::  LP(100)
      type(StationParameters)  ::  SP(999)
      type(FaultParameters)    ::  FP(999)
      integer*4  ::  irank, nsize, master
      integer*4  ::  ierror, istatus(MPI_STATUS_SIZE)
      integer*4  ::  jstart, jend, nstarty, nendy, iLayer, j

      irank = GP%irank; nsize = GP%nsizeTotal; master = GP%master
      if(irank.eq.master) then
          write(*,*) 'computing parameters (e.g., Coriolis) ...'
          write(*,*)
      endif

      do iLayer=1, GP%NumLayers
      if(irank.lt.LP(iLayer)%nsize) then
          !////// Coriolis parameters //////!
          if(GP%CoordinatesType.eq.0) then
              nstarty = LP(iLayer)%PartitionInfo(irank+1,3)
              nendy = LP(iLayer)%PartitionInfo(irank+1,4)
              jstart = MAX(1,nstarty)
              jend = MIN(LP(iLayer)%NY,nendy)
              do j=jstart, jend
                  LP(iLayer)%CPX(j)=2.0*GP%Omega*SIN(LP(iLayer)%Y(j)*GP%PI/180.0)*LP(iLayer)%dt
              enddo
              jend = MIN(LP(iLayer)%NY-1,nendy);
              do j=jstart, jend
                  LP(iLayer)%CPY(j)=2.0*GP%Omega*SIN((LP(iLayer)%Y(j)+0.5*LP(iLayer)%dy)*GP%PI/180.0)*LP(iLayer)%dt
              enddo
          endif

          !//// mass/moment equations parameters //////!
          if(GP%CoordinatesType.eq.0) then
              nstarty = LP(iLayer)%PartitionInfo(irank+1,3)
              nendy = LP(iLayer)%PartitionInfo(irank+1,4)
              jstart = MAX(1,nstarty-1)
              jend = MIN(LP(iLayer)%NY,nendy+1)
              do j=jstart, jend
                  LP(iLayer)%CSY(j) = COS(LP(iLayer)%Y(j)*GP%PI/180.0)
                  LP(iLayer)%CS1(j) = LP(iLayer)%dt/GP%R_EARTH/LP(iLayer)%dx/GP%PI*180.0/LP(iLayer)%CSY(j)
                  LP(iLayer)%CS2(j) = LP(iLayer)%dt/GP%R_EARTH/LP(iLayer)%dy/GP%PI*180.0/LP(iLayer)%CSY(j)
                  LP(iLayer)%CS3(j) = LP(iLayer)%CS1(j)*GP%GRAV*0.5
                  LP(iLayer)%CS4    = LP(iLayer)%dt*GP%GRAV*0.5/GP%R_EARTH/LP(iLayer)%dy/GP%PI*180.0
              enddo
          else
              LP(iLayer)%CC1 = LP(iLayer)%dt/LP(iLayer)%dx
              LP(iLayer)%CC2 = LP(iLayer)%dt/LP(iLayer)%dy
              LP(iLayer)%CC3 = LP(iLayer)%CC1*GP%GRAV*0.5
              LP(iLayer)%CC4 = LP(iLayer)%CC2*GP%GRAV*0.5
          endif
      endif
      enddo

      end subroutine computeParameters



      subroutine removeStationFiles(GP, LP)

      use VariableDefination
      implicit NONE
      type(GlobalParameters)   ::  GP
      type(LayerParameters)    ::  LP(100)
      integer*4                ::  i
      character(999)           ::  s

      do i=1, GP%NumStations
          write(s,'(a7,i4.4,a4)') 'Station',i,'.dat'
          open(1000+i,file=s,form='unformatted')
          write(1000+i) GP%TotalTimeSteps+1, GP%NumFaults, GP%ComputeGreen
          close(1000+i)
      enddo
      if(GP%SaveFlux.eq.1) then
          do i=1, GP%NumStations
              write(s,'(a7,i4.4,a6)') 'Station',i,'_M.dat'
              open(1000+i,file=s,form='unformatted')
              write(1000+i) GP%TotalTimeSteps+1, GP%NumFaults, GP%ComputeGreen
              close(1000+i)
          enddo
          do i=1, GP%NumStations
              write(s,'(a7,i4.4,a6)') 'Station',i,'_N.dat'
              open(1000+i,file=s,form='unformatted')
              write(1000+i) GP%TotalTimeSteps+1, GP%NumFaults, GP%ComputeGreen
              close(1000+i)
          enddo
      endif

      end subroutine removeStationFiles



      subroutine modifyFaultParameters(GP, FP, iCal)

      use VariableDefination
      implicit NONE
      type(GlobalParameters)   ::  GP
      type(FaultParameters)    ::  FP(999)
      integer*4                ::  iCal, iFault
      integer*4                ::  irank, nsize, master

      irank = GP%irank; nsize = GP%nsizeTotal; master = GP%master
      if(irank.eq.master) then
          write(*,*)
          write(*,*) '***************************************************'
          write(*,'(a,i4,a,i4)') ' Green function for subfault ',iCal,' of ',GP%NumFaults
          write(*,*) '***************************************************'
          write(*,*)
      endif
      do iFault=1, GP%Numfaults
          if(iFault.eq.iCal) then
              FP(iFault)%T0 = 0
              FP(iFault)%NT = 0
              FP(iFault)%Slip = 1.0
          else
              FP(iFault)%T0 = GP%TotalTime*2
              FP(iFault)%NT = GP%TotalTimeSteps*2
          endif
      enddo

      end subroutine modifyFaultParameters



      subroutine getInitialCondition(GP, LP, SP, FP, iCal, iTimeStep, LocalData, LocalDataLength)

      use VariableDefination
      implicit NONE
      include "mpif.h"
      type(GlobalParameters)   ::  GP
      type(LayerParameters)    ::  LP(100)
      type(StationParameters)  ::  SP(999)
      type(FaultParameters)    ::  FP(999)
      integer*4  ::  iCal, iTimeStep
      integer*4  ::  LocalDataLength
      real*8     ::  LocalData(LocalDataLength)
      integer*4  ::  irank, nsize, master
      integer*4  ::  ierror, istatus(MPI_STATUS_SIZE)
      integer*4  ::  iLayer, iLayerLevel, iNode, iFault, iFlag, i, j
      integer*4  ::  istartx, iendx, istarty, iendy, nRowBoundary
      real*8     ::  x, y, h, u1, u2
      real*8     ::  depth,xlength,xwidth
      real*8     ::  slip,strike,dip,rake,x0,y0

      irank = GP%irank; nsize = GP%nsizeTotal; master = GP%master

      nRowBoundary = GP%nRowBoundary
      if(GP%InitialConditionType.eq.0.and.iTimeStep.eq.0) then
          if(irank.eq.master) then
              write(*,*) 'reading initial water elevation from file ',TRIM(ADJUSTL(GP%InitialConditionFileName)),' ...'
              write(*,*)
              if(GP%InitialConditionFileFormat.eq.1) then
                  call readInitialConditionNetCDF(GP, LP)
              else
                  open(23,file=TRIM(ADJUSTL(GP%InitialConditionFileName)),status='old',form='formatted')
                  do j=1, LP(1)%NY
                      do i=1, LP(1)%NX
                          read(23,*) x,y,LP(1)%H(1,i,j)
                          if(LP(1)%Z(i,j).le.GP%WaterDepthLimit) LP(1)%H(1,i,j)=0.0
                      enddo
                  enddo
                  close(23)
              endif
              do iLayerLevel=2, GP%NumLayerLevels
                  do iLayer=1, GP%NumLayers
                      if(LP(iLayer)%Level.eq.iLayerLevel) then
                          do j=1, LP(iLayer)%NY
                              do i=1, LP(iLayer)%NX
                                  x = LP(iLayer)%X(i); y = LP(iLayer)%Y(j); iFlag = -1
                                  call interpData(GP, LP, LP(iLayer)%Parent, iFlag, x, y, h)
                                  LP(iLayer)%H(1,i,j) = h
                              enddo
                          enddo
                      endif
                  enddo
              enddo
          endif
          if(irank.eq.master) then
              do iLayer=1, GP%NumLayers
                  do iNode=0, LP(iLayer)%nsize-1
                      if(iNode.ne.master) then
                          istartx = MAX(1,LP(iLayer)%PartitionInfo(iNode+1,1)-nRowBoundary)
                          iendx   = MIN(LP(iLayer)%NX,LP(iLayer)%PartitionInfo(iNode+1,2)+nRowBoundary)
                          istarty = MAX(1,LP(iLayer)%PartitionInfo(iNode+1,3)-nRowBoundary)
                          iendy   = MIN(LP(iLayer)%NY,LP(iLayer)%PartitionInfo(iNode+1,4)+nRowBoundary)
                          do i=istartx, iendx
                              do j=istarty, iendy
                                  LocalData(i-istartx+1+(j-istarty)*(iendx-istartx+1)) = LP(iLayer)%H(1,i,j)
                              enddo
                          enddo
                          call MPI_SEND(LocalData,(iendx-istartx+1)*(iendy-istarty+1), &
                              MPI_DOUBLE_PRECISION,iNode,iLayer,MPI_COMM_WORLD,ierror)
                      endif
                  enddo
              enddo
          else
              do iLayer=1, GP%NumLayers
                  if(irank.lt.LP(iLayer)%nsize) then
                      istartx = MAX(1,LP(iLayer)%PartitionInfo(irank+1,1)-nRowBoundary)
                      iendx   = MIN(LP(iLayer)%NX,LP(iLayer)%PartitionInfo(irank+1,2)+nRowBoundary)
                      istarty = MAX(1,LP(iLayer)%PartitionInfo(irank+1,3)-nRowBoundary)
                      iendy   = MIN(LP(iLayer)%NY,LP(iLayer)%PartitionInfo(irank+1,4)+nRowBoundary)
                      call MPI_RECV(LocalData,(iendx-istartx+1)*(iendy-istarty+1),  &
                          MPI_DOUBLE_PRECISION,master,iLayer,MPI_COMM_WORLD,istatus,ierror)
                      do i=istartx, iendx
                          do j=istarty, iendy
                              LP(iLayer)%H(1,i,j) = LocalData(i-istartx+1+(j-istarty)*(iendx-istartx+1))
                          enddo
                      enddo
                  endif
              enddo
          endif
      endif

      if(GP%InitialConditionType.eq.1) then
          if(irank.eq.master.and.iTimeStep.eq.0) then
              write(*,*) 'calculating initial water elevation from Okada''s model...'
              write(*,*)
          endif
          do iFault=1, GP%NumFaults
              if(FP(iFault)%NT.eq.iTimeStep) then
                  if(irank.eq.master.and.GP%ComputeGreen.ne.1) then
                      write(*,'(a,i5,a,i5)') '        sub fault ',iFault,'    of',GP%Numfaults
                  endif
                  do iLayer=1, GP%NumLayers
                  if(irank.lt.LP(iLayer)%nsize) then
                      istartx = MAX(1,LP(iLayer)%PartitionInfo(irank+1,1)-nRowBoundary)
                      iendx   = MIN(LP(iLayer)%NX,LP(iLayer)%PartitionInfo(irank+1,2)+nRowBoundary)
                      istarty = MAX(1,LP(iLayer)%PartitionInfo(irank+1,3)-nRowBoundary)
                      iendy   = MIN(LP(iLayer)%NY,LP(iLayer)%PartitionInfo(irank+1,4)+nRowBoundary)
                      do i=istartx,iendx
                          do j=istarty,iendy
                              if(LP(iLayer)%Z(i,j).gt.GP%WaterDepthLimit) then
                                  x = LP(iLayer)%X(i); y = LP(iLayer)%Y(j);
                                  depth = FP(iFault)%Depth; xlength = FP(ifault)%Length
                                  xwidth = FP(iFault)%Width; strike = FP(iFault)%Strike
                                  dip = FP(iFault)%Dip; slip = FP(iFault)%Slip
                                  rake = FP(iFault)%Rake; x0 = FP(iFault)%X0; y0 = FP(iFault)%Y0
                                  if(GP%CoordinatesType.eq.0) then
                                      x = GP%R_Earth*COS(y0/180.0*GP%PI)*(x-x0)/180.0*GP%PI
                                      y = GP%R_Earth*(y-y0)/180.0*GP%PI
                                      x0 = 0.0
                                      y0 = 0.0
                                  endif
                                  call okada1985(x,y,u1,u2,h,depth,xlength,xwidth,slip,strike,dip,rake,x0,y0)
                                  LP(iLayer)%H(1,i,j) = LP(iLayer)%H(1,i,j)+h
                              endif
                          enddo
                      enddo
                  endif
                  enddo
              endif
          enddo
          if(irank.eq.master.and.iTimeStep.eq.0) write(*,*)
      endif

      end subroutine getInitialCondition



      subroutine interpData(GP, LP, iLayer, iFlag, x, y, val)
      !! iFlag = -1,-2,-3: interpolate H,M,N(1,i,j)
      !! iFlag = 1,2,3: interpolate H,M,N(2,i,j)

      use VariableDefination
      implicit NONE
      type(GlobalParameters)   ::  GP
      type(LayerParameters)    ::  LP(100)
      integer*4                ::  iLayer, iFlag
      real*8                   ::  x, y, val
      integer*4                ::  ii, jj, ii2, jj2
      integer*4                ::  i1, i2, i3, i4
      real*8                   ::  x1, y1, x2, y2, invx, invy
      real*8                   ::  z1, z2, z3, z4

      if(iFlag.eq.-1) then
          ii = FLOOR((x-LP(iLayer)%X(1))/LP(iLayer)%dx)+1
          jj = FLOOR((y-LP(iLayer)%Y(1))/LP(iLayer)%dy)+1
          ii2 = MIN(ii+1, LP(iLayer)%NX)
          jj2 = MIN(jj+1, LP(iLayer)%NY)
          z1 = LP(iLayer)%H(1,ii,jj)
          z2 = LP(iLayer)%H(1,ii2,jj)
          z3 = LP(iLayer)%H(1,ii,jj2)
          z4 = LP(iLayer)%H(1,ii2,jj2)
          x1 = LP(iLayer)%X(ii); y1 = LP(iLayer)%Y(jj)
          x2 = LP(iLayer)%X(ii2); y2 = LP(iLayer)%Y(jj2)
      else if(iFlag.eq.1) then
          ii = FLOOR((x-LP(iLayer)%X(1))/LP(iLayer)%dx)+1
          jj = FLOOR((y-LP(iLayer)%Y(1))/LP(iLayer)%dy)+1
          ii2 = MIN(ii+1, LP(iLayer)%NX)
          jj2 = MIN(jj+1, LP(iLayer)%NY)
          z1 = LP(iLayer)%H(2,ii,jj)
          z2 = LP(iLayer)%H(2,ii2,jj)
          z3 = LP(iLayer)%H(2,ii,jj2)
          z4 = LP(iLayer)%H(2,ii2,jj2)
          x1 = LP(iLayer)%X(ii); y1 = LP(iLayer)%Y(jj)
          x2 = LP(iLayer)%X(ii2); y2 = LP(iLayer)%Y(jj2)
      else if(iFlag.eq.-2) then
          ii = FLOOR((x-LP(iLayer)%X(1)-LP(iLayer)%dx*0.5)/LP(iLayer)%dx)+1
          jj = FLOOR((y-LP(iLayer)%Y(1))/LP(iLayer)%dy)+1
          ii2 = MIN(ii+1, LP(iLayer)%NX-1)
          jj2 = MIN(jj+1, LP(iLayer)%NY)
          if(ii.eq.0) then
              ii = 1; ii2 = 1
          endif
          z1 = LP(iLayer)%M(1,ii,jj)
          z2 = LP(iLayer)%M(1,ii2,jj)
          z3 = LP(iLayer)%M(1,ii,jj2)
          z4 = LP(iLayer)%M(1,ii2,jj2)
          x1 = LP(iLayer)%X(ii)+LP(iLayer)%dx*0.5; y1 = LP(iLayer)%Y(jj)
          x2 = LP(iLayer)%X(ii2)+LP(iLayer)%dx*0.5; y2 = LP(iLayer)%Y(jj2)
      else if(iFlag.eq.2) then
          ii = FLOOR((x-LP(iLayer)%X(1)-LP(iLayer)%dx*0.5)/LP(iLayer)%dx)+1
          jj = FLOOR((y-LP(iLayer)%Y(1))/LP(iLayer)%dy)+1
          ii2 = MIN(ii+1, LP(iLayer)%NX-1)
          jj2 = MIN(jj+1, LP(iLayer)%NY)
          if(ii.eq.0) then
              ii = 1; ii2 = 1
          endif
          z1 = LP(iLayer)%M(2,ii,jj)
          z2 = LP(iLayer)%M(2,ii2,jj)
          z3 = LP(iLayer)%M(2,ii,jj2)
          z4 = LP(iLayer)%M(2,ii2,jj2)
          x1 = LP(iLayer)%X(ii)+LP(iLayer)%dx*0.5; y1 = LP(iLayer)%Y(jj)
          x2 = LP(iLayer)%X(ii2)+LP(iLayer)%dx*0.5; y2 = LP(iLayer)%Y(jj2)
      else if(iFlag.eq.-3) then
          ii = FLOOR((x-LP(iLayer)%X(1))/LP(iLayer)%dx)+1
          jj = FLOOR((y-LP(iLayer)%Y(1)-LP(iLayer)%dy*0.5)/LP(iLayer)%dy)+1
          ii2 = MIN(ii+1, LP(iLayer)%NX)
          jj2 = MIN(jj+1, LP(iLayer)%NY-1)
          if(jj.eq.0) then
              jj = 1; jj2 = 1
          endif
          z1 = LP(iLayer)%N(1,ii,jj)
          z2 = LP(iLayer)%N(1,ii2,jj)
          z3 = LP(iLayer)%N(1,ii,jj2)
          z4 = LP(iLayer)%N(1,ii2,jj2)
          x1 = LP(iLayer)%X(ii); y1 = LP(iLayer)%Y(jj)+LP(iLayer)%dy/2
          x2 = LP(iLayer)%X(ii2); y2 = LP(iLayer)%Y(jj2)+LP(iLayer)%dy/2
      else if(iFlag.eq.3) then
          ii = FLOOR((x-LP(iLayer)%X(1))/LP(iLayer)%dx)+1
          jj = FLOOR((y-LP(iLayer)%Y(1)-LP(iLayer)%dy*0.5)/LP(iLayer)%dy)+1
          ii2 = MIN(ii+1, LP(iLayer)%NX)
          jj2 = MIN(jj+1, LP(iLayer)%NY-1)
          if(jj.eq.0) then
              jj = 1; jj2 = 1
          endif
          z1 = LP(iLayer)%N(2,ii,jj)
          z2 = LP(iLayer)%N(2,ii2,jj)
          z3 = LP(iLayer)%N(2,ii,jj2)
          z4 = LP(iLayer)%N(2,ii2,jj2)
          x1 = LP(iLayer)%X(ii); y1 = LP(iLayer)%Y(jj)+LP(iLayer)%dy*0.5
          x2 = LP(iLayer)%X(ii2); y2 = LP(iLayer)%Y(jj2)+LP(iLayer)%dy*0.5
      endif
      i1 = 1; i2 = 1; i3 = 1; i4 = 1
      if(LP(iLayer)%Z(ii,jj).le.GP%WaterDepthLimit) i1 = 0
      if(LP(iLayer)%Z(ii2,jj).le.GP%WaterDepthLimit) i2 = 0
      if(LP(iLayer)%Z(ii,jj2).le.GP%WaterDepthLimit) i3 = 0
      if(LP(iLayer)%Z(ii2,jj2).le.GP%WaterDepthLimit) i4 = 0
      if(i1.eq.1.and.i2.eq.0.and.i3.eq.0.and.i4.eq.0) then
          ii2 = ii; jj2 = jj
      elseif(i1.eq.0.and.i2.eq.1.and.i3.eq.0.and.i4.eq.0) then
          ii = ii2; jj2 = jj; z1 = z2
      elseif(i1.eq.0.and.i2.eq.0.and.i3.eq.1.and.i4.eq.0) then
          ii2 = ii; jj = jj2; z1 = z3
      elseif(i1.eq.0.and.i2.eq.0.and.i3.eq.0.and.i4.eq.1) then
          ii = ii2; jj = jj2; z1 = z4
      elseif(i1.eq.1.and.i2.eq.1.and.i3.eq.0.and.i4.eq.0) then
          jj2 = jj;
      elseif(i1.eq.1.and.i2.eq.0.and.i3.eq.1.and.i4.eq.0) then
          ii2 = ii;
      elseif(i1.eq.0.and.i2.eq.1.and.i3.eq.0.and.i4.eq.1) then
          ii = ii2; z1 = z2; z3 = z4
      elseif(i1.eq.0.and.i2.eq.0.and.i3.eq.1.and.i4.eq.1) then
          jj = jj2; z1 = z3; z2 = z4
      elseif((i1.eq.0.or.i2.eq.0).and.i3.eq.1.and.i4.eq.1) then
          jj = jj2; z1 = z3; z2 = z4
      elseif(i1.eq.1.and.i2.eq.1.and.(i3.eq.0.or.i4.eq.0)) then
          jj2 = jj
      endif

      if(ii.ne.ii2) invx = 1.0/(x2-x1)
      if(jj.ne.jj2) invy = 1.0/(y2-y1)
      if(ii.eq.ii2.or.jj.eq.jj2) then
          if(ii.eq.ii2.and.jj.eq.jj2) then
              val = z1
          else if(ii.eq.ii2) then
              val = -z1*(y-y2)*invy+z3*(y-y1)*invy
          else
              val = -z1*(x-x2)*invx+z2*(x-x1)*invx
          endif
      else
          val = z1*(x2-x)*(y2-y)+z2*(x-x1)*(y2-y)+z3*(x2-x)*(y-y1)+z4*(x-x1)*(y-y1)
          val = val*invx*invy
      endif

      end subroutine interpData



      subroutine screenOutput(GP, LP, iTimeStep)

      use VariableDefination
      implicit NONE
      include "mpif.h"
      type(GlobalParameters)   ::  GP
      type(LayerParameters)    ::  LP(100)
      integer*4  ::  iTimeStep
      real*8     ::  CPUTime
      real*8     ::  ETA
      integer*4  ::  nh, nm, ns, i, j, iNode
      integer*4  ::  ierror, istatus(MPI_STATUS_SIZE)
      real*8     ::  LocalData(999*10), tt, calt, comt, writ

      if(GP%irank.eq.GP%master) then
      if(iTimeStep.eq.0) then
          write(*,*)
          write(*,*) '================= OUTPUT RESULTS =================='
          write(*,*) '  TimeStep   Seconds   Percentage   ETA(hh:mm:ss)'
          write(*,'(i11,f10.1,7x,f4.1,a2)',advance='no') &
              iTimeStep,iTimeStep*GP%dt,iTimeStep*1.0/GP%TotalTimeSteps*100,' %'
      elseif(iTimeStep.le.GP%TotalTimeSteps) then
          call CPU_TIME(CPUTime)
          if((iTimeStep.eq.10).or. &
              (iTimeStep.gt.10.and.CPUTime-GP%CPUTimeInitial-GP%CPUTime(GP%master+1,1).ge.15.0).or.&
              (iTimeStep.gt.10.and.MOD(iTimeStep,MAX(1,GP%TotalTimeSteps/20)).eq.0)) then
              ETA = (GP%TotalTimeSteps-iTimeStep)*(CPUTime-GP%CPUTimeInitialCalculation)/iTimeStep
              nh = FLOOR(ETA/3600)
              nm = FLOOR((ETA-nh*3600)/60)
              ns = INT(ETA-nh*3600-nm*60)
              if(iTimeStep.eq.10) then
                  write(*,'(i9.2,a1,i2.2,a1,i2.2)') nh,':',nm,':',ns
              else
                write(*,'(i11,f10.1,6x,f5.1,a2,i9.2,a1,i2.2,a1,i2.2)') &
                  iTimeStep,iTimeStep*GP%dt,iTimeStep*1.0/GP%TotalTimeSteps*100,'%',&
                  nh,':',nm,':',ns
              endif
              GP%CPUTime(GP%master+1,1) = CPUTime-GP%CPUTimeInitial
          endif
      endif
      endif

      if(iTimeStep.eq.GP%TotalTimeSteps+100) then
          if(GP%irank.ne.GP%master) then
              do i=1, GP%nsizeTotal
              do j=1, 5
                  LocalData((i-1)+j) = GP%CPUTime(GP%irank+1,j)
              enddo
              enddo
              call MPI_SEND(LocalData,GP%nsizeTotal*5,MPI_DOUBLE_PRECISION, &
                  GP%master,2015,MPI_COMM_WORLD,ierror)
          else
              tt   = GP%CPUTime(GP%irank+1,1); calt = GP%CPUTime(GP%irank+1,3);
              comt = GP%CPUTime(GP%irank+1,4); writ = GP%CPUTime(GP%irank+1,5);
              do iNode=0, GP%nsizeTotal-1
              if(iNode.ne.GP%master) then
                  call MPI_RECV(LocalData,GP%nsizeTotal*5,MPI_DOUBLE_PRECISION, &
                      iNode,2015,MPI_COMM_WORLD,istatus,ierror)
                  do i=1, GP%nsizeTotal
                  do j=1, 5
                      GP%CPUTime(iNode+1,j) = LocalData((i-1)+j)
                  enddo
                  enddo
              endif
              enddo
              do iNode=0, GP%nsizeTotal-1
                  tt   = MAX(  tt,GP%CPUTime(iNode+1,1)); calt = MAX(calt,GP%CPUTime(iNode+1,3));
                  comt = MAX(comt,GP%CPUTime(iNode+1,4)); writ = MAX(writ,GP%CPUTime(iNode+1,5));
              enddo
              write(*,*) '==================================================='
              write(*,*)
              write(*,*) 'pcomcot calculation complete.'
              write(*,*)
              write(*,'(a)')              ' Max Time On Node : '
              write(*,'(a,f7.1,f10.1,a)') ' Total Time       : ',  tt,100.0d0,' %'
              write(*,'(a,f7.1,f10.1,a)') ' Calculation      : ',calt,calt/tt*100,' %'
              write(*,'(a,f7.1,f10.1,a)') ' Communication    : ',comt,comt/tt*100,' %'
              write(*,'(a,f7.1,f10.1,a)') ' Write to File    : ',writ,writ/tt*100,' %'
              write(*,'(a,f7.1,f10.1,a)') ' Other            : ', &
                  tt-calt-comt-writ,(tt-calt-comt-writ)/tt*100,' %'
              write(*,*)
          endif
      endif

      end subroutine screenOutput



      subroutine bcastComputeDomainBoundaryValue(GP, LP, iLayer, iCal, iTimeStep, LocalData, LocalDataLength)

      use VariableDefination
      implicit NONE
      include "mpif.h"
      type(GlobalParameters)   ::  GP
      type(LayerParameters)    ::  LP(100)
      integer*4  ::  iLayer, iCal, iTimeStep
      integer*4  ::  LocalDataLength
      real*8     ::  LocalData(LocalDataLength)
      integer*4  ::  irank, nsize, master
      integer*4  ::  ierror, istatus(MPI_STATUS_SIZE)
      integer*4  ::  iBC, i, j
      integer*4  ::  nstartx, nendx, nstarty, nendy
      real*8     ::  CPUTime1, CPUTime2

      irank = GP%irank; nsize = GP%nsizeTotal; master = GP%master
      call CPU_TIME(CPUTime1)

      do iBC=1, LP(iLayer)%BoundarySendRecvCount
          if(irank.eq.LP(iLayer)%BoundarySendRecv(iBC, 1)) then
              nstartx = LP(iLayer)%BoundarySendRecv(iBC,3)
              nendx   = LP(iLayer)%BoundarySendRecv(iBC,4)
              nstarty = LP(iLayer)%BoundarySendRecv(iBC,5)
              nendy   = LP(iLayer)%BoundarySendRecv(iBC,6)
              do j=nstarty, nendy
                  do i=nstartx, nendx
                      LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1) = LP(iLayer)%H(2,i,j)
                  enddo
              enddo
              call MPI_SEND(LocalData,(nendx-nstartx+1)*(nendy-nstarty+1),MPI_DOUBLE_PRECISION, &
                  LP(iLayer)%BoundarySendRecv(iBC,2),iLayer,MPI_COMM_WORLD,ierror)
          elseif(irank.eq.LP(iLayer)%BoundarySendRecv(iBC,2)) then
              nstartx = LP(iLayer)%BoundarySendRecv(iBC,3)
              nendx   = LP(iLayer)%BoundarySendRecv(iBC,4)
              nstarty = LP(iLayer)%BoundarySendRecv(iBC,5)
              nendy   = LP(iLayer)%BoundarySendRecv(iBC,6)
              call MPI_RECV(LocalData,(nendx-nstartx+1)*(nendy-nstarty+1),MPI_DOUBLE_PRECISION, &
                  LP(iLayer)%BoundarySendRecv(iBC,1),iLayer,MPI_COMM_WORLD,istatus,ierror)
              do j=nstarty, nendy
                  do i=nstartx, nendx
                      LP(iLayer)%H(2,i,j) = LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1)
                  enddo
              enddo
          endif
      enddo
      call CPU_TIME(CPUTime2); GP%CPUTime(irank+1,4) = GP%CPUTime(irank+1,4)+CPUTime2-CPUTime1

      end subroutine bcastComputeDomainBoundaryValue



      subroutine bcastComputeDomainBoundaryFlux(GP, LP, iLayer, iCal, iTimeStep, LocalData, LocalDataLength)

      use VariableDefination
      implicit NONE
      include "mpif.h"
      type(GlobalParameters)   ::  GP
      type(LayerParameters)    ::  LP(100)
      integer*4  ::  iLayer, iCal, iTimeStep
      integer*4  ::  LocalDataLength
      real*8     ::  LocalData(LocalDataLength)
      integer*4  ::  irank, nsize, master
      integer*4  ::  ierror, istatus(MPI_STATUS_SIZE)
      integer*4  ::  iBC, i, j
      integer*4  ::  nstartx, nendx, nstarty, nendy
      real*8     ::  CPUTime1, CPUTime2

      irank = GP%irank; nsize = GP%nsizeTotal; master = GP%master
      call CPU_TIME(CPUTime1)

      do iBC=1, LP(iLayer)%BoundarySendRecvCount
          if(irank.eq.LP(iLayer)%BoundarySendRecv(iBC, 1)) then
              nstartx = LP(iLayer)%BoundarySendRecv(iBC, 7)
              nendx   = LP(iLayer)%BoundarySendRecv(iBC, 8)
              nstarty = LP(iLayer)%BoundarySendRecv(iBC, 9)
              nendy   = LP(iLayer)%BoundarySendRecv(iBC,10)
              do j=nstarty, nendy
                  do i=nstartx, nendx
                      LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1) = LP(iLayer)%M(2,i,j)
                  enddo
              enddo
              call MPI_SEND(LocalData,(nendx-nstartx+1)*(nendy-nstarty+1),MPI_DOUBLE_PRECISION, &
                  LP(iLayer)%BoundarySendRecv(iBC,2),iLayer,MPI_COMM_WORLD,ierror)
              nstartx = LP(iLayer)%BoundarySendRecv(iBC,11)
              nendx   = LP(iLayer)%BoundarySendRecv(iBC,12)
              nstarty = LP(iLayer)%BoundarySendRecv(iBC,13)
              nendy   = LP(iLayer)%BoundarySendRecv(iBC,14)
              do j=nstarty, nendy
                  do i=nstartx, nendx
                      LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1) = LP(iLayer)%N(2,i,j)
                  enddo
              enddo
              call MPI_SEND(LocalData,(nendx-nstartx+1)*(nendy-nstarty+1),MPI_DOUBLE_PRECISION, &
                  LP(iLayer)%BoundarySendRecv(iBC,2),iLayer,MPI_COMM_WORLD,ierror)
          elseif(irank.eq.LP(iLayer)%BoundarySendRecv(iBC,2)) then
              nstartx = LP(iLayer)%BoundarySendRecv(iBC, 7)
              nendx   = LP(iLayer)%BoundarySendRecv(iBC, 8)
              nstarty = LP(iLayer)%BoundarySendRecv(iBC, 9)
              nendy   = LP(iLayer)%BoundarySendRecv(iBC,10)
              call MPI_RECV(LocalData,(nendx-nstartx+1)*(nendy-nstarty+1),MPI_DOUBLE_PRECISION, &
                  LP(iLayer)%BoundarySendRecv(iBC,1),iLayer,MPI_COMM_WORLD,istatus,ierror)
              do j=nstarty, nendy
                  do i=nstartx, nendx
                      LP(iLayer)%M(2,i,j) = LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1)
                  enddo
              enddo
              nstartx = LP(iLayer)%BoundarySendRecv(iBC,11)
              nendx   = LP(iLayer)%BoundarySendRecv(iBC,12)
              nstarty = LP(iLayer)%BoundarySendRecv(iBC,13)
              nendy   = LP(iLayer)%BoundarySendRecv(iBC,14)
              call MPI_RECV(LocalData,(nendx-nstartx+1)*(nendy-nstarty+1),MPI_DOUBLE_PRECISION, &
                  LP(iLayer)%BoundarySendRecv(iBC,1),iLayer,MPI_COMM_WORLD,istatus,ierror)
              do j=nstarty, nendy
                  do i=nstartx, nendx
                      LP(iLayer)%N(2,i,j) = LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1)
                  enddo
              enddo
          endif
      enddo
      call CPU_TIME(CPUTime2); GP%CPUTime(irank+1,4) = GP%CPUTime(irank+1,4)+CPUTime2-CPUTime1

      end subroutine bcastComputeDomainBoundaryFlux



      subroutine updateComputedResults(GP, LP, iLayer)

      use VariableDefination
      implicit NONE
      type(GlobalParameters)   ::  GP
      type(LayerParameters)    ::  LP(100)
      integer*4  ::  iLayer
      integer*4  ::  irank, nstartx, nendx, nstarty, nendy
      integer*4  ::  istartx, iendx, istarty, iendy
      integer*4  ::  i, j

      irank = GP%irank
      if(irank.lt.LP(iLayer)%nsize) then
          istartx = MAX(1,LP(iLayer)%PartitionInfo(irank+1,1)-1)
          iendx   = MIN(LP(iLayer)%NX,LP(iLayer)%PartitionInfo(irank+1,2)+1)
          istarty = MAX(1,LP(iLayer)%PartitionInfo(irank+1,3)-1)
          iendy   = MIN(LP(iLayer)%NY,LP(iLayer)%PartitionInfo(irank+1,4)+1)
          do j=istarty, iendy
              do i=istartx, iendx
                  LP(iLayer)%H(1,i,j) = LP(iLayer)%H(2,i,j)
                  if(i.ne.LP(iLayer)%NX) LP(iLayer)%M(1,i,j) = LP(iLayer)%M(2,i,j)
                  if(j.ne.LP(iLayer)%NY) LP(iLayer)%N(1,i,j) = LP(iLayer)%N(2,i,j)
              enddo
          enddo
      endif

      end subroutine updateComputedResults



      subroutine calculateStationData(GP, LP, SP, iCal, iTimeStep)

      use VariableDefination
      implicit NONE
      type(GlobalParameters)   ::  GP
      type(LayerParameters)    ::  LP(100)
      type(StationParameters)  ::  SP(999)
      integer*4  ::  iCal, iTimeStep, irank
      integer*4  ::  iSta, iFlag
      real*8     ::  h

      irank = GP%irank
      do iSta=1, GP%NumStations
          if(irank.eq.SP(iSta)%nNode) then
              iFlag = -1
              call interpData(GP, LP, SP(iSta)%nLayer, iFlag, SP(iSta)%X, SP(iSta)%Y, h)
              SP(iSta)%H(iTimeStep+1) = h

              iFlag = -2
              call interpData(GP, LP, SP(iSta)%nLayer, iFlag, SP(iSta)%X, SP(iSta)%Y, h)
              SP(iSta)%M(iTimeStep+1) = h

              iFlag = -3
              call interpData(GP, LP, SP(iSta)%nLayer, iFlag, SP(iSta)%X, SP(iSta)%Y, h)
              SP(iSta)%N(iTimeStep+1) = h
          endif
      enddo

      end subroutine calculateStationData



      subroutine saveStationData(GP, LP, SP, iCal)

      use VariableDefination
      implicit NONE
      include "mpif.h"
      type(GlobalParameters)    ::  GP
      type(LayerParameters)     ::  LP(100)
      type(StationParameters)   ::  SP(999)
      integer*4        ::  iCal
      integer*4        ::  irank, nsize, master
      integer*4        ::  ierror, istatus(MPI_STATUS_SIZE)
      integer*4        ::  iSta, i, j
      character(999)   ::  s

      irank = GP%irank; nsize = GP%nsizeTotal; master = GP%master
      do iSta=1, GP%NumStations
          if(irank.eq.SP(iSta)%nNode.and.irank.ne.master) then
              call MPI_SEND(SP(iSta)%H, GP%TotalTimeSteps+1, &
                  MPI_DOUBLE_PRECISION, master, 2015, MPI_COMM_WORLD, ierror)
              call MPI_SEND(SP(iSta)%M, GP%TotalTimeSteps+1, &
                  MPI_DOUBLE_PRECISION, master, 2015, MPI_COMM_WORLD, ierror)
              call MPI_SEND(SP(iSta)%N, GP%TotalTimeSteps+1, &
                  MPI_DOUBLE_PRECISION, master, 2015, MPI_COMM_WORLD, ierror)
          endif
      enddo
      do iSta=1, GP%NumStations
          if(irank.eq.master.and.irank.ne.SP(iSta)%nNode) then
              call MPI_RECV(SP(iSta)%H, GP%TotalTimeSteps+1, &
                  MPI_DOUBLE_PRECISION, SP(iSta)%nNode, 2015, MPI_COMM_WORLD, istatus, ierror)
              call MPI_RECV(SP(iSta)%M, GP%TotalTimeSteps+1, &
                  MPI_DOUBLE_PRECISION, SP(iSta)%nNode, 2015, MPI_COMM_WORLD, istatus, ierror)
              call MPI_RECV(SP(iSta)%N, GP%TotalTimeSteps+1, &
                  MPI_DOUBLE_PRECISION, SP(iSta)%nNode, 2015, MPI_COMM_WORLD, istatus, ierror)
          endif
      enddo

      if(irank.eq.master) then
      do i=1, GP%NumStations
          write(s,'(a7,i4.4,a4)') 'Station',i,'.dat'
          open(1000+i,file=s,form='unformatted',access='append')
          if(iCal.eq.1) then
              write(1000+i) (GP%t(j), j=1,GP%TotalTimeSteps+1)
          endif
          write(1000+i) (SP(i)%H(j), j=1,GP%TotalTimeSteps+1)
          close(1000+i)

          if(GP%SaveFlux.eq.1) then

          write(s,'(a7,i4.4,a6)') 'Station',i,'_M.dat'
          open(1000+i,file=s,form='unformatted',access='append')
          if(iCal.eq.1) then
              write(1000+i) (GP%t(j), j=1,GP%TotalTimeSteps+1)
          endif
          write(1000+i) (SP(i)%M(j), j=1,GP%TotalTimeSteps+1)
          close(1000+i)

          write(s,'(a7,i4.4,a6)') 'Station',i,'_N.dat'
          open(1000+i,file=s,form='unformatted',access='append')
          if(iCal.eq.1) then
              write(1000+i) (GP%t(j), j=1,GP%TotalTimeSteps+1)
          endif
          write(1000+i) (SP(i)%N(j), j=1,GP%TotalTimeSteps+1)
          close(1000+i)

          endif
      enddo
      endif

      end subroutine saveStationData



      subroutine saveSnapshot(GP, LP, iCal, iTimeStep, LocalData, LocalDataLength)

      use VariableDefination
      implicit NONE
      include "mpif.h"
      type(GlobalParameters)   ::  GP
      type(LayerParameters)    ::  LP(100)
      integer*4  ::  iCal, iTimeStep
      integer*4  ::  LocalDataLength
      real*8     ::  LocalData(LocalDataLength)
      integer*4  ::  irank, nsize, master
      integer*4  ::  ierror, istatus(MPI_STATUS_SIZE)
      integer*4  ::  iLayer, iNode, i, j
      integer*4  ::  nstartx, nendx, nstarty, nendy
      character(999)  ::  s
      real*8     ::  CPUTime1, CPUTime2

      irank = GP%irank; nsize = GP%nsizeTotal; master = GP%master
      call CPU_TIME(CPUTime1)

      do iLayer=1, GP%NumLayers
          if(irank.lt.LP(iLayer)%nsize.and.irank.ne.master) then
              nstartx = LP(iLayer)%PartitionInfo(irank+1,1)
              nendx   = LP(iLayer)%PartitionInfo(irank+1,2)
              nstarty = LP(iLayer)%PartitionInfo(irank+1,3)
              nendy   = LP(iLayer)%PartitionInfo(irank+1,4)
              if(iTimeStep.eq.0) then
                  do j=nstarty, nendy
                  do i=nstartx, nendx
                      LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1) = LP(iLayer)%H(1,i,j)
                  enddo
                  enddo
              else
                  do j=nstarty, nendy
                  do i=nstartx, nendx
                      LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1) = LP(iLayer)%H(2,i,j)
                  enddo
                  enddo
              endif
              call MPI_SEND(LocalData, (nendx-nstartx+1)*(nendy-nstarty+1), &
                  MPI_DOUBLE_PRECISION,master,2015,MPI_COMM_WORLD,ierror)
          endif
          if(irank.eq.master) then
              do iNode=0, LP(iLayer)%nsize-1
                  if(iNode.ne.master) then
                      nstartx = LP(iLayer)%PartitionInfo(iNode+1,1)
                      nendx   = LP(iLayer)%PartitionInfo(iNode+1,2)
                      nstarty = LP(iLayer)%PartitionInfo(iNode+1,3)
                      nendy   = LP(iLayer)%PartitionInfo(iNode+1,4)
                      call MPI_RECV(LocalData,(nendx-nstartx+1)*(nendy-nstarty+1), &
                          MPI_DOUBLE_PRECISION,iNode,2015,MPI_COMM_WORLD,istatus,ierror)
                      if(iTimeStep.eq.0) then
                          do j=nstarty, nendy
                          do i=nstartx, nendx
                              LP(iLayer)%H(1,i,j) = LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1)
                          enddo
                          enddo
                      else
                          do j=nstarty, nendy
                          do i=nstartx, nendx
                              LP(iLayer)%H(2,i,j) = LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1)
                          enddo
                          enddo
                      endif
                  endif
              enddo
              write(s,'(a2,i2.2,a1,i6.6,a4)') 'z_',iLayer,'_',iTimeStep,'.dat'
              open(22,file=s,form='unformatted')
              write(22) LP(iLayer)%NX, LP(iLayer)%NY
              if(iTimeStep.eq.0) then
                  do j=1, LP(iLayer)%NY
                      write(22) (LP(iLayer)%H(1,i,j), i=1,LP(iLayer)%NX)
                  enddo
              else
                  do j=1, LP(iLayer)%NY
                      write(22) (LP(iLayer)%H(2,i,j), i=1,LP(iLayer)%NX)
                  enddo
              endif
              close(22)
          endif
      enddo


      if(GP%SaveFlux.eq.1) then

      do iLayer=1, GP%NumLayers
          if(irank.lt.LP(iLayer)%nsize.and.irank.ne.master) then
              nstartx = LP(iLayer)%PartitionInfo(irank+1,1)
              nendx   = MIN(LP(iLayer)%PartitionInfo(irank+1,2),LP(iLayer)%NX-1)
              nstarty = LP(iLayer)%PartitionInfo(irank+1,3)
              nendy   = LP(iLayer)%PartitionInfo(irank+1,4)
              if(iTimeStep.eq.0) then
                  do j=nstarty, nendy
                  do i=nstartx, nendx
                      LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1) = LP(iLayer)%M(1,i,j)
                  enddo
                  enddo
              else
                  do j=nstarty, nendy
                  do i=nstartx, nendx
                      LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1) = LP(iLayer)%M(2,i,j)
                  enddo
                  enddo
              endif
              call MPI_SEND(LocalData, (nendx-nstartx+1)*(nendy-nstarty+1),&
                  MPI_DOUBLE_PRECISION,master,2015,MPI_COMM_WORLD,ierror)
          endif
          if(irank.eq.master) then
              do iNode=0, LP(iLayer)%nsize-1
                  if(iNode.ne.master) then
                      nstartx = LP(iLayer)%PartitionInfo(iNode+1,1)
                      nendx   = MIN(LP(iLayer)%PartitionInfo(iNode+1,2),LP(iLayer)%NX-1)
                      nstarty = LP(iLayer)%PartitionInfo(iNode+1,3)
                      nendy   = LP(iLayer)%PartitionInfo(iNode+1,4)
                      call MPI_RECV(LocalData,(nendx-nstartx+1)*(nendy-nstarty+1), &
                          MPI_DOUBLE_PRECISION,iNode,2015,MPI_COMM_WORLD,istatus,ierror)
                      if(iTimeStep.eq.0) then
                          do j=nstarty, nendy
                          do i=nstartx, nendx
                              LP(iLayer)%M(1,i,j) = LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1)
                          enddo
                          enddo
                      else
                          do j=nstarty, nendy
                          do i=nstartx, nendx
                              LP(iLayer)%M(2,i,j) = LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1)
                          enddo
                          enddo
                      endif
                  endif
              enddo
              write(s,'(a2,i2.2,a1,i6.6,a4)') 'M_',iLayer,'_',iTimeStep,'.dat'
              open(22,file=s,form='unformatted')
              write(22) LP(iLayer)%NX-1, LP(iLayer)%NY
              if(iTimeStep.eq.0) then
                  do j=1, LP(iLayer)%NY
                      write(22) (LP(iLayer)%M(1,i,j), i=1,LP(iLayer)%NX-1)
                  enddo
              else
                  do j=1, LP(iLayer)%NY
                      write(22) (LP(iLayer)%M(2,i,j), i=1,LP(iLayer)%NX-1)
                  enddo
              endif
              close(22)
          endif
      enddo

      do iLayer=1, GP%NumLayers
          if(irank.lt.LP(iLayer)%nsize.and.irank.ne.master) then
              nstartx = LP(iLayer)%PartitionInfo(irank+1,1)
              nendx   = LP(iLayer)%PartitionInfo(irank+1,2)
              nstarty = LP(iLayer)%PartitionInfo(irank+1,3)
              nendy   = MIN(LP(iLayer)%PartitionInfo(irank+1,4),LP(iLayer)%NY-1)
              if(iTimeStep.eq.0) then
                  do j=nstarty, nendy
                  do i=nstartx, nendx
                      LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1) = LP(iLayer)%N(1,i,j)
                  enddo
                  enddo
              else
                  do j=nstarty, nendy
                  do i=nstartx, nendx
                      LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1) = LP(iLayer)%N(2,i,j)
                  enddo
                  enddo
              endif
              call MPI_SEND(LocalData, (nendx-nstartx+1)*(nendy-nstarty+1),&
                  MPI_DOUBLE_PRECISION,master,2015,MPI_COMM_WORLD,ierror)
          endif
          if(irank.eq.master) then
              do iNode=0, LP(iLayer)%nsize-1
                  if(iNode.ne.master) then
                      nstartx = LP(iLayer)%PartitionInfo(iNode+1,1)
                      nendx   = LP(iLayer)%PartitionInfo(iNode+1,2)
                      nstarty = LP(iLayer)%PartitionInfo(iNode+1,3)
                      nendy   = MIN(LP(iLayer)%PartitionInfo(iNode+1,4),LP(iLayer)%NY-1)
                      call MPI_RECV(LocalData,(nendx-nstartx+1)*(nendy-nstarty+1), &
                          MPI_DOUBLE_PRECISION,iNode,2015,MPI_COMM_WORLD,istatus,ierror)
                      if(iTimeStep.eq.0) then
                          do j=nstarty, nendy
                          do i=nstartx, nendx
                              LP(iLayer)%N(1,i,j) = LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1)
                          enddo
                          enddo
                      else
                          do j=nstarty, nendy
                          do i=nstartx, nendx
                              LP(iLayer)%N(2,i,j) = LocalData((i-nstartx)*(nendy-nstarty+1)+j-nstarty+1)
                          enddo
                          enddo
                      endif
                  endif
              enddo
              write(s,'(a2,i2.2,a1,i6.6,a4)') 'N_',iLayer,'_',iTimeStep,'.dat'
              open(22,file=s,form='unformatted')
              write(22) LP(iLayer)%NX, LP(iLayer)%NY-1
              if(iTimeStep.eq.0) then
                  do j=1, LP(iLayer)%NY-1
                      write(22) (LP(iLayer)%N(1,i,j), i=1,LP(iLayer)%NX)
                  enddo
              else
                  do j=1, LP(iLayer)%NY-1
                      write(22) (LP(iLayer)%N(2,i,j), i=1,LP(iLayer)%NX)
                  enddo
              endif
              close(22)
          endif
      enddo

      endif
      call CPU_TIME(CPUTime2); GP%CPUTime(irank+1,5) = GP%CPUTime(irank+1,5)+CPUTime2-CPUTime1

      end subroutine saveSnapshot



      subroutine getLayerBoundaryValueFromParent(GP, LP, iLayer, iCal, iTimeStep, LocalData, LocalDataLength)

      use VariableDefination
      implicit NONE
      include "mpif.h"
      type(GlobalParameters)   ::  GP
      type(LayerParameters)    ::  LP(100)
      integer*4  ::  iLayer, iCal, iTimeStep
      integer*4  ::  LocalDataLength
      real*8     ::  LocalData(LocalDataLength)
      integer*4  ::  irank, nsize, master
      integer*4  ::  ierror, istatus(MPI_STATUS_SIZE)
      integer*4  ::  pLayer
      integer*4  ::  iPC, iNodeFrom, iNodeTo, iBoundary, nstart, nend
      integer*4  ::  istart, iend, jstart, jend, i, j, iFlag
      real*8     ::  x, y, h
      real*8     ::  CPUTime1, CPUTime2

      irank = GP%irank; nsize = GP%nsizeTotal; master = GP%master
      pLayer = LP(iLayer)%Parent   ! interpolate data from pLayer, parent layer of ilayer
      LP(iLayer)%hxmin = 0.0d0; LP(iLayer)%hxmax = 0.0d0;
      LP(iLayer)%hymin = 0.0d0; LP(iLayer)%hymax = 0.0d0;
      call CPU_TIME(CPUTime1)

      do iPC=1, LP(iLayer)%ParentToChildSendRecvCount
          iNodeFrom  = LP(iLayer)%ParentToChildSendRecv(iPC,1)
          iNodeTo    = LP(iLayer)%ParentToChildSendRecv(iPC,2)
          iBoundary  = LP(iLayer)%ParentToChildSendRecv(iPC,3)
          nstart     = LP(iLayer)%ParentToChildSendRecv(iPC,4)
          nend       = LP(iLayer)%ParentToChildSendRecv(iPC,5)
          if(iBoundary.eq.1) then
              istart = 1; iend = istart; jstart = nstart; jend = nend
          elseif(iBoundary.eq.2) then
              istart = LP(iLayer)%NX; iend = istart; jstart = nstart; jend = nend
          elseif(iBoundary.eq.3) then
              istart = nstart; iend = nend; jstart = 1; jend = jstart
          elseif(iBoundary.eq.4) then
              istart = nstart; iend = nend; jstart = LP(iLayer)%NY; jend = jstart
          endif

          if(irank.eq.iNodeTo) then
              do j=jstart, jend
              do i=istart, iend
                  if(iBoundary.eq.1) LP(iLayer)%hxmin(1,j) = LP(iLayer)%H(1,i,j)
                  if(iBoundary.eq.2) LP(iLayer)%hxmax(1,j) = LP(iLayer)%H(1,i,j)
                  if(iBoundary.eq.3) LP(iLayer)%hymin(1,i) = LP(iLayer)%H(1,i,j)
                  if(iBoundary.eq.4) LP(iLayer)%hymax(1,i) = LP(iLayer)%H(1,i,j)
              enddo
              enddo
          endif

          if(irank.eq.iNodeFrom) then
              do j=jstart, jend
              do i=istart, iend
                  x = LP(iLayer)%X(i); y = LP(iLayer)%Y(j)
                  iFlag = 1
                  call interpData(GP, LP, pLayer, iFlag, x, y, h)
                  if(iBoundary.eq.1) LP(iLayer)%hxmin(2,j) = h
                  if(iBoundary.eq.2) LP(iLayer)%hxmax(2,j) = h
                  if(iBoundary.eq.3) LP(iLayer)%hymin(2,i) = h
                  if(iBoundary.eq.4) LP(iLayer)%hymax(2,i) = h
              enddo
              enddo
          endif

          if(irank.eq.iNodeFrom.and.iNodeFrom.ne.iNodeTo) then
              do i=nstart, nend
                  if(iBoundary.eq.1) then
                      LocalData(i-nstart+1) = LP(iLayer)%hxmin(2,i)
                  elseif(iBoundary.eq.2) then
                      LocalData(i-nstart+1) = LP(iLayer)%hxmax(2,i)
                  elseif(iBoundary.eq.3) then
                      LocalData(i-nstart+1) = LP(iLayer)%hymin(2,i)
                  elseif(iBoundary.eq.4) then
                      LocalData(i-nstart+1) = LP(iLayer)%hymax(2,i)
                  endif
              enddo
              call MPI_SEND(LocalData,nend-nstart+1,MPI_DOUBLE_PRECISION,iNodeTo,2015,MPI_COMM_WORLD,ierror)
          endif

          if(irank.eq.iNodeTo.and.iNodeFrom.ne.iNodeTo) then
              call MPI_RECV(LocalData,nend-nstart+1,MPI_DOUBLE_PRECISION, &
                  iNodeFrom,2015,MPI_COMM_WORLD,istatus,ierror)
              do i=nstart, nend
                  if(iBoundary.eq.1) then
                      LP(iLayer)%hxmin(2,i) = LocalData(i-nstart+1)
                  elseif(iBoundary.eq.2) then
                      LP(iLayer)%hxmax(2,i) = LocalData(i-nstart+1)
                  elseif(iBoundary.eq.3) then
                      LP(iLayer)%hymin(2,i) = LocalData(i-nstart+1)
                  elseif(iBoundary.eq.4) then
                      LP(iLayer)%hymax(2,i) = LocalData(i-nstart+1)
                  endif
              enddo
          endif

      enddo
      call CPU_TIME(CPUTime2); GP%CPUTime(irank+1,4) = GP%CPUTime(irank+1,4)+CPUTime2-CPUTime1

      end subroutine getLayerBoundaryValueFromParent



      subroutine getLayerBoundaryValueAtFineTimeStep(GP, LP, iLayer, iStep)

      use VariableDefination
      implicit NONE
      type(GlobalParameters)   ::  GP
      type(LayerParameters)    ::  LP(100)
      integer*4  ::  iLayer, iStep
      integer*4  ::  nSteps, irank, i, j
      integer*4  ::  nstartx, nendx, nstarty, nendy
      real*8     ::  t1, t2, h1, h2, t, deltat

      irank = GP%irank
      if(irank.lt.LP(iLayer)%nsize) then

      nSteps = LP(iLayer)%nStepsPerTimeStep
      t1 = 0.0d0; t2 = GP%dt; deltat = 1.0d0/(t2-t1)
      t = (iStep-1)*LP(iLayer)%dt+LP(iLayer)%dt*LP(iLayer)%dtratio
      do j=1, LP(iLayer)%NY
      do i=1, LP(iLayer)%NX
          if(j.eq.1.or.j.eq.LP(iLayer)%NY.or.i.eq.1.or.i.eq.LP(iLayer)%NX) then
              nstartx = LP(iLayer)%PartitionInfo(irank+1,1)
              nendx   = LP(iLayer)%PartitionInfo(irank+1,2)
              nstarty = LP(iLayer)%PartitionInfo(irank+1,3)
              nendy   = LP(iLayer)%PartitionInfo(irank+1,4)
              if(i+GP%nRowBoundary.ge.nstartx.and.i-GP%nRowBoundary.le.nendx.and. &
                  j+GP%nRowBoundary.ge.nstarty.and.j-GP%nRowBoundary.le.nendy) then
              if(LP(iLayer)%Z(i,j).gt.GP%WaterDepthLimit) then
                  if(i.eq.1) then
                      h1 = LP(iLayer)%hxmin(1,j); h2 = LP(iLayer)%hxmin(2,j)
                  else if(i.eq.LP(iLayer)%NX) then
                      h1 = LP(iLayer)%hxmax(1,j); h2 = LP(iLayer)%hxmax(2,j)
                  else if(j.eq.1) then
                      h1 = LP(iLayer)%hymin(1,i); h2 = LP(iLayer)%hymin(2,i)
                  else
                      h1 = LP(iLayer)%hymax(1,i); h2 = LP(iLayer)%hymax(2,i)
                  endif
                  LP(iLayer)%H(2,i,j) = (h2-h1)*deltat*(t-t1)+h1
              else
                  LP(iLayer)%H(2,i,j) = 0.0d0
              endif
              endif
          endif
      enddo
      enddo

      endif

      end subroutine getLayerBoundaryValueAtFineTimeStep



      subroutine interpLayerValueToGlobalTime(GP, LP, iLayer, iStep)

      use VariableDefination
      implicit NONE
      type(GlobalParameters)   ::  GP
      type(LayerParameters)    ::  LP(100)
      integer*4  ::  iLayer, iStep
      integer*4  ::  i, j, irank, nSteps
      integer*4  ::  nstartx, nendx, nstarty, nendy
      real*8     ::  t1, t2, h1, h2, t, m1, m2, n1, n2

      irank = GP%irank
      if(irank.lt.LP(iLayer)%nsize) then

      nSteps = LP(iLayer)%nStepsPerTimeStep
      t1 = (nSteps-1)*LP(iLayer)%dt; t2 = t1+LP(iLayer)%dt; t = GP%dt
      do j=1, LP(iLayer)%NY
          do i=1, LP(iLayer)%NX
              nstartx = MAX(LP(iLayer)%PartitionInfo(irank+1,1)-1,1)
              nendx   = MIN(LP(iLayer)%PartitionInfo(irank+1,2)+1,LP(iLayer)%NX)
              nstarty = MAX(LP(iLayer)%PartitionInfo(irank+1,3)-1,1)
              nendy   = MAX(LP(iLayer)%PartitionInfo(irank+1,4)+1,LP(iLayer)%NY)
              if(i.ge.nstartx.and.i.le.nendx.and.j.ge.nstarty.and.j.le.nendy) then
                  h1 = LP(iLayer)%H(1,i,j); h2 = LP(iLayer)%H(2,i,j)
                  LP(iLayer)%H(2,i,j) = (h2-h1)/(t2-t1)*(t-t1)+h1
                  if(i.ne.LP(iLayer)%NX) then
                      m1 = LP(iLayer)%M(1,i,j); m2 = LP(iLayer)%M(2,i,j)
                      LP(iLayer)%M(2,i,j) = (m2-m1)/(t2-t1)*(t-t1)+m1
                  endif
                  if(j.ne.LP(iLayer)%NY) then
                      n1 = LP(iLayer)%N(1,i,j); n2 = LP(iLayer)%N(2,i,j)
                      LP(iLayer)%N(2,i,j) = (n2-n1)/(t2-t1)*(t-t1)+n1
                  endif
              endif
          enddo
      enddo

      endif
      
      end subroutine interpLayerValueToGlobalTime



      subroutine feedbackToParentLayer(GP, LP, iLayer, iCal, iTimeStep, LocalData, LocalDataLength)

      use VariableDefination
      implicit NONE
      include "mpif.h"
      type(GlobalParameters)   ::  GP
      type(LayerParameters)    ::  LP(100)
      integer*4  ::  iLayer, iCal, iTimeStep
      integer*4  ::  LocalDataLength
      real*8     ::  LocalData(LocalDataLength)
      integer*4  ::  irank, nsize, master
      integer*4  ::  ierror, istatus(MPI_STATUS_SIZE)
      integer*4  ::  pLayer
      integer*4  ::  iPC, iNodeFrom, iNodeTo, iHMN, istart, iend, jstart, jend
      integer*4  ::  i, j, iFlag
      real*8     ::  x, y, h
      real*8     ::  CPUTime1, CPUTime2

      irank = GP%irank; nsize = GP%nsizeTotal; master = GP%master
      pLayer = LP(iLayer)%Parent   ! interpolate data at iLayer, feedback to pLayer
      call CPU_TIME(CPUTime1)

      do iPC=1, LP(iLayer)%ChildToParentSendRecvCount
          iNodeFrom  = LP(iLayer)%ChildToParentSendRecv(iPC,1)
          iNodeTo    = LP(iLayer)%ChildToParentSendRecv(iPC,2)
          iHMN       = LP(iLayer)%ChildToParentSendRecv(iPC,3)
          istart     = LP(iLayer)%ChildToParentSendRecv(iPC,4)
          iend       = LP(iLayer)%ChildToParentSendRecv(iPC,5)
          jstart     = LP(iLayer)%ChildToParentSendRecv(iPC,6)
          jend       = LP(iLayer)%ChildToParentSendRecv(iPC,7)

          if(irank.eq.iNodeFrom) then
              do j=jstart, jend
              do i=istart, iend
              if(iNodeFrom.eq.iNodeTo.and.LP(pLayer)%Z(i,j).le.GP%WaterDepthLimit) then
                  if(iHMN.eq.1) then
                      LP(pLayer)%H(2,i,j) = 0.0d0
                  elseif(iHMN.eq.2) then
                      LP(pLayer)%M(2,i,j) = 0.0d0
                  else if(iHMN.eq.3) then
                      LP(pLayer)%N(2,i,j) = 0.0d0
                  endif
              else
                  x = LP(pLayer)%X(i); y = LP(pLayer)%Y(j)
                  if(iHMN.eq.2) x = x+0.5*LP(pLayer)%dx
                  if(iHMN.eq.3) y = y+0.5*LP(pLayer)%dy
                  iFlag = iHMN
                  call interpData(GP, LP, iLayer, iFlag, x, y, h)
                  if(iHMN.eq.1) then
                      LP(pLayer)%H(2,i,j) = h
                  elseif(iHMN.eq.2) then
                      LP(pLayer)%M(2,i,j) = h
                  else if(iHMN.eq.3) then
                      LP(pLayer)%N(2,i,j) = h
                  endif
              endif
              enddo
              enddo
          endif

          if(irank.eq.iNodeFrom.and.iNodeFrom.ne.iNodeTo) then
              do j=jstart, jend
              do i=istart, iend
                  if(iHMN.eq.1) then
                      LocalData((i-istart)*(jend-jstart+1)+j-jstart+1) = LP(pLayer)%H(2,i,j)
                  elseif(iHMN.eq.2) then
                      LocalData((i-istart)*(jend-jstart+1)+j-jstart+1) = LP(pLayer)%M(2,i,j)
                  else if(iHMN.eq.3) then
                      LocalData((i-istart)*(jend-jstart+1)+j-jstart+1) = LP(pLayer)%N(2,i,j)
                  endif
              enddo
              enddo
              call MPI_SEND(LocalData,(iend-istart+1)*(jend-jstart+1), &
                  MPI_DOUBLE_PRECISION,iNodeTo,2015,MPI_COMM_WORLD,ierror)
          endif

          if(irank.eq.iNodeTo.and.iNodeFrom.ne.iNodeTo) then
              call MPI_RECV(LocalData,(iend-istart+1)*(jend-jstart+1), &
                  MPI_DOUBLE_PRECISION,iNodeFrom,2015,MPI_COMM_WORLD,istatus,ierror)
              do j=jstart, jend
              do i=istart, iend
              if(LP(pLayer)%Z(i,j).gt.GP%WaterDepthLimit) then
                  if(iHMN.eq.1) then
                      LP(pLayer)%H(2,i,j) = LocalData((i-istart)*(jend-jstart+1)+j-jstart+1)
                  elseif(iHMN.eq.2) then
                      LP(pLayer)%M(2,i,j) = LocalData((i-istart)*(jend-jstart+1)+j-jstart+1)
                  else if(iHMN.eq.3) then
                      LP(pLayer)%N(2,i,j) = LocalData((i-istart)*(jend-jstart+1)+j-jstart+1)
                  endif
              endif
              enddo
              enddo
          endif

      enddo
      call CPU_TIME(CPUTime2); GP%CPUTime(irank+1,4) = GP%CPUTime(irank+1,4)+CPUTime2-CPUTime1

      end subroutine feedbackToParentLayer
