!************************************************************************!
!                                                                        !
!         Develped by Chao An, Cornell University (2012~2015)            !
!         Report bugs to: ca298@cornell.edu                              !
!         Last updated: 2015/7/21                                        !
!                                                                        !
!************************************************************************!

      program pcomcot
      use VariableDefination
      
      implicit NONE
      include 'mpif.h'
      integer*4  ::  irank, nsize, master
      integer*4  ::  ierror, istatus(MPI_STATUS_SIZE)
      
      type(GlobalParameters)     ::  GP
      type(LayerParameters)      ::  LP(100)
      type(StationParameters)    ::  SP(999)
      type(FaultParameters)      ::  FP(999)
      integer*4                  ::  iLayer, iLayerLevel, iSta, iCal
      integer*4                  ::  iTimeStep, iStep
      real*8                     ::  CPUTime1, CPUTime2
      integer*4                  ::  LocalDataLength
      real*8, allocatable        ::  LocalData(:)
      character*999              ::  s

      !****** time reversal ******!
      real*8                     ::  TR_R, tmp
      integer*4                  ::  TR_I, TR_J, TR_MinI, TR_MaxI
      integer*4                  ::  TR_MinJ, TR_MaxJ

      call MPI_INIT(ierror)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nsize, ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD, irank, ierror)
      call CPU_TIME(GP%CPUTimeInitial); GP%CPUTime(irank+1,1) = 0.0d0
      master = 0; GP%irank=irank; GP%nsizeTotal=nsize; GP%master=master

      !########## read data on master node ##########!
      if(irank.eq.master) then
          call readConfig(GP)
          call checkFiles(GP)
          call getBathymetryDataSize(GP, LP)
          call determineLayerDependency(GP, LP)
          do iLayer=1, GP%NumLayers
              ALLOCATE(LP(iLayer)%X(LP(iLayer)%NX))
              ALLOCATE(LP(iLayer)%Y(LP(iLayer)%NY))
              ALLOCATE(LP(iLayer)%Z(LP(iLayer)%NX,LP(iLayer)%NY))
          enddo
          call readBathymetry(GP, LP)
! time reversal; remove land in bathymetry !
! do TR_J = 1,LP(1)%NY
! do TR_I = 2,LP(1)%NX
! if(LP(1)%Z(TR_I,TR_J).le.0.0d0) then
! LP(1)%Z(TR_I,TR_J) = LP(1)%Z(TR_I-1,TR_J)+10.0d0
! endif
! enddo
! enddo

! time reversal; remove land in bathymetry !
          call cflCheck(GP, LP)
          call readStations(GP, LP, SP)
          if(GP%InitialConditionType.eq.1) then
              call readFaultParameters(GP, FP)
          endif
          call partitionDomain(GP, LP, SP)
      endif

      !#### broadcast dimensions, allocate memoery ####!
      !#### bacst common info of time reversal ####!
      call bcastCommonInfo(GP, LP, SP, FP)
      ALLOCATE(GP%t(GP%TotalTimeSteps+1))
      do iLayer=1, GP%NumLayers
          if(irank.ne.master) then
              ALLOCATE(LP(iLayer)%X(LP(iLayer)%NX))
              ALLOCATE(LP(iLayer)%Y(LP(iLayer)%NY))
              ALLOCATE(LP(iLayer)%Z(LP(iLayer)%NX,LP(iLayer)%NY))
          endif
          ALLOCATE(LP(iLayer)%H(2, LP(iLayer)%NX, LP(iLayer)%NY))
          ALLOCATE(LP(iLayer)%M(2, LP(iLayer)%NX-1, LP(iLayer)%NY))
          ALLOCATE(LP(iLayer)%N(2, LP(iLayer)%NX, LP(iLayer)%NY-1))
          ALLOCATE(LP(iLayer)%hxmin(2, LP(iLayer)%NY))
          ALLOCATE(LP(iLayer)%hxmax(2, LP(iLayer)%NY))
          ALLOCATE(LP(iLayer)%hymin(2, LP(iLayer)%NX))
          ALLOCATE(LP(iLayer)%hymax(2, LP(iLayer)%NX))
          if(GP%CoordinatesType.eq.0) then
              ALLOCATE(LP(iLayer)%CPX(LP(iLayer)%NY))
              ALLOCATE(LP(iLayer)%CPY(LP(iLayer)%NY-1))
              ALLOCATE(LP(iLayer)%CSY(LP(iLayer)%NY))
              ALLOCATE(LP(iLayer)%CS1(LP(iLayer)%NY))
              ALLOCATE(LP(iLayer)%CS2(LP(iLayer)%NY))
              ALLOCATE(LP(iLayer)%CS3(LP(iLayer)%NY))
          endif
      enddo
      do iSta=1, GP%NumStations
          ALLOCATE(SP(iSta)%H(GP%TotalTimeSteps+1))
          ALLOCATE(SP(iSta)%M(GP%TotalTimeSteps+1))
          ALLOCATE(SP(iSta)%N(GP%TotalTimeSteps+1))
      enddo
      LocalDataLength = (GP%MaxNX+10)*(GP%MaxNY+10)
      ALLOCATE(LocalData(LocalDataLength))
      
      !#### broadcast, compute arrarys after memories allocated ####!
      call bcastBathymetry(GP, LP, SP, FP, LocalData, LocalDataLength)
      call computeParameters(GP, LP, SP, FP)
      if(irank.eq.master) call removeStationFiles(GP, SP)
      call CPU_TIME(GP%CPUTime(irank+1,2))
      GP%CPUTime(irank+1,2) = GP%CPUTime(irank+1,2)-GP%CPUTimeInitial

      ! time reversal; for code simplicity, data for stations are calcuated on every comouting node !
      !***** read station data for time reversal ******!
      if(irank.eq.master) then
          open(1,file='i_TimeReversal',form='formatted')
          read(1,'(a)') SP(1)%TR_DataPath
          read(1,*) SP(1)%TR_NR
          read(1,*) SP(1)%TR_PointSourceType
          read(1,*) SP(1)%TR_AmpCorrection
          read(1,*) SP(1)%TR_EpiX,SP(1)%TR_EpiY
          if(GP%CoordinatesType.eq.0) then
              if((GP%StartEastWest.eq.1).and.(SP(1)%TR_EPiX.lt.0)) &  
                       SP(1)%TR_EpiX = SP(1)%TR_EpiX+360
              if((GP%StartEastWest.eq.2).and.(SP(1)%TR_EPiX.gt.180))  &
                       SP(1)%TR_EpiX = SP(1)%TR_EpiX-360
          endif
          read(1,*) (SP(iSta)%TR_Weight, iSta=1,GP%NumStations)
          close(1)
          do iSta = 1,GP%NumStations
              SP(iSta)%TR_NR = SP(1)%TR_NR
              SP(iSta)%TR_PointSourceType = SP(1)%TR_PointSourceType
              SP(iSta)%TR_AmpCorrection = SP(1)%TR_AmpCorrection
              SP(iSta)%TR_EpiX = SP(1)%TR_EpiX
              SP(iSta)%TR_EpiY = SP(1)%TR_EpiY
          enddo
          do iSta = 1,GP%NumStations
              write(s,'(a7,i4.4,a4)') 'Station',iSta,'.dat'
              open(1,file=TRIM(ADJUSTL(SP(1)%TR_DataPath))  &
                     //TRIM(ADJUSTL(s)),form='unformatted')
              read(1) iCal,iCal,iCal
              read(1) (SP(iSta)%H(iCal), iCal=1,GP%TotalTimeSteps+1) !Time
              read(1) (SP(iSta)%H(iCal), iCal=1,GP%TotalTimeSteps+1) !Water Height
              SP(iSta)%TR_HMax = 0.0
              do iCal = 1,GP%TotalTimeSteps+1
               SP(iSta)%TR_HMax=MAX(SP(iSta)%TR_HMax,  &
                     ABS(SP(iSta)%H(iCal)))
              enddo
              if(SP(iSta)%TR_AmpCorrection.gt.0.and.SP(iSta)%TR_Weight.gt.0.0) then
                  do iCal=1,GP%TotalTimeSteps+1 !Amplitude correction
                      if(GP%CoordinatesType.eq.1) then
                          TR_R = SQRT((SP(iSta)%X-SP(iSta)%TR_EpiX)**2 &
                                   +(SP(iSta)%Y-SP(iSta)%TR_EpiY)**2)
                      else
                          call distazOnSphere(SP(iSta)%Y,SP(iSta)%X,  &
                                   SP(iSta)%TR_EpiY,SP(iSta)%TR_EpiX, &
                              TR_R,tmp,tmp)
                      endif
                      if(SP(iSta)%TR_AmpCorrection.eq.1) then
                          SP(iSta)%H(iCal) = SP(iSta)%H(iCal)  &
                             *SQRT(2*TR_R)*SP(iSta)%TR_Weight !Amplitude correction place
                      elseif(SP(iSta)%TR_AmpCorrection.eq.2) then
                          iLayer = SP(iSta)%nLayer
                          SP(iSta)%IX = NINT((SP(iSta)%X  &
                                   -LP(iLayer)%X(1))/LP(iLayer)%dx)+1
                          SP(iSta)%IY = NINT((SP(iSta)%Y-   &
                                   LP(iLayer)%Y(1))/LP(iLayer)%dy)+1
                          SP(iSta)%H(iCal) = SP(iSta)%H(iCal)*   &
                      SQRT(1*TR_R)*SP(iSta)%TR_Weight/SP(iSta)%TR_HMax &
                              /LP(iLayer)%Z(SP(iSta)%IX,SP(iSta)%IY)
                      endif
                  enddo
              endif
              close(1)
          enddo
      endif
      do iSta = 1,GP%NumStations
          call MPI_BCAST(SP(iSta)%TR_NR,1,MPI_INTEGER,master,  &
                MPI_COMM_WORLD,ierror)
          call MPI_BCAST(SP(iSta)%TR_PointSourceType,1,MPI_INTEGER,  &
                master,MPI_COMM_WORLD,ierror)
          call MPI_BCAST(SP(iSta)%TR_AmpCorrection,1,MPI_INTEGER,  &
                master,MPI_COMM_WORLD,ierror)
          call MPI_BCAST(SP(iSta)%TR_EpiX,1,MPI_DOUBLE_PRECISION,  &
                master,MPI_COMM_WORLD,ierror)
          call MPI_BCAST(SP(iSta)%TR_EpiY,1,MPI_DOUBLE_PRECISION,  &
                master,MPI_COMM_WORLD,ierror)
          call MPI_BCAST(SP(iSta)%TR_Weight,1,MPI_DOUBLE_PRECISION,  &
                master,MPI_COMM_WORLD,ierror)
          call MPI_BCAST(SP(iSta)%TR_HMax,1,MPI_DOUBLE_PRECISION,  &
                master,MPI_COMM_WORLD,ierror)
          call MPI_BCAST(SP(iSta)%H,GP%TotalTimeSteps+1,  &
                MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
      enddo
      do iSta=1, GP%NumStations
          iLayer = SP(iSta)%nLayer
          SP(iSta)%IX=NINT((SP(iSta)%X-LP(iLayer)%X(1))/LP(iLayer)%dx)+1
          SP(iSta)%IY=NINT((SP(iSta)%Y-LP(iLayer)%Y(1))/LP(iLayer)%dy)+1
          TR_R = SP(iSta)%TR_NR*LP(iLayer)%dt*  & 
                  SQRT(9.8*LP(iLayer)%Z(SP(iSta)%IX,SP(iSta)%IY))
          if(GP%CoordinatesType.eq.0) TR_R=TR_R/6371000.0*180.0/3.1415927
          SP(iSta)%TR_NPoints = 0
          TR_MinI = MAX(SP(iSta)%IX-(NINT(TR_R/LP(iLayer)%dx)+1),1)
          TR_MaxI = MIN(SP(iSta)%IX+(NINT(TR_R/LP(iLayer)%dx)+1),  &
                      LP(iLayer)%NX)
          TR_MinJ = MAX(SP(iSta)%IY-(NINT(TR_R/LP(iLayer)%dy)+1),1)
          TR_MaxJ = MIN(SP(iSta)%IY+(NINT(TR_R/LP(iLayer)%dy)+1),  &
                      LP(iLayer)%NY)
          do TR_I = TR_MinI,TR_MaxI
          do TR_J = TR_MinJ,TR_MaxJ
              if(SP(iSta)%TR_PointSourceType.eq.1) then
                  if((LP(iLayer)%X(TR_I)-LP(iLayer)%X(SP(iSta)%IX))**2 &
                    +(LP(iLayer)%Y(TR_J)-LP(iLayer)%Y(SP(iSta)%IY))**2 &
                     .le. TR_R**2) &
                      SP(iSta)%TR_NPoints = SP(iSta)%TR_NPoints+1
              elseif(SP(iSta)%TR_PointSourceType.eq.2) then
                  if((LP(iLayer)%X(TR_I)-LP(iLayer)%X(SP(iSta)%IX))**2 &
                    +(LP(iLayer)%Y(TR_J)-LP(iLayer)%Y(SP(iSta)%IY))**2 &
                     .le. 4*TR_R**2) &
                      SP(iSta)%TR_NPoints = SP(iSta)%TR_NPoints+1
              endif
          enddo
          enddo
          ALLOCATE(SP(iSta)%TR_Points(SP(iSta)%TR_NPoints,2))
          ALLOCATE(SP(iSta)%TR_PointsValue(SP(iSta)%TR_NPoints,   &
                  GP%TotalTimeSteps+1))
          SP(iSta)%TR_NPoints = 0
          do TR_I = TR_MinI,TR_MaxI
          do TR_J = TR_MinJ,TR_MaxJ
              if(SP(iSta)%TR_PointSourceType.eq.1) then
                  if((LP(iLayer)%X(TR_I)-LP(iLayer)%X(SP(iSta)%IX))**2 &
                   +(LP(iLayer)%Y(TR_J)-LP(iLayer)%Y(SP(iSta)%IY))**2 &
                   .le. TR_R**2) then
                      SP(iSta)%TR_NPoints = SP(iSta)%TR_NPoints+1
                      SP(iSta)%TR_Points(SP(iSta)%TR_NPoints,1) = TR_I
                      SP(iSta)%TR_Points(SP(iSta)%TR_NPoints,2) = TR_J
                      do iTimeStep = 1,GP%TotalTimeSteps+1
                          SP(iSta)%TR_PointsValue(  &
                            SP(iSta)%TR_NPoints,iTimeStep) = &
                            SP(iSta)%H(GP%TotalTimeSteps+2-iTimeStep)
                      enddo
                  endif
              elseif(SP(iSta)%TR_PointSourceType.eq.2) then
                  if((LP(iLayer)%X(TR_I)-LP(iLayer)%X(SP(iSta)%IX))**2 &
                    +(LP(iLayer)%Y(TR_J)-LP(iLayer)%Y(SP(iSta)%IY))**2 &
                    .le. 4*TR_R**2) then
                      SP(iSta)%TR_NPoints = SP(iSta)%TR_NPoints+1
                      SP(iSta)%TR_Points(SP(iSta)%TR_NPoints,1) = TR_I
                      SP(iSta)%TR_Points(SP(iSta)%TR_NPoints,2) = TR_J
                      do iTimeStep = 1,GP%TotalTimeSteps+1
                          SP(iSta)%TR_PointsValue(  &
                            SP(iSta)%TR_NPoints,iTimeStep) = &
                            SP(iSta)%H(GP%TotalTimeSteps+2-iTimeStep)* &
                            EXP(-((LP(iLayer)%X(TR_I)-  &
                            LP(iLayer)%X(SP(iSta)%IX))**2 &
                            +(LP(iLayer)%Y(TR_J)-   &
                            LP(iLayer)%Y(SP(iSta)%IY))**2)/TR_R**2)
                      enddo
                  endif
              endif
          enddo
          enddo
      enddo
      ! time reversal !

      !############ start computing ############!
      do iCal=1, GP%nCalculations

       if(GP%ComputeGreen.eq.1) call modifyFaultParameters(GP, FP, iCal)
       do iLayer=1, GP%NumLayers
          LP(iLayer)%H = 0.0d0; LP(iLayer)%M = 0.0d0; LP(iLayer)%N = 0.0d0
       enddo
       do iTimeStep=0, GP%TotalTimeSteps
          GP%t(iTimeStep+1) = iTimeStep*GP%dt
          call getInitialCondition(GP, LP, SP, FP, iCal, iTimeStep, &
                   LocalData, LocalDataLength)
          if(iTimeStep.eq.1) call CPU_TIME(GP%CPUTimeInitialCalculation)
          call screenOutput(GP, LP, iTimeStep)
! time reversal only works for single layer !
          do iSta = 1,GP%NumStations
           if(SP(iSta)%TR_Weight.gt.0.0.and.SP(iSta)%nNode.eq.irank) then
            do TR_I = 1,SP(iSta)%TR_NPoints
              LP(1)%H(1,SP(iSta)%TR_Points(TR_I,1),SP(iSta)%TR_Points( &
                    TR_I,2)) = SP(iSta)%TR_PointsValue(TR_I,iTimeStep+1)
              LP(1)%H(2,SP(iSta)%TR_Points(TR_I,1),SP(iSta)%TR_Points( &
                    TR_I,2)) = SP(iSta)%TR_PointsValue(TR_I,iTimeStep+1)
            enddo
           endif
          enddo
! time reversal !!

          if(iTimeStep.ne.0) then
           do iLayerLevel=1, GP%NumLayerLevels
            do iLayer=1, GP%NumLayers
             if(LP(iLayer)%Level.eq.iLayerLevel) then

              if(iLayerLevel.gt.1) then
                  call getLayerBoundaryValueFromParent(GP, LP, iLayer, &
                          iCal, iTimeStep, LocalData, LocalDataLength)
              endif
              do iStep=1, LP(iLayer)%nStepsPerTimeStep
                  if(iLayerLevel.gt.1) then
                      call getLayerBoundaryValueAtFineTimeStep(GP, LP, &
                            iLayer, iStep)
                  endif
                  call mass(GP, LP, iLayer, iStep)
! time reversal only works for single layer !
                  do iSta = 1,GP%NumStations
                   if(SP(iSta)%TR_Weight.gt.0.0.and.SP(iSta)%nNode.eq. &
                     irank) then
                    do TR_I = 1,SP(iSta)%TR_NPoints
                      LP(iLayer)%H(2,SP(iSta)%TR_Points(TR_I,1), &
                          SP(iSta)%TR_Points(TR_I,2)) = &
                          SP(iSta)%TR_PointsValue(TR_I,iTimeStep+1)
                    enddo
                   endif
                  enddo
! time reversal !!
                  call bcastComputeDomainBoundaryValue(GP, LP, iLayer, &
                        iCal, iTimeStep, LocalData, LocalDataLength)
                  call moment(GP, LP, iLayer, iStep)
                  call bcastComputeDomainBoundaryFlux(GP, LP, iLayer, &
                        iCal, iTimeStep, LocalData, LocalDataLength)
                  call updateComputedResults(GP, LP, iLayer)
              enddo

             endif
            enddo
           enddo
          endif

          if(iTimeStep.ne.0.and.GP%FeedbackToParentLayer.eq.1) then
           do iLayerLevel=GP%NumLayerLevels, 2, -1
            do iLayer=1, GP%NumLayers
             if(LP(iLayer)%Level.eq.iLayerLevel) then
              call feedbackToParentLayer(GP, LP, iLayer, iCal, &
                   iTimeStep, LocalData, LocalDataLength)
              call updateComputedResults(GP, LP, LP(iLayer)%Parent)
             endif
            enddo
           enddo
          endif

          call calculateStationData(GP, LP, SP, iCal, iTimeStep)
          if(iTimeStep .gt. GP%DTSaveSTART-0.01) then
            if(MOD(int(iTimeStep-GP%DTSaveSTART),GP%NDTSaveData).eq.0) then
              call saveSnapshot(GP, LP, iCal, iTimeStep, & 
                    LocalData, LocalDataLength)
            endif
          endif
       enddo !do loop for iTimeStep

       call saveStationData(GP, LP, SP, iCal)
      enddo !do loop for iCal

      call CPU_TIME(GP%CPUTIME(irank+1,1))
      GP%CPUTime(irank+1,1) = GP%CPUTime(irank+1,1)-GP%CPUTimeInitial
      call screenOutput(GP, LP, GP%TotalTimeSteps+100)

      call MPI_FINALIZE(ierror)

      end program pcomcot



      subroutine distazOnSphere(lat1, lon1, lat2, lon2, dis, az, baz)
      implicit NONE
      real*8   ::  lat1, lon1, lat2, lon2, dis, az, baz
      real*8   ::  phi1, phi2, theta1, theta2, sind, cosd
      real*8   ::  DegToRad

      DegToRad = 3.14159265359d0/180.0
      phi1 = lat1*DegToRad; phi2 = lat2*DegToRad
      theta1 = lon1*DegToRad; theta2 = lon2*DegToRad
      sind = SQRT((COS(theta1)*COS(phi1)-COS(theta2)*COS(phi2))**2 + &
             (SIN(theta1)*COS(phi1)-SIN(theta2)*COS(phi2))**2 + &
             (SIN(phi1)-SIN(phi2))**2) * &
             SQRT((COS(theta1)*COS(phi1)+COS(theta2)*COS(phi2))**2 + &
             (SIN(theta1)*COS(phi1)+SIN(theta2)*COS(phi2))**2 + &
             (SIN(phi1)+SIN(phi2))**2)*0.5
      cosd = SIN(phi1)*SIN(phi2)+COS(phi1)*COS(phi2)*COS(theta1-theta2)
      dis = ATAN2(sind, cosd)/DegToRad
      az = ATAN2((COS(phi2)*SIN(theta2-theta1))/sind, &
          (SIN(phi2)-SIN(phi1)*cosd)/COS(phi1)/sind)/DegToRad
      if(az.lt.0.0) az = az+360.0
      baz = ATAN2(COS(phi1)*SIN(theta2-theta1)/sind, &
          (SIN(phi1)-SIN(phi2)*cosd)/COS(phi2)/sind)/DegToRad
      baz = 360-baz
      if(baz.ge.360.0) baz = baz-360.0

      end subroutine distazOnSphere
