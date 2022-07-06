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
      integer*4                  ::  iLayer, iLayerLevel, iSta, iCal, iTimeStep, iStep
      real*8                     ::  CPUTime1, CPUTime2
      integer*4                  ::  LocalDataLength
      real*8, allocatable        ::  LocalData(:)

      call MPI_INIT(ierror)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nsize, ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD, irank, ierror)
      call CPU_TIME(GP%CPUTimeInitial); GP%CPUTime(irank+1,1) = 0.0d0
      master = 0; GP%irank = irank; GP%nsizeTotal = nsize; GP%master = master

      !########## read data on master node ##########!
      if(irank.eq.master) then
          call readConfig(GP)
          call checkFiles(GP)
          call getBathymetryDataSize(GP, LP)
          call determineLayerDependency(GP, LP)
          do iLayer = 1,GP%NumLayers
              ALLOCATE(LP(iLayer)%X(LP(iLayer)%NX))
              ALLOCATE(LP(iLayer)%Y(LP(iLayer)%NY))
              ALLOCATE(LP(iLayer)%Z(LP(iLayer)%NX,LP(iLayer)%NY))
          enddo
          call readBathymetry(GP, LP)
          call cflCheck(GP, LP)
          call readStations(GP, LP, SP)
          if(GP%InitialConditionType.eq.1) then
              call readFaultParameters(GP, FP)
          endif
          call partitionDomain(GP, LP, SP)
      endif

      !#### broadcast dimensions, allocate memoery ####!
      call bcastCommonInfo(GP, LP, SP, FP)
      ALLOCATE(GP%t(GP%TotalTimeSteps+1))
      do iLayer = 1,GP%NumLayers
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
      do iSta = 1,GP%NumStations
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

      !############ start computing ############!
      do iCal = 1,GP%nCalculations

      if(GP%ComputeGreen.eq.1) call modifyFaultParameters(GP, FP, iCal)
      do iLayer = 1,GP%NumLayers
          LP(iLayer)%H = 0.0d0; LP(iLayer)%M = 0.0d0; LP(iLayer)%N = 0.0d0
      enddo
      do iTimeStep=0, GP%TotalTimeSteps
          GP%t(iTimeStep+1) = iTimeStep*GP%dt
          call getInitialCondition(GP, LP, SP, FP, iCal, iTimeStep, LocalData, LocalDataLength)
          if(iTimeStep.eq.1) call CPU_TIME(GP%CPUTimeInitialCalculation)
          call screenOutput(GP, LP, iTimeStep)

          if(iTimeStep.ne.0) then
          do iLayerLevel = 1,GP%NumLayerLevels
          do iLayer = 1,GP%NumLayers
          if(LP(iLayer)%Level.eq.iLayerLevel) then

              if(iLayerLevel.gt.1) then
                  call getLayerBoundaryValueFromParent(GP, LP, iLayer, iCal, iTimeStep, LocalData, LocalDataLength)
              endif
              do iStep = 1,LP(iLayer)%nStepsPerTimeStep
                  if(iLayerLevel.gt.1) then
                      call getLayerBoundaryValueAtFineTimeStep(GP, LP, iLayer, iStep)
                  endif
                  call mass(GP, LP, iLayer, iStep)
                  call bcastComputeDomainBoundaryValue(GP, LP, iLayer, iCal, iTimeStep, LocalData, LocalDataLength)
                  call moment(GP, LP, iLayer, iStep)
                  call bcastComputeDomainBoundaryFlux(GP, LP, iLayer, iCal, iTimeStep, LocalData, LocalDataLength)
                  call updateComputedResults(GP, LP, iLayer)
              enddo

          endif
          enddo
          enddo
          endif

          if(iTimeStep.ne.0.and.GP%FeedbackToParentLayer.eq.1) then
          do iLayerLevel = GP%NumLayerLevels,2,-1
          do iLayer = 1,GP%NumLayers
          if(LP(iLayer)%Level.eq.iLayerLevel) then
              call feedbackToParentLayer(GP, LP, iLayer, iCal, iTimeStep, LocalData, LocalDataLength)
              call updateComputedResults(GP, LP, LP(iLayer)%Parent)
          endif
          enddo
          enddo
          endif

          call calculateStationData(GP, LP, SP, iCal, iTimeStep)
          if(MOD(iTimeStep,GP%NDTSaveData).eq.0.or.iTimeStep.eq.GP%TotalTimeSteps) then
              call saveSnapshot(GP, LP, iCal, iTimeStep, LocalData, LocalDataLength)
          endif
      enddo !do loop for iTimeStep

      call saveStationData(GP, LP, SP, iCal)
      enddo !do loop for iCal

      call CPU_TIME(GP%CPUTIME(irank+1,1))
      GP%CPUTime(irank+1,1) = GP%CPUTime(irank+1,1)-GP%CPUTimeInitial
      call screenOutput(GP, LP, GP%TotalTimeSteps+100)

      call MPI_FINALIZE(ierror)

      end program pcomcot

