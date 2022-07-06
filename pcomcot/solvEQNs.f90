      
      subroutine mass(GP, LP, iLayer, iStep)

      use VariableDefination
      implicit NONE
      type(GlobalParameters)    ::  GP
      type(LayerParameters)     ::  LP(100)
      integer*4  ::  iLayer, iStep, irank
      integer*4  ::  nstartx, nendx, nstarty, nendy
      integer*4  ::  istart, iend, jstart, jend
      integer*4  ::  i, j, iBoundary, HasThisBoundary
      real*8     ::  m1, m2, n1, n2, tmp
      real*8     ::  CPUTime1, CPUTime2

      irank = GP%irank
      call CPU_TIME(CPUTime1)

      if(irank.lt.LP(iLayer)%nsize) then

      nstartx = LP(iLayer)%PartitionInfo(irank+1,1)
      nendx   = LP(iLayer)%PartitionInfo(irank+1,2)
      nstarty = LP(iLayer)%PartitionInfo(irank+1,3)
      nendy   = LP(iLayer)%PartitionInfo(irank+1,4)
      istart = MAX(2,nstartx); iend = MIN(LP(iLayer)%NX-1,nendx)
      jstart = MAX(2,nstarty); jend = MIN(LP(iLayer)%NY-1,nendy)
      do j = jstart,jend
      do i = istart,iend
          if(LP(iLayer)%Z(i,j).gt.GP%WaterDepthLimit) then
              m1 = LP(iLayer)%M(1, i-1, j);  m2 = LP(iLayer)%M(1, i, j)
              n1 = LP(iLayer)%N(1, i, j-1);  n2 = LP(iLayer)%N(1, i, j)
              if(LP(iLayer)%Z(i-1,j).le.GP%WaterDepthLimit) m1 = -m2
              if(LP(iLayer)%Z(i+1,j).le.GP%WaterDepthLimit) m2 = -m1
              if(LP(iLayer)%Z(i,j-1).le.GP%WaterDepthLimit) n1 = -n2
              if(LP(iLayer)%Z(i,j+1).le.GP%WaterDepthLimit) n2 = -n1
              if(GP%CoordinatesType.eq.1) then
                  if(iStep.ne.LP(iLayer)%nStepsPerTimeStep) then
                      LP(iLayer)%H(2,i,j) = LP(iLayer)%H(1,i,j) -  &
                          LP(iLayer)%CC1*(m2-m1) - LP(iLayer)%CC2*(n2-n1)
                  else
                      LP(iLayer)%H(2,i,j) = LP(iLayer)%H(1,i,j) -  &
                          LP(iLayer)%CC1*LP(iLayer)%dtratio*(m2-m1) - LP(iLayer)%CC2*LP(iLayer)%dtratio*(n2-n1)
                  endif
                  if(GP%CoordinatesType.eq.1) then
                      !nonlinear terms
                  endif
              else
                  if(iStep.ne.LP(iLayer)%nStepsPerTimeStep) then
                      LP(iLayer)%H(2,i,j) = LP(iLayer)%H(1,i,j) - LP(iLayer)%CS1(j)*(m2-m1) - &
                          LP(iLayer)%CS2(j)*(n2*LP(iLayer)%CSY(j)-n1*LP(iLayer)%CSY(j-1))
                  else
                      LP(iLayer)%H(2,i,j) = LP(iLayer)%H(1,i,j) - LP(iLayer)%CS1(j)*LP(iLayer)%dtratio*(m2-m1) - &
                          LP(iLayer)%CS2(j)*LP(iLayer)%dtratio*(n2*LP(iLayer)%CSY(j)-n1*LP(iLayer)%CSY(j-1))
                  endif
                  if(GP%CoordinatesType.eq.1) then
                      !nonlinear terms
                  endif
              endif
          else
              LP(iLayer)%H(2,i,j) = 0.0d0
          endif
      enddo
      enddo

      !//// boundaries ////!
      if(LP(iLayer)%Level.eq.1) then
          do iBoundary = 1,4
              HasThisBoundary = 0
              if(iBoundary.eq.1.and.nstartx.eq.1) HasThisBoundary = 1
              if(iBoundary.eq.2.and.nendx.eq.LP(iLayer)%NX) HasThisBoundary = 1
              if(iBoundary.eq.3.and.nstarty.eq.1) HasThisBoundary = 1
              if(iBoundary.eq.4.and.nendy.eq.LP(iLayer)%NY) HasThisBoundary = 1
              if(HasThisBoundary.eq.1) then
                  if(iBoundary.eq.1) then
                      istart = 1; iend = 1; jstart = nstarty; jend = nendy
                  else if(iBoundary.eq.2) then
                      istart = LP(iLayer)%NX; iend = LP(iLayer)%NX; jstart = nstarty; jend = nendy
                  else if(iBoundary.eq.3) then
                      istart = MAX(2,nstartx); iend = MIN(LP(iLayer)%NX-1,nendx); jstart = 1; jend = 1
                  else
                      istart = MAX(2,nstartx); iend = MIN(LP(iLayer)%NX-1,nendx);
                      jstart = LP(iLayer)%NY;  jend = LP(iLayer)%NY
                  endif
                  do j = jstart,jend
                  do i = istart,iend
                      if(LP(iLayer)%Z(i,j).gt.GP%WaterDepthLimit) then
                          if(GP%BoundaryConditionType.eq.1) then ! Wall boundary
                              if(iBoundary.eq.1.or.iBoundary.eq.2) then
                                  if(iBoundary.eq.1) then
                                      m2 = LP(iLayer)%M(1,i,j); m1 = -m2
                                  else
                                      m1 = LP(iLayer)%M(1,i-1,j);  m2 = -m1
                                  endif
                                  if(j.ne.1) n1 = LP(iLayer)%N(1,i,j-1)
                                  if(j.ne.LP(iLayer)%NY) n2 = LP(iLayer)%N(1,i,j)
                                  if(j.eq.1) n1 = -n2
                                  if(j.eq.LP(iLayer)%NY) n2 = -n1
                                  if(j.ne.1.and.j.ne.LP(iLayer)%NY) then
                                      if(LP(iLayer)%Z(i,j-1).le.GP%WaterDepthLimit) n1 = -n2
                                      if(LP(iLayer)%Z(i,j+1).le.GP%WaterDepthLimit) n2 = -n1
                                  endif
                              else if(iBoundary.eq.3.or.iBoundary.eq.4) then
                                  if(iBoundary.eq.3) then
                                      n2 = LP(iLayer)%N(1,i,j); n1 = -n2
                                  else
                                      n1 = LP(iLayer)%N(1,i,j-1); n2 = -n1
                                  endif
                                  m1 = LP(iLayer)%M(1,i-1,j)
                                  m2 = LP(iLayer)%M(1,i,j)
                                  if(LP(iLayer)%Z(i-1,j).le.GP%WaterDepthLimit) m1 = -m2
                                  if(LP(iLayer)%Z(i+1,j).le.GP%WaterDepthLimit) m2 = -m1
                              endif
                              if(GP%CoordinatesType.eq.1) then
                                  if(iStep.ne.LP(iLayer)%nStepsPerTimeStep) then
                                      LP(iLayer)%H(2,i,j) = LP(iLayer)%H(1,i,j) -  &
                                          LP(iLayer)%CC1*(m2-m1) - LP(iLayer)%CC2*(n2-n1)
                                  else
                                      LP(iLayer)%H(2,i,j) = LP(iLayer)%H(1,i,j) -  &
                                          LP(iLayer)%CC1*LP(iLayer)%dtratio*(m2-m1) - &
                                          LP(iLayer)%CC2*LP(iLayer)%dtratio*(n2-n1)
                                  endif
                              else
                                  if(j.ne.1) tmp = LP(iLayer)%CSY(j-1)
                                  if(j.eq.1) tmp = COS((LP(iLayer)%Y(j)-LP(iLayer)%dy)*GP%PI/180.0)
                                  if(iStep.ne.LP(iLayer)%nStepsPerTimeStep) then
                                      LP(iLayer)%H(2,i,j) = LP(iLayer)%H(1,i,j) -  &
                                          LP(iLayer)%CS1(j)*(m2-m1) - LP(iLayer)%CS2(j)*(n2*LP(iLayer)%CSY(j)-n1*tmp)
                                  else
                                      LP(iLayer)%H(2,i,j) = LP(iLayer)%H(1,i,j) - &
                                          LP(iLayer)%CS1(j)*LP(iLayer)%dtratio*(m2-m1) - &
                                          LP(iLayer)%CS2(j)*(n2*LP(iLayer)%CSY(j)-n1*tmp)
                                  endif
                              endif
                          else ! Other boundary conditions
                          endif
                      else
                          LP(iLayer)%H(2,i,j) = 0.0d0
                      endif
                  enddo
                  enddo
              endif
          enddo
      endif !if: boundaries

      endif
      call CPU_TIME(CPUTime2); GP%CPUTime(irank+1,3) = GP%CPUTime(irank+1,3)+CPUTime2-CPUTime1

      end subroutine mass



      subroutine moment(GP, LP, iLayer, iStep)

      use VariableDefination
      implicit NONE
      type(GlobalParameters)    ::  GP
      type(LayerParameters)     ::  LP(100)
      integer*4  ::  iLayer, iStep, irank
      integer*4  ::  nstartx, nendx, nstarty, nendy
      integer*4  ::  istart, iend, jstart, jend
      integer*4  ::  i, j
      real*8     ::  m1, m2, n1, n2, tmp, flux
      real*8     ::  CPUTime1, CPUTime2

      irank = GP%irank
      call CPU_TIME(CPUTime1)

      if(irank.lt.LP(iLayer)%nsize) then

      nstartx = LP(iLayer)%PartitionInfo(irank+1,1)
      nendx   = LP(iLayer)%PartitionInfo(irank+1,2)
      nstarty = LP(iLayer)%PartitionInfo(irank+1,3)
      nendy   = LP(iLayer)%PartitionInfo(irank+1,4)
      istart = MAX(1,nstartx); iend = MIN(LP(iLayer)%NX-1,nendx)
      jstart = MAX(1,nstarty); jend = MIN(LP(iLayer)%NY,nendy)
      do j= jstart,jend
      do i = istart,iend
          if(LP(iLayer)%Z(i,j).gt.GP%WaterDepthLimit) then
              if(GP%CoordinatesType.eq.1) then
                  if(iStep.ne.LP(iLayer)%nStepsPerTimeStep) then
                      LP(iLayer)%M(2,i,j) = LP(iLayer)%M(1,i,j) -                     &
                          LP(iLayer)%CC3*(LP(iLayer)%Z(i,j)+LP(iLayer)%Z(i+1,j))*     &
                          (LP(iLayer)%H(2,i+1,j)-LP(iLayer)%H(2,i,j))
                  else
                      LP(iLayer)%M(2,i,j) = LP(iLayer)%M(1,i,j) - &
                          LP(iLayer)%CC3*LP(iLayer)%dtratio*(LP(iLayer)%Z(i,j)+LP(iLayer)%Z(i+1,j))* &
                          (LP(iLayer)%H(2,i+1,j)-LP(iLayer)%H(2,i,j))
                  endif
              else
                  if(iStep.ne.LP(iLayer)%nStepsPerTimeStep) then
                      LP(iLayer)%M(2,i,j) = LP(iLayer)%M(1,i,j) -                     &
                          LP(iLayer)%CS3(j)*(LP(iLayer)%Z(i,j)+LP(iLayer)%Z(i+1,j))*  &
                          (LP(iLayer)%H(2,i+1,j)-LP(iLayer)%H(2,i,j))
                  else
                      LP(iLayer)%M(2,i,j) = LP(iLayer)%M(1,i,j) - &
                          LP(iLayer)%CS3(j)*LP(iLayer)%dtratio*(LP(iLayer)%Z(i,j)+LP(iLayer)%Z(i+1,j))* &
                          (LP(iLayer)%H(2,i+1,j)-LP(iLayer)%H(2,i,j))
                  endif
                  !/// Coriolis force  ///!
                  if(j.eq.1) then
                      flux = 0.5*(LP(iLayer)%N(1,i,j)+LP(iLayer)%N(1,i+1,j))
                  else if(j.eq.LP(iLayer)%NY) then
                      flux = 0.5*(LP(iLayer)%N(1,i,j-1)+LP(iLayer)%N(1,i+1,j-1))
                  else
                      flux=0.25*(LP(iLayer)%N(1,i,j)+LP(iLayer)%N(1,i+1,j)+ &
                          LP(iLayer)%N(1,i,j-1)+LP(iLayer)%N(1,i+1,j-1))
                  endif
                  if(iStep.ne.LP(iLayer)%nStepsPerTimeStep) then
                      LP(iLayer)%M(2,i,j) = LP(iLayer)%M(2,i,j)+LP(iLayer)%CPX(j)*flux
                  else
                      LP(iLayer)%M(2,i,j) = LP(iLayer)%M(2,i,j)+LP(iLayer)%CPX(j)*LP(iLayer)%dtratio*flux
                  endif
              endif
          else
              LP(iLayer)%M(2,i,j) = 0.0d0
          endif
      enddo
      enddo
      istart = MAX(1,nstartx); iend = MIN(LP(iLayer)%NX,nendx)
      jstart = MAX(1,nstarty); jend = MIN(LP(iLayer)%NY-1,nendy)
      do j = jstart,jend
      do i = istart,iend
          if(LP(iLayer)%Z(i,j).gt.GP%WaterDepthLimit) then
              if(GP%CoordinatesType.eq.1) then
                  if(iStep.ne.LP(iLayer)%nStepsPerTimeStep) then
                      LP(iLayer)%N(2,i,j) = LP(iLayer)%N(1,i,j) -                     &
                          LP(iLayer)%CC4*(LP(iLayer)%Z(i,j)+LP(iLayer)%Z(i,j+1))*     &
                          (LP(iLayer)%H(2,i,j+1)-LP(iLayer)%H(2,i,j))
                  else
                      LP(iLayer)%N(2,i,j) = LP(iLayer)%N(1,i,j) -                     &
                          LP(iLayer)%CC4*LP(iLayer)%dtratio*(LP(iLayer)%Z(i,j)+LP(iLayer)%Z(i,j+1))*     &
                          (LP(iLayer)%H(2,i,j+1)-LP(iLayer)%H(2,i,j))
                  endif
              else
                  if(iStep.ne.LP(iLayer)%nStepsPerTimeStep) then
                      LP(iLayer)%N(2,i,j) = LP(iLayer)%N(1,i,j) -                     &
                          LP(iLayer)%CS4*(LP(iLayer)%Z(i,j)+LP(iLayer)%Z(i,j+1))*  &
                          (LP(iLayer)%H(2,i,j+1)-LP(iLayer)%H(2,i,j))
                  else
                      LP(iLayer)%N(2,i,j) = LP(iLayer)%N(1,i,j) -                     &
                          LP(iLayer)%CS4*LP(iLayer)%dtratio*(LP(iLayer)%Z(i,j)+LP(iLayer)%Z(i,j+1))*  &
                          (LP(iLayer)%H(2,i,j+1)-LP(iLayer)%H(2,i,j))
                  endif
                  !/// Coriolis force  ///!
                  if(i.eq.1) then
                      flux = 0.5*(LP(iLayer)%M(1,i,j)+LP(iLayer)%M(1,i,j+1))
                  else if(i.eq.LP(iLayer)%NX) then
                      flux = 0.5*(LP(iLayer)%M(1,i-1,j)+LP(iLayer)%M(1,i-1,j+1))
                  else
                      flux = 0.25*(LP(iLayer)%M(1,i,j)+LP(iLayer)%M(1,i,j+1)+ &
                          LP(iLayer)%M(1,i-1,j)+LP(iLayer)%M(1,i-1,j+1))
                  endif
                  if(iStep.ne.LP(iLayer)%nStepsPerTimeStep) then
                      LP(iLayer)%N(2,i,j) = LP(iLayer)%N(2,i,j)-LP(iLayer)%CPY(j)*flux
                  else
                      LP(iLayer)%N(2,i,j) = LP(iLayer)%N(2,i,j)-LP(iLayer)%CPY(j)*LP(iLayer)%dtratio*flux
                  endif
              endif
          else
              LP(iLayer)%N(2,i,j) = 0.0d0
          endif
      enddo
      enddo

      endif
      call CPU_TIME(CPUTime2); GP%CPUTime(irank+1,3) = GP%CPUTime(irank+1,3)+CPUTime2-CPUTime1

      end subroutine moment

