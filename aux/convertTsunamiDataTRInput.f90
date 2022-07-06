!!    Input:     ./Stations.ctl, ./station.dat2 ./cut_sta_para.ctl
!!    Output:    Stationxxxx.dat (bin)

      implicit NONE
      integer*4               ::  MaxDataLength
      parameter(MaxDataLength=100000)
      real*8                  ::  t(MaxDataLength), f(MaxDataLength)
      real*8                  ::  t2(MaxDataLength), f2(MaxDataLength)
      real*8                  ::  tmin, tmax, dt, ttmp
      real*8                  ::  CutTime
      character*999           ::  StaName, stmp
      integer*4               ::  ios, i, ndat, ndat2, nsta

      open(398,file='cut_sta_para.ctl',form='formatted')
      read(398,*),tmin,tmax,dt
      close(398)

      ndat2 = 0
      do
          ttmp = tmin+ndat2*dt
          if(ttmp.gt.tmax) exit
          ndat2 = ndat2+1
          t2(ndat2) = ttmp
      enddo

      open(1,file='./Stations.ctl',form='formatted')
      !open(2,file='./TsunamiData/StationsMaxTime.ctl',form='formatted')
      nsta = 0
      do
          read(1,'(a)',iostat=ios) stmp
          if(ios.ne.0) exit

          nsta = nsta+1
          !read(2,*) CutTime
          !CutTime = CutTime*60.0

          stmp = ADJUSTR(stmp)
          i = INDEX(stmp,' ',back=.true.)
          StaName = TRIM(stmp(i+1:))//'.dat2'

          open(3,file='./'//TRIM(StaName),form='formatted')
          ndat = 0
          do
              read(3,'(a)',iostat=ios) stmp
              if(ios.ne.0) exit
              ndat = ndat+1
              read(stmp,*) t(ndat),f(ndat)
              !if(t(ndat).gt.CutTime) f(ndat) = 0.0d0
          enddo
          close(3)

          call interpData(t,f,MaxDataLength,ndat,t2,f2,MaxDataLength,ndat2)

          i = 0
          write(stmp,'(a,i4.4,a)') 'Station',nsta,'.dat'
          open(4,file=TRIM(stmp),form='unformatted')
          write(4) ndat2,i,i
          write(4) (t2(i), i=1,ndat2)
          write(4) (f2(i), i=1,ndat2)
          close(4)

          write(*,*) TRIM(StaName)//'    -->    '//TRIM(stmp)

      enddo
      !close(2)
      close(1)


      end



      subroutine interpData(t,f,mlen,n,t2,f2,mlen2,n2)
      implicit NONE
      real*8     ::  t(mlen),f(mlen),t2(mlen2),f2(mlen2)
      integer*4  ::  mlen,n,mlen2,n2
      
      integer*4  ::  i, j, jleft, jright

      do i=1,n2
          if(t2(i).le.t(1)) then
              f2(i) = f(1)
          elseif(t2(i).ge.t(n)) then
              f2(i) = f(n)
          else
              j = n/2; jleft = 1; jright = n;
              do
                  if((t2(i).ge.t(j)).and.(t2(i).le.t(j+1))) exit
                  if(t2(i).lt.t(j)) then
                      jright = j; j = (j+jleft)/2
                  else
                      jleft = j; j = (j+jright)/2
                  endif
              enddo
              f2(i) = (f(j+1)-f(j))/(t(j+1)-t(j))*(t2(i)-t(j))+f(j)
          endif
      enddo

      end subroutine interpData

      subroutine distazOnSphere(lat1, lon1, lat2, lon2, dis, az, baz)
      implicit NONE
      real*4   ::  lat1, lon1, lat2, lon2, dis, az, baz
      real*4   ::  phi1, phi2, theta1, theta2, sind, cosd
      real*4   ::  DegToRad

      DegToRad = 3.14159265359/180.0
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
