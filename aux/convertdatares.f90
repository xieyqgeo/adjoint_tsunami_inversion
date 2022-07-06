!!    Input:     ./TsunamiData/Stations.ctl, ./TsunamiData/TsunamiData.dat
!!    Output:    Stationxxxx.dat

      implicit NONE
      integer*4               ::  MaxDataLength
      parameter(MaxDataLength=50000)
      real*8                  ::  t(MaxDataLength), f(MaxDataLength)
      real*8                  ::  t2(MaxDataLength), f2(MaxDataLength)
      real*8                  ::  tmin, tmax, dt, ttmp
      real*8                  ::  CutTime
      character*999           ::  StaName, stmp
      real*8                  ::  ndataa, nstaa
      integer*4               ::  i, ista, ndata, nsta

      open(3,file='data_res_matlab.bin',access='stream',form='unformatted')
      read(3) nstaa
      read(3) ndataa
      read(3) dt
      nsta=int(nstaa); ndata=int(ndataa)
      write(*,*),nsta,ndata,dt
      ndata = ndata + 1
      do i=1,ndata
          t2(i) = dt * (i - 1)
      enddo
      do ista=1,nsta
          f2(1) = 0
          read(3) (f2(i),i=2,ndata)
          i = 0
          write(stmp,'(a,i4.4,a)') 'Station',ista,'.dat'
          open(4,file=TRIM(stmp),form='unformatted')
          write(4) ndata,i,i
          write(4) (t2(i), i=1,ndata)
          write(4) (f2(i), i=1,ndata)
          close(4)

          !write(*,*) '-->'//TRIM(stmp)

      enddo
      close(3)


      end



