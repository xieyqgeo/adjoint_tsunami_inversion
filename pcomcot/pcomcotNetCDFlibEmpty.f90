
      subroutine getBathymetryDataSizeNetCDF(GP, LP, iLayer)
      
      use VariableDefination
      implicit NONE
      include "mpif.h"
      type(GlobalParameters)   ::  GP
      type(LayerParameters)    ::  LP(GP%NumLayers)
      integer*4                ::  iLayer, ierror

      write(*,*) 'ERROR: pcomcot compiled without netcdf lib. can''t read .nf files.'
      call MPI_ABORT(MPI_COMM_WORLD, ierror)

      end subroutine getBathymetryDataSizeNetCDF



      subroutine readBathymetryNetCDF(GP, LP, iLayer)
      
      use VariableDefination
      implicit NONE
      include "mpif.h"
      type(GlobalParameters)   ::  GP
      type(LayerParameters)    ::  LP(GP%NumLayers)
      integer*4                ::  iLayer, ierror

      write(*,*) 'ERROR: pcomcot compiled without netcdf lib. can''t read .nf files.'
      call MPI_ABORT(MPI_COMM_WORLD, ierror)

      end subroutine readBathymetryNetCDF



      subroutine readInitialConditionNetCDF(GP, LP)
      
      use VariableDefination
      implicit NONE
      include "mpif.h"
      type(GlobalParameters)   ::  GP
      type(LayerParameters)    ::  LP(GP%NumLayers)
      integer*4                ::  ierror

      write(*,*) 'ERROR: pcomcot compiled without netcdf lib. can''t read .nf files.'
      call MPI_ABORT(MPI_COMM_WORLD, ierror)

      end subroutine readInitialConditionNetCDF
