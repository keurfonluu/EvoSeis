!=======================================================================
! EvoWave
!-----------------------------------------------------------------------
! EvoWave is a program that generates synthetic seismograms.
!
! Created by
!     Keurfon Luu <keurfon.luu@mines-paristech.fr>
!     MINES ParisTech - Centre de Géosciences
!     PSL - Research University
!
! Last updated
!     2016-12-14 16:51
!=======================================================================

program evowave

  use omp_lib
  use forlab, only: IPRE, RPRE
  use evoseis, only: Network, Cluster, Vel3D

  implicit none

  ! Input parameters
  !==================


  ! Local variables
  !=================

  print *, "Il est vivaaaaant!"

  print *
  stop
end program evowave
