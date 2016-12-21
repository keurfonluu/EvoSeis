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
!     2016-12-21 15:11
!=======================================================================

program evowave

  use omp_lib
  use forlab, only: IPRE, RPRE, disp, loadtxt, savebin
  use evoseis, only: Network, Source, Vel3D, Wavelet, Waveform

  implicit none

  ! Input parameters
  !==================
  integer(kind = IPRE) :: form, unit
  real(kind = RPRE) :: qp, qs, rho, length, fs, duration, f0, par

  ! Local variables
  !=================
  integer(kind = IPRE) :: nrcv, k
  real(kind = RPRE) :: M0
  real(kind = RPRE), dimension(:), allocatable :: MT
  real(kind = RPRE), dimension(:,:), allocatable :: stations
  real(kind = RPRE), dimension(:,:,:), allocatable :: velp3d, vels3d, csgather
  type(Vel3D) :: velp, vels
  type(Network) :: ntwrk
  type(Source) :: src
  type(Wavelet) :: wave
  type(Waveform) :: seis

  ! Initialize velocity models
  !============================
  allocate(velp3d(1, 1, 1), vels3d(1, 1, 1))
  velp3d = 3000.
  vels3d = 2000.
  qp = 25.
  qs = 15.
  rho = 2300.

  velp = Vel3D(velp3d, 1., 1., 1., qp)
  vels = Vel3D(vels3d, 1., 1., 1., qs)

  ! Initialize network
  !====================
  stations = loadtxt("./example_waveform/line.txt", 3)
  nrcv = size(stations, 1)
  ntwrk = Network(stations(:,1), stations(:,2), stations(:,3))

  ! Initialize source
  !===================
  MT = [ 1., 0., 0., 1., 0., 1. ] / 3.
  M0 = 1e9
  src = Source(550., 0., 550., seismic_moment = M0, moment_tensor = MT)

  ! Initialize wavelet
  !====================
  fs = 1000.
  length = 2.048
  form = 0
  duration = 0.005
  f0 = 100.

  select case(form)
  case(0)
    par = duration
  case(1)
    par = f0
  end select

  wave = Wavelet(length, fs, form)
  call wave % model(par)

  ! Initialize Waveform
  !=====================
  unit = 1
  seis = Waveform(length, fs, unit)

  ! Common shot gather
  !====================
  allocate(csgather(int(fs * length), nrcv, 3))
  do k = 1, nrcv
    call seis % model(velp, vels, rho, src, ntwrk % stations(k), wave, par)
    call seis % lqt2zne()
    ! call seis % zne2lqt()
    csgather(:,k,1) = seis % comp1
    csgather(:,k,2) = seis % comp2
    csgather(:,k,3) = seis % comp3
  end do

  ! Save
  !======
  call system("rm -rf ./example_waveform/comp*.bin")
  call savebin("./example_waveform/comp1.bin", csgather(:,:,1))
  call savebin("./example_waveform/comp2.bin", csgather(:,:,2))
  call savebin("./example_waveform/comp3.bin", csgather(:,:,3))

  print *
  stop
end program evowave
