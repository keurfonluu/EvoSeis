!=======================================================================
! EvoLoc
!-----------------------------------------------------------------------
! EvoLoc is a program to locate seismic events given arrival times and
! velocity models.
!
! Created by
!     Keurfon Luu <keurfon.luu@mines-paristech.fr>
!     MINES ParisTech - Centre de Géosciences
!     PSL - Research University
!
! Last updated
!     2016-11-27 23:09
!=======================================================================

program evoloc
  
  use omp_lib
  use forlab, only: IPRE, RPRE, CLEN, File, rng, loadtxt, loadbin, zeros, &
    horzcat, tic, toc, num2str
  use evoseis, only: Network, Cluster, Vel3D
  
  implicit none
  
  ! Input parameters
  !==================
  integer(kind = IPRE) :: nz, nx, ny, year, month, day, zcoord, &
    do_bayesh, coord_unit, pdf, num_threads, nsweep
  integer(kind = IPRE) :: method, np, itermax
  real(kind = RPRE) :: dz, dx, dy, spmin, spmax, ssmin, ssmax, &
    west, south, east, north
  real(kind = RPRE) :: omega, phi, F, CR
  integer(kind = IPRE), dimension(:), allocatable :: emitted_signal
  character(len = CLEN) :: velp_filename, vels_filename, stat_filename, &
    tobs_filename, sigma_filename, output_dirname
  character(len = :), allocatable :: python, input_filename
  
  ! Local variables
  !=================
  integer(kind = IPRE) :: i, j, k, nstat, nsh = 0, nobs, nsrc
  real(kind = 8) :: toc_time
  real(kind = RPRE) :: zs = 0, par(2)
  real(kind = RPRE), dimension(:,:), allocatable :: stations, tobs, sigma, &
    sigmap, sigmas, tobsp, tobss
  real(kind = RPRE), dimension(:,:,:), allocatable :: tmp, velp3d, vels3d
  type(Network) :: ntwrk
  type(Cluster) :: clust
  type(Vel3D) :: velp, vels
  
  ! Read command line arguments
  !=============================
  call command_arguments()

  ! Open GUI
  !==========
  call system(python // " -W ignore bin/gui/evolocGUI.py")
  
  ! Check emitted signal
  !======================
  emitted_signal = loadbin(".evotmp/signal_emit", 4, 1)
  if ( emitted_signal(1) .eq. 0 ) then
    call system("rm -rf .evotmp")   ! Remove .evotmp directory
    stop
  end if
  
  ! Read input parameters
  !=======================
  input_filename = ".evotmp/parameters"
  call read_input()
  call system("rm -rf .evotmp")     ! Remove .evotmp directory
  
  ! Initialize OpenMP
  !===================
  call omp_set_num_threads(num_threads)   ! OpenMP number of threads

  ! Initialize random number seed
  !===============================
  call rng()
  
  ! Create output directory
  !=========================
  call system("rm -rf " // trim(output_dirname))
  call system("mkdir -p " // trim(output_dirname))
  
  ! Stations file
  !===============
  stations = loadtxt(stat_filename, 3)
  nstat = size(stations, 1)
  
  ! Observed arrival times file
  !=============================
  tobs = loadtxt(tobs_filename, 2)
  nobs = size(tobs, 1)
  nsrc = nobs / nstat
  
  ! Observed uncertainties file
  !=============================
  if ( do_bayesh .eq. 0 ) then
    sigma = loadtxt(sigma_filename, 2)
  else
    stop "Hierarchical bayesian formulation not implemented yet!"
  end if
  
  ! Create Network and Cluster objects (array)
  !============================================
  ntwrk = Network(stations(:,1), stations(:,2), stations(:,3), tobs, sigma)
  clust = Cluster(nsrc, EID = [ ( i, i = 1, nsrc ) ])
  
  ! Apply vertical shift (if in elevation)
  !========================================
  if ( zcoord .eq. 1 ) call ntwrk % zshift(zs, nsh, nz, dz)
  
  ! Velocity files
  !================
  tmp = loadbin(velp_filename, 4, nz-nsh-1, nx-1, ny-1)
  velp3d = zeros(nz-1, nx-1, ny-1)
  do i = 1, nsh
    velp3d(i,:,:) = tmp(1,:,:)
  end do
  velp3d(nsh+1:,:,:) = tmp
  
  tmp = loadbin(vels_filename, 4, nz-nsh-1, nx-1, ny-1)
  vels3d = zeros(nz-1, nx-1, ny-1)
  do i = 1, nsh
    vels3d(i,:,:) = tmp(1,:,:)
  end do
  vels3d(nsh+1:,:,:) = tmp
  
  ! P and S velocity models
  !=========================
  velp = Vel3D(velp3d, dx, dy, dz)
  vels = Vel3D(vels3d, dx, dy, dz)
  
  ! Deallocate useless arrays
  !===========================
  deallocate(tmp, stations, tobs)
  
  ! Convert degrees to UTM
  !========================
  if (coord_unit .eq. 1) call ntwrk % convert_to_utm(south, west)
  
  ! Use an Eikonal solver to compute traveltime grids
  !===================================================
  print *
  
  write(*, fmt = "(A)", advance = "no") " Computing traveltime grids ..."
  call tic()
  call ntwrk % compute_ttgrids(velp, vels)
  call toc(toc_time)
  write(*,*) "done (" // num2str(toc_time) // " seconds)"
  
  ! Set optimization method parameters
  !====================================
  select case(method)
    case(1)
      par = [ omega, phi ]
    case(2)
      par = [ CR, F ]
  end select
  
  ! Locate the sources within the cluster
  !=======================================
  write(*, fmt = "(A)", advance = "no") " Locating events ..."
  call tic()
  call clust % locate(velp, vels, ntwrk, method, np, itermax, par)
  call toc(toc_time)
  write(*,*) "done (" // num2str(toc_time) // " seconds)"
  
  ! Undo vertical shift
  !=====================
  if ( zs .gt. 0.0d0 ) clust % sources(:) % z = clust % sources(:) % z - zs
  
  ! Convert locations in UTM to latitude / longitude
  !==================================================
  if ( coord_unit .eq. 1 ) call clust % convert_to_latlon(south, west)
  
  ! Save location results
  !=======================
  call clust % savehypo(trim(output_dirname) // "location_results.txt", &
                        coord_unit, year, month, day)
  print *, "Location results saved in " // trim(output_dirname) // "location_results.txt"
  
  print *
  stop
  
contains

  !---------------------------------------------------------------------
  ! Subroutine read_input
  !---------------------------------------------------------------------
  subroutine read_input() 
    type(File) :: input_file

    input_file = File(999, input_filename)
    call input_file%open()
    read(input_file%unit,*); read(input_file%unit,*) nz, nx, ny
    read(input_file%unit,*); read(input_file%unit,*) dz, dx, dy
    read(input_file%unit,*); read(input_file%unit,*) velp_filename
    read(input_file%unit,*); read(input_file%unit,*) vels_filename
    read(input_file%unit,*); read(input_file%unit,*) stat_filename
    read(input_file%unit,*); read(input_file%unit,*) tobs_filename
    read(input_file%unit,*); read(input_file%unit,*) year, month, day
    read(input_file%unit,*); read(input_file%unit,*) zcoord
    read(input_file%unit,*); read(input_file%unit,*) do_bayesh
    read(input_file%unit,*); read(input_file%unit,*) sigma_filename
    read(input_file%unit,*); read(input_file%unit,*) spmin, spmax
    read(input_file%unit,*); read(input_file%unit,*) ssmin, ssmax
    read(input_file%unit,*); read(input_file%unit,*) coord_unit
    read(input_file%unit,*); read(input_file%unit,*) west, south, east, north
    read(input_file%unit,*); read(input_file%unit,*) method
    read(input_file%unit,*); read(input_file%unit,*) np
    read(input_file%unit,*); read(input_file%unit,*) itermax
    read(input_file%unit,*); read(input_file%unit,*) omega, phi
    read(input_file%unit,*); read(input_file%unit,*) CR, F
    read(input_file%unit,*); read(input_file%unit,*) output_dirname
    read(input_file%unit,*); read(input_file%unit,*) pdf
    read(input_file%unit,*); read(input_file%unit,*) num_threads
    call input_file%close()
    return
  end subroutine read_input
  
  !---------------------------------------------------------------------
  ! Subroutine command_arguments
  !---------------------------------------------------------------------
  subroutine command_arguments() 
    integer(kind = IPRE) :: i, nargin
    character(len = CLEN) :: argin
    character(len = :), allocatable :: argname, argval
    
    nsweep = 1
    python = "python"
    nargin = command_argument_count()
    do i = 1, nargin
      call get_command_argument(i, argin)
      call split_argument(argin, argname, argval)
      if (argname .eq. "py") then
        python = argval
      elseif (argname .eq. "nsweep") then
        read(argval, *) nsweep
      else
        print *, "Error: Unknown argument '" // trim(argname) // "'."
        stop
      end if
    end do
    return
  end subroutine
  
  !---------------------------------------------------------------------
  ! Subroutine split_argument
  !---------------------------------------------------------------------
  subroutine split_argument(argin, argname, argval)
    character(len = *), intent(in) :: argin
    character(len = :), allocatable, intent(out) :: argname, argval
    integer(kind = IPRE) :: idx
    
    idx = index(argin, "=")
    if (idx .ne. 0) then
      argname = trim(argin(:idx-1))
      argval = trim(argin(idx+1:))
    else
      print *, "Error: Missing '=' in argument '" // trim(argin) // "'."
      stop
    end if
    return
  end subroutine split_argument
  
end program evoloc
