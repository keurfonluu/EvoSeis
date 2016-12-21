!=======================================================================
! EvoSeis
!-----------------------------------------------------------------------
! EvoSeis provides a packages of functions for earthquake seismology
! using Evolutionary Algorithms as optimization methods.
!
! Created by
!     Keurfon Luu <keurfon.luu@mines-paristech.fr>
!     MINES ParisTech - Centre de Géosciences
!     PSL - Research University
!
! Last updated
!     2016-11-27 23:11
!
! Notes
!-----------------------------------------------------------------------
! All attributes are public to ease implementation: Fortran does not
! allow (yet?) accessing encapsulated attributes through methods.
!
! TODO
!-----------------------------------------------------------------------
! Add magnitude to object Source.
!=======================================================================

module evoseis

  use omp_lib
  use forlab, only: IPRE, RPRE, File, horzcat, interp3, linspace, &
    deg2utm, utm2deg, datenum, datestr
  use fftpack, only: fft, ifft
  use optimizers, only: optifit, optifunc, cpso, de

!=======================================================================
! Object Location
!=======================================================================

  type Location
    ! Coordinates (in meters)
    real(kind = RPRE), public :: x, y, z
  end type Location

  interface Location
    module procedure init_Location
  end interface Location
  public :: Location
  private :: init_Location

!=======================================================================
! Object Hypocenter
!=======================================================================

  type, extends(Location) :: Hypocenter
    ! Event ID
    integer(kind = IPRE), public :: eventid
    ! Hypocenter origin time (in seconds)
    real(kind = 8), public :: origin_time
  contains
    procedure, private, pass :: locate_hypo
    generic, public :: locate => locate_hypo
  end type Hypocenter

  interface Hypocenter
    module procedure init_Hypocenter
  end interface Hypocenter
  public :: Hypocenter
  private :: init_Hypocenter

!=======================================================================
! Object Source
!=======================================================================

  type, extends(Hypocenter) :: Source
    ! Moment tensor
    real(kind = RPRE), dimension(:), allocatable, public :: moment_tensor
  end type Source

  interface Source
    module procedure init_Source
  end interface Source
  public :: Source
  private :: init_Source

!=======================================================================
! Object Cluster
!=======================================================================

  type Cluster
    ! Sources
    type(Source), dimension(:), allocatable, public :: sources
    ! Number of sources
    integer(kind = IPRE), public :: nsrc
  contains
    procedure, private, pass :: locate_cluster
    generic, public :: locate => locate_cluster
    procedure, public, pass :: convert_to_latlon, savehypo
  end type Cluster

  interface Cluster
    module procedure init_Cluster
  end interface Cluster
  public :: Cluster
  private :: init_Cluster

!=======================================================================
! Object Station
!=======================================================================

  type, extends(Location) :: Station
    ! Observed arrival time (in seconds)
    real(kind = RPRE), dimension(:,:), allocatable, public :: arrivals
    ! Observed uncertainties (in seconds)
    real(kind = RPRE), dimension(:,:), allocatable, public :: uncertainties
    ! Computed traveltime grid (in seconds)
    real(kind = RPRE), dimension(:,:,:,:), allocatable, public :: ttgrid
  contains
    procedure, public, pass :: compute_ttgrid, get_ttime
  end type Station

  interface Station
    module procedure init_Station
  end interface Station
  public :: Station
  private :: init_Station

!=======================================================================
! Object Network
!=======================================================================

  type Network
    ! Stations
    type(Station), dimension(:), allocatable, public :: stations
    ! Number of stations
    integer(kind = IPRE), public :: nstat
    ! Computed traveltime grids (in seconds)
    real(kind = RPRE), dimension(:,:,:,:,:), allocatable, public :: ttgrids
  contains
    procedure, public, pass :: compute_ttgrids, zshift, convert_to_utm
  end type Network

  interface Network
    module procedure init_Network
  end interface Network
  public :: Network
  private :: init_Network

!=======================================================================
! Object Vel3D
!=======================================================================

  type Vel3D
    ! Earth velocity model (in m/s)
    real(kind = RPRE), dimension(:,:,:), pointer, public :: model => null()
    ! Grid size
    integer(kind = IPRE), public :: nx, ny, nz
    ! Mesh size (in meters)
    real(kind = RPRE), public :: dx, dy, dz
  end type Vel3D

  interface Vel3D
    module procedure init_Vel3D
  end interface Vel3D
  public :: Vel3D
  private :: init_Vel3D

contains

!=======================================================================
! Object Location methods
!=======================================================================

  type(Location) function init_Location(x, y, z)
    real(kind = RPRE), intent(in) :: x, y, z

    init_Location % x = x
    init_Location % y = y
    init_Location % z = z
    return
  end function init_Location

!=======================================================================
! Object Hypocenter methods
!=======================================================================

  type(Hypocenter) function init_Hypocenter(x, y, z, eventid, origin_time)
    integer(kind = IPRE), intent(in), optional :: eventid
    real(kind = RPRE), intent(in), optional :: x, y, z
    real(kind = 8), intent(in), optional :: origin_time

    if (present(x)) init_Hypocenter % x = x
    if (present(y)) init_Hypocenter % y = y
    if (present(z)) init_Hypocenter % z = z
    if (present(eventid)) init_Hypocenter % eventid = eventid
    if (present(origin_time)) init_Hypocenter % origin_time = origin_time
    return
  end function init_Hypocenter

  subroutine locate_hypo(self, velp, vels, ntwrk, method, np, itermax, par, fit)

    ! Input arguments
    !=================
    class(Hypocenter), intent(inout) :: self
    type(Vel3D), intent(in), target :: velp, vels
    type(Network), intent(in), target :: ntwrk
    integer(kind = IPRE), intent(in), optional :: method, np, itermax
    real(kind = RPRE), intent(in), optional :: par(2)
    real(kind = RPRE), intent(inout), optional :: fit

    ! Local variables
    !=================
    integer(kind = IPRE) :: i, opt_method, opt_np, opt_itermax, niter, flag
    real(kind = RPRE) :: opt_par(2), tobs_sum, tobs_mean
    real(kind = RPRE) :: xl(3), xu(3), opt_fit, w, phi1, phi2, eps1, eps2, eps3
    real(kind = 8) :: t0
    real(kind = RPRE), dimension(:), allocatable :: tcalcp, tcalcs
    real(kind = RPRE), dimension(:), allocatable :: loc
    procedure(optifunc), pointer :: costfunc => null()

    ! Global variables to transfer to the cost function
    !===================================================
    integer(kind = IPRE), save :: nstat
    real(kind = RPRE), save :: sigma_sum
    real(kind = RPRE), dimension(:), allocatable, save :: tobsp, tobss, &
      tobsp_tild, tobss_tild, sigmap, sigmas
    logical, dimension(:), allocatable, save :: maskp, masks
    type(Vel3D), pointer, save :: velp_ptr, vels_ptr
    type(Network), pointer, save :: ntwrk_ptr

    ! Optional arguments default values
    !===================================
    opt_method = 2              ! Method (1: CPSO, 2: DE)
    opt_np = 20                 ! Number of individuals
    opt_itermax = 100           ! Maximum number of iterations
    if (present(method)) opt_method = method
    if (present(np)) opt_np = np
    if (present(itermax)) opt_itermax = itermax

    ! Retrieve arrival times associated to the current hypocenter
    !=============================================================
    nstat = size(ntwrk % stations)
    tobsp = [ ( ntwrk % stations(i) % arrivals(self % eventid,1), i = 1, nstat ) ]
    tobss = [ ( ntwrk % stations(i) % arrivals(self % eventid,2), i = 1, nstat ) ]
    sigmap = [ ( ntwrk % stations(i) % uncertainties(self % eventid,1), i = 1, nstat ) ]
    sigmas = [ ( ntwrk % stations(i) % uncertainties(self % eventid,2), i = 1, nstat ) ]
    maskp = tobsp .gt. 0.0d0
    masks = tobss .gt. 0.0d0

    ! Remove the weighted mean from the observed arrival times
    !==========================================================
    allocate(tobsp_tild(nstat), tobss_tild(nstat))
    tobs_sum = sum(tobsp/sigmap**2, mask = maskp) &
               + sum(tobss/sigmas**2, mask = masks)
    sigma_sum = sum(1/sigmap**2, mask = maskp) &
                + sum(1/sigmas**2, mask = masks)
    tobs_mean = tobs_sum/sigma_sum
    tobsp_tild = merge(tobsp - tobs_mean, real(-5d-3, RPRE), maskp)
    tobss_tild = merge(tobss - tobs_mean, real(-5d-3, RPRE), masks)

    ! Optimization method parameters
    !================================
    costfunc => bayfunc
    xl = 0.0d0
    xu = [ velp % nx * velp % dx, &
           velp % ny * velp % dy, &
           velp % nz * velp % dz ]
    eps1 = 1.0d-4
    eps2 = -1.0d+8
    eps3 = 1.0d-4

    ! Set pointers (not to allocate too much)
    !=========================================
    velp_ptr => velp
    vels_ptr => vels
    ntwrk_ptr => ntwrk

    ! Locate
    !========
    select case(method)

      ! Competitive PSO
      !=================
      case(1)
        opt_par = [ 0.7, 1.5 ]
        if (present(par)) opt_par = par

        loc = cpso(costfunc, xl, xu, np = opt_np, itermax = opt_itermax, &
                   w = opt_par(1), phi1 = opt_par(2), phi2 = opt_par(2), &
                   alpha = real(1.25, RPRE), eps1 = eps1, eps2 = eps2, eps3 = eps3, &
                   fit = opt_fit, niter = niter, flag = flag, snap = .false.)

      ! Differential Evolution
      !========================
      case(2)
        opt_par = [ 0.1, 0.5 ]
        if (present(par)) opt_par = par

        loc = de(costfunc, xl, xu, np = opt_np, itermax = opt_itermax, &
                 CR = opt_par(1), F = opt_par(2), eps1 = eps1, eps2 = eps2, &
                 fit = opt_fit, niter = niter, flag = flag, snap = .false.)

    end select

    ! Compute origin time
    !=====================
    allocate(tcalcp(nstat), tcalcs(nstat))
    !$omp parallel default(shared)
    !$omp do schedule(runtime)
    do k = 1, nstat
      tcalcp(k) = ntwrk % stations(k) % get_ttime(velp, loc(1), loc(2), loc(3), 1)
      tcalcs(k) = ntwrk % stations(k) % get_ttime(vels, loc(1), loc(2), loc(3), 2)
    end do
    !$omp end parallel

    t0 = sum((dble(tobsp) - dble(tcalcp)) / dble(sigmap)**2, maskp) &
         + sum((dble(tobss) - dble(tcalcs)) / dble(sigmas)**2, masks)
    t0 = t0 / ( sum(1.0d0 / dble(sigmap)**2, maskp) &
                + sum(1.0d0 / dble(sigmas)**2, masks) )

    ! Update hypocenter
    !===================
    self % x = loc(1)
    self % y = loc(2)
    self % z = loc(3)
    self % origin_time = t0

    if (present(fit)) fit = opt_fit
    deallocate(tobsp_tild, tobss_tild)
    return

  contains

  !---------------------------------------------------------------------
  ! Function bayfunc
  !---------------------------------------------------------------------
    type(optifit) function bayfunc(x)
      real(kind = RPRE), dimension(:), intent(in) :: x
      integer(kind = IPRE) :: k
      real(kind = RPRE) :: tcalc_sum, tcalc_mean
      real(kind = RPRE), dimension(:), allocatable :: tcalcp, tcalcs, &
        tcalcp_tild, tcalcs_tild

      ! Compute traveltimes for each station
      !======================================
      allocate(tcalcp(nstat), tcalcs(nstat))
      !$omp parallel default(shared)
      !$omp do schedule(runtime)
      do k = 1, nstat
        tcalcp(k) = ntwrk_ptr % stations(k) % get_ttime(velp_ptr, x(1), x(2), x(3), 1)
        tcalcs(k) = ntwrk_ptr % stations(k) % get_ttime(vels_ptr, x(1), x(2), x(3), 2)
      end do
      !$omp end parallel

      ! Remove the weighted mean
      !==========================
      allocate(tcalcp_tild(nstat), tcalcs_tild(nstat))
      tcalc_sum = sum(tcalcp/sigmap**2, mask = maskp) &
                  + sum(tcalcs/sigmas**2, mask = masks)
      tcalc_mean = tcalc_sum/sigma_sum
      tcalcp_tild = merge(tcalcp - tcalc_mean, real(-5d-3, RPRE), maskp)
      tcalcs_tild = merge(tcalcs - tcalc_mean, real(-5d-3, RPRE), masks)

      ! Compute cost function
      !=======================
      bayfunc % val = sum((tobsp_tild - tcalcp_tild)**2/sigmap**2, mask = maskp) &
                      + sum((tobss_tild - tcalcs_tild)**2/sigmas**2, mask = masks)
      bayfunc % val = 0.5d0 * bayfunc % val
      return
    end function bayfunc

  !---------------------------------------------------------------------
  ! Function edtfunc
  !---------------------------------------------------------------------
  ! http://alomax.free.fr/nlloc/edt/edt.html
  !---------------------------------------------------------------------
    type(optifit) function edtfunc(x)
      real(kind = RPRE), dimension(:), intent(in) :: x
      integer(kind = IPRE) :: k, l
      real(kind = RPRE) :: tcalc_sum, tcalc_mean
      real(kind = RPRE), dimension(:), allocatable :: tcalcp, tcalcs, &
        tcalcp_tild, tcalcs_tild, sigmap_edt, sigmas_edt

      ! Compute traveltimes for each station
      !======================================
      allocate(tcalcp(nstat), tcalcs(nstat))
      !$omp parallel default(shared)
      !$omp do schedule(runtime)
      do k = 1, nstat
        tcalcp(k) = ntwrk_ptr % stations(k) % get_ttime(velp_ptr, x(1), x(2), x(3), 1)
        tcalcs(k) = ntwrk_ptr % stations(k) % get_ttime(vels_ptr, x(1), x(2), x(3), 2)
      end do
      !$omp end parallel

      ! Compute cost function
      !=======================
      edtfunc % val = 0.0d0
      do k = 1, nstat-1
        if (maskp(k)) then
          do l = k+1, nstat
            if (maskp(l)) then
              edtfunc % val = edtfunc % val &
                              + exp( - ( ( tobsp(k) - tobsp(l) ) - ( tcalcp(k) - tcalcp(l) ) )**2 &
                                     / ( sigmap(k)**2 + sigmap(l)**2 ) ) &
                              / sqrt( sigmap(k)**2 + sigmap(l)**2 )
            end if
          end do
        end if
        if (masks(k)) then
          do l = k+1, nstat
            if (masks(l)) then
              edtfunc % val = edtfunc % val &
                              + exp( - ( ( tobss(k) - tobss(l) ) - ( tcalcs(k) - tcalcs(l) ) )**2 &
                                     / ( sigmas(k)**2 + sigmas(l)**2 ) ) &
                              / sqrt( sigmas(k)**2 + sigmas(l)**2 )
            end if
          end do
        end if
      end do
      edtfunc % val = - edtfunc % val
      return
    end function edtfunc

  end subroutine locate_hypo

!=======================================================================
! Object Source methods
!=======================================================================

  type(Source) function init_Source(x, y, z, eventid, origin_time, moment_tensor)
    integer(kind = IPRE), intent(in), optional :: eventid
    real(kind = RPRE), intent(in), optional :: x, y, z
    real(kind = 8), intent(in), optional :: origin_time
    real(kind = RPRE), dimension(:), intent(in), optional :: moment_tensor

    if (present(x)) init_Source % x = x
    if (present(y)) init_Source % y = y
    if (present(z)) init_Source % z = z
    if (present(eventid)) init_Source % eventid = eventid
    if (present(origin_time)) init_Source % origin_time = origin_time
    if (present(moment_tensor)) init_Source % moment_tensor = moment_tensor
    return
  end function init_Source

!=======================================================================
! Object Cluster methods
!=======================================================================

  type(Cluster) function init_Cluster(nsrc, X, Y, Z, EID, T0, MT)
    integer(kind = IPRE), intent(in) :: nsrc
    real(kind = RPRE), dimension(:), intent(in), optional :: X, Y, Z, T0
    integer(kind = IPRE), dimension(:), intent(in), optional :: EID
    real(kind = RPRE), dimension(:,:), intent(in), optional :: MT
    integer(kind = IPRE) :: j

    init_Cluster % nsrc = nsrc
    allocate(init_Cluster % sources(nsrc))
    if (present(X)) init_Cluster % sources(:) % x = X
    if (present(Y)) init_Cluster % sources(:) % y = Y
    if (present(Z)) init_Cluster % sources(:) % z = Z
    if (present(EID)) init_Cluster % sources(:) % eventid = EID
    if (present(T0)) init_Cluster % sources(:) % origin_time = T0
    if (present(MT)) then
      do j = 1, nsrc
        init_Cluster % sources(j) % moment_tensor(:) = MT(j,:)
      end do
    end if
    return
  end function init_Cluster

  subroutine locate_cluster(self, velp, vels, ntwrk, method, np, itermax, par, fit)

    ! Input arguments
    !=================
    class(Cluster), intent(inout) :: self
    type(Vel3D), intent(in) :: velp, vels
    type(Network), intent(in) :: ntwrk
    integer(kind = IPRE), intent(in), optional :: method, np, itermax
    real(kind = RPRE), intent(in), optional :: par(2)
    real(kind = RPRE), intent(inout), optional :: fit

    ! Local variables
    !=================
    integer(kind = IPRE) :: opt_method, opt_np, opt_itermax
    real(kind = RPRE) :: opt_fit, opt_par(2)
    integer(kind = IPRE) :: j

    ! Optional arguments default values
    !===================================
    opt_method = 2              ! Method (1: CPSO, 2: DE)
    opt_np = 20                 ! Number of individuals
    opt_itermax = 100           ! Maximum number of iterations
    if (present(method)) opt_method = method
    if (present(np)) opt_np = np
    if (present(itermax)) opt_itermax = itermax

    select case(opt_method)
      case(1)
        opt_par = [ 0.7, 1.5 ]
        if (present(par)) opt_par = par
      case(2)
        opt_par = [ 0.1, 0.5 ]
        if (present(par)) opt_par = par
    end select

    do j = 1, self % nsrc
      call self % sources(j) % locate(velp, vels, ntwrk, &
                                      opt_method, opt_np, opt_itermax, &
                                      opt_par, opt_fit)
    end do

    if (present(fit)) fit = opt_fit
    return
  end subroutine locate_cluster

  subroutine convert_to_latlon(self, south, west)
    class(Cluster), intent(inout) :: self
    real(kind = RPRE), intent(in) :: south, west
    integer(kind = IPRE) :: zn
    real(kind = RPRE) :: east0, north0
    character(len = 1) :: zl
    integer(kind = IPRE), dimension(:), allocatable :: loc_zn
    real(kind = RPRE), dimension(:), allocatable :: loc_lat, loc_lon
    character(len = 1), dimension(:), allocatable :: loc_zl

    call deg2utm(south, west, east0, north0, zn, zl)
    allocate(loc_zn(self % nsrc), loc_zl(self % nsrc))
    loc_zn = zn
    loc_zl = zl
    self % sources(:) % x = self % sources(:) % x + east0
    self % sources(:) % y = self % sources(:) % y + north0
    call utm2deg( self % sources(:) % x, self % sources(:) % y, &
                  loc_zn, loc_zl, loc_lat, loc_lon )
    self % sources(:) % x = loc_lat
    self % sources(:) % y = loc_lon
    return
  end subroutine convert_to_latlon

  subroutine savehypo(self, outfile, coord_unit, year, month, day)
    class(Cluster), intent(in) :: self
    character(len = *), intent(in) :: outfile
    integer(kind = IPRE), intent(in) :: coord_unit, year, month, day
    integer(kind = IPRE) :: j
    type(File) :: fout

    fout = File(999, trim(outfile))
    call fout % open()
    select case(coord_unit)
      case(0)
        if ((year .eq. 0) .and. (month .eq. 0) .and. (day .eq. 0)) then
          write(fout % unit, *), "Time (s)                    ", &
                                 "X (m)            ", &
                                 "Y (m)            ", &
                                 "Z (m)"
          do j = 1, self % nsrc
            write(fout % unit, *), self % sources(j) % origin_time, &
                                   self % sources(j) % x, &
                                   self % sources(j) % y, &
                                   self % sources(j) % z
          end do
        else
          datemin = datenum(year, month, day)
          write(fout % unit, *), "Date                          ", &
                                  "X (m)           ", &
                                  "Y (m)           ", &
                                  "Z (m)"
          do j = 1, self % nsrc
            write(fout % unit, *), datestr(datemin + self % sources(j) % origin_time &
                                                   / ( 24.0d0 * 60.0d0 * 60.0d0 )), &
                                   self % sources(j) % x, &
                                   self % sources(j) % y, &
                                   self % sources(j) % z
          end do
        end if
      case(1)
        if ((year .eq. 0) .and. (month .eq. 0) .and. (day .eq. 0)) then
          write(fout % unit, *), "Time (s)                    ", &
                                 "Latitude (°)     ", &
                                 "Longitude (°)    ", &
                                 "Depth (m)"
          do j = 1, self % nsrc
            write(fout % unit, *), self % sources(j) % origin_time, &
                                   self % sources(j) % x, &
                                   self % sources(j) % y, &
                                   self % sources(j) % z
          end do
        else
          datemin = datenum(year, month, day)
          write(fout % unit, *), "Date                          ", &
                                 "Latitude (°)     ", &
                                 "Longitude (°)    ", &
                                 "Depth (m)"
          do j = 1, self % nsrc
            write(fout % unit, *), datestr(datemin + self % sources(j) % origin_time &
                                                   / ( 24.0d0 * 60.0d0 * 60.0d0 )), &
                                   self % sources(j) % x, &
                                   self % sources(j) % y, &
                                   self % sources(j) % z
          end do
        end if
    end select
    call fout % close()
    return
  end subroutine savehypo

!=======================================================================
! Object Station methods
!=======================================================================

  type(Station) function init_Station(x, y, z, tobs, sigma, ttgrid)
    real(kind = RPRE), intent(in) :: x, y, z
    real(kind = RPRE), dimension(:,:), intent(in) :: tobs
    real(kind = RPRE), dimension(:,:), intent(in), optional :: sigma
    real(kind = RPRE), dimension(:,:,:,:), intent(in), optional :: ttgrid

    init_Station % x = x
    init_Station % y = y
    init_Station % z = z
    init_Station % arrivals = tobs
    if (present(sigma)) init_Station % uncertainties = sigma
    if (present(ttgrid)) init_Station % ttgrid = ttgrid(:,:,:,:)
    return
  end function init_Station

  subroutine compute_ttgrid(self, velp, vels, nsweep)
    class(Station), intent(inout) :: self
    type(Vel3D), intent(in) :: velp, vels
    integer(kind = IPRE), intent(in), optional :: nsweep
    integer(kind = IPRE) :: opt_nsweep, nz, nx, ny
    real(kind = 4), dimension(:,:,:), allocatable :: tmp
    real(kind = RPRE), dimension(:,:,:,:), allocatable :: tt

    opt_nsweep = 1
    if (present(nsweep)) opt_nsweep = nsweep

    nz = velp % nz + 1
    nx = velp % nx + 1
    ny = velp % ny + 1

    allocate(tmp(nz, nx, ny), tt(nz, nx, ny, 2))
    call FTeik3D_2(sngl(velp % model), tmp, nz, nx, ny, &
                   sngl(self % z), sngl(self % x), sngl(self % y), &
                   sngl(velp % dz), sngl(velp % dx), sngl(velp % dy), &
                   opt_nsweep, 5.)
    tt(:,:,:,1) = tmp
    call FTeik3D_2(sngl(vels % model), tmp, nz, nx, ny, &
                   sngl(self % z), sngl(self % x), sngl(self % y), &
                   sngl(vels % dz), sngl(vels % dx), sngl(vels % dy), &
                   opt_nsweep, 5.)
    tt(:,:,:,2) = tmp
    self % ttgrid = tt(:,:,:,:)
    return
  end subroutine compute_ttgrid

  function get_ttime(self, vel, xq, yq, zq, phase) result(ttime)

    ! Input arguments
    !=================
    real(kind = RPRE) :: ttime
    class(Station), intent(in) :: self
    type(Vel3D), intent(in) :: vel
    real(kind = RPRE), intent(in) :: xq, yq, zq
    integer(kind = IPRE), intent(in) :: phase

    ! Local variables
    !=================
    integer(kind = IPRE) :: i, nz, nx, ny
    real(kind = RPRE), dimension(:), allocatable :: az, ax, ay

    integer(kind = IPRE) :: zi(1), xi(1), yi(1), iz(8), ix(8), iy(8)
    real(kind = RPRE) :: z(2), x(2), y(2), ttr(2,2,2), v(2,2,2), d, dq

    nz = vel % nz + 1
    nx = vel % nx + 1
    ny = vel % ny + 1
    az = linspace(0, (nz-1)*vel % dz, nz)
    ax = linspace(0, (nx-1)*vel % dx, nx)
    ay = linspace(0, (ny-1)*vel % dy, ny)

    zi = min( minloc(zq - az, mask = zq .ge. az), nz-2 )
    xi = min( minloc(xq - ax, mask = xq .ge. ax), nx-2 )
    yi = min( minloc(yq - ay, mask = yq .ge. ay), ny-2 )
    ! zi = minloc(zq - az, mask = zq .ge. az)
    ! xi = minloc(xq - ax, mask = xq .ge. ax)
    ! yi = minloc(yq - ay, mask = yq .ge. ay)
    z = [ az(zi(1)), az(zi(1) + 1) ]
    x = [ ax(xi(1)), ax(xi(1) + 1) ]
    y = [ ay(yi(1)), ay(yi(1) + 1) ]
    iz = [ 1, 2, 1, 2, 1, 2, 1, 2 ]
    ix = [ 1, 1, 2, 2, 1, 1, 2, 2 ]
    iy = [ 1, 1, 1, 1, 2, 2, 2, 2 ]

    ttr = self % ttgrid(zi(1):zi(1) + 1, xi(1):xi(1) + 1, yi(1):yi(1) + 1, phase)
    dq = sqrt((self % x - xq)**2 &
              + (self % y - yq)**2 &
              + (self % z - zq)**2)
    do i = 1, 8
      d = sqrt((self % x - x(ix(i)))**2 &
                + (self % y - y(iy(i)))**2 &
                + (self % z - z(iz(i)))**2)
      v(iz(i), ix(i), iy(i)) = d / ttr(iz(i), ix(i), iy(i))
    end do
    ttime = dq / interp3(z, x, y, v, zq, xq, yq)
    return
  end function get_ttime

!=======================================================================
! Object Network methods
!=======================================================================

  type(Network) function init_Network(X, Y, Z, tobs, sigma, ttgrids)
    real(kind = RPRE), dimension(:), intent(in) :: X, Y, Z
    real(kind = RPRE), dimension(:,:), intent(in) :: tobs
    real(kind = RPRE), dimension(:,:), intent(in), optional :: sigma
    real(kind = RPRE), dimension(:,:,:,:,:), intent(in), optional :: ttgrids
    integer(kind = IPRE) :: k, nstat, nobs, nev
    real(kind = RPRE), dimension(:,:), allocatable :: tobsp, tobss, &
      sigmap, sigmas

    nstat = size(X)
    nobs = size(tobs, 1)
    nev = nobs / nstat
    init_Network % nstat = nstat
    allocate(init_Network % stations(nstat))

    tobsp = reshape(tobs(:,1), [ nstat, nev ])
    tobss = reshape(tobs(:,2), [ nstat, nev ])

    do k = 1, nstat
      init_Network % stations(k) = Station(X(k), Y(k), Z(k), &
                                           horzcat(tobsp(k,:), tobss(k,:)))
      if (present(sigma)) then
        sigmap = reshape(sigma(:,1), [ nstat, nev ])
        sigmas = reshape(sigma(:,2), [ nstat, nev ])
        init_Network % stations(k) % uncertainties = horzcat(sigmap(k,:), sigmas(k,:))
      end if
      if (present(ttgrids)) then
        init_Network % stations(k) % ttgrid = ttgrids(:,:,:,:,k)
      end if
    end do
    return
  end function init_Network

  subroutine compute_ttgrids(self, velp, vels)
    class(Network), intent(inout) :: self
    type(Vel3D), intent(in) :: velp, vels
    integer(kind = IPRE) :: k

    !$omp parallel default(shared)
    !$omp do schedule(runtime)
    do k = 1, self % nstat
      call self % stations(k) % compute_ttgrid(velp, vels)
    end do
    !$omp end parallel
    return
  end subroutine compute_ttgrids

  subroutine zshift(self, zs, nsh, nz, dz)
    class(Network), intent(inout) :: self
    real(kind = RPRE), intent(out) :: zs
    integer(kind = IPRE), intent(inout) :: nz, nsh
    real(kind = RPRE), intent(in) :: dz

    nsh = ceiling( maxval(self % stations(:) % z) / dz )
    nz = nz + nsh
    zs = nsh * dz
    if ( zs .gt. 0.0d0 ) self % stations(:) % z = -self % stations(:) % z + zs
    return
  end subroutine zshift

  subroutine convert_to_utm(self, south, west)
    class(Network), intent(inout) :: self
    real(kind = RPRE), intent(in) :: south, west
    integer(kind = IPRE) :: zn
    real(kind = RPRE) :: east0, north0
    character(len = 1) :: zl
    integer(kind = IPRE), dimension(:), allocatable :: stat_zn
    real(kind = RPRE), dimension(:), allocatable :: stat_east, stat_north
    character(len = 1), dimension(:), allocatable :: stat_zl

    call deg2utm(south, west, east0, north0, zn, zl)
    call deg2utm(self % stations(:) % x, self % stations(:) % y, &
                 stat_east, stat_north, stat_zn, stat_zl)
    self % stations(:) % x = stat_east - east0
    self % stations(:) % y = stat_north - north0
    return
  end subroutine convert_to_utm

!=======================================================================
! Object Vel3D methods
!=======================================================================

  type(Vel3D) function init_Vel3D(model, dx, dy, dz)
    real(kind = RPRE), dimension(:,:,:), target, intent(in) :: model
    real(kind = RPRE), intent(in) :: dx, dy, dz

    init_Vel3D % model => model
    init_Vel3D % dx = dx
    init_Vel3D % dy = dy
    init_Vel3D % dz = dz
    init_Vel3D % nx = size(model,2)
    init_Vel3D % ny = size(model,3)
    init_Vel3D % nz = size(model,1)
    return
  end function init_Vel3D

end module evoseis
