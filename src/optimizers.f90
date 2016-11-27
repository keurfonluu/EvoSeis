!=======================================================================
! Optimizers
!-----------------------------------------------------------------------
! Optimizers provides functions to solve non-linear optimization
! problems.
!
! Created by
!     Keurfon Luu <keurfon.luu@mines-paristech.fr>
!     MINES ParisTech - Centre de Géosciences
!     PSL - Research University
!
! Last updated
!     2016-11-04 14:35
!
! Functions
!-----------------------------------------------------------------------
!   -   cmaes
!   -   cpso
!   -   de
!
! Subroutines
!-----------------------------------------------------------------------
!   -   export_parameters
!   -   snap_positions
!=======================================================================

module optimizers

  use forlab

  implicit none

!=======================================================================
! Type optifit
!=======================================================================
  type optifit
    real(kind = RPRE) :: val
    real(kind = RPRE), dimension(:), allocatable :: extra
  end type optifit

!=======================================================================
! Function type optifunc
!=======================================================================
  abstract interface
    type(optifit) function optifunc(x)
      use forlab, only : RPRE
      import
      real(kind = RPRE), dimension(:), intent(in) :: x
    end function optifunc
  end interface

contains

!=======================================================================
! cmaes
!-----------------------------------------------------------------------
! cmaes computes the nd parameters that minimize the input function
! using a Covariance Matrix Adaptation - Evolution Strategy (CMA-ES).
!
! Syntax
!-----------------------------------------------------------------------
! x = cmaes(fitness, xl, xu)
! x = cmaes(fitness, xl, xu, [options = ])
!
! Description
!-----------------------------------------------------------------------
! x = cmaes(fitness, xl, xu) returns a nd vector with the parameters
! that minimizes the function fitness, given the lower and upper
! boundary xl and xu.
!
! Inputs
!-----------------------------------------------------------------------
! fitness             Function to be minimized (real)
! xl                  Lower boundary (real nd array)
! xu                  Upper boundary (real nd array)
!
! Options
!-----------------------------------------------------------------------
! Xi                  Initial mean vector
! lambda = 4+3log(nd) Number of offsprings
! mu = lambda/2       Number of parents
! sigma = 1.0d0       Step size
! itermax = 100       Maximum number of iterations
! eps3 = 1.0d-8       Minimum fitness precision
! fit                 Global best position fitness
! niter               Output number of iterations
! neval               Output number of fitness evaluations
! flag                Output reason of termination
!                       -   0   Fitness is lower than precision eps,
!                       -   1   NoEffectAxis,
!                       -   2   NoEffectCoord,
!                       -   3   ConditionCov,
!                       -   4   EqualFunValues,
!                       -   5   TolXUp,
!                       -   6   TolFun,
!                       -   7   TolX,
!                       -   8   Maximum number of iterations is reached.
! extra               Additional outputs from fitness function
! snap = .false.      Save offsprings positions after each iteration
! snapdir = "snap/"   Snapshot directory
!
! Note
!-----------------------------------------------------------------------
! The code has been translated from the Matlab code described in:
! "The CMA Evolution Strategy: A Tutorial" by Nikolaus Hansen.
! Additional features have been added following CMA-ES codes webpage:
! https://www.lri.fr/~hansen/cmaes_inmatlab.html
!
! In an inversion problem with parameters of different scales, a
! practical hint is to scale the parameters so that they all lie between
! 0 and 10 for instance (therefore, each parameter has the same
! sensitivity to sigma). Rescaling has to be done in the fitness
! function:
! scaling = 10.0d0 / xu
! x = cmaes(fitness, xl*scaling, xu*scaling, sigma = 2.)
! x = x / scaling
!=======================================================================

  function cmaes(fitness, xl, xu, Xi, lambda, mu, sigma, itermax, &
               eps, fit, niter, neval, flag, extra, &
               snap, snapdir)

    ! Input arguments
    !=================
    real(kind = RPRE), dimension(:), allocatable :: cmaes
    procedure(optifunc) :: fitness
    real(kind = RPRE), dimension(:), intent(in) :: xl, xu
    real(kind = RPRE), dimension(:), intent(in), optional :: Xi
    integer(kind = IPRE), intent(in), optional :: lambda, mu, itermax
    real(kind = RPRE), intent(in), optional :: sigma, eps
    real(kind = RPRE), intent(inout), optional :: fit
    integer(kind = IPRE), intent(inout), optional :: niter, flag, neval
    real(kind = RPRE), dimension(:,:,:), allocatable, intent(inout), optional :: extra
    logical, intent(in), optional :: snap
    character(len = *), intent(in), optional :: snapdir

    ! Local variables
    !=================
    integer(kind = IPRE) :: opt_lambda, opt_mu, opt_itermax, nd, i, k, &
      iter, counteval, eigeneval, ilim, ne
    real(kind = RPRE) :: opt_sigma, opt_eps, mueff, cc, cs, c1, cmu, &
      damps, chind, hsig, perc(2), delta, insigma
    integer(kind = IPRE), dimension(:), allocatable :: arindex
    real(kind = RPRE), dimension(:), allocatable :: xmean, &
      xold, weights, pc, ps, D, arfitness, arbestfitness, &
      bnd_weights, bnd_scale, dfithist, tx
    real(kind = RPRE), dimension(:,:), allocatable :: B, C, invsqrtC, &
      arx, arxvalid, artmp, pctmp
    real(kind = RPRE), dimension(:,:,:), allocatable :: opt_extra
    logical :: opt_snap, validfitval = .false., iniphase = .true., &
      converge = .false.
    logical, dimension(:), allocatable :: ti, idx
    character(len = :), allocatable :: opt_snapdir
    type(optifit) :: tmpfit

    ! Population initial positions
    !==============================
    nd = size(xl)               ! Number of dimensions
    if (present(Xi)) then
      xmean = Xi
    else
      xmean = xl + randu(nd)*(xu - xl)
    end if

    ! Optional arguments default values
    !===================================
    ! Number of offsprings
    if (present(lambda)) then
      opt_lambda = lambda
    else
      opt_lambda = 4 + 3 * log(real(nd))
    end if

    ! Number of parents
    if (present(mu)) then
      opt_mu = mu
    else
      opt_mu = opt_lambda / 2
    end if

    opt_itermax = 100           ! Maximum number of iterations
    opt_sigma = 1.0d0           ! Step size
    opt_eps = 1.0d-8            ! Minimum fitness precision
    opt_snap = .false.          ! Save offsprings positions after each iteration
    opt_snapdir = "snap/"       ! Snapshot directory
    if (present(itermax)) opt_itermax = itermax
    if (present(sigma)) opt_sigma = sigma
    if (present(eps)) opt_eps = eps
    if (present(snap)) opt_snap = snap
    if (present(snapdir)) opt_snapdir = snapdir

    ! Strategy parameter setting: Selection
    !=======================================
    weights = log(opt_mu+0.5d0) - log(linspace(1, opt_mu, opt_mu))
    weights = weights / sum(weights)
    mueff = sum(weights)**2 / sum(weights**2)

    ! Strategy parameter setting: Adaptation
    !========================================
    cc = ( 4.0d0 + mueff / nd ) / ( nd + 4.0d0 + 2.0d0 * mueff / nd )
    cs = ( mueff + 2.0d0 ) / ( nd + mueff + 5.0d0 )
    c1 = 2.0d0 / ( (nd + 1.3d0)**2 + mueff )
    cmu = min( 1.0d0-c1, 2.0d0 * ( mueff - 2.0d0 + 1.0d0 / mueff ) / ( (nd+2.0d0)**2 + mueff ) )
    damps = 1.0d0 + 2.0d0 * max( 0.0d0, sqrt( ( mueff - 1.0d0 ) / ( nd + 1.0d0 ) ) - 1.0d0 ) + cs

    ! Initialize dynamic (internal) strategy parameters and constants
    !=================================================================
    pc = zeros(nd)
    ps = zeros(nd)
    B = eye(nd, nd)
    D = ones(nd)
    C = eye(nd, nd)
    invsqrtC = eye(nd, nd)
    chind = sqrt(real(nd, RPRE)) * ( 1.0d0 - 1.0d0 / ( 4.0d0*nd ) + 1.0d0 / ( 21.0d0*nd**2 ) )

    ! Initialize boundaries weights
    !===============================
    bnd_weights = zeros(nd)
    bnd_scale = zeros(nd)
    dfithist = [ real(1., RPRE) ]

    ! Save offsprings positions
    !===========================
    if ( opt_snap ) call system("rm -rf " // opt_snapdir)

    ! (opt_mu, opt_lambda)-CMA-ES
    !=============================
    iter = 0
    counteval = 0
    eigeneval = 0
    arx = zeros(nd, opt_lambda)
    arxvalid = zeros(nd, opt_lambda)
    arfitness = zeros(opt_lambda)
    arbestfitness = zeros(opt_itermax)
    ilim = 10 + 30*nd/opt_lambda
    insigma = opt_sigma
    validfitval = .false.
    iniphase = .true.
    converge = .false.

    do while ( .not. converge )
      iter = iter + 1

      ! Generate lambda offspring
      !===========================
      do k = 1, opt_lambda
        arx(:,k) = xmean + opt_sigma * matmul(B, D*randn(nd))
        arxvalid(:,k) = arx(:,k)
        arxvalid(:,k) = merge(arxvalid(:,k), xl, arxvalid(:,k) .ge. xl)
        arxvalid(:,k) = merge(arxvalid(:,k), xu, arxvalid(:,k) .le. xu)
      end do

      ! Evaluate fitness
      !==================
      arfitness = zeros(opt_lambda)
      do k = 1, opt_lambda
        tmpfit = fitness(arxvalid(:,k))
        arfitness(k) = tmpfit % val
        if (present(extra)) then
          if (.not. allocated(opt_extra)) then
            if (.not. allocated(tmpfit % extra)) then
              print *, "Error: in cmaes, extra requested but not defined in fitness."
              stop
            else
              ne = size(tmpfit % extra)
              opt_extra = zeros(ne, opt_lambda, opt_itermax)
            end if
          end if
          opt_extra(:,k,iter) = tmpfit % extra
        end if
        counteval = counteval + 1
      end do

      ! Handle boundaries by penalizing fitness
      !=========================================
      ! Get delta fitness values
      perc = prctile(arfitness, [ 25, 75 ])
      delta = ( perc(2) - perc(1) ) / real(nd, RPRE) / mean(diag(C)) / opt_sigma**2

      ! Catch non-sensible values
      if ( delta .eq. 0.0d0 ) then
        delta = minval( dfithist( find( dfithist .gt. 0.0d0 ) ) )
      elseif ( .not. validfitval ) then
        dfithist = zeros(0)
        validfitval = .true.
      end if

      ! Store delta fitness values
      if ( size(dfithist) .lt. 20 + (3*nd)/opt_lambda ) then
        dfithist = [ dfithist, delta ]
      else
        dfithist = [ dfithist(2:size(dfithist)), delta ]
      end if

      ! Corrected mean
      ti = xmean .lt. xl .or. xmean .gt. xu
      tx = xmean
      tx = merge(tx, xl, tx .ge. xl)
      tx = merge(tx, xu, tx .le. xu)

      ! Set initial weights
      if ( iniphase ) then
        if ( any(ti) ) then
          bnd_weights = 2.0002d0 * median(dfithist)
          if ( validfitval .and. iter .gt. 2 ) iniphase = .false.
        end if
      end if

      if ( any(ti) ) then
        tx = xmean - tx
        idx = ti &
              .and. abs(tx) .gt. 3.0d0 * max( 1.0d0, sqrt(real(nd, RPRE))/mueff ) &
                           * opt_sigma * sqrt(diag(C))
        idx = idx .and. signum(tx) .eq. signum(xmean - xold)
        where ( idx )
          bnd_weights = bnd_weights * 1.2d0**( min( 1.0d0, mueff/10.0d0/real(nd, RPRE) ) )
        end where
      end if

      ! Calculate scaling biased to unity, product is one
      bnd_scale = exp( 0.9d0 * ( log(diag(C)) - mean(log(diag(C))) ) )

      ! Assigned penalized fitness
      arfitness = arfitness &
                  + matmul(bnd_weights / bnd_scale, (arxvalid - arx)**2)

      ! Sort by fitness and compute weighted mean into xmean
      !======================================================
      arindex = argsort(arfitness)
      xold = xmean
      xmean = matmul(arx(:,arindex(:opt_mu)), weights)

      ! Save best fitness
      !===================
      arbestfitness(iter) = arfitness(arindex(1))

      ! Cumulation: Update evolution paths
      !====================================
      ps = ( 1.0d0 - cs ) * ps &
           + sqrt( cs * ( 2.0d0 - cs ) * mueff ) * matmul(invsqrtC, (xmean - xold)) / opt_sigma
      if ( norm(ps) / sqrt( 1.0d0 - ( 1.0d0 - cs )**(2.0d0*counteval/opt_lambda) ) / chind &
           .lt. 1.4d0 + 2.0d0 / ( nd + 1.0d0 ) ) then
        hsig = 1.0d0
      else
        hsig = 0.0d0
      end if
      pc = ( 1.0d0 - cc ) * pc &
           + hsig * sqrt( cc * ( 2.0d0 - cc ) * mueff ) * (xmean - xold) / opt_sigma

      ! Adapt covariance matrix C
      !===========================
      artmp = ( arx(:,arindex(:opt_mu)) - repmat(xold, opt_mu) ) / opt_sigma
      pctmp = reshape( pc, shape = [ size(pc), 1 ], order = [ 1, 2 ] )
      C = ( 1.0d0 - c1 - cmu ) * C &
          + c1 * ( matmul(pctmp, transpose(pctmp)) &
                   + ( 1.0d0 - hsig ) * cc * ( 2.0d0 - cc ) * C ) &
          + cmu * matmul(matmul(artmp, diag(weights)), transpose(artmp))

      ! Adapt step size sigma
      !=======================
      opt_sigma = opt_sigma * exp( ( cs / damps ) * ( norm(ps) / chind - 1.0d0 ) )

      ! Diagonalization of C
      !======================
      if ( counteval - eigeneval .gt. opt_lambda / ( c1 + cmu ) / nd / 10.0d0 ) then
        eigeneval = counteval
        C = triu(C) + transpose(triu(C, 1))
        call eig(C, B, D)
        D = sqrt(D)
        invsqrtC = matmul(matmul(B, diag(real(1.0d0/D, RPRE))), transpose(B))
      end if

      ! Stop if fitness is lower than precision eps
      !=============================================
      if ( arfitness(arindex(1)) .le. opt_eps ) then
        converge = .true.
        if (present(flag)) flag = 0
      end if

      ! NoEffectAxis: stop if numerical precision problem
      !===================================================
      i = floor( mod( real(iter), real(nd) ) ) + 1
      if ( all( 0.1d0 * opt_sigma * B(:,i) * D(i) .lt. 1.0d-10 ) ) then
        converge = .true.
        if (present(flag)) flag = 1
      end if

      ! NoEffectCoord: stop if too low coordinate axis deviations
      !===========================================================
      if ( any( 0.2d0 * opt_sigma * sqrt(diag(C)) .lt. 1.0d-10 ) ) then
        converge = .true.
        if (present(flag)) flag = 2
      end if

      ! ConditionCov: stop if the condition number exceeds 10e14
      !==========================================================
      if ( maxval(D) .gt. 1.0d7 * minval(D) ) then
        converge = .true.
        if (present(flag)) flag = 3
      end if

      ! EqualFunValues: stop if the range of fitness values is zero
      !=============================================================
      if ( iter .ge. ilim ) then
        if ( maxval(arbestfitness(iter-ilim+1:iter)) &
             - minval(arbestfitness(iter-ilim+1:iter)) .lt. 1.0d-10 ) then
          converge = .true.
          if (present(flag)) flag = 4
        end if
      end if

      ! TolXUp: stop if x-changes larger than 1.0d3 times initial sigma
      !=================================================================
      if ( any( opt_sigma * sqrt(diag(C)) .gt. 1.0d3 * insigma ) ) then
        converge = .true.
        if (present(flag)) flag = 5
      end if

      ! TolFun: stop if fun-changes smaller than 1.0d-12
      !==================================================
      if ( iter .gt. 2 &
           .and. maxval( [ arfitness, arbestfitness ] ) &
                 - minval( [ arfitness, arbestfitness ] ) .lt. 1.0d-12 ) then
        converge = .true.
        if (present(flag)) flag = 6
      end if

      ! TolX: stop if x-changes smaller than 1.0d-11 times initial sigma
      !==================================================================
      if ( all( opt_sigma * ( max( abs(pc), sqrt(diag(C)) ) ) .lt. 1.0d-11 * insigma ) ) then
        converge = .true.
        if (present(flag)) flag = 7
      end if

      ! Stop if maximum iteration is reached
      !======================================
      if ( iter .ge. opt_itermax ) then
        converge = .true.
        if (present(flag)) flag = 8
      end if

      if ( .not. converge .and. opt_snap ) then
        call snap_positions(arx, arfitness, iter, opt_snapdir)
      end if
    end do

    cmaes = arx(:,arindex(1))
    if (present(fit)) fit = arfitness(arindex(1))
    if (present(niter)) niter = iter
    if (present(neval)) neval = counteval
    if (present(extra)) extra = opt_extra(:,:,:iter)
    if ( opt_snap ) then
      call export_parameters(opt_lambda, nd, iter, counteval, xl, xu, opt_snapdir)
      call snap_positions(arx, arfitness, iter, opt_snapdir)
    end if
    return
  end function cmaes

!=======================================================================
! cpso
!-----------------------------------------------------------------------
! cpso computes the nd parameters that minimize the input function using
! a Competitive Particle Swarm Optimization (Luu et al., 2016).
!
! Syntax
!-----------------------------------------------------------------------
! x = cpso(fitness, xl, xu)
! x = cpso(fitness, xl, xu, [options = ])
!
! Description
!-----------------------------------------------------------------------
! x = cpso(fitness, xl, xu) returns a nd vector with the parameters that
! minimizes the function fitness, given the lower and upper boundary xl
! and xu.
!
! Inputs
!-----------------------------------------------------------------------
! fitness             Function to be minimized (real)
! xl                  Lower boundary (real nd array)
! xu                  Upper boundary (real nd array)
!
! Options
!-----------------------------------------------------------------------
! Xi                  Particles initial positions
! np = 100            Number of particles
! itermax = 100       Maximum number of iterations
! w = 0.7d0           Particles inertia
! l = 0.1d0           Velocity clamping
! phi1 = 1.5d0        Cognition parameter
! phi2 = 1.5d0        Global sociability parameter
! alpha = 1.25d0      Competitivity parameter
! delta = f(nd)       Swarm maximum radius
! eps1 = 1.0d-8       Minimum change in best position
! eps2 = 1.0d-8       Minimum fitness precision
! eps3 = 1.0d-8       Minimum change in global best position fitness
! fit                 Global best position fitness
! niter               Output number of iterations
! neval               Output number of fitness evaluations
! flag                Output reason of termination
!                       -   0   Best position changes less than eps1,
!                       -   1   Maximum number of iterations is reached,
!                       -   2   Best position fitness changes less than
!                               eps3.
! extra               Additional outputs from fitness function
! snap = .false.      Save particles positions after each iteration
! snapdir = "snap/"   Snapshot directory
!
! Notes
!-----------------------------------------------------------------------
! Set alpha = 0 for classical PSO.
!=======================================================================

  function cpso(fitness, xl, xu, Xi, np, itermax, w, l, alpha, delta, &
                phi1, phi2, eps1, eps2, eps3, fit, niter, neval, flag, extra, &
                snap, snapdir)

    ! Input arguments
    !=================
    real(kind = RPRE), dimension(:), allocatable :: cpso
    procedure(optifunc) :: fitness
    real(kind = RPRE), dimension(:), intent(in) :: xl, xu
    real(kind = RPRE), dimension(:,:), intent(in), optional :: Xi
    integer(kind = IPRE), intent(in), optional :: np, itermax
    real(kind = RPRE), intent(in), optional :: w, l, alpha, delta, phi1, phi2, &
      eps1, eps2, eps3
    real(kind = RPRE), intent(inout), optional :: fit
    integer(kind = IPRE), intent(inout), optional :: niter, neval, flag
    real(kind = RPRE), dimension(:,:,:), allocatable, intent(inout), optional :: extra
    logical, intent(in), optional :: snap
    character(len = *), intent(in), optional :: snapdir

    ! Local variables
    !=================
    integer(kind = IPRE) :: opt_np, opt_itermax, nd, i, counteval, &
      iter, gbidx(1), ne
    real(kind = RPRE) :: opt_w, opt_l, opt_alpha, opt_delta, opt_phi1, opt_phi2, &
      opt_eps1, opt_eps2, opt_eps3, pfit, gfit
    integer(kind = IPRE), dimension(:), allocatable :: idx
    real(kind = RPRE), dimension(:), allocatable :: pbestfit, gbest, vu, vl
    real(kind = RPRE), dimension(:,:), allocatable :: X, V, pbest, r1, r2
    real(kind = RPRE), dimension(:,:,:), allocatable :: opt_extra
    logical :: opt_snap, converge, reset
    character(len = :), allocatable :: opt_snapdir
    type(optifit) :: tmpfit

    ! Optional arguments default values
    !===================================
    opt_np = 100                ! Number of particles
    opt_itermax = 100           ! Maximum number of iterations
    opt_w = 0.7d0               ! Particles inertia
    opt_l = 0.1d0               ! Velocity clamping
    opt_phi1 = 1.5d0            ! Cognition parameter
    opt_phi2 = 1.5d0            ! Global sociability parameter
    opt_alpha = 1.25d0          ! Competitivity parameter
    opt_eps1 = 1.0d-8           ! Minimum change in best position
    opt_eps2 = 1.0d-8           ! Minimum fitness precision
    opt_eps3 = 1.0d-8           ! Minimum change in global best position fitness
    opt_snap = .false.          ! Save particles positions after each iteration
    opt_snapdir = "snap/"       ! Snapshot directory
    if (present(np)) opt_np = np
    if (present(itermax)) opt_itermax = itermax
    if (present(w)) opt_w = w
    if (present(l)) opt_l = l
    if (present(phi1)) opt_phi1 = phi1
    if (present(phi2)) opt_phi2 = phi2
    if (present(alpha)) opt_alpha = alpha
    if (present(eps1)) opt_eps1 = eps1
    if (present(eps2)) opt_eps2 = eps2
    if (present(eps3)) opt_eps3 = eps3
    if (present(snap)) opt_snap = snap
    if (present(snapdir)) opt_snapdir = snapdir

    ! Swarm maximum radius
    if (present(delta)) then
      opt_delta = delta
    else
      opt_delta = 0.08686d0 * log( 1.0d0 + 0.004d0 * opt_np )
    end if

    ! Particles initial positions
    !=============================
    nd = size(xl)               ! Number of dimensions
    if ( present(Xi) ) then
      X = Xi
      opt_np = size(Xi, 2)
      if ( present(np) .and. (np .ne. opt_np) ) then
        print *, "Warning: in cpso, np differs from size(Xi, 2), np " &
          // "set to size(Xi, 2) (" // num2str(opt_np) // ")."
      end if
    else
      X = repmat(xl, opt_np) + randu(nd, opt_np)*repmat(xu - xl, opt_np)
    end if

    ! Initialize swarm
    !==================
    V = randu(nd, opt_np)       ! Particles velocities
    pbest = X                   ! Particles personal best positions
    pbestfit = zeros(opt_np)    ! Particles personal best positions fitness
    vu = opt_l*(xu - xl)        ! Velocity upper boundary
    vl = -vu                    ! Velocity lower boundary
    counteval = 0               ! Fitness evaluations counter

    ! Initialize particle velocity
    !==============================
    do i = 1, opt_np
      V(:,i) = vl + V(:,i)*(vu - vl)
    end do

    ! Compute fitness
    !=================
    pbestfit = zeros(opt_np)
    do i = 1, opt_np
      tmpfit = fitness(X(:,i))
      pbestfit(i) = tmpfit % val
      if (present(extra)) then
        if (.not. allocated(opt_extra)) then
          if (.not. allocated(tmpfit % extra)) then
            print *, "Error: in cpso, extra requested but not defined in fitness."
            stop
          else
            ne = size(tmpfit % extra)
            opt_extra = zeros(ne, opt_np, opt_itermax)
          end if
        end if
        opt_extra(:,i,1) = tmpfit % extra
      end if
      counteval = counteval + 1
    end do

    ! Initialize global best
    !========================
    gbidx = minloc(pbestfit)
    gfit = pbestfit(gbidx(1))
    gbest = pbest(:,gbidx(1))

    ! Iterate until one of the termination criterion is satisfied
    !=============================================================
    iter = 1
    reset = .false.
    converge = .false.

    if ( opt_snap ) then
      call system("rm -rf " // opt_snapdir)
      call snap_positions(X, pbestfit, iter, opt_snapdir)
    end if

    do while ( .not. converge )
      iter = iter + 1
      r1 = randu(nd, opt_np)
      r2 = randu(nd, opt_np)
      
      do i = 1, opt_np
      
        ! Update particle velocity
        !==========================
        V(:,i) = opt_w*V(:,i) + opt_phi1*r1(:,i)*(pbest(:,i) - X(:,i)) &
                              + opt_phi2*r2(:,i)*(gbest - X(:,i))

        ! Update particle position
        !==========================
        X(:,i) = X(:,i) + V(:,i)
        X(:,i) = merge(X(:,i), xl, X(:,i) .ge. xl)
        X(:,i) = merge(X(:,i), xu, X(:,i) .le. xu)
        tmpfit = fitness(X(:,i))
        pfit = tmpfit % val
        if (present(extra)) opt_extra(:,i,iter) = tmpfit % extra
        counteval = counteval + 1
        
        ! Update particle best position
        !===============================
        if ( pfit .lt. pbestfit(i) ) then
          pbest(:,i) = X(:,i)
          pbestfit(i) = pfit
        end if
        
        ! Update global best position
        !=============================
        if ( pfit .lt. gfit ) then
        
          ! Stop if global best position changes less than eps1
          if ( (norm(gbest - X(:,i)) .le. opt_eps1) &
               .and. (pfit .le. opt_eps2) &
               .and. (opt_np .gt. 1) ) then
            converge = .true.
            cpso = X(:,i)
            if (present(fit)) fit = pfit
            if (present(flag)) flag = 0
            
          ! Stop if best fitness changes less than eps3
          elseif ( (abs(gfit - pfit) .le. opt_eps3) &
                   .and. (pfit .le. opt_eps2) ) then
            converge = .true.
            cpso = X(:,i)
            if (present(fit)) fit = pfit
            if (present(flag)) flag = 2
          
          ! Otherwise, update global best position
          else
            gbest = X(:,i)
            gfit = pfit
          end if
          
        end if
        
        if ( converge ) exit
      end do
      
      ! Stop if maximum iteration is reached
      !======================================
      if ( iter .ge. opt_itermax ) then
        converge = .true.
        cpso = gbest
        if (present(fit)) fit = gfit
        if (present(flag)) flag = 1
      end if
      
      if ( .not. converge ) then
        if ( opt_alpha .gt. 0.0d0 ) call update_reset()
        if ( opt_snap ) call snap_positions(X, pbestfit, iter, opt_snapdir)
      end if
    end do
    
    if (present(niter)) niter = iter
    if (present(neval)) neval = counteval
    if (present(extra)) extra = opt_extra(:,:,:iter)
    if ( opt_snap ) then
      call export_parameters(opt_np, nd, iter, counteval, xl, xu, opt_snapdir)
      call snap_positions(X, pbestfit, iter, opt_snapdir)
    end if

    return
  contains

    !-------------------------------------------------------------------
    ! Subroutine update_reset
    !-------------------------------------------------------------------
    subroutine update_reset()
      integer(kind = IPRE) :: i, nw
      real(kind = RPRE) :: inorm, ls
      real(kind = RPRE), dimension(:), allocatable :: dist, newfit
      integer(kind = IPRE), dimension(:), allocatable :: idx
      real(kind = RPRE), dimension(:,:), allocatable :: newextra

      ! Evaluate swarm size
      !=====================
      dist = zeros(opt_np)
      do i = 1, opt_np
        dist(i) = norm(X(:,i) - gbest)
      end do
      dist = dist / norm(xu - xl)

      ! Restart particles if swarm size lower than threshold
      !======================================================
      if ( maxval(dist) .lt. opt_delta ) then
        reset = .true.

        ! Rank particles
        !================
        inorm = real(iter) / real(opt_itermax)
        ls = -1.0d0 / 0.09d0
        nw = ( opt_np - 1 ) * 1.0d0 / ( 1.0d0 + exp( -ls * ( inorm - opt_alpha + 0.5d0 ) ) )
        idx = argsort(pbestfit, 2)
        idx = idx(:nw)

        ! Reset "worst" particles velocities
        !====================================
        ! Reset position, velocity and personal best
        do i = 1, nw
          V(:,idx(i)) = vl + randu(nd)*(vu - vl)
          X(:,idx(i)) = xl + randu(nd)*(xu - xl)
          pbest(:,idx(i)) = X(:,idx(i))
        end do

        ! Reset personal best fitness
        newfit = zeros(nw)
        if (present(extra)) newextra = zeros(ne, nw)
        do i = 1, nw
          tmpfit = fitness(pbest(:,idx(i)))
          newfit(i) = tmpfit % val
          if (present(extra)) newextra(:,i) = tmpfit % extra
          counteval = counteval + 1
        end do
        pbestfit(idx) = newfit
        if (present(extra)) opt_extra(:,idx,iter) = newextra

      else
        reset = .false.
      end if
      return
    end subroutine update_reset

  end function cpso

!=======================================================================
! de
!-----------------------------------------------------------------------
! de computes the nd parameters that minimize the input function using
! Differential Evolution.
!
! Syntax
!-----------------------------------------------------------------------
! x = de(fitness, xl, xu)
! x = de(fitness, xl, xu, [options = ])
!
! Description
!-----------------------------------------------------------------------
! x = de(fitness, xl, xu) returns a nd vector with the parameters that
! minimizes the function fitness, given the lower and upper boundary xl
! and xu.
!
! Inputs
!-----------------------------------------------------------------------
! fitness             Function to be minimized (real)
! xl                  Lower boundary (real nd array)
! xu                  Upper boundary (real nd array)
!
! Options
!-----------------------------------------------------------------------
! Xi                  Population initial positions
! np = 100            Number of individuals
! itermax = 100       Maximum number of iterations
! F = 1.0d0           Differential weight
! CR = 0.5d0          Crossover probability
! eps1 = 1.0d-8       Minimum change in best position
! eps2 = 1.0d-8       Minimum fitness precision
! fit                 Best individual position fitness
! niter               Output number of iterations
! neval               Output number of fitness evaluations
! flag                Output reason of termination
!                       -   0   Best position changes less than eps1,
!                       -   1   Maximum number of iterations is reached.
! extra               Additional outputs from fitness function
! snap = .false.      Save individuals positions after each iteration
! snapdir = "snap/"   Snapshot directory
!=======================================================================

  function de(fitness, xl, xu, Xi, np, itermax, F, CR, eps1, eps2, &
              fit, niter, neval, flag, extra, snap, snapdir)

    ! Input arguments
    !=================
    real(kind = RPRE), dimension(:), allocatable :: de
    procedure(optifunc) :: fitness
    real(kind = RPRE), dimension(:), intent(in) :: xl, xu
    real(kind = RPRE), dimension(:,:), intent(in), optional :: Xi
    integer(kind = IPRE), intent(in), optional :: np, itermax
    real(kind = RPRE), intent(in), optional :: F, CR, eps1, eps2
    real(kind = RPRE), intent(inout), optional :: fit
    integer(kind = IPRE), intent(inout), optional :: niter, neval, flag
    real(kind = RPRE), dimension(:,:,:), allocatable, intent(inout), optional :: extra
    logical, intent(in), optional :: snap
    character(len = *), intent(in), optional :: snapdir

    ! Local variables
    !=================
    integer(kind = IPRE) :: opt_np, opt_itermax, nd, i, irand, &
      counteval, iter, gbidx(1), ne
    real(kind = RPRE) :: opt_F, opt_CR, opt_eps1, opt_eps2, gfit
    integer(kind = IPRE), dimension(:), allocatable :: idx
    real(kind = RPRE), dimension(:), allocatable :: pfit, &
      pbestfit, gbest, x1, x2, x3
    real(kind = RPRE), dimension(:,:), allocatable :: X, V, U, r1
    real(kind = RPRE), dimension(:,:,:), allocatable :: opt_extra
    logical :: opt_snap, converge
    character(len = :), allocatable :: opt_snapdir
    type(optifit) :: tmpfit

    ! Optional arguments default values
    !===================================
    opt_np = 100                ! Number of individuals
    opt_itermax = 100           ! Maximum number of iterations
    opt_F = 1.0d0               ! Differential weight
    opt_CR = 0.5d0              ! Crossover probability
    opt_eps1 = 1.0d-8           ! Minimum change in best position
    opt_eps2 = 1.0d-8           ! Minimum fitness precision
    opt_snap = .false.          ! Save particles positions after each iteration
    opt_snapdir = "snap/"       ! Snapshot directory
    if (present(np)) opt_np = np
    if (present(itermax)) opt_itermax = itermax
    if (present(F)) opt_F = F
    if (present(CR)) opt_CR = CR
    if (present(eps1)) opt_eps1 = eps1
    if (present(eps2)) opt_eps2 = eps2
    if (present(snap)) opt_snap = snap
    if (present(snapdir)) opt_snapdir = snapdir

    ! Population initial positions
    !==============================
    nd = size(xl)               ! Number of dimensions
    if ( present(Xi) ) then
      X = Xi
      opt_np = size(Xi, 2)
      if ( present(np) .and. (np .ne. opt_np) ) then
        print *, "Warning: in de, np differs from size(Xi, 2), np " &
          // "set to size(Xi, 2) (" // num2str(opt_np) // ")."
      end if
    else
      X = repmat(xl, opt_np) + randu(nd, opt_np)*repmat(xu - xl, opt_np)
    end if

    ! Initialize population
    !=======================
    V = zeros(nd, opt_np)       ! Donor vectors
    U = zeros(nd, opt_np)       ! Trial vectors
    pfit = zeros(opt_np)        ! Individuals current fitness
    counteval = 0               ! Fitness evaluations counter

    ! Compute fitness
    !=================
    pfit = zeros(opt_np)
    do i = 1, opt_np
      tmpfit = fitness(X(:,i))
      pfit(i) = tmpfit % val
      if (present(extra)) then
        if (.not. allocated(opt_extra)) then
          if (.not. allocated(tmpfit % extra)) then
            print *, "Error: in de, extra requested but not defined in fitness."
            stop
          else
            ne = size(tmpfit % extra)
            opt_extra = zeros(ne, opt_np, opt_itermax)
          end if
        end if
        opt_extra(:,i,1) = tmpfit % extra
      end if
      counteval = counteval + 1
    end do
    pbestfit = pfit

    ! Initialize best individual
    !============================
    gbidx = minloc(pbestfit)
    gfit = pbestfit(gbidx(1))
    gbest = X(:,gbidx(1))

    ! Iterate until one of the termination criterion is satisfied
    !=============================================================
    iter = 1
    converge = .false.
    
    if ( opt_snap ) then
      call system("rm -rf " // opt_snapdir)
      call snap_positions(X, pbestfit, iter, opt_snapdir)
    end if

    do while ( .not. converge )
      iter = iter + 1
      r1 = randu(nd, opt_np)

      ! Mutation
      !==========
      do i = 1, opt_np
        idx = randperm(opt_np, 4)
        idx = idx(find(idx .ne. i))
        x1 = X(:,idx(1))
        x2 = X(:,idx(2))
        x3 = X(:,idx(3))
        V(:,i) = x1 + opt_F*(x2 - x3)
        V(:,i) = merge(V(:,i), xl, V(:,i) .ge. xl)
        V(:,i) = merge(V(:,i), xu, V(:,i) .le. xu)
      end do

      ! Recombination
      !===============
      irand = randi(nd)
      do i = 1, opt_np
        U(:,i) = merge( V(:,i), X(:,i), r1(:,i) .le. opt_CR &
                                        .or. [ ( i, i = 1, nd ) ] .eq. irand )
      end do

      ! Compute fitness
      !=================
      pfit = zeros(opt_np)
      do i = 1, opt_np
        tmpfit = fitness(U(:,i))
        pfit(i) = tmpfit % val
        if (present(extra)) opt_extra(:,i,iter) = tmpfit % extra
        counteval = counteval + 1
      end do

      ! Selection
      !===========
      idx = find( pfit .lt. pbestfit )
      pbestfit(idx) = pfit(idx)
      X(:,idx) = U(:,idx)

      ! Update best individual
      !========================
      gbidx = minloc(pbestfit)

      ! Stop if best individual position changes less than eps1
      !=========================================================
      if ( (norm(gbest - X(:,gbidx(1))) .le. opt_eps1) &
           .and. (pbestfit(gbidx(1)) .le. opt_eps2) &
           .and. (opt_np .gt. 1) ) then
        converge = .true.
        de = X(:,gbidx(1))
        if (present(fit)) fit = pbestfit(gbidx(1))
        if (present(flag)) flag = 0

      ! Stop if maximum iteration is reached
      !======================================
      elseif ( iter .ge. opt_itermax ) then
        converge = .true.
        de = X(:,gbidx(1))
        if (present(fit)) fit = pbestfit(gbidx(1))
        if (present(flag)) flag = 1

      ! Otherwise, update best individual
      !===================================
      else
        gbest = X(:,gbidx(1))
        gfit = pbestfit(gbidx(1))
      end if

      if ( .not. converge .and. opt_snap ) then
        call snap_positions(X, pbestfit, iter, opt_snapdir)
      end if
    end do
    
    if (present(niter)) niter = iter
    if (present(neval)) neval = counteval
    if (present(extra)) extra = opt_extra(:,:,:iter)
    if ( opt_snap ) then
      call export_parameters(opt_np, nd, iter, counteval, xl, xu, opt_snapdir)
      call snap_positions(X, pbestfit, iter, opt_snapdir)
    end if

    return
  end function de

!=======================================================================
! export_parameters
!=======================================================================

  subroutine export_parameters(np, nd, niter, neval, xl, xu, snapdir)
    integer(kind = IPRE), intent(in) :: np, nd, niter, neval
    real(kind = RPRE), dimension(:), intent(in) :: xl, xu
    character(len = :), allocatable, intent(in) :: snapdir
    integer(kind = IPRE) :: i
    type(File) :: parameters

    parameters = File(999, snapdir // "parameters.txt")
    call parameters%open()
    write(parameters%unit, *) "# Number of individuals np"
    write(parameters%unit, *) num2str(np)
    write(parameters%unit, *) "# Number of dimensions nd"
    write(parameters%unit, *) num2str(nd)
    write(parameters%unit, *) "# Number of iterations niter"
    write(parameters%unit, *) num2str(niter)
    write(parameters%unit, *) "# Number of iterations neval"
    write(parameters%unit, *) num2str(neval)
    write(parameters%unit, *) "# Lower boundaries xl"
    write(parameters%unit, *) ( xl(i), i = 1, nd )
    write(parameters%unit, *) "# Upper boundaries xu"
    write(parameters%unit, *) ( xu(i), i = 1, nd )
    call parameters%close()
    return
  end subroutine export_parameters

!=======================================================================
! snap_positions
!=======================================================================

  subroutine snap_positions(x, fit, iter, snapdir)
    real(kind = RPRE), dimension(:,:), intent(in) :: x
    real(kind = RPRE), dimension(:), intent(in) :: fit
    integer(kind = IPRE), intent(in) :: iter
    character(len = *), intent(in) :: snapdir

    ! Create snap directory
    !=======================
    call system("mkdir -p " // snapdir)

    ! Save individuals positions and fitness
    !========================================
    call savebin(snapdir // "xiter" // num2str(iter) // ".bin", x)
    call savebin(snapdir // "fititer" // num2str(iter) // ".bin", fit)
    return
  end subroutine snap_positions

!=======================================================================
! Test functions
!=======================================================================

  !---------------------------------------------------------------------
  ! Function Ackley
  !---------------------------------------------------------------------
  type(optifit) function ackley(x)
    real(kind = RPRE), dimension(:), intent(in) :: x
    integer(kind = IPRE) :: nd
    real(kind = RPRE) :: sum1, sum2, e = 2.7182818284590451d0

    nd = size(x)
    sum1 = sqrt( 1.0d0 / nd * sum( x**2 ) )
    sum2 = 1.0d0 / nd * sum( cos( 2.0d0 * pi * x ) )
    ackley % val = 20.0d0 + e - 20.0d0 * exp( -0.2d0 * sum1 ) - exp(sum2)
    return
  end function ackley

  !---------------------------------------------------------------------
  ! Function Griewank
  !---------------------------------------------------------------------
  type(optifit) function griewank(x)
    real(kind = RPRE), dimension(:), intent(in) :: x
    integer(kind = IPRE) :: nd
    real(kind = RPRE) :: sum1, prod1

    nd = size(x)
    sum1 = sum( x**2 ) / 4000.0d0
    prod1 = product( cos( x / sqrt( linspace(1, nd, nd) ) ) )
    griewank % val = 1.0d0 + sum1 - prod1
    return
  end function griewank

  !---------------------------------------------------------------------
  ! Function Quartic
  !---------------------------------------------------------------------
  type(optifit) function quartic(x)
    real(kind = RPRE), dimension(:), intent(in) :: x
    integer(kind = IPRE) :: nd

    nd = size(x)
    quartic % val = sum( linspace(1, nd, nd) * x**4 )
    return
  end function quartic

  !---------------------------------------------------------------------
  ! Function Quartic with noise
  !---------------------------------------------------------------------
  type(optifit) function quartic_noise(x)
    real(kind = RPRE), dimension(:), intent(in) :: x
    type(optifit) :: tmpfit

    tmpfit = quartic(x)
    quartic_noise % val = tmpfit % val + randu()
    return
  end function quartic_noise

  !---------------------------------------------------------------------
  ! Function Rastrigin
  !---------------------------------------------------------------------
  type(optifit) function rastrigin(x)
    real(kind = RPRE), dimension(:), intent(in) :: x
    integer(kind = IPRE) :: nd
    real(kind = RPRE) :: sum1

    nd = size(x)
    sum1 = sum( x**2 - 10.0d0 * cos( 2.0d0 * pi * x ) )
    rastrigin % val = 10.0d0 * nd + sum1
    return
  end function rastrigin

  !---------------------------------------------------------------------
  ! Function Rosenbrock
  !---------------------------------------------------------------------
  type(optifit) function rosenbrock(x)
    real(kind = RPRE), dimension(:), intent(in) :: x
    integer(kind = IPRE) :: nd
    real(kind = RPRE) :: sum1, sum2

    nd = size(x)
    sum1 = sum( ( x(2:) - x(:nd-1)**2 )**2 )
    sum2 = sum( ( 1.0d0 - x(:nd-1) )**2 )
    rosenbrock % val = 100.0d0 * sum1 + sum2
    return
  end function rosenbrock

  !---------------------------------------------------------------------
  ! Function Schwefel
  !---------------------------------------------------------------------
  type(optifit) function schwefel(x)
    real(kind = RPRE), dimension(:), intent(in) :: x
    integer(kind = IPRE) :: nd
    real(kind = RPRE) :: sum1

    nd = size(x)
    sum1 = sum( x * sin( sqrt( abs(x) ) ) )
    schwefel % val = 418.9829d0 * nd - sum1
    return
  end function schwefel

  !---------------------------------------------------------------------
  ! Function Sphere
  !---------------------------------------------------------------------
  type(optifit) function sphere(x)
    real(kind = RPRE), dimension(:), intent(in) :: x

    sphere % val = sum( x**2 )
    return
  end function sphere

  !---------------------------------------------------------------------
  ! Function Styblinski-Tang
  !---------------------------------------------------------------------
  type(optifit) function styblinski_tang(x)
    real(kind = RPRE), dimension(:), intent(in) :: x
    real(kind = RPRE) :: sum1

    sum1 = sum( x**4 - 16.0d0 * x**2 + 5.0d0 * x )
    styblinski_tang % val = sum1 / 2.0d0
    return
  end function styblinski_tang

end module optimizers
