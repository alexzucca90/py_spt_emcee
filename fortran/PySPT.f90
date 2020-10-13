!==============================================================
! Module to use SPT pol 500d 2019 B mode data.
! Based on SPT publicly available likelihood
! at www. ......
! and cliklike.f90
!
! Version: 0.1
!
! Written by Alex Zucca, azucca@sfu.ca
!
!   The code is slightly modified, and it's used in older versions
!   of CosmoMC
!
!
!==============================================================

module PySPT_mod

    implicit none

contains

    !------------------------------------------------------------------
    !> modified version of the SPT likelihood
    function SPT_LnLike(dls_lcdm, dls_pmf, dls_poiss, dls_galdust, &
                            &Add, PoissonLevels, MapBcal150, MapBcal90, &
                            &BeamFactors, beam_err, spec, windows, cov) result(LnLike)

        implicit none
        !> IO parameters
        double precision, intent(in) :: dls_lcdm(2998)
        double precision, intent(in) :: dls_pmf(2998)
        double precision, intent(in) :: dls_poiss(2998)
        double precision, intent(in) :: dls_galdust(2998)
        double precision, intent(in) :: Add
        double precision, intent(in) :: PoissonLevels(3)
        double precision, intent(in) :: MapBcal150, MapBcal90
        double precision, intent(in) :: BeamFactors(7)
        double precision, intent(in) :: beam_err(21,7)
        double precision, intent(in) :: spec(7,3)
        double precision, intent(in) :: windows(2998, 21)
        double precision, intent(inout) :: cov(21, 21)
        double precision :: LnLike

        !> internal parameters (from SPT likelihood)
        integer, parameter :: nband = 3
        integer, parameter :: nall = 21
        integer, parameter :: nfreq = 2
        integer, parameter :: nbin = 7
        integer, parameter :: spt_windows_lmax = 3000
        integer, parameter :: spt_windows_lmin = 3
        double precision :: CalFactors(3)
        double precision, parameter :: meanAdd = 0.0094d0
        double precision, parameter :: sigmaAdd = 0.0021d0

        integer :: i, j, k
        double precision, dimension(1:21) :: deltacb, tmp2cb, BeamFac
        double precision, dimension(1:7) :: tmpcb
        double precision, dimension(1) :: junk, detcov
        double precision, dimension(2998) :: dl_fgs, Dls_dust150ghz
        integer, parameter :: N_BEAM_EXPECTED = 7
        integer, parameter :: neff = 21
        integer, parameter :: n_beam_terms = 7
        double precision :: PriorLnLike
        double precision :: InvCalCov(2,2)
        double precision :: y1, y2
        double precision :: effFreqs(2, 3)
        double precision :: eff_dust_freqs(2)

        !print*, 'dls_lcdm fortran:', dls_lcdm(1), dls_lcdm(2998)
        !print*, 'dls_pmf fortran:', dls_pmf(1), dls_pmf(2998)

        !> set up eff_dust_freq
        eff_dust_freqs(1) = 96.2
        eff_dust_freqs(2) = 149.5


        !> set up the calibration inverse covariance
        InvCalCov(1,1) = 367.669099005
        InvCalCov(1,2) = -367.669099005
        InvCalCov(2,1) = InvCalCov(1,2)
        InvCalCov(2,2) = 38911.2850797

        k=0
        !> set up efffreq
        do i=1,nfreq
            do j=i,nfreq
                k=k+1
                effFreqs(1,k) = eff_dust_freqs(i)
                effFreqs(2,k) = eff_dust_freqs(j)
            enddo
        enddo

        !> prepare the calibration factors
        CalFactors(1) = MapBcal150*MapBcal150
        CalFactors(2) = MapBcal150*MapBcal90
        CalFactors(3) = MapBcal90*MapBcal90

        !> multiply dust emission by constant
        Dls_dust150ghz(:) = Add * dls_galdust(:)

        tmpcb(:)=0
        do k=1,nband !150x150, 95x150, 95x95ghz
            !First get model foreground spectrum (in Dl).
            !Note all the foregrounds are recorded in Dl at l=3000, so we
            !divide by d3000 to get to a normalized Cl spectrum.
            !Start with Poisson power
            dl_fgs(:) = PoissonLevels(k) * dls_poiss(:)

            !Add galactic dust
            dl_fgs(:) = dl_fgs(:) + Dls_dust150ghz(:) * dustFreqScalingFrom150GHz(effFreqs(:,k))

            !Now add model CMB: LCDM + PMF
            dl_fgs(:) = dl_fgs(:) + dls_lcdm(:) + dls_pmf(:)

            !Now bin into bandpowers with the window functions.
            call dgemv('T',spt_windows_lmax-spt_windows_lmin+1,nbin,1.d0,&
            windows(:,1+nbin*(k-1):nbin+nbin*(k-1)),spt_windows_lmax-spt_windows_lmin+1,&
            dl_fgs,1,0,tmpcb,1)

            !scale theory spectrum by calibration:
            tmpcb = tmpcb(:) / CalFactors(k)

            tmp2cb(1+(k-1)*nbin:nbin+(k-1)*nbin) = tmpcb(:)
        end do

        BeamFac = 1.d0

        !> adjust with the beam factors
        do i=1,N_BEAM_EXPECTED
            BeamFac = BeamFac(:) * (1.d0 + beam_err(:,i) * BeamFactors(i))
        enddo

        deltacb(:) = tmp2cb(:) * BeamFac(:)
        !deltacb(:) = tmp2cb(:)

        !> subtract the measurements
        do k=1,nband
            deltacb(1+(k-1)*nbin:nbin+(k-1)*nbin) = deltacb(1+(k-1)*nbin:nbin+(k-1)*nbin) - spec(:,k)
        enddo

        LnLike =  Matrix_GaussianLogLikeChol(nall, cov, deltacb)

        !> add the priors
        !beam prior
        PriorLnLike =  0.5 * sum(BeamFactors**2)

        !> add the dust prior
        PriorLnLike = PriorLnLike + 0.5*((Add - meanAdd) /sigmaAdd)**2

        !> add the calibration prior
        y1 = log(MapBcal90)
        y2 = log(MapBcal150)
        !hardwired 2x2 matrix multiply:
        ! 0.5 yt C^-1 y
        PriorLnLike = PriorLnLike + 0.5*( InvCalCov(1,1) *y1*y1 + &
            & 2 * InvCalCov(1,2) * y1*y2 + &
            & InvCalCov(2,2) * y2*y2 )

        LnLike = LnLike + PriorLnLike

    end function SPT_LnLike

    !-------------------------------------------------------------------
    !
    !   OTHER FUNCTIONS FOR SPT LIKELIHOOD, keep them in f2py
    !
    !-------------------------------------------------------------------

    !-------------------------------------------------------------------
    !> nu,nu0 in GHz
    !
    !> dBdT is proportional to derivative of planck function
    !  but is normalized so its equal to 1 at nu0
    !
    function dBdT(nu,nu0)

        double precision :: x, x0, dBdT, dBdT0, nu, nu0

        x0 = nu0/56.78
        dBdT0 = x0**4 * exp(x0) / (exp(x0)-1)**2

        x = nu/56.78
        dBdT = x**4 * exp(x) / (exp(x)-1)**2 / dbdT0

    end function dBdT


    !-------------------------------------------------------------------
    !> proportional to the Planck function normalized to 1 at nu0
    function Bnu(nu,nu0,T)
        double precision :: Bnu, nu,nu0,T
        !h/k
        !4.799237 x 10-11 s K
        !expect GHz
        ! so 4.799237e-2 K/GHz
        double precision, parameter :: hk = 4.799237e-2

        Bnu = (nu/nu0)**3
        Bnu = Bnu * (exp( hk*nu0/T)-1) / (exp( hk*nu/T)-1)

    end function Bnu


    !-------------------------------------------------------------------
    !> dust scaling by frequency from 150x150Ghz
    function dustFreqScalingFrom150GHz(effFreqs)

        double precision, dimension(2), intent(in) :: effFreqs
        double precision dustFreqScalingFrom150GHz
        double precision, parameter :: beta = 1.59, Tdust=19.6 !Kelvin

        dustFreqScalingFrom150GHz = ((effFreqs(1)*effFreqs(2))/(150.d0*150.d0))**beta
        dustFreqScalingFrom150GHz = dustFreqScalingFrom150GHz * &
        Bnu(effFreqs(1),150.d0,Tdust) * Bnu(effFreqs(2),150.d0,Tdust)
        dustFreqScalingFrom150GHz = dustFreqScalingFrom150GHz / &
        dBdT(effFreqs(1),150.d0)/ dBdT(effFreqs(2),150.d0)

    end function dustFreqScalingFrom150GHz


    !-------------------------------------------------------------------
    !> generate t
    function power_law(index,fr0,freqPair)
        double precision, intent(in) :: index,fr0,freqPair(2)
        double precision ::   fri,frj
        double precision :: power_law
        fri=freqPair(1)
        frj=freqPair(2)
        power_law = 1.d0/dBdT(fri,fr0)/dBdT(frj,fr0)*&
        (fri/fr0*frj/fr0)**(index)
    end function power_law

    !-------------------------------------------------------------------
    !> Returns -Log Likelihood for Gaussian: (d^T Cov^{-1} d + log|Cov|)/2
    !  Cov is already the cholesky factorization of the covariance matrix
    !
    function Matrix_GaussianLogLikeChol(n, Cov, d) result(LogLike)
        implicit none
        integer, intent(in) :: n
        double precision, intent(inout):: Cov(n,n)
        double precision, intent(in):: d(n)
        double precision :: tmp(n)
        double precision :: LogLike
        integer info,i

        LogLike = 0
        !Log Det term:
        do i=1, n
            !Don't divide det term by 2 since we're dealing with the Cholesky-factored cov,
            !not the whole matrix.
            LogLike = LogLike  + log(Cov(i,i))
        end do

        !Solve for Cov^{-1}d [could use faster symmetric method]
        tmp = d

        !> adjust to single precision
        call DPOTRS('L', n, 1, Cov, n, tmp, n, INFO )
        !call SPOTRS('L', n, 1, Cov, n, tmp, n, INFO )
        if (INFO/=0) then
            print*, 'Matrix_GaussianLogLikeCholDouble: error in solving for cov^{-1}d'
            stop
        end if

        !Add together
        LogLike = LogLike + dot_product(tmp,d)/2.d0

    end function Matrix_GaussianLogLikeChol

    !-------------------------------------------------------------------
    !> Cholesky decomposition of a matrix
    subroutine Matrix_Cholesky(n, M)
        !Note upper triangular is not zeroed
        integer, intent(in) :: n
        double precision, intent(inout):: M(n,n)
        integer info

        !> call LAPACK subroutine
        call dpotrf ('L', n, M, n, info)

        if (info/=0) then
            print*, 'Matrix_Cholesky: not positive definite '
            stop
        end if

    end subroutine Matrix_Cholesky


end module PySPT_mod
