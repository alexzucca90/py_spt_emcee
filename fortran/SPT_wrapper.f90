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
!==============================================================



!------------------------------------------------------------------
!> wrapper for Likelihood
subroutine SPT_LnLike_wrapper(dls_lcdm, dls_pmf, dls_poiss, Dls_galdust, &
                                &Add, PoissonLevels, MapBcal150, MapBcal90, &
                                &BeamFactors, beam_err, spec, windows, cov, LnLike)
    use PySPT_mod
    implicit none

    double precision, intent(in) :: dls_lcdm(2998)
    double precision, intent(in) :: dls_pmf(2998)
    double precision, intent(in) :: dls_poiss(2998)
    double precision, intent(in) :: Dls_galdust(2998)
    double precision, intent(in) :: Add
    double precision, intent(in) :: PoissonLevels(3)
    double precision, intent(in) :: MapBcal150, MapBcal90
    double precision, intent(in) :: BeamFactors(7)
    double precision, intent(in) :: beam_err(21,7)
    double precision, intent(in) :: spec(7,3)
    double precision, intent(in) :: windows(2998, 21)
    double precision, intent(inout) :: cov(21, 21)
    double precision, intent(out) :: LnLike

    !f2py intent(in) :: dls_lcdm
    !f2py intent(in) :: dls_pmf
    !f2py intent(in) :: dls_poiss
    !f2py intent(in) :: Dls_galdust
    !f2py intent(in) :: Add
    !f2py intent(in) :: PoissonLevels
    !f2py intent(in) :: MapBcal150, MapBcal90
    !f2py intent(in) :: BeamFactors
    !f2py intent(in) :: beam_err
    !f2py intent(in) :: spec
    !f2py intent(in) :: windows
    !f2py intent(inout) :: cov
    !f2py intent(out) :: LnLike


    LnLike = SPT_LnLike(dls_lcdm, dls_pmf, dls_poiss, Dls_galdust, &
        &Add, PoissonLevels, MapBcal150, MapBcal90, &
        &BeamFactors, beam_err, spec, windows, cov)

end subroutine SPT_LnLike_wrapper

!------------------------------------------------------------------
!> wrapper for Cholesky Matrix decomposition
subroutine Matrix_Cholesky_wrapper(n, M, Mchol)
    use PySPT_mod
    implicit none
    integer, intent(in) :: n
    double precision, intent(inout) :: M(n,n)
    double precision, intent(out) :: Mchol(n, n)

    !f2py intent(in) :: M
    !f2py intent(hide), depend(M) :: n=shape(M,1)
    !f2py intent(out) :: Mchol

    call Matrix_Cholesky(n, M)
    Mchol = M

end subroutine Matrix_Cholesky_wrapper


