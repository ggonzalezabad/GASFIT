MODULE GASFIT_MODULE

! General control input variables
  CHARACTER(256) :: fitin, fitout, specin, specout, inputline, general_line
  LOGICAL :: wrt_scr, write_fit, write_spec, mirror, autodiff, if_residuals
  INTEGER*4 :: npoints, nfirstfit
  REAL*8 :: delchi, provar, automult
  
! Parameters
  CHARACTER(3), PARAMETER :: ad1_str='ad1', ad2_str='ad2', ble_str='ble', &
       scp_str='Scp', bsp_str='Bsp', alb_str='ALB', hwe_str='HWE', &
       shi_str='SHI', sqe_str='SQE', sha_str='SHA', asy_str='ASY', &
       us1_str='us1', us2_str='us2'

! High resolution solar spectrum variables
  REAL*8, ALLOCATABLE, DIMENSION(:) :: kppos, kpspec, kppos_ss, kpspec_gauss
  INTEGER*4 :: nkppos

! Solar fitting control & variables (read or derived from input control file)
  CHARACTER(256) :: solar_line
  CHARACTER(30), ALLOCATABLE, DIMENSION(:) :: sun_par_names
  CHARACTER(3), ALLOCATABLE, DIMENSION(:) :: sun_par_str
  LOGICAL :: iterate_sun, weight_sun
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: if_var_sun
  INTEGER*4 :: n_solar_pars, ll_sun, lu_sun, nvar_sun, &
       nalb, nhwe, nshi, nsqe, nsha, nasy
  INTEGER*4, ALLOCATABLE, DIMENSION(:) :: list_sun
  REAL*8 :: div_sun
  REAL*8, ALLOCATABLE, DIMENSION(:) :: var_sun, diffsun, init_sun, var_sun_factor

! Solar spectrum
  REAL*8, ALLOCATABLE, DIMENSION(:) :: pos_sun, spec_sun, sig_sun

! Database of reference cross sections
  REAL*8, ALLOCATABLE, DIMENSION(:,:) :: database

! Radiance fitting control & variables (read from input control file)
  CHARACTER(256) :: radiance_line
  CHARACTER(30), ALLOCATABLE, DIMENSION(:) :: par_names
  CHARACTER(3), ALLOCATABLE, DIMENSION(:) :: par_str
  LOGICAL :: iterate_rad, weight_rad, renorm, update_pars
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: if_varied
  INTEGER*4 :: ll_rad, lu_rad, nreport, cldmax, npars, nvaried, &
       nralb, nrhwe, nrshi, nrsqe, nrsha, nrasy
  REAL*8 :: div_rad, phase, szamax, szamin, latmax, latmin, report_mult
  REAL*8, ALLOCATABLE, DIMENSION(:) :: initial, var, diff, var_factor

! Variables for fit. To be allocated depending on the number of
  ! spectral pixels and parameters before each call to specfit
  REAL*8, ALLOCATABLE, DIMENSION(:,:) :: correl, covar
  INTEGER*4 :: iteration
  REAL*8 :: rms, chisq

CONTAINS

  SUBROUTINE GASFIT()
    ! Fits an orbit of satellite-measured Earth radiance spectra, or a sequence of
    ! ground-based atmospheric spectra.
    !
    ! This version of the gas fitting code includes the option to do BOAS
    ! fitting, that is to fit radiances directly, or to fit a high-pass filtered
    ! logarithm of (radiance / irradiance). BOAS normally gives an improvement of a
    ! factor of 2-3 in the fitting statistics, with average parameter values very
    ! close to the same.
    
    ! Please communicate any corrections, additions and improvements:
    
    ! Kelly Chance
    ! Atomic and Molecular Physics Division
    ! Harvard-Smithsonian Center for Astrophysics
    ! 60 Garden Street
    ! Cambridge, MA 02138, USA
    ! phone +1-617-495-7389
    ! e-mail kchance@cfa.harvard.edu
    ! http://cfa-www.harvard.edu/~kchance/
    
    ! December 30, 2010

    IMPLICIT NONE
    LOGICAL :: iprovar
    REAL*8, ALLOCATABLE, DIMENSION(:) :: fit, temp

    INTEGER*4, ALLOCATABLE, DIMENSION(:) :: list_rad
    REAL*8, ALLOCATABLE, DIMENSION(:,:) :: underspec
    REAL*8, ALLOCATABLE, DIMENSION(:) :: pos_rad, spec_rad, sig_rad, &
         residual

    REAL*8 :: asum, avg, &
         davg,  dgas, dhw1e, dshap, dasym, drelavg, dshift, dshiftavg, &
         dshiftrad, gas, hw1e, shap, asym, remult, rmsavg, &
         shift, sigsum, squeeze, ssum, ilat, ilon, isza, isaa, ivza, ivaa, &
         lat, lon, sza, saa, vza, vaa
    
    INTEGER*4 :: i, ipix, icld, iyear, imonth, iday, ihour, imin, isec, &
         j, nfirst, ngas, nrads, npix, npixfit, cld, year, month, day, &
         hour, minu, sec
    
    write (*,'(5x, a)') 'enter fitting input file.'
    read (*, '(a)') fitin
    write (*, '(a)') 'Reading file... '//fitin
    open (unit = 21, file = fitin, status='old')
    
    ! Open fitting output file.
    if (wrt_scr) write (*,'(5x, a)') 'enter fitting output file.'
    read (*, '(a)') fitout
    write (*, '(a)') 'Opening file... '//fitout
    open (unit = 22, file = fitout, status='unknown')
    
    ! Read general parameters:
    !
    ! inputline and general_line are comments that may be echoed in the output.
    !
    ! npoints is the number of spectral points in the irradiance and each radiance.
    !
    ! delchi and provar are fitting convergence criteria: a change in chi-squared
    ! and a relative change in each fitted parameter. Automult factor to calculate
    ! fitting variable variations in optimization.
    !
    ! If (wrt_scr) enter a more verbose manner of running.
    !
    ! If (write_fit) write conditions and diagnostics to an output file. this is
    ! not the same as writing out the fitting results, which is done automatically.
    !
    ! If (write_spec) write measured spectrum out???
    !
    ! If (mirror) mirror the input file onto the output file, with updated
    ! parameters. Useful for tuning fitting parameters, particularly when
    ! optimizing by using one spectrum.
    !
    ! If (autodiff) ignore the input parameter differences used for finite
    ! differencing to calculate jacobian matrix elements and use abs (automult *
    ! parameter) instead. Experience shows that this is most often a reasonable
    ! approach to fitting, provided no initial parameter guesses are zero or very
    ! close to it.
    !
    ! If (if_residuals) write out the averaged fitting residual at the end of the
    ! spectral fitting results. Used to generate the "common mode" correction for
    ! systematic fitting residuals initially implemented at SAO.
    !
    read (21, '(a)') inputline
    read (21, '(a)') general_line
    read (21, *) npoints, delchi, provar, automult
    read (21, *) wrt_scr, write_fit, write_spec, mirror, autodiff, &
         if_residuals
    
    ! Read solar fitting parameters:
    !
    ! solar_line is a comment that may be echoed in the output.
    !
    ! if (iterate_sun) perform wavelength calibration and slit width calibration on
    ! the satellite irradiance spectrum for this orbit. Then, perform a wavelength
    ! calibration for one radiance spectrum, selected by nfirstfit (read in below).
    !
    ! f (weight_sun) make the fit a weighted fit, with highly-weighted wavelength
    ! selected by ll_sun and lu_sun. This is used to select the active fitting
    ! region, weighted as 1, whereas the inactive region is de-weighted by (1.d6)**2
    !
    ! n_solar_pars is the number of solar fitting parameters to read, normally 12.
    ! this separates the solar parameters from the radiance fitting controls and
    ! parameters which come afterward. Parameters are grouped in fours, a legacy
    ! from when my fitting code was developed to analyze laboratory spectrum. each
    ! line in a spectrum required (1) a gaussian (doppler) width; (2) a lorentzian
    ! (pressure broadening and/or lifetime) width; (3) an intensity; and (4) a
    ! position in frequency (usually MHz, GHz, THz, or cm-1) space. I find that
    ! maintaining the grouping by fours is convenient for atmospheric spectral
    ! analysis, even though it means some potential parameters are unused.
    !
    ! ll_sun and lu_sun select the active fitting region. They are integers which
    ! select the fitting location in detector pixel space from the 1-npoints
    ! wavelength in the spectrum.
    !
    ! div_sun divides the irradiance spectrum by a (usually large) number to avoid
    ! potential arithmetical overflow problems.
    !
    ! var_sun are the solar fitting variables, normally describing a cubic baseline
    ! polynomial, a cubic overall scaling polynomial, an "albedo" (intensity
    ! multiplier), slit width (hw1/e), spectral wavelength shift and spectral
    ! (accordion) squeeze. Normally, parameter 5 (the zeroth order scaling term)
    ! and parameter 9, the intensity multiplier are not varied simultaneously since
    ! they are very highly correlated.
    !
    ! If (if_var_sun (i)) vary parameter (i) in the fitting.
    !
    ! diffsun (i) is the difference in parameter (i) used to calculate the jacobian
    ! matrix element for parameters (i). If autodiff is true, diffsun (i) is
    ! ignored and 0.001 * abs (var_sun (i)) is used instead.
    !
    read (21, '(a)') solar_line
    read (21, *) iterate_sun, weight_sun, n_solar_pars, ll_sun, lu_sun, div_sun

    ! Now that I know the number of solar fitting variables allocate
    ALLOCATE(var_sun(1:n_solar_pars), diffsun(1:n_solar_pars), if_var_sun(1:n_solar_pars), &
         init_sun(1:n_solar_pars), list_sun(1:n_solar_pars), sun_par_names(1:n_solar_pars), &
         sun_par_str(1:n_solar_pars), var_sun_factor(1:n_solar_pars))
    
    ! Read variables
    do i = 1, n_solar_pars
       read (21, *) j, var_sun(i), if_var_sun(i), diffsun(i), &
            sun_par_str(i), var_sun_factor(i), sun_par_names(i)
    end do
    
    ! Read radiance fitting parameters:
    !
    ! radiance_line is a comment that may be echoed in the output.
    !
    ! If (iterate_rad) fit the radiances.
    !
    ! If (renorm) renormalize radiances to their weighted average before fitting.
    ! Numerical overflow insurance. Good for BOAS.
    !
    ! If (weight_rad) ) make the fit a weighted fit, as with weight_sun.
    !
    ! If (update_pars) use the parameters from fitting radiance n as the initial
    ! guesses for fitting radiance n+1.
    !
    ! ll_rad and lu_rad select the active fitting region.
    !
    ! div_rad divides the radiance spectra by a (usually large) number to avoid
    ! potential arithmetical overflow problems.
    !
    ! phase is a (non-fitted) parameter controlling the calculation of the
    ! correction for spectral undersampling. See subroutine undersample and
    ! Applied Optics 44, 1296-1304, 2005 for more detail.
    !
    ! nreport is the parameter number outputted to the fitting output file.
    ! report_mult is a multiplier for its value and fitting uncertainty. report_mult
    ! is required because the reference spectra are scaled (below) to avoid
    ! potential numerical underflow issues. The fitted parameter values must be
    ! scaled by the  same amount
    !
    ! Ground pixels (scenes) with solar zenith angle (sza) <= szamax, (sza) => szamin,
    ! latitudes <= latmax, and latitudes >= latmin are analyzed.
    !
    ! Maximum cloud index clasification (cldmax). Only process pixels with
    ! (cloud index) <= cldmax.
    !
    ! nfirstfit is the number of the radiance (counting consecutively along the
    ! orbit, as given by ipix) to fit for wavelength in order to establish a
    ! standard wavelength grid and sample the reference spectra.
    !
    read (21, '(a)') radiance_line
    read (21, *) iterate_rad, renorm, weight_rad, update_pars, ll_rad, lu_rad, &
         div_rad, phase, nreport, report_mult, npars
    read (21, *) szamax, szamin, latmax, latmin, cldmax, nfirstfit
    
    ! Allocate radiance fitting parameter variables
    ALLOCATE(var(1:npars), if_varied(1:npars),diff(1:npars), initial(1:npars), &
         par_names(1:npars), var_factor(1:npars), par_str(1:npars),list_rad(1:npars))
    
    ! Read variables
    DO i = 1, npars
       read (21, *) j, var(i), if_varied(i), diff(i), par_str(i), var_factor(i), par_names(i)
    END DO
    
    ! Automatic difference-taking.
    if (autodiff) then
       do i = 1, npars
          if (if_varied (i)) diff (i) = dabs (automult * var (i))
       end do
    end if

    do i = 1, npars
       IF (par_str(i) .EQ. alb_str) nralb = i 
       IF (par_str(i) .EQ. hwe_str) nrhwe = i 
       IF (par_str(i) .EQ. shi_str) nrshi = i 
       IF (par_str(i) .EQ. sqe_str) nrsqe = i 
       IF (par_str(i) .EQ. sha_str) nrsha = i 
       IF (par_str(i) .EQ. asy_str) nrasy = i 
    end do
    
    ! Order the coefficients for mrqmin and save the initial values of the
    ! coefficients.
    nvaried = 0
    do i = 1, npars
       initial (i) = var (i)
       if (if_varied (i)) then
          nvaried = nvaried + 1
          list_rad (nvaried) = i
          if (i .eq. nreport) ngas = nvaried
       end if
    end do
    if (wrt_scr) write (*, *) 'nvaried  =', nvaried
    if (wrt_scr) write (*, *) 'Max. sza =', szamax
    if (wrt_scr) write (*, *) 'Min. sza =', szamin
    if (wrt_scr) write (*, *) 'Max. lat =', latmax
    if (wrt_scr) write (*, *) 'Min. lat =', latmin
    if (wrt_scr) write (*, *) 'Max. cld =', cldmax
    if (wrt_scr) write (*, *) 'nhw1e    =', nrhwe
    if (wrt_scr) write (*, *) 'nshift   =', nrshi
    if (wrt_scr) write (*, *) 'nsqueeze =', nrsqe
    if (wrt_scr) write (*, *) 'nshape   =', nrsha
    if (wrt_scr) write (*, *) 'nasym    =', nrasy
    
    ! Open irradiance and radiance level 1 file.
    if (wrt_scr) write (*,'(5x, a)') 'enter spectrum input file.'
    read (*, '(a)') specin
    write (*, '(a)') 'Opening file... '//specin
    open (unit = 23, file = specin, status='old')
    
    ! Open spectrum output file.
    if (write_spec) then
       if (wrt_scr) write (*,'(5x, a)') 'enter spectrum output file.'
       read (*, '(a)') specout
       if (wrt_scr) write (*,'(a)') 'Opening file... '//TRIM(specout)
       open (unit = 24, file = specout, status='unknown')
    end if
    
    ! Read the solar spectrum. A 2-column (position and irradiance value) spectrum
    ! of npoints points, with no header, is expected.
    ALLOCATE(pos_sun(1:npoints), spec_sun(1:npoints), sig_sun(1:npoints))
    do i = 1, npoints
       read (23, *) pos_sun (i), spec_sun (i)
       ! divide to avoid overflow issues.
       spec_sun (i) = spec_sun (i) / div_sun
       ! select fitting window by the use of weighting.
       sig_sun(i) = 1.d30
       if (i .ge. ll_sun .and. i .le. lu_sun) sig_sun(i) = 1.d0
    end do
    
    ! Write out the input line and other input information.
    IF (write_fit) WRITE (22, '(a)') inputline
    IF (write_fit) CALL write_input()
    
    ! Proceed with solar fit for slit function and wavelength calibration
    ! Perform solar wavelength calibration and slit width fitting:
    ! Calculate avg here, for use in calculated spectra.
    if (.not. weight_sun) then
       avg = (pos_sun (npoints) + pos_sun (1)) / 2.
    else
       asum = 0.
       ssum = 0.
       do i = 1, npoints
          asum = asum + pos_sun (i) / (sig_sun (i)**2)
          ssum = ssum + 1. / (sig_sun (i)**2)
       end do
       avg = asum / ssum
    end if
    if (wrt_scr) write (*, *) ' avg = ', avg
    
    ! Calculate and iterate on the irradiance spectrum. The flag for convergence
    ! from lack of change in the variables, iprovar, is set to .false. before
    ! iteration begins.
    iprovar = .false.
    if (iterate_sun) then
    
       ! Use strings to assign specific fitting parameter indices
       do i = 1, n_solar_pars
          IF (sun_par_str(i) .EQ. alb_str) nalb = i 
          IF (sun_par_str(i) .EQ. hwe_str) nhwe = i 
          IF (sun_par_str(i) .EQ. shi_str) nshi = i 
          IF (sun_par_str(i) .EQ. sqe_str) nsqe = i 
          IF (sun_par_str(i) .EQ. sha_str) nsha = i 
          IF (sun_par_str(i) .EQ. asy_str) nasy = i 
       end do
       
       ! Automatic difference-taking
       if (autodiff) then
          do i = 1, n_solar_pars
             if (if_var_sun(i)) diffsun(i) = ABS(automult * var_sun(i))
          end do
       end if
       
       ! Order the parameters for mrqmin (housekeeping for the fitting process) and
       ! save the initial values of the parameters, in case they are required at some
       ! later stage. At present, just used in writing the fitting output file.
       nvar_sun = 0
       do i = 1, n_solar_pars
          init_sun (i) = var_sun (i)
          if (if_var_sun (i)) then
             nvar_sun = nvar_sun + 1
             list_sun (nvar_sun) = i
          end if
       end do
       if (wrt_scr) write (*, *) 'nvar_sun =', nvar_sun
       
       ! Allocate variables neeed for fitting (specfit, mrqcof & mrqmin)
       ALLOCATE(correl(1:n_solar_pars,1:n_solar_pars),covar(1:n_solar_pars,1:n_solar_pars), &
            fit(1:npoints),database(npars+4,1:npoints))
       correl = 0.0d0; covar=0.0d0; fit=0.0d0; database = 0.0d0

       ! Perform fit
       call specfit (npoints, n_solar_pars, nvar_sun, list_sun(1:n_solar_pars), &
            avg, spec_sun(1:npoints), pos_sun(1:npoints), sig_sun(1:npoints), &
            fit(1:npoints), var_sun(1:n_solar_pars), diffsun(1:n_solar_pars), &
            var_sun_factor(1:n_solar_pars), sun_par_str(1:n_solar_pars), &
            iprovar, 1)

       ! Shift and squeeze solar spectrum.
       shift = var_sun(nshi); squeeze = var_sun(nsqe)
       do i = 1, npoints
          pos_sun(i) = (pos_sun(i) - var_sun(nshi)) / (1.d0 + var_sun(nsqe))
       end do
       dshift = rms * sqrt (covar (nshi, nshi) * float (npoints) / &
            float (npoints - nvar_sun))
       write (*, '(a, 1pe11.3)') 'solar wavelength calibration: rms = ', rms
       write (*, '(a, 1p2e14.6)') 'irrad: shift, 1 sigma = ', - var_sun (nshi), dshift
       write (*, *) 'nvar_sun = ', nvar_sun
              
       ! Set slit function parameters for later use in undersampling correction. 
       ! HW1E means the gaussian half-width at 1/e of the maximum intensity.
       hw1e = var_sun(nhwe); shap = var_sun(nsha); asym = var_sun(nasy)
       dhw1e = rms * sqrt (covar (nhwe, nhwe) * float (npoints) / &
            float (npoints - nvar_sun))
       dshap = rms * sqrt (covar (nsha, nsha) * float (npoints) / &
            float (npoints - nvar_sun))
       dasym = rms * sqrt (covar (nasy, nasy) * float (npoints) / &
            float (npoints - nvar_sun))
       write (*, '(a, 1p2e14.6)') 'irrad: hw1e, 1 sigma = ', hw1e, dhw1e
       write (*, '(a, 1p2e14.6)') 'irrad: shap, 1 sigma = ', shap, dshap
       write (*, '(a, 1p2e14.6)') 'irrad: asym, 1 sigma = ', asym, dasym

       ! Deallocate variables needed for fitting
       DEALLOCATE(correl, covar, fit)  

    end if ! on iterate_sun
    write (*, *) 'finished with iterate_sun'

    ! Finish building database of reference spectra
    ! Calculate the undersampled spectrum.
    ALLOCATE(underspec(1:2,1:npoints))
    call undersample (pos_sun(1:npoints), npoints, underspec, phase)

    ! Save solar spectrum and undersampling spectrum in to database
    ! Sun
    database(1,1:npoints) = spec_sun(1:npoints)
    ! Under sampling
    DO i = 1, npars
       IF (par_names(i) .EQ. us1_str) THEN
          database(i,1:npoints)  = underspec(1,1:npoints)
       ELSE IF (par_names(i) .EQ. us2_str) THEN
          database(i,1:npoints)  = underspec(2,1:npoints)
       END IF
    END DO
    DEALLOCATE(underspec)

    ! Calculate the splined fitting database. Note that the undersampled
    ! spectrum has just been done. Finish filling in reference database array.
    call dataspline (pos_sun(1:npoints), npoints, hw1e, shap, asym)

    ! Deallocate variables
    DEALLOCATE(var_sun, var_sun_factor, sun_par_str, sun_par_names, diffsun, if_var_sun, &
         init_sun, list_sun, pos_sun, spec_sun, kppos, kpspec, kppos_ss, &
         kpspec_gauss)
    
    ! ----------------------------------------------------------------- !
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
    ! Fit the radiances (loop over to only keep one spectrum in memory) !
    ! ----------------------------------------------------------------- !
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
    ! Initialize number of radiance pixels & number of radiance pixel fitted
    npix = 0; npixfit = 0

    ! Initialize several diagnostics.
    ! rmsavg: average fitting rms.
    ! davg: average fitting uncertainty for reported parameter.
    ! drelavg: average relative fitting uncertainty for reported parameter.
    ! dshiftavg:  average fitting uncertainty for spectral shift parameter.
    rmsavg = 0.d0
    davg = 0.d0
    drelavg = 0.d0
    dshiftavg = 0.d0

    ! Allocate general fitting variables
    ALLOCATE(pos_rad(1:npoints),spec_rad(1:npoints),sig_rad(1:npoints))

    ! Add the weighting, in the single array.
    do i = 1, npoints
       sig_rad(i) = 1.d30
       if (i .ge. ll_rad .and. i .le. lu_rad) sig_rad(i) = 1.d0
    end do

    ! Radiance loop
    radfit: DO 
       ! Read the measured spectra.
       !
       ! ipix: The number of the radiance, counting consecutively along the orbit.
       !
       ! iyear (Year data was collected), imonth (Month data was collected), iday 
       ! (Day data was collected), ihour (Hour data wav collected), imin 
       ! (Minuted data was collected) and isec (Second data was collected)
       !
       ! ilat (Latitude data was collected), ilon
       ! 1 is normally selected when we extract data from a GOME el1 file. Because of
       ! longer integration times when the light conditions are low, sometimes our
       ! selected band has not been read out, and we need to skip over these ground
       ! pixels.
       !
       ! isub: Subpixel counter: 0 = east, 1 = center, 2 = west, 3 = flyback.
       !
       ! sza1, saa1, sza2, saa2, sza3, saa3: Solar zenith and azimuth angles for
       ! east/west edges and center of ground pixel.
       !
       ! zs, re: Satellite height and earth radius for each ground pixel measurement.
       !
       ! lat1, lon1, lat2, lon2, lat3, lon3, lat4, lon4,  latc, lonc: Latitudes and
       ! longitudes for 4 corners and center of ground pixel.
       !
       READ (23, *, end = 50) ipix, icld
       READ (23, *) iyear, imonth, iday, ihour, imin, isec
       READ (23, *) ilat, ilon
       READ (23, *) isza, isaa, ivza, ivaa
       npix = npix+1
       
       ! Decide whether to process. If so, update the appropriate arrays
       IF ( (icld .LE. cldmax) .AND. (isza .LE. szamax) .AND. (isza .GE. szamin) .AND. &
            (ilat .LE. latmax) .AND. (ilat .GE. latmin) ) THEN
          if (wrt_scr) write(*,'(a,I5)') 'Processing pixel # ',npix

          ! If the pixel is to be processed read the radiance and save it in the right
          ! variables
          do j = 1, npoints
             read (23, *) pos_rad(j), spec_rad(j)
             spec_rad(j) = spec_rad(j) / div_rad
          end do
          npixfit = npixfit + 1
       
          ! Renormalize, if requested
          if (renorm) then
             remult = 0.d0
             if (weight_rad) then
                sigsum = 0.d0
                do j = 1, npoints
                   remult = remult + spec_rad(j) / (sig_rad(j) * sig_rad(j))
                   sigsum = sigsum + 1.d0 / (sig_rad(j) * sig_rad(j))
                end do
                remult = remult / sigsum
             else
                do j = 1, npoints
                   remult = remult + spec_rad(j)
                end do
                remult = remult / npoints
             end if
             do j = 1, npoints
                spec_rad(j) = spec_rad(j) / remult
             end do
          end if
       
       else 

          if (wrt_scr) write(*,'(a,I5)') 'Skipping pixel # ',npix
          ! Dummy read non-processed radiance
          do j = 1, npoints
             read (23, *)
          end do
          ! Cycle this skipped pixel
          cycle
       end if     

    END DO radfit
50  continue
    DEALLOCATE(var, var_factor, par_str, par_names, diff, if_varied, initial, &
         pos_rad, spec_rad, sig_rad)

    print*, shift, squeeze
    stop

!!$       
!!$       if (.not. weight_rad) then
!!$          avg = (pos (npoints) + pos (1)) / 2.
!!$       else
!!$          asum = 0.
!!$          ssum = 0.
!!$          do j = 1, npoints
!!$             asum = asum + pos (j) / (sig (j)**2)
!!$             ssum = ssum + 1. / (sig (j)**2)
!!$          end do
!!$          avg = asum / ssum
!!$       end if
!!$       iprovar = .false.
!!$       if (iterate_rad) then
!!$          if (.not. update_pars) then
!!$             if (wrt_scr) write (*, *) ' updating parameters'
!!$             do j = 1, npars
!!$                var (j) = initial (j)
!!$             end do
!!$          end if
!!$          
!!$          !   BOAS fitting now done here.
!!$    call specfit (npoints, npars, nvaried, list_rad, iteration, avg, &
!!$      spec, pos, sig, fit, var, chisq, delchi, diff, correl, covar, rms, &
!!$      provar, iprovar, wrt_scr, 4)

          !   write out radiance fit: keep here as possible future diagnostic.
          !   do j = ll_rad, lu_rad
          !     write (24, '(f11.6, 1p3e13.5)') pos (j), spec (j), fit (j), &
          !     spec (j) - fit (j)
          !   end do
          
          !   Write out fitting results.
!!$          gas = report_mult * var (nreport)
!!$          dgas = report_mult * rms * sqrt (covar (ngas, ngas) * float (npoints) / &
!!$               float (npoints - nvaried))
!!$          write (22, '(i5, 1p3e12.4)') npix(i), gas, dgas, rms
!!$          write (22, '(4f7.2)') sza(i), saa(i), vza(i), vaa(i)
!!$          write (22, '(10f7.2)') lat(i), lon(i)
!!$          dshiftrad = rms * sqrt (covar (nshiftrad, nshiftrad) * float (npoints) / &
!!$               float (npoints - nvaried))
!!$          dshiftavg = dshiftavg + dshiftrad
!!$          rmsavg = rmsavg + rms
!!$          davg = davg + dgas
!!$          drelavg = drelavg + dabs (dgas / gas)
!!$          
!!$          ! Write out extra diagnostics.
!!$          if (write_fit) &
!!$               call write_output (nvaried, list_rad, iteration, write_fit, write_spec, &
!!$               iterate_rad, pos, spec, fit, var, initial, covar, chisq, delchi, correl, &
!!$               rms, npoints, weight_rad, provar, iprovar)
!!$          
!!$          ! else (TBD later) calculate sun and radiance spectrum, no iteration.
!!$       end if ! on iterate_rad
!!$       if (if_residuals) then
!!$          do j = 1, npoints
!!$             residual (j) = residual (j) + spec (j) - fit (j)
!!$          end do
!!$       end if
!!$    end do ! on nrads.
!!$    
!!$    ! Write out the average fitting statistics to screen.
!!$    write (*, '(a33, 1pe13.5)') '                       Avg rms = ', rmsavg / nrads
!!$    write (*, '(a33, 1pe13.5)') '                     Avg dgas = ', davg / nrads
!!$    write (*, '(a33, 1pe13.5)') ' Avg relative gas uncertainty = ', drelavg / nrads
!!$    write (*, '(a33, 1pe13.5)') ' Avg radiance shift uncertainty = ', dshiftavg / &
!!$         nrads
!!$    
!!$    ! Write out the fitting statistics and average residuals.
!!$    if (if_residuals) then
!!$       write (22, '(a28, 1pe13.5)') '                  Avg rms = ', rmsavg / nrads
!!$       write (22, '(a28, 1pe13.5)') '                Avg dgas = ', davg / nrads
!!$       write (22, '(a28, 1pe13.5)') ' Avg rel gas uncertainty = ', drelavg / nrads
!!$       do i = 1, npoints
!!$          write (22, '(f10.5, 1pe13.5)') pos (i), residual (i) / nrads
!!$       end do
!!$    end if
!!$    
!!$    ! If (mirror), mirror the input file, with updated parameters, onto the output
!!$    ! fitting file.
!!$    if (mirror) then
!!$       write (22, *)
!!$       write (22, '(a)') inputline
!!$       write (22, '(a)') general_line
!!$       write (22, '(i5, 1p2e12.4, 0pf8.4)') npoints, delchi, provar, automult
!!$       write (22, *) wrt_scr, write_fit, write_spec, mirror, autodiff, &
!!$            if_residuals
!!$       write (22, '(a)') solar_line
!!$       write (22, '(2l2, 3i4, 1pe9.2)') iterate_sun, weight_sun, n_solar_pars, &
!!$            ll_sun, lu_sun, div_sun
!!$       do i = 1, n_solar_pars
!!$          write (22, '(1pe15.7, l4, 1pe14.5)') var_sun (i), if_var_sun (i), &
!!$               diffsun (i)
!!$       end do
!!$       write (22, '(a)') radiance_line
!!$       write (22, '(4l2, 2i4, 1pe9.2, 0pf5.2)') iterate_rad, renorm, weight_rad, &
!!$            update_pars, ll_rad, lu_rad, div_rad, phase
!!$       write (22, '(4f7.2, i6)') szamax, szamin, latmax, latmin, nfirstfit
!!$       do i = 1, npars
!!$          write (22, '(1pe15.7, l4, 1pe14.5)') var (i), if_varied (i), diff (i)
!!$       end do
!!$    end if
    
    close (unit = 21)
    close (unit = 22)
    close (unit = 23)
    close (unit = 24)
    
  END SUBROUTINE GASFIT

  SUBROUTINE specfit (np, npars, nvaried, lista, avg, &
       spec, pos, sig, fit_spec, var, diff, fac, str, iprovar, ntype)
    ! Driver for Levenberg-Marquardt nonlinear least squares minimization.
    
    IMPLICIT NONE
    
    ! Input variables
    INTEGER*4, INTENT(IN) :: np, npars, nvaried, ntype
    REAL*8, INTENT(IN) :: avg
    REAL*8, INTENT(IN), DIMENSION(1:np) :: spec, pos, sig
    REAL*8, INTENT(IN), DIMENSION(1:npars) :: diff, fac
    CHARACTER(3), INTENT(IN), DIMENSION(1:npars) :: str
    
    ! Output variables
    REAL*8, INTENT(OUT), DIMENSION(1:np) :: fit_spec
    
    ! Modified variables
    REAL*8, INTENT(INOUT), DIMENSION(1:npars) :: var
    LOGICAL, INTENT(INOUT) :: iprovar
    INTEGER*4, INTENT(INOUT), DIMENSION(1:npars) :: lista

    ! Local variables
    REAL*8 :: alamda, ochisq, difftest, prop, rsum
    REAL*8, DIMENSION(1:npars) :: var0
    INTEGER*4 :: itest, i, j    
    external funcs
    
    alamda = - 1
    ! First call to minimization routine
    call mrqmin (pos(1:np), spec(1:np), sig(1:np), np, var(1:npars), &
         fac(1:npars), str(1:npars), npars, lista(1:npars), nvaried, &
         funcs, alamda, avg, diff(1:npars), ntype)
    iteration = 1
    itest = 0

    iter_loop: DO
       if (wrt_scr) then
          write (22,'(/1x, a, i2, t18, a, 1pe10.4, t44, a, e9.2)') 'iteration #', &
               iteration, 'chi-squared: ', chisq, 'alamda:', alamda
          write (22,'(1x, a)') 'variables:'
          do i = 1, npars
             write (22,'(A3,3x, E10.2)') str(i), var(i)
          end do
          write (*,'(/1x, a, i2, t18, a, 1pe10.4, t44, a, e9.2)') 'iteration #', &
               iteration, 'chi-squared: ', chisq, 'alamda:', alamda
       end if
      
       iteration = iteration + 1
       ochisq = chisq
       do i = 1, nvaried
          var0 (lista (i)) = var (lista (i))
       end do
       ! Iteration over the minimization
       call mrqmin (pos(1:np), spec(1:np), sig(1:np), np, var(1:npars), &
            fac(1:npars), str(1:npars), npars, lista(1:npars), nvaried, &
            funcs, alamda, avg, diff(1:npars), ntype)

       do i = 1, nvaried
          prop = abs (var0 (lista (i)) - var (lista (i)))
          difftest = provar * diff (lista (i))
          ! Check iteration conditions to finish loop
          if (prop .gt. difftest) then
             if (chisq .gt. ochisq) itest = 0
             if (abs(ochisq - chisq) .lt. delchi) itest = itest + 1
             if (itest .lt. 1) cycle iter_loop !go to 2
             iprovar = .true.
             exit iter_loop
          end if
       end do
    END DO iter_loop

    if (wrt_scr) write (*, '(1x, a)') 'iprovar = .true.'
    alamda = 0.0    
    ! Final call to minimization subroutine
    call mrqmin (pos(1:np), spec(1:np), sig(1:np), np, var(1:npars), &
         fac(1:npars), str(1:npars), npars, lista(1:npars), nvaried, &
         funcs, alamda, avg, diff(1:npars), ntype)
    if (wrt_scr) then
       write (22,'(/1x, a, i2, t18, a, e10.4, t43, a, e9.2)') 'iteration #', &
            iteration, 'chi-squared: ', chisq, 'alamda:', alamda
       write (*,'(/1x, a, i2, t18, a, e10.4, t43, a, e9.2)') 'iteration #', &
            iteration, 'chi-squared: ', chisq, 'alamda:', alamda
       write (*,'(1x, a)') 'variables:'
       do i = 1, npars
          write (*,'(A3,3x, E10.2)') str(i), var(i)
       end do
    end if
    
    ! Calculate the correlation matrix.
    do i = 1, nvaried
       correl (i, i) = 1.
       do j = 1, i
          if (i .ne. j) correl (i, j) = covar (i, j) / sqrt (covar (i, i) * &
               covar (j, j))
       end do
    end do
    
    ! Calculate the final spectrum.
    call spectrum (np, npars, avg, pos(1:np), fit_spec(1:np), &
    var(1:npars), fac(1:npars), str(1:npars), ntype)
    
    ! Calculate the rms of the fit.
    rsum = 0.
    do i = 1, np
       rsum = rsum + ((fit_spec (i) - spec (i)) / sig (i))**2
    end do
    rms = sqrt (rsum / np)

    return
  END SUBROUTINE specfit

  subroutine spectrum (npoints, npars, avg, pos, fit, var, fac, str, ntype)
    
    ! Spectrum calculation for both fitting and non-fitting cases. For fitting, it
    ! is called with 4 different qualifiers ("ntypes") for purposes of determining
    ! the calculation type and initializing the fitting reference database.
    ! ntype = 1: Solar wavelength calibration. Open kpnospec, the high spectral
    ! resolution irradiance reference spectrum.
    ! ntype /= 1: Later radiance fits (refspec already filled).
    IMPLICIT NONE

    ! Input variables
    INTEGER*4, INTENT(IN) :: npoints, npars, ntype
    REAL*8, INTENT(IN) :: avg
    REAL*8, INTENT(IN), DIMENSION(1:npars) :: var, fac
    REAL*8, INTENT(IN), DIMENSION(1:npoints) :: pos
    CHARACTER(3), INTENT(IN), DIMENSION(1:npars) :: str

    ! Output variables
    REAL*8, INTENT(OUT), DIMENSION(1:npoints) :: fit

    ! Local variables
    INTEGER*4, PARAMETER :: maxpts = 10000, maxkpno = 20000
    INTEGER*4 :: i, jmax, j
    REAL*8 :: hw, sh, as, sun_avg
    REAL*8, DIMENSION(1:npoints) :: del, scaling, baseline, &
         sunpos, sunpos_ss, sunspec, sunspec_ss
    REAL*8, ALLOCATABLE, DIMENSION(:) :: d2sun
    CHARACTER*120 :: kpnospec
    LOGICAL :: first_sun=.true.

    if (ntype .eq. 1) then
       if (first_sun) then
          first_sun = .false.
          if (wrt_scr) write (*, '(5x,a)') 'Reading high resolution solar spectrum...'          
          read (*, '(a)') kpnospec
          if (wrt_scr) write (*, '(5x,a)') '... from file '//TRIM(kpnospec)
          ! Get total number of points in the spectrum file (so we can allocated space for them)
          open (unit = 25, file = kpnospec, status='old')
          nkppos=0
          do
             read(25,*,END=11)
             nkppos=nkppos+1
          end do
11        rewind(unit=25)
          ! Allocate solar variables
          ALLOCATE(kppos(1:nkppos),kpspec(1:nkppos),kppos_ss(1:nkppos),kpspec_gauss(1:nkppos))
          ! Read solar spectrum
          do i = 1, nkppos     
             read (25, *) kppos(i), kpspec(i)
          end do
          ! Normalize solar spectrum
          sun_avg = SUM(kpspec(1:nkppos)) / REAL(nkppos,KIND=8)
          kpspec(1:nkppos) = kpspec(1:nkppos) / sun_avg
          close (unit = 25)
       end if ! on first_sun
    end if

    ! Calculate the spectrum: First do the shift and squeeze. Shift var (nshi),
    ! squeeze by 1 + var (nsqe); do in absolute sense, to make it easy to back-convert
    ! radiance data.
    do i = 1, nkppos
       if (ntype .eq. 1) then
          kppos_ss(i) = kppos(i) * (1.d0 + var(nsqe)) + var(nshi)
       else
          sunpos_ss(i) = sunpos(i) * (1.d0 + var(nsqe)) + var(nshi)
       end if
    end do
    
    ! Broadening and re-sampling of solar spectrum.
    if (ntype .eq. 1) then
       ! Case for wavelength fitting of irradiance and radiance.
       ! Broaden the solar reference by the hw1e value with shape
       ! factor sh and asymmetric factor as.
       hw = var(nhwe); sh = var(nsha); as = var(nasy)
       call super_gauss (nkppos, kppos_ss, kpspec, hw, sh, as, kpspec_gauss)
       ! Re-sample the solar reference spectrum to the radiance grid.
       ALLOCATE(d2sun(1:nkppos))
       call spline (kppos_ss, kpspec_gauss, nkppos, d2sun)
       call splint (kppos_ss, kpspec_gauss, d2sun, nkppos, npoints, pos, sunspec_ss)
       DEALLOCATE(d2sun)
    else
       ! Re-sample the solar reference spectrum to the radiance grid.
       ALLOCATE(d2sun(1:npoints))
       call spline (sunpos_ss, sunspec, npoints, d2sun)
       call splint (sunpos_ss, sunspec, d2sun, npoints, npoints, pos, sunspec_ss)
       DEALLOCATE(d2sun)
    end if

    ! Add contributions from each  variable
    ! First add albedo contribution
    fit(1:npoints) = var(nalb) * sunspec_ss(1:npoints)
    ! Now loop over the rest of parameters to add BOAS, and polynomial contributions
    ! Initial add-on contributions
    do i = 1, npars
       IF (str(i) .EQ. ad1_str) fit(1:npoints) = fit(1:npoints) + var(i) * database(i,1:npoints)
    end do
    ! Beer's law contributions
    do i = 1, npars
       IF (str(i) .EQ. ble_str) fit(1:npoints) = fit(1:npoints) * dexp(-var(i) * database(i,1:npoints))
    end do
    ! Final add-on contributions
    do i = 1, npars
       IF (str(i) .EQ. ad2_str) fit(1:npoints) = fit(1:npoints) + var(i) * database(i,1:npoints)
    end do
    ! Polynomials
    del(1:npoints) = pos(1:npoints) - avg
    ! Scaling polynomial
    scaling(1:npoints) = 0.0d0
    do i = 1, npars
       IF (str(i) .EQ. scp_str) scaling(1:npoints) = scaling(1:npoints) + var(i) * del(1:npoints)**fac(i)
    end do
    fit(1:npoints) = fit(1:npoints) * scaling(1:npoints)
    ! Baseline polynomial
    baseline(1:npoints) = 0.0d0
    do i = 1, npars
       IF (str(i) .EQ. bsp_str) baseline(1:npoints) = baseline(1:npoints) + var(i) * del(1:npoints)**fac(i)
    end do
    fit(1:npoints) = fit(1:npoints) + baseline(1:npoints)
    
    return
  end subroutine spectrum
  
  subroutine mrqmin (x, y, sig, ndata, a, fac, str, ma, lista, mfit, &
       funcs, alamda, avg, diff, ntype)
    
    ! Levenberg-Marquardt minimization adapted from Numerical Recipes.
    IMPLICIT NONE
    ! Input variables
    REAL*8, INTENT(IN), DIMENSION(1:ndata) :: x, y, sig
    REAL*8, INTENT(IN), DIMENSION(1:ma) :: diff, fac
    REAL*8, INTENT(IN) :: avg
    INTEGER*4, INTENT(IN) :: ndata, ma, mfit, ntype
    CHARACTER(3), INTENT(IN), DIMENSION(1:ma) :: str
    
    ! Modified variables
    REAL*8, INTENT(INOUT), DIMENSION(1:ma) :: a
    REAL*8, INTENT(INOUT) :: alamda
    INTEGER*4, INTENT(INOUT), DIMENSION(1:ma) :: lista
    
    ! Local variables
    REAL*8 :: ochisq, even
    REAL*8, ALLOCATABLE, DIMENSION(:) :: atry, da
    REAL*8, DIMENSION(1:ma,1:ma) :: identity, inverse
    REAL*8, ALLOCATABLE, DIMENSION(:) :: beta
    REAL*8, ALLOCATABLE, DIMENSION(:,:) :: alpha, covar0
    LOGICAL :: first_call = .TRUE.
    INTEGER*4 :: j,k,kk,ihit
    INTEGER*4, DIMENSION(1:mfit) :: index
    
    external funcs
    save ochisq,atry,da,beta,first_call,alpha

    if (first_call) then
       first_call = .false.
       ALLOCATE(atry(1:ma),da(1:ma),beta(1:ma),alpha(1:ma,1:ma))
    end if

    if (alamda .lt. 0.) then
       ! Just for security we check that lista has right size
       kk = mfit + 1
       do j = 1, ma
          ihit = 0
          do k = 1, mfit
             if (lista (k) .eq. j) ihit = ihit + 1
          end do
          if (ihit .eq. 0) then
             lista (kk) = j
             kk = kk + 1
          else if (ihit .gt. 1) then
             print*, 'improper permutation in lista'
          end if
       end do
       if (kk .ne. (ma + 1)) print*, 'improper permutation in lista'
       alamda = 0.001
       ! Initial fit
       call mrqcof (x, y, sig, ndata, a, fac, str, ma, lista, mfit, alpha, &
            beta, avg, diff, ntype)
       ! Save chisq for later use
       ochisq = chisq
       ! Initialize new fitting parameters
       do j = 1, ma
          atry (j) = a (j)
       end do

    end if

    ! Save covar and beta
    do j = 1, mfit
       do k = 1, mfit
          covar(j,k) = alpha(j,k)
       end do
       covar(j,j) = alpha(j,j) * (1.0+alamda)
       da(j) = beta(j)
    end do

    ! Solve matrix
    call ludcmp (covar, mfit, ma, index, even)
    call lubksb (covar, mfit, ma, index, da)    

    if (alamda .eq. 0.) then
       identity = 0.d0
       do j = 1, ma
          identity (j, j) = 1.d0
       end do
       ! Calculate inverse of curvature matrix to obtain covariance matrix.
       inverse = identity
       do j = 1, mfit
          !   call lubksb (covar, mfit, ma, index, inverse (1, j))
          call lubksb (covar, mfit, ma, index, inverse (:, j))
          !   Luke Valin correction 31mar2014; answers are exactly the same
       end do
       ALLOCATE(covar0(1:ma,1:ma))
       covar0 = inverse
       covar  = identity
       do j = 1, mfit
          do k = 1, mfit
             ! Save covar
             covar (lista(j), lista(k)) = covar0(j,k)
          end do
       end do
       DEALLOCATE(atry,da,beta,alpha,covar0)
       return       
    end if

    ! Work out new parameters based in da
    do j = 1, mfit
       atry (lista (j)) = a (lista (j)) + da (j)
    end do

    ! Try out new parameters (atry)
    call mrqcof (x, y, sig, ndata, atry, fac, str, ma, lista, mfit, covar, da, &
         avg, diff, ntype)

    ! If the new parameters are more succesful
    if (chisq .lt. ochisq) then
       alamda = 0.1 * alamda
       ! Save chisq
       ochisq = chisq
       do j = 1, mfit
          do k = 1, mfit
             ! Save alpha
             alpha (j, k) = covar (j, k)
          end do
          ! Save beta
          beta (j) = da (j)
          a (lista (j)) = atry (lista (j))
       end do
    else
       alamda = 10. * alamda
       chisq = ochisq
    end if
    
    return
  end subroutine mrqmin

  subroutine mrqcof (x, y, sig, ndata, a, fac, str, ma, lista, mfit, alpha, beta, &
       avg, diff, ntype)
    
    ! Calculates chi-squared gradient and curvature matrices for the Levenberg-
    ! Marquardt minimization. From Numerical Recipes.
    IMPLICIT NONE

    ! Input variables
    INTEGER*4, INTENT(IN) :: ndata, ma, mfit, ntype
    INTEGER*4, INTENT(IN), DIMENSION(1:ma) :: lista
    REAL*8, INTENT(IN) :: avg
    REAL*8, INTENT(IN), DIMENSION(1:ndata) :: x, y, sig
    REAL*8, INTENT(IN), DIMENSION(1:ma) :: a, diff, fac
    CHARACTER(3), INTENT(IN), DIMENSION(1:ma) :: str

    ! Modified variables
    REAL*8, INTENT(OUT), DIMENSION(1:ma) :: beta
    REAL*8, INTENT(OUT), DIMENSION(1:ma,1:ma) :: alpha

    ! Local variabes
    REAL*8, DIMENSION(1:ma,1:ndata) :: dyda
    REAL*8, DIMENSION(1:ndata) :: ymod
    REAL*8 :: sig2i, dy, wt
    INTEGER*4 :: j, k, i
    external funcs

    do j = 1, mfit
       do k = 1, j
          alpha (j, k) = 0.
       end do
       beta (j) = 0.
    end do
    chisq = 0.

    call funcs (x, ndata, a, fac, str, ymod, dyda, ma, lista, mfit, avg, diff, &
         ntype)
    do i = 1, ndata
       sig2i = 1./ (sig (i) * sig (i))
       dy = y (i) - ymod (i)
       do j = 1, mfit
          wt = dyda (lista (j), i) * sig2i
          do k = 1, j
             alpha (j, k) = alpha (j, k) + wt * dyda (lista (k), i)
          end do
          beta (j) = beta (j) + dy * wt
       end do
       chisq = chisq + dy * dy * sig2i
    end do
    do j = 2, mfit
       do k = 1, j - 1
          alpha (k, j) = alpha (j, k)
       end do
    end do
    return
  end subroutine mrqcof

  SUBROUTINE write_input()

    IMPLICIT NONE
    INTEGER*4 :: i
    
    ! Summarizes the inputs for later reference.
    write (22, '(1x, t2, a)') 'BOAS fitting'
    write (22, '(1x, t2, i3, t6, a, t29, i3, t33, a)') npars, &
         'total parameters used,', nvaried, 'varied.'
    write (22, '(1x, t2, i3, t6, a, t29, i3, t33, a)') n_solar_pars, &
         'solar parameters used,', nvar_sun, 'varied.'
    write (22, '(1x, t2, a, t13, i4, t18, a)') 'spectra of', npoints, &
         'points were calculated.'
    write (22, *)
    if (iterate_rad) then
       write (22, *) 'radiance iteration is enabled.'
    else
       write (22, *) 'radiance iteration is disabled.'
    end if
    if (iterate_sun) then
       write (22, *) 'solar iteration is enabled.'
    else
       write (22, *) 'solar iteration is disabled.'
    end if
    if (weight_rad) then
       write (22, *) 'radiance weighting is enabled.'
       write (22, *) 'll_rad = ', ll_rad, 'lu_rad = ', lu_rad
    else
       write (22, *) 'radiance weighting is disabled.'
    end if
    if (weight_sun) then
       write (22, *) 'solar weighting is enabled.'
       write (22, *) 'll_sun = ', ll_sun, 'lu_sun = ', lu_sun
    else
       write (22, *) 'solar weighting is disabled.'
    end if
    if (autodiff) then
       write (22, *) 'autodiff enabled; finite differences taken'
       write (22, '(2x,a,1pE9.3,a)') 'with ', automult, ' times initial variable.'
    end if
    write (22, '(1x, t2, a, t51, 1pe9.3)') &
         'convergence for relative parameter change set at ', provar
    write (22, *) 'times the differentiation parameter for each variable.'
    write (22, *) ' renorm = ', renorm
    write (22, *) ' update_pars = ', update_pars
    write (22, *) ' wrt_scr = ', wrt_scr
    write (22, '(a, f7.2)') '  szamax =', szamax
    write (22, '(a, f7.2)') '  szamin =', szamin
    write (22, '(a, I3)') '  cldmax =', cldmax
    write (22, '(a, f7.2)') '  latmax =', latmax
    write (22, '(a, f7.2)') '  latmin =', latmin
    write (22, '(a, i6)') '  nfirst =', nfirstfit
    write (22, '(a, f8.4)') '  phase =', phase
    
    write (22, *)
    write (22,'(1x, a)') 'Solar Fitting'
    write (22,'(T3,A9,T33,A13,T50,A6)') 'Parameter','Initial value','Varied'
    DO i = 1, n_solar_pars
       write (22,'(T4,A,T32,1pE13.5,T52,L)') TRIM(sun_par_names(i)), init_sun(i), if_var_sun(i)
    END DO
    
    write (22, *)
    write (22,'(1x, a)') 'Radiance Fitting'
    write (22,'(T3,A9,T33,A13,T50,A6)') 'Parameter','Initial value','Varied'
    DO i = 1, npars
       write (22,'(T4,A,T32,1pE13.5,T52,L)') TRIM(par_names(i)), initial(i), if_varied(i)
    END DO
    write (22, *)
    
    return
  end subroutine write_input

  subroutine write_output (nvaried, lista, iteration, write_fit, write_spec, &
       iterate, pos, spec, fit, var, initial, covar, chisq, delchi, correl, rms, &
       npoints, weight, provar, iprovar)
    
    ! Detailed fitting diagnostics. Best when only fitting a very few spectra!

    implicit real*8 (a - h, o - z)
    parameter (mmax = 64)
    parameter (maxpts = 7000)
    logical write_fit, write_spec, iterate, weight, iprovar
    real*8 initial (mmax)
    dimension pos (maxpts), spec (maxpts), fit (maxpts)
    dimension lista (mmax)
    dimension covar (mmax, mmax), correl (mmax, mmax)
    dimension var (mmax)
    save
    
    write (*, *) 'writing fitting output file'
    if (iterate) then
       write (22, '(/,1x,t2,a,t28,i3,t32,a)') 'convergence reached after', &
            iteration, ' iterations.'
       if (iprovar) then
          write (22, '(1x,t2,a)') &
               'convergence criterion is a change in every variable by less than or'
          write (22, '(1x,t2,a,t11,1pe9.3,t21,a)') &
               'equal to ', provar, 'times the differential parameters.'
          write (22, '(1x,t2,a,t26,1pe10.4,t36,a,//)') 'final value of chisq is', &
               chisq, '.'
       else
          write (22, '(1x,t2,a,t58,1pe10.4,t68,a)') &
               'convergence criterion is a change in chisq by less than', delchi, &
               ' for two'
          write (22, '(1x,t2,a,t50,1pe10.4,t60,a)') &
               'successive iterations.  final value of chisq is', chisq, '.'
       end if
       write (22, '(/,1x,t2,a,t8,1pe10.4,t18,a,/)') 'rms = ', rms, '.'
       write (22, '(1x,t32,a,/)') 'final parameters'
       write (22, '(1x,t2,a)') &
            'parameter    initial value    final value      std. deviation'
       do i = 1, nvaried
          j = lista (i)
          if (weight) then
             write (22, '(1x,t5,i2,t15,1pe14.7,t32,e14.7,t49,e14.7)') j, initial (j), &
                  var (j), rms * sqrt (covar (i, i) * float (npoints) / &
                  float (npoints - nvaried))
          else
             write (22, '(1x,t5,i2,t15,1pe14.7,t32,e14.7,t49,e14.7)') j, initial (j), &
                  var (j), rms * sqrt (covar (i, i) * float (npoints) / &
                  float (npoints - nvaried))
          end if
       end do
       write (22, '(/,1x,t31,a)') 'correlation matrix'
       do i = 1, nvaried
          write (22, '(/,1x,12(2x,f8.5))') (correl (i, j), j = 1, i)
       end do
    end if
    
    if (write_spec) then
       if (write_fit) write (*, *) 'writing spectrum output file'
       do i = 1, npoints
          write (24, '(f10.5, 1p3e13.5)') pos (i), spec (i), fit (i), spec (i) - &
               fit (i)
       end do
    end if
    
    return
  end subroutine write_output

  subroutine super_gauss (np, x, y, hw, sh, as, yc)
    
    ! Convolves input spectrum with Super Gaussian slit function of specified hw1e, !
    ! shape factor and asymmetric factor (half-width at 1/e intensity).
    ! For evenly-spaced spectra!!!
    IMPLICIT NONE
    
    ! Input variables
    INTEGER*4, INTENT(IN) :: np !Number of points
    REAL*8, INTENT(IN) :: hw, sh, as !HW1E, SHAPE, ASYM
    REAL*8, INTENT(IN), DIMENSION(1:np) :: x, y ! Wavelengths, spectrum
    
    ! Output variables
    REAL*8, INTENT(OUT), DIMENSION(1:np) :: yc ! Convolved spectrum
    
    ! Local variables
    INTEGER*4, PARAMETER :: nslit = 1000 ! Number of point in the slit function
    INTEGER*4 :: i, j, ns, nlo, nhi, nhalf, ss
    REAL*8, DIMENSION(1:nslit) :: slit
    REAL*8, DIMENSION(1:3*np) :: ytemp
    REAL*8, DIMENSION(1:np) :: ycc
    REAL*8 :: delpos, slitsum, wvl
    
    ! If there is no HW1E then no convolution to be done
    if (hw .eq. 0) then
       write (*, *) ' no gaussian convolution applied.'
       yc = y
       return
    end if
    
    ! ------------------------------------------------------------------
    ! Temporary input spectrum mirroed at the ends to ensure convolution
    ! ------------------------------------------------------------------
    ytemp(np+1:2*np) = y(1:np)
    do i = 1, np
       ytemp(np+1-i) = y(i)
       ytemp(2*np+i) = y(np+1-i)
    end DO
    
    ! Work out spacing of input spectrum
    delpos = x(2) - x(1)
    
    ! Apply slit function convolution.
    ! Calculate slit function values out to 0.001 times x0 value, normalize so that
    ! sum = 1.
    slitsum = 1.d0
    nhalf = nslit / 2
    slit = 0.d0
    slit(nhalf) = 1.d0
    getslit: do i = 1, nhalf-1
       ! Right branch
       wvl = delpos * i
       slit(nhalf+i) = EXP(-(ABS( wvl / ( hw + sign(REAL(1,KIND=8),wvl) * as ) ) )**sh)
       ! Left branch
       slit(nhalf-i) = EXP(-(ABS( wvl / ( hw + sign(REAL(1,KIND=8),-wvl) * as ) ) )**sh)
       ns = i
       if (slit (nhalf+i) / slit(nhalf) .le. 0.001) exit getslit
    end do getslit
    slitsum = SUM(slit(nhalf-ns:nhalf+ns))
    ! Normalize
    do i = nhalf-ns, nhalf+ns
       slit (i) = slit (i) / slitsum
    end do
    
    ! Convolve spectrum, reflect at endpoints.
    do i = 1, np
       
       ! ------------------------------------------------------
       ! Starting index of spectrum contributing to convolution
       ! ------------------------------------------------------
       ss = np+i-ns
       ! ----------------
       ! Safe convolution
       ! ----------------
       yc(i) = dot_product(slit(nhalf-ns:nhalf+ns), ytemp(ss:ss+2*ns))
    end do
    return
  end subroutine super_gauss
  
  subroutine super_gauss_uneven (np, x, y, hw, sh, as, yc)
    
    ! Convolves input spectrum with Super Gaussian slit function of specified hw1e,
    ! shape factor and asymmetric factor (half-width at 1/e intensity).
    ! For non evenly-spaced spectra!!!
    IMPLICIT NONE
     
    ! Input variables
    INTEGER*4, INTENT(IN) :: np !Number of points
    REAL*8, INTENT(IN) :: hw, sh, as !HW1E, SHAPE, ASYM
    REAL*8, INTENT(IN), DIMENSION(1:np) :: x, y ! Wavelengths, spectrum
    
    ! Output variables
    REAL*8, INTENT(OUT), DIMENSION(1:np) :: yc ! Convolved spectrum
    
    ! Local variables
    INTEGER*4, PARAMETER :: nslit = 1000 ! Number of point in the slit function
    INTEGER*4 :: i, j, ns, nlo, nhi, nhalf, ss
    REAL*8, DIMENSION(1:nslit) :: slit
    REAL*8, DIMENSION(1:3*np) :: ytemp, xtemp
    REAL*8, DIMENSION(1:np) :: ycc
    REAL*8 :: delpos, slitsum, wvl1, wvl2
    
    ! If there is no HW1E then no convolution to be done
    if (hw .eq. 0) then
       write (*, *) ' no gaussian convolution applied.'
       yc = y
       return
    end if
    
    ! ------------------------------------------------------------------
    ! Temporary input spectrum mirroed at the ends to ensure convolution
    ! ------------------------------------------------------------------
    ytemp(np+1:2*np) = y(1:np)
    xtemp(np+1:2*np) = x(1:np)
    do i = 1, np-1
       ytemp(np+1-i) = y(i)
       ytemp(2*np+i) = y(np-i)
       delpos = xtemp(np+i)-xtemp(np+i+1)
       xtemp(np+1-i) = xtemp(np+2-i)+delpos
       delpos = xtemp(2*np+1-i)-xtemp(2*np-i)
       xtemp(2*np+i) = xtemp(2*np+i-1)+delpos
    end DO
    
    ! Convolve spectrum, reflect at endpoints.
    do i = 1, np
       
       ! For each point work out the slit function
       ! Calculate slit function values out to 0.001 times x0 value, normalize so that
       ! sum = 1.
       nhalf = nslit / 2
       slit = 0.d0
       slit(nhalf) = 1.d0
       getslit: do j = 1, nhalf-1
          ! Right branch
          wvl1 = x(i)-xtemp(np+i+j)
          slit(nhalf+j) = EXP(-(ABS( wvl1 / ( hw + sign(REAL(1,KIND=8),wvl1) * as ) ) )**sh)
          ! Left branch
          wvl2 = x(i)-xtemp(np+i-j)
          slit(nhalf-j) = EXP(-(ABS( wvl2 / ( hw + sign(REAL(1,KIND=8),wvl2) * as ) ) )**sh)
          ns = j
          if (slit (nhalf+j) / slit(nhalf) .le. 0.001) exit getslit
       end do getslit
       slitsum = SUM(slit(nhalf-ns:nhalf+ns))
       ! Normalize
       do j = nhalf-ns, nhalf+ns
          slit (j) = slit (j) / slitsum
       end do
       
       ! ------------------------------------------------------
       ! Starting index of spectrum contributing to convolution
       ! ------------------------------------------------------
       ss = np+i-ns
       ! ----------------
       ! Safe convolution
       ! ----------------
       yc(i) = dot_product(slit(nhalf-ns:nhalf+ns), ytemp(ss:ss+2*ns))
    end do
    return
  end subroutine super_gauss_uneven
  
  subroutine spline(x, y, n, y2)
    ! Cubic spline derivative calculation from Numerical Recipes. Modified to always
    ! use "natural" boundary conditions (2nd derivatives = 0 at boundaries).
    IMPLICIT NONE
    
    ! Input variables
    INTEGER*4, INTENT(IN) :: n
    REAL*8, INTENT(IN), DIMENSION(1:n) :: x, y
    
    ! Output variables
    REAL*8, INTENT(INOUT), DIMENSION(1:n) :: y2
    
    ! Local variables
    REAL*8, DIMENSION(1:n) :: u
    REAL*8 :: sig, p, qn, un
    INTEGER*4 :: i, k
    
    y2 (1) = 0.
    u (1) = 0.
    do i = 2, n - 1
       sig = (x (i) - x (i - 1)) / (x (i + 1) - x (i - 1))
       p = sig * y2 (i - 1) + 2.
       y2 (i) = (sig - 1.) / p
       u (i) = (6. * ((y (i + 1) - y (i)) / (x (i + 1) - x (i)) - (y (i) - &
            y (i - 1)) / (x (i) - x (i - 1))) / (x (i + 1) - x (i - 1)) - sig * &
            u (i - 1)) / p
    end do
    qn = 0.
    un = 0.
    y2 (n) = (un - qn * u (n - 1)) / (qn * y2 (n - 1) + 1.)
    do k = n - 1, 1, -1
       y2 (k) = y2 (k) * y2 (k + 1) + u (k)
    end do
    
    return
  end subroutine spline
  
  subroutine splint (xa, ya, y2a, n, np, x, y)
    ! Cubic spline interpolation from Numerical Recipes, using results of
    ! subroutine spline.
    
    IMPLICIT NONE
    ! Input variables
    INTEGER*4, INTENT(IN) :: n, np
    REAL*8, INTENT(IN), DIMENSION(1:n) :: xa, ya, y2a
    REAL*8, INTENT(IN), DIMENSION(1:np) :: x
    
    ! Ouput variables
    REAL*8, INTENT(OUT), DIMENSION(1:np) :: y
    
    ! Local variables
    INTEGER*4 :: i, k, klo, khi
    REAL*8 :: h, a, b
    
    do i = 1, np
       klo = 1
       khi = n
       do while (khi - klo .gt. 1)
          k = (khi + klo) / 2
          if (xa (k) .gt. x (i)) then
             khi = k
          else
             klo = k
          end if
       end do
       h = xa (khi) - xa (klo)
       if (h .eq. 0.) print*, 'bad xa input.'
       a = (xa (khi) - x (i)) / h
       b = (x (i) - xa (klo)) / h
       y (i) = a * ya (klo) + b * ya (khi) + ((a**3 - a) * y2a (klo) + (b**3 - b) * &
            y2a (khi)) * (h**2) / 6.
    end do
    
    return
  end subroutine splint
  
  subroutine undersample (wav, nw, underspec, fraction)
    
    ! Convolves input spectrum with Gaussian slit function of specified HW1e (half-
    ! width at 1/e of maximum intensity), and samples at a particular input phase to
    ! give the undersampling spectrum, that is, the correction for not Nyquist
    ! sampling the spectra (see Applied Optics 44, 1296-1304, 2005 for more detail).
    ! This version calculates both phases of the undersampling spectrum, phase1 -
    ! i.e., underspec (1, i) - being the more common in radiance spectra. "KPNO"
    ! data means the high spectral resolution solar reference spectrum (now, most
    ! usually, a portion of the SAO2010 Solar Irradiance Reference Spectrum,
    ! available at http://www.cfa.harvard.edu/atmosphere/).
    IMPLICIT NONE
    
    ! Input variables
    INTEGER*4, INTENT(in) :: nw
    REAL*8, INTENT(in) :: fraction
    REAL*8, INTENT(in), DIMENSION(1:nw) :: wav
    
    ! Output variable
    REAL*8, INTENT(out), DIMENSION(2,1:nw) :: underspec
    
    ! Local variables
    INTEGER*4 :: i
    REAL*8, DIMENSION(1:nw) :: resample, wav2, over, under
    REAL*8, ALLOCATABLE, DIMENSION(:) :: d2spec, d2res
    
    ! The high resolution spectrum is already in memory (nkppos, kppos, kpspec
    ! and convolved kpspec_gauss)
    
    ! Allocate variables for interpolation
    ALLOCATE(d2spec(1:nkppos),d2res(1:nw))
    
    ! Phase1 calculation: Calculate spline derivatives for KPNO data.
    call spline (kppos, kpspec_gauss, nkppos, d2spec)
    ! Calculate solar spectrum at radiance positions
    call splint (kppos, kpspec_gauss, d2spec, nkppos, nw, wav, resample)
    ! Calculate spline derivatives for resampled data.
    call spline (wav, resample, nw, d2res)
    
    ! Calculate solar spectrum at radiance + phase positions, original and
    ! resampled.
    wav2 (1) = wav (1)
    do i = 2, nw
       wav2 (i) = (1.d0 - fraction) * wav (i - 1) + fraction * wav (i)
    end do
    
    call splint (kppos, kpspec_gauss, d2spec, nkppos, nw, wav2, over)
    call splint (wav, resample, d2res, nw, nw, wav2, under)
    do i = 1, nw
       underspec (1, i) = over (i) - under (i)
       resample (i) = over (i)
    end do
    
    ! Phase2 calculation: Calculate spline derivatives for KPNO data.
    ! Calculate spline derivatives for resampled data.
    call spline (wav2, resample, nw, d2res)
    ! Calculate solar spectrum at radiance positions, original and resampled.
    call splint (kppos, kpspec_gauss, d2spec, nkppos, nw, wav, over)
    call splint (wav2, resample, d2res, nw, nw, wav, under)
    ! Calculate undersample spectrum
    do i = 1, nw
       underspec (2, i) = over (i) - under (i)
    end do
    
    ! Deallocate interpolation variables
    DEALLOCATE(d2spec,d2res)
    
    return
  end subroutine undersample
  
  subroutine dataspline (wav, nwav, hw, sh, as)

    ! Opens reference spectrum files.
    ! Convolves reference spectra with supplied slit width if requested.
    ! Samples reference spectra to standard grid determined from wavelength
    ! calibration of a radiance spectrum.
    ! Loads spectra into the database file.
    
    ! Local variables
    IMPLICIT NONE
    
    ! Input variables
    INTEGER*4, INTENT(IN) :: nwav
    REAL*8, INTENT(IN) :: hw, sh, as
    REAL*8, INTENT(IN), DIMENSION(1:nwav) :: wav
    
    ! Local variables
    INTEGER*4 :: i, j, np
    REAL*8, ALLOCATABLE, DIMENSION(:) :: x, y,yc,yt,der2spl
    
    ! Loop over fitting variables
    do i = 1, npars
       ! Only process those parameters with cross-section files
       If ( (par_str(i)   .eq. ad1_str .or. par_str(i) .eq. ad2_str .or. par_str(i) .eq. ble_str) .and. &
            (par_names(i) .ne. us1_str .and. par_names(i) .ne. us2_str) ) THEN
          ! Open file
          OPEN(unit=11, file=TRIM(par_names(i)), status='old')
          ! Find number of points
          np=0
          do
             read(11,*,END=11)
             np=np+1
          end do
11        rewind(unit=11)
          ! Log message
          if (wrt_scr) write(*,'(A,I5,A)') 'Reading file... '//TRIM(par_names(i))//'  with ',np,' points'
          ! Allocate variable to hold wavlengths and spectrum
          ALLOCATE(x(1:np),y(1:np),yc(1:np),der2spl(1:np),yt(1:nwav))
          ! Read wavelengths and spectrum
          DO j = 1, np
             READ(11,*) x(j), y(j)
          END DO
          ! Convolve and spline spectrum
          call super_gauss_uneven(np,x,y,hw,sh,as,yc)
          call spline (x, yc, np, der2spl)
          call splint (x, yc, der2spl, np, nwav, wav, yt)
          
          ! Save cross section in reference spectrum
          database (i,1:nwav) = yt(1:nwav)
          
          ! Deallocate variables
          DEALLOCATE(x,y,yc,der2spl,yt)
          CLOSE(11)
       END If
    end do
    
    return
  end subroutine dataspline

  subroutine ludcmp (a, n, np, indx, d)
    
    ! LU decomposition, from Numerical Recipes.
    IMPLICIT NONE
    
    ! Input variables
    INTEGER*4, INTENT(IN) :: n, np
    
    ! Modified variables
    REAL*8, INTENT(INOUT), DIMENSION(1:np,1:np) :: a
    REAL*8, INTENT(INOUT) :: d
    INTEGER*4, INTENT(INOUT), DIMENSION(1:np) :: indx
    
    ! Local variables
    REAL*8, PARAMETER :: tiny = 1.0d-20
    REAL*8, DIMENSION(1:n) :: vv 
    REAL*8 :: aamax, sum, dum
    INTEGER*4 :: i, j, k, imax
    d = 1.
    
    do i = 1, n
       aamax = 0.d0
       do j = 1, n
          if (abs (a (i, j)) .gt. aamax) aamax = abs (a (i, j))
       end do
       if (aamax .eq. 0.) print*, 'singular matrix.'
       vv (i) = 1.d0 / aamax
    end do
    do j = 1, n
       if (j .gt. 1) then
          do i = 1, j - 1
             sum = a (i, j)
             if (i .gt. 1)then
                do k = 1, i - 1
                   sum = sum - a (i, k) * a (k, j)
                end do
                a (i, j) = sum
             end if
          end do
       end if
       aamax = 0.d0
       do i = j, n
          sum = a (i, j)
          if (j .gt. 1)then
             do k = 1, j - 1
                sum = sum - a (i, k) * a (k, j)
             end do
             a (i, j) = sum
          end if
          dum = vv (i) * abs (sum)
          if (dum .ge. aamax) then
             imax = i
             aamax = dum
          end if
       end do
       if (j .ne. imax) then
          do k = 1, n
             dum = a (imax, k)
             a (imax, k) = a (j, k)
             a (j, k) = dum
          end do
          d = -d
          vv (imax) = vv (j)
       end if
       indx (j) = imax
       if (j .ne. n) then
          if (a (j, j) .eq. 0.) a (j, j) = tiny
          dum = 1.d0 / a (j, j)
          do i = j + 1, n
             a (i, j) = a (i, j) * dum
          end do
       end if
    end do
    if (a (n, n) .eq. 0.) a (n, n) = tiny
    
    return
  end subroutine ludcmp
  
  subroutine lubksb (a, n, np, indx, b)

    ! LU forward- and back-substitution, from Numerical Recipes.
    IMPLICIT NONE
    
    ! Input variables
    INTEGER*4, INTENT(IN) :: n, np
    
    ! Modified variables
    REAL*8, INTENT(INOUT), DIMENSION(1:np,1:np) :: a
    REAL*8, INTENT(INOUT), DIMENSION(1:n) :: b
    INTEGER*4, INTENT(INOUT), DIMENSION(1:n) :: indx
    
    ! Local variables
    REAL*8 :: sum
    INTEGER*4 :: ii, i, ll, j
    
    ii = 0
    do i = 1, n
       ll = indx (i)
       sum = b (ll)
       b (ll) = b (i)
       if (ii .ne. 0) then
          do j = ii, i - 1
             sum = sum-a (i, j) * b (j)
          end do
       else if (sum .ne. 0.) then
          ii = i
       end if
       b (i) = sum
    end do
    do i = n, 1, -1
       sum = b (i)
       if (i .lt. n) then
          do j = i + 1, n
             sum = sum - a (i, j) * b (j)
          end do
       end if
       b (i) = sum / a (i, i)
    end do
    
    return
  end subroutine lubksb
  
END MODULE GASFIT_MODULE

SUBROUTINE funcs (pos, npoints, vars, fac, str, ymod, dyda, npars, lista, nvaried, &
     avg, diff, ntype)

  USE gasfit_module, ONLY: spectrum

  ! Calculates the spectrum and its derivatives for mrqcof.
  IMPLICIT NONE
  
  ! Input variables
  INTEGER*4, INTENT(IN) :: npoints, npars, nvaried, ntype
  REAL*8 :: avg
  REAL*8, INTENT(IN), DIMENSION(1:npoints) :: pos
  REAL*8, INTENT(IN), DIMENSION(1:npars) :: diff, vars, fac
  INTEGER*4, INTENT(IN), DIMENSION(1:npars) :: lista
  CHARACTER(3), INTENT(IN), DIMENSION(1:npars) :: str
  
  ! Modified variables
  REAL*8, INTENT(INOUT), DIMENSION(1:npoints) :: ymod
  REAL*8, INTENT(INOUT), DIMENSION(1:npars,1:npoints) :: dyda
  
  ! Local variables
  INTEGER*4 :: i, j
  REAL*8 :: var0
  REAL*8, DIMENSION(1:npoints) :: dyplus
  REAL*8, DIMENSION(1:npars) :: vars0

  ! Calculate the spectrum.
  call spectrum (npoints, npars, avg, pos, ymod, vars, fac, str, ntype)

  ! Calculate the derivatives by finite differences. dyplus is y (var + finite
  ! difference in variable)
  vars0 = vars
  do i = 1, nvaried
     var0 = vars0 (lista (i))
     vars0 (lista (i)) = vars0 (lista (i)) + diff (lista (i))
     call spectrum (npoints, npars, avg, pos, dyplus, vars0, fac, str, ntype)
     do j = 1, npoints
        dyda (lista (i), j) = (dyplus (j) - ymod (j)) / diff (lista (i))
     end do
     vars0 (lista (i)) = var0
  end do
  
  return
END SUBROUTINE funcs
