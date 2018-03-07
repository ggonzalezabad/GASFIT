PROGRAM GASBETA

! Fits an orbit of satellite-measured Earth radiance spectra, or a sequence of
! ground-based atmospheric spectra.
!
! This version of the gas fitting code includes the option to do BOAS or DOAS
! fitting, that is to fit radiances directly, or to fit a high-pass filtered
! logarithm of (radiance / irradiance). BOAS normally gives an improvement of a
! factor of 2-3 in the fitting statistics, with average parameter values very
! close to the same. DOAS usually runs somewhat faster.

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
CHARACTER(120):: fitin, fitout, specin, specout, &
     inputline, general_line, solar_line, radiance_line
INTEGER*4, PARAMETER :: mmax = 64
INTEGER*4, PARAMETER :: maxpts = 7000
INTEGER*4, PARAMETER :: maxpix = 2000
REAL*8, PARAMETER ::pi = 3.14159265358979d0
LOGICAL :: wrt_scr, smooth, write_fit, write_spec, mirror, autodiff, doas, &
  if_residuals, iterate_sun, weight_sun, if_var_sun (mmax), iterate_rad, &
  weight_rad, renorm, update_pars, if_varied (mmax), iprovar
REAL*8, DIMENSION(mmax) :: initial, init_sun, var_sun, diffsun, var, diff
REAL*8, DIMENSION(maxpts) :: pos_sun, spec_sun, sig_sun, spline_sun, &
     deriv2sun, pos, spec, sig, fit, temp, residual
REAL*8, DIMENSION(2,maxpts) :: underspec
REAL*8, DIMENSION(11,maxpts) :: database
REAL*8, DIMENSION(maxpix,maxpts) :: pos_rad, spec_rad
REAL*8, DIMENSION(mmax,mmax) :: correl, covar
REAL*8 :: latmax, latmin, asum, automult, avg, chisq, &
     davg, delchi, dgas, dhw1e, div_rad, div_sun, drelavg, dshift, dshiftavg, &
     dshiftrad, gas, hw1e, phase, provar, remult, report_mult, rms, rmsavg, &
     shift, sigsum, squeeze, ssum, ilat, ilon, isza, isaa, ivza, ivaa, &
     szamax, szamin
INTEGER*4 :: i, ipix, icld, iyear, imonth, iday, ihour, imin, isec, &
     iteration, j, ll_rad, ll_sun, lu_rad, lu_sun, &
     n_solar_pars, nfirst,nfirstfit, ngas, nhw1e, npars, npoints, nrads, nreport, &
     nshift, nshiftrad, nvar_sun, nvaried, cldmax
INTEGER*4, DIMENSION(mmax) :: list_sun, list_rad
INTEGER*4, DIMENSION(maxpix) :: npix, cld, year, month, day, hour, minu, sec
REAL*8, DIMENSION(maxpix) :: lat, lon, sza, saa, vza, vaa

! Initialize the residual spectrum here.
residual = 0.d0

write (*,'(5x, a)') 'enter fitting input file.'
read (*, '(a)') fitin
open (unit = 21, file = fitin, status='old')

! Read general parameters:
!
! inputline and general_line are comments that may be echoed in the output.
!
! npoints is the number of spectral points in the irradiance and each radiance.
!
! delchi and provar are fitting convergence criteria: a change in chi-squared
! and a relative change in each fitted parameter.
!
! If (wrt_scr) enter a more verbose manner of running.
!
! If (smooth) smooth the spectra and reference spectra before fitting. not
! recommended, but implemented for comparison with results of colleagues who
! smooth rather than correct for spectral undersampling.
!
! If (write_fit) write conditions and diagnostics to an output file. this is
! not the same as writing out the fitting results, which is done automatically.
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
! If (doas) then do a DOAS (differential optical absorption spectroscopy) fit
! rather than a BOAS (basic optical absorption spectroscopy) fit. Useful to
! compare results and processing speed. DOAS was developed for ground-based
! measurements where the baseline could not reliably be fitted. For satellite
! measurements it is normally not necessary nor an improvement to the fitting
! results.
!
! If (if_residuals) write out the averaged fitting residual at the end of the
! spectral fitting results. Used to generate the "common mode" correction for
! systematic fitting residuals initially implemented at SAO.
!
read (21, '(a)') inputline
read (21, '(a)') general_line
read (21, *) npoints, delchi, provar, automult
read (21, *) wrt_scr, smooth, write_fit, write_spec, mirror, autodiff, doas, &
  if_residuals

! Read solar fitting parameters:
!
! solar_line is a comment that may be echoed in the output.
!
! f (iterate_sun) perform wavelength calibration and slit width calibration on
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
do i = 1, n_solar_pars
  read (21, *) var_sun (i), if_var_sun (i), diffsun (i)
end do

! Automatic difference-taking
if (autodiff) then
  do i = 1, n_solar_pars
    if (if_var_sun (i)) diffsun (i) = dabs (automult * var_sun (i))
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
  if (i .eq. 10) nhw1e = nvar_sun
  if (i .eq. 11) nshift = nvar_sun
end do
if (nvar_sun .gt. mmax) then
  write (22, *) ' maximum number of solar parameters exceeded.'
  go to 1000
end if
if (wrt_scr) write (*, *) 'nvar_sun =', nvar_sun

! Read radiance fitting parameters:
!
! radiance_line is a comment that may be echoed in the output.
!
! If (iterate_rad) fit the radiances.
!
! If (renorm) renormalize radiances to their weighted average before fitting.
! Numerical overflow insurance. Good for BOAS, don't normally use with DOAS.
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
! scaled by the  same amount (if (.not. doas)). If (doas), report_mult should be
! the negative of the scalar for the reported parameter.
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
  div_rad, phase, nreport, report_mult
read (21, *) szamax, szamin, latmax, latmin, cldmax, nfirstfit
i = 1
20 read (21, *, end = 30) var (i), if_varied (i), diff (i)
  i = i + 1
  go to 20
30 npars = i - 1

! Automatic difference-taking.
if (autodiff) then
  do i = 1, npars
    if (if_varied (i)) diff (i) = dabs (automult * var (i))
  end do
end if

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
  if (i .eq. 11) nshiftrad = nvaried
end do
if (nvaried .gt. mmax) then
  write (22, *) ' maximum number of radiance parameters exceeded.'
  go to 1000
end if
if (wrt_scr) write (*, *) 'nvaried =', nvaried
if (wrt_scr) write (*, *) 'Max. sza =', szamax
if (wrt_scr) write (*, *) 'Min. sza =', szamin
if (wrt_scr) write (*, *) 'Max. lat =', latmax
if (wrt_scr) write (*, *) 'Min. lat =', latmin
if (wrt_scr) write (*, *) 'Max. cld =', cldmax

! Open fitting output file.
if (wrt_scr) write (*,'(5x, a)') 'enter fitting output file.'
read (*, '(a)') fitout
open (unit = 22, file = fitout, status='unknown')

! Open irradiance and radiance level 1 file.
if (wrt_scr) write (*,'(5x, a)') 'enter spectrum input file.'
read (*, '(a)') specin
open (unit = 23, file = specin, status='old')

! Open spectrum output file.
if (write_spec) then
  if (wrt_scr) write (*,'(5x, a)') 'enter spectrum output file.'
  read (*, '(a)') specout
  open (unit = 24, file = specout, status='unknown')
end if

! Read the solar spectrum. A 2-column (position and irradiance value) spectrum
! of npoints points, with no header, is expected.
do i = 1, npoints
  read (23, *) pos_sun (i), spec_sun (i)
! divide to avoid overflow issues.
  spec_sun (i) = spec_sun (i) / div_sun
! select fitting window by the use of weighting.
  sig_sun (i) = 1.d6
  if (i .ge. ll_sun .and. i .le. lu_sun) sig_sun (i) = 1.d0
end do

! Initialize several diagnostics.
! rmsavg: average fitting rms.
! davg: average fitting uncertainty for reported parameter.
! drelavg: average relative fitting uncertainty for reported parameter.
! dshiftavg:  average fitting uncertainty for spectral shift parameter.
! nfirst is the number among fitted radiances corresponding to the ground pixel
! nfirstfit.
rmsavg = 0.d0
davg = 0.d0
drelavg = 0.d0
dshiftavg = 0.d0
nfirst = 0

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
i = 1
40 READ (23, *, end = 50) ipix, icld
  READ (23, *) iyear, imonth, iday, ihour, imin, isec
  READ (23, *) ilat, ilon
  READ (23, *) isza, isaa, ivza, ivaa

! Decide whether to process. If so, update the appropriate arrays and the
! counter for radiances being fitted.
  IF ( (icld .LE. cldmax) .AND. (isza .LE. szamax) .AND. (isza .GE. szamin) .AND. &
       (ilat .LE. latmax) .AND. (ilat .GE. latmin) ) THEN
     npix(i) = ipix
     cld(i) = icld
     year(i) = iyear
     month(i) = imonth
     day(i) = iday
     hour(i) = ihour
     minu(i) = imin
     sec(i) = isec
     sza(i) = isza
     saa(i) = isaa
     vza(i) = ivza
     vaa(i) = ivaa
     lat(i) = ilat
     lon(i) = ilon
     !   Read radiance spectrum (positions and radiances) for each ground pixel.
     do j = 1, npoints
        read (23, *) pos_rad (i, j), spec_rad (i, j)
        spec_rad (i, j) = spec_rad (i, j) / div_rad ! usually 1.d13 or 1.d14
     end do
     !   Is this the pixel selected for full wavelength calibration?
     if (ipix .eq. nfirstfit) nfirst = i
     i = i + 1
     !   Ignore spectra not selected for processing.
  else 
     do j = 1, npoints
        read (23, *)
     end do
  end if
  go to 40
50 nrads = i - 1
  if (nrads .eq. 0) then
     write (*, *) 'nrads = 0'
     stop
  end if

! if there is no spectrum corresponding to nfirstfit, select nfirst halfway
! along the orbit.
if (nfirst .eq. 0) nfirst = CEILING(FLOAT(nrads)/2.0)

! Add the weighting, in the single array.
do i = 1, npoints
  sig(i) = 1.d6
  if (i .ge. ll_rad .and. i .le. lu_rad) sig(i) = 1.d0
end do

! Write out the input line and other input information.
IF (write_fit) WRITE (22, '(a)') inputline
IF (write_fit) CALL write_input (npoints, npars, n_solar_pars, nvaried, &
  nvar_sun, iterate_rad, iterate_sun, weight_rad, weight_sun, if_varied, &
  if_var_sun, initial, init_sun, nfirst, update_pars, provar, autodiff, &
  ll_rad, lu_rad, ll_sun, lu_sun, phase, wrt_scr, smooth, renorm, szamax, &
  szamin, cldmax, latmax, latmin, doas)

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
  call specfit (npoints, n_solar_pars, nvar_sun, list_sun, &
  iteration, smooth, avg, spec_sun, pos_sun, sig_sun, fit, var_sun, chisq, &
  delchi, diffsun, correl, covar, rms, provar, iprovar, wrt_scr, 1, database, &
  doas)
print*, iprovar
stop
! Shift and squeeze solar spectrum.
  do i = 1, npoints
    pos_sun (i) = (pos_sun (i) - var_sun (11)) / (1.d0 + var_sun (12))
  end do
  dshift = rms * sqrt (covar (nshift, nshift) * float (npoints) / &
    float (npoints - nvar_sun))
  write (*, '(a, 1pe11.3)') 'solar wavelength calibration: rms = ', rms
  write (*, '(a, 1p2e14.6)') 'irrad: shift, 1 sigma = ', - var_sun (11), dshift
  write (*, *) 'nvar_sun = ', nvar_sun

! Freeze variation of slit width for radiance fitting is not currently
! implemented.
! if_var_sun (10) = .false.

! Set width for later use in undersampling correction. hw1e means the gaussian
! half-width at 1/e of the maximum intensity.
  hw1e = var_sun (10)
  dhw1e = rms * sqrt (covar (nhw1e, nhw1e) * float (npoints) / &
    float (npoints - nvar_sun))
  write (*, '(a, 1p2e14.6)') 'irrad: hw1e, 1 sigma = ', hw1e, dhw1e

! Order the coefficients for mrqmin and save the initial values of the
! coefficients (uncomment if slit variation freeze is to be implemented).
! nvar_sun = 0
! do i = 1, n_solar_pars
!   init_sun (i) = var_sun (i)
!   if (if_var_sun (i)) then
!     nvar_sun = nvar_sun + 1
!     list_sun (nvar_sun) = i
!   end if
! end do

! Select and wavelength calibrate selected radiance spectrum.
  write (*, *) ' nfirst = ', nfirst
  do i = 1, npoints
    pos (i) = pos_rad (nfirst, i)
    spec (i) = spec_rad (nfirst, i)
  end do
  if (.not. weight_sun) then
    avg = (pos (npoints) + pos (1)) / 2.
  else
    asum = 0.
    ssum = 0.
    do i = 1, npoints
      asum = asum + pos (i) / (sig (i)**2)
      ssum = ssum + 1. / (sig (i)**2)
    end do
    avg = asum / ssum
  end if
  iprovar = .false.
  call specfit (npoints, n_solar_pars, nvar_sun, list_sun, iteration, smooth, &
    avg, spec, pos, sig, fit, var_sun, chisq, delchi, diffsun, correl, covar, &
    rms, provar, iprovar, wrt_scr, 2, database, doas)
  dshift = rms * sqrt (covar (10, 10) * float (npoints) / float (npoints - &
    nvar_sun))
  write (*, '(a, 1p2e14.6)') 'rad: shift, 1 sigma = ', - var_sun (11), dshift

! Save parameters, to shift and squeeze radiance spectra. Remember SIGN of shift
! determined this way, because it is minus the shift required to wavelength
! correct the satellite radiance/irradiance spectra!
  shift = var_sun (11)
  squeeze = 1.d0 + var_sun (12)

! Apply them to current positions, for undersampling and data-splining
! calculations.
  do i = 1, npoints
    pos (i) = (pos (i) - shift) / squeeze
  end do
! else (later) calculate sun and radiance spectrum
end if ! on iterate_sun
write (*, *) 'finished with iterate_sun'

! Calculate the undersampled spectrum.
call undersample (pos, npoints, underspec, hw1e, phase, wrt_scr)

! Calculate the splined fitting database. Note that the undersampled
! spectrum has just been done. Finish filling in reference database array.
call dataspline (pos, npoints, database, wrt_scr)

! Spline irradiance spectrum onto radiance grid.
call spline (pos_sun, spec_sun, npoints, deriv2sun)
call splint (pos_sun, spec_sun, deriv2sun, npoints, npoints, pos, spline_sun)
do i = 1, npoints
  database (1, i) = spline_sun (i)
  database (9, i) = underspec (1, i)
  database (10, i) = underspec (2, i)
! Set up Ring spectrum for high pass filtering of divided spectrum (afterward,
! add solar spectrum back in order to divide by solar spectrum with altered
! wavelength calibration in subroutine spectrum) - DOAS fitting only.
  if (doas) then
    database (2, i) = database (2, i) / spline_sun (i)
  end if
end do

! For the DOAS case, high-pass filter the reference spectra.
if (doas) then
  call subtract_cubic (pos, npoints, database, ll_rad, lu_rad)
  do i = 1, npoints
    database (2, i) = database (2, i) * spline_sun (i)
  end do
end if

! Load database using nfirst wavelength scale.
call spectrum (npoints, npars, smooth, avg, pos, fit, var, wrt_scr, 3, &
  database, doas)

! Fit radiance spectra that meet the input file selection requirements.
do i = 1, nrads
  if (mod (i, 100) .eq. 0) write (*, *) ' nrad = ', i
  do j = 1, npoints
    pos (j) =  (pos_rad (i, j) - shift) / squeeze
    spec (j) =  spec_rad (i, j)
  end do

! Renormalize, if requested (normally don't do for DOAS since we are dividing
! by the solar spectrum just below, and logarithmizing).
  if (renorm) then
    remult = 0.d0
    if (weight_rad) then
      sigsum = 0.d0
      do j = 1, npoints
        remult = remult + spec (j) / (sig (j) * sig (j))
        sigsum = sigsum + 1.d0 / (sig (j) * sig (j))
      end do
      remult = remult / sigsum
    else
      do j = 1, npoints
        remult = remult + spec (j)
      end do
      remult = remult / npoints
    end if
    do j = 1, npoints
      spec (j) = spec (j) / remult
    end do
  end if

! High pass filtering for DOAS. First, take log (rad/irrad), then filter by
! subtracting a cubic, then re-add the log (irradiance). This way we are fitting
! the log (rad) with the proper filtering having been done. The spectrum
! subroutine will start by subtracting the irradiance spectrum, where the proper
! shifting and squeezing can take place. In this high-pass filtering, we ignore
! the small extra wavelength calibration change, as for Ring effect, above.
  if (doas) then
    do j = 1, npoints
      spec (j) = dlog (spec (j) / spline_sun (j))
    end do
    call subtract_cubic_meas (pos, npoints, spec, ll_rad, lu_rad)
    do j = 1, npoints
      spec (j) = spec (j) + dlog (spline_sun (j))
    end do
  end if

! Apply smoothing (1/16,1/4,3/8,1/4,1/16); 2/98 UHE/IFE recommendation: This is
! done by several colleagues as an alternative to undersampling correction. It
! is included here in case it is needed for comparison, but is not normally
! needed or used (or desired).
  if (smooth) then
    do j = 1, npoints
      temp (j) = spec (j)
    end do
    do j = 3, npoints - 2
      spec (j) = 0.375d0 * temp (j) + 0.25d0 * (temp (j + 1) + temp (j - 1)) + &
        0.0625d0 * (temp (j + 2) + temp (j - 2))
    end do
  end if
  if (.not. weight_rad) then
    avg = (pos (npoints) + pos (1)) / 2.
  else
    asum = 0.
    ssum = 0.
    do j = 1, npoints
      asum = asum + pos (j) / (sig (j)**2)
      ssum = ssum + 1. / (sig (j)**2)
    end do
    avg = asum / ssum
  end if
  iprovar = .false.
  if (iterate_rad) then
    if (.not. update_pars) then
      if (wrt_scr) write (*, *) ' updating parameters'
      do j = 1, npars
        var (j) = initial (j)
      end do
    end if

!   BOAS or DOAS fitting now done here.
    call specfit (npoints, npars, nvaried, list_rad, iteration, smooth, avg, &
      spec, pos, sig, fit, var, chisq, delchi, diff, correl, covar, rms, &
      provar, iprovar, wrt_scr, 4, database, doas)

!   write out radiance fit: keep here as possible future diagnostic.
!   do j = ll_rad, lu_rad
!     write (24, '(f11.6, 1p3e13.5)') pos (j), spec (j), fit (j), &
!     spec (j) - fit (j)
!   end do

!   Write out fitting results.
    gas = report_mult * var (nreport)
    dgas = report_mult * rms * sqrt (covar (ngas, ngas) * float (npoints) / &
      float (npoints - nvaried))
    write (22, '(i5, 1p3e12.4)') npix(i), gas, dgas, rms
    write (22, '(4f7.2)') sza(i), saa(i), vza(i), vaa(i)
    write (22, '(10f7.2)') lat(i), lon(i)
    dshiftrad = rms * sqrt (covar (nshiftrad, nshiftrad) * float (npoints) / &
      float (npoints - nvaried))
    dshiftavg = dshiftavg + dshiftrad
    rmsavg = rmsavg + rms
    davg = davg + dgas
    drelavg = drelavg + dabs (dgas / gas)

! Write out extra diagnostics.
if (write_fit) &
  call write_output (nvaried, list_rad, iteration, write_fit, write_spec, &
  iterate_rad, pos, spec, fit, var, initial, covar, chisq, delchi, correl, &
  rms, npoints, weight_rad, provar, iprovar)

! else (TBD later) calculate sun and radiance spectrum, no iteration.
  end if ! on iterate_rad
  if (if_residuals) then
    do j = 1, npoints
      residual (j) = residual (j) + spec (j) - fit (j)
    end do
  end if
end do ! on nrads.

! Write out the average fitting statistics to screen.
write (*, '(a33, 1pe13.5)') '                       Avg rms = ', rmsavg / nrads
write (*, '(a33, 1pe13.5)') '                     Avg dgas = ', davg / nrads
write (*, '(a33, 1pe13.5)') ' Avg relative gas uncertainty = ', drelavg / nrads
write (*, '(a33, 1pe13.5)') ' Avg radiance shift uncertainty = ', dshiftavg / &
  nrads

! Write out the fitting statistics and average residuals.
if (if_residuals) then
  write (22, '(a28, 1pe13.5)') '                  Avg rms = ', rmsavg / nrads
  write (22, '(a28, 1pe13.5)') '                Avg dgas = ', davg / nrads
  write (22, '(a28, 1pe13.5)') ' Avg rel gas uncertainty = ', drelavg / nrads
  do i = 1, npoints
    write (22, '(f10.5, 1pe13.5)') pos (i), residual (i) / nrads
  end do
end if

! If (mirror), mirror the input file, with updated parameters, onto the output
! fitting file.
if (mirror) then
  write (22, *)
  write (22, '(a)') inputline
  write (22, '(a)') general_line
  write (22, '(i5, 1p2e12.4, 0pf8.4)') npoints, delchi, provar, automult
  write (22, *) wrt_scr, smooth, write_fit, write_spec, mirror, autodiff, &
    doas, if_residuals
  write (22, '(a)') solar_line
  write (22, '(2l2, 3i4, 1pe9.2)') iterate_sun, weight_sun, n_solar_pars, &
    ll_sun, lu_sun, div_sun
  do i = 1, n_solar_pars
    write (22, '(1pe15.7, l4, 1pe14.5)') var_sun (i), if_var_sun (i), &
    diffsun (i)
  end do
  write (22, '(a)') radiance_line
  write (22, '(4l2, 2i4, 1pe9.2, 0pf5.2)') iterate_rad, renorm, weight_rad, &
    update_pars, ll_rad, lu_rad, div_rad, phase
  write (22, '(4f7.2, i6)') szamax, szamin, latmax, latmin, nfirstfit
  do i = 1, npars
    write (22, '(1pe15.7, l4, 1pe14.5)') var (i), if_varied (i), diff (i)
  end do
end if

1000 close (unit = 21)
close (unit = 22)
close (unit = 23)
close (unit = 24)

stop
end program
!
subroutine specfit (npoints, npars, nvaried, lista, iteration, smooth, avg, &
  spec, pos, sig, fit, var, chisq, delchi, diff, correl, covar, rms, provar, &
  iprovar, wrt_scr, ntype, database, doas)

! Driver for Levenberg-Marquardt nonlinear least squares minimization.

implicit real*8 (a - h, o - z)
external funcs
parameter (mmax = 64)
parameter (maxpts = 7000)
logical iprovar, wrt_scr, smooth, doas
dimension spec (maxpts), pos (maxpts), sig (maxpts), fit (maxpts)
dimension lista (mmax)
dimension covar (mmax, mmax), alpha (mmax, mmax), correl (mmax, mmax)
dimension var (mmax), diff (mmax)
dimension var0 (mmax)
dimension database (11, maxpts)
save

alamda = - 1
call mrqmin (pos, spec, sig, npoints, var, npars, lista, nvaried, covar, &
  alpha, chisq, funcs, alamda, delchi, smooth, avg, diff, wrt_scr, ntype, &
  database, doas)

iteration = 1
itest = 0
1 if (wrt_scr) then
  write (22,'(/1x, a, i2, t18, a, 1pe10.4, t44, a, e9.2)') 'iteration #', &
    iteration, 'chi-squared: ', chisq, 'alamda:', alamda
  write (22,'(1x, a)') 'variables:'
  write (22,'(1x, 1p6e13.4)') (var (lista (i)), i = 1, nvaried)
  write (*,'(/1x, a, i2, t18, a, 1pe10.4, t44, a, e9.2)') 'iteration #', &
    iteration, 'chi-squared: ', chisq, 'alamda:', alamda
  write (*,'(1x, a)') 'variables:'
  write (*,'(1x, 1p6e12.4)') (var (lista (i)), i = 1, nvaried)
end if
iteration = iteration + 1
ochisq = chisq
do i = 1, nvaried
  var0 (lista (i)) = var (lista (i))
end do
call mrqmin (pos, spec, sig, npoints, var, npars, lista, nvaried, covar, &
  alpha, chisq, funcs, alamda, delchi, smooth, avg, diff, wrt_scr, ntype, &
  database, doas)
do i = 1, nvaried
  prop = abs (var0 (lista (i)) - var (lista (i)))
  difftest = provar * diff (lista (i))
  if (prop .gt. difftest) then
    go to 2
  end if
end do
iprovar = .true.
if (wrt_scr) write (*, '(1x, a)') 'iprovar = .true.'
go to 3
2 if (chisq .gt. ochisq) then
  itest = 0
else if (abs (ochisq - chisq) .lt. delchi) then
  itest = itest + 1
end if
if (itest .lt. 2) then
  go to 1
end if
3 alamda = 0.0
call mrqmin (pos, spec, sig, npoints, var, npars, lista, nvaried, covar, &
  alpha, chisq, funcs, alamda, delchi, smooth, avg, diff, wrt_scr, ntype, &
  database, doas)
if (wrt_scr) then
  write (22,'(/1x, a, i2, t18, a, e10.4, t43, a, e9.2)') 'iteration #', &
    iteration, 'chi-squared:', chisq, 'alamda:', alamda
  write (*,'(/1x, a, i2, t18, a, e10.4, t43, a, e9.2)') 'iteration #', &
    iteration, 'chi-squared:', chisq, 'alamda:', alamda
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
call spectrum (npoints, npars, smooth, avg, pos, fit, var, wrt_scr, ntype, &
  database, doas)

! Calculate the rms of the fit.
rsum = 0.
do i = 1, npoints
  rsum = rsum + ((fit (i) - spec (i)) / sig (i))**2
end do
rms = sqrt (rsum / npoints)

return
end subroutine specfit
!
subroutine spectrum (npoints, npars, smooth, avg, pos, fit, var, wrt_scr, &
  ntype, database, doas)

! Spectrum calculation for both fitting and non-fitting cases. For fitting, it
! is called with 4 different qualifiers ("ntypes") for purposes of determining
! the calculation type and initializing the fitting reference database.
! ntype = 1: Solar wavelength calibration. Open kpnospec, the high spectral
! resolution irradiance reference spectrum.
! ntype = 2: Radiance wavelength calibration (kpnospec already open).
! ntype = 3: First radiance fit (fill refspec from database).
! ntype = 4: Later radiance fits (refspec already filled).

implicit real*8 (a - h, o - z)
parameter (mmax = 64)
parameter (maxpts = 7000)
parameter (maxkpno = 7000)
dimension pos (maxpts), fit (maxpts)
dimension sunpos (maxpts), refspec (maxpts, 10)
dimension database (11, maxpts)
dimension sunpos_ss (maxpts), sunspec (maxpts)
dimension sunspec_ss (maxpts)
dimension suntemp (maxpts), reftemp (maxpts, 10)
dimension d2sun (maxpts)
dimension var (mmax)
real*8 kppos (maxkpno), kppos_ss (maxkpno), kpspec (maxkpno), kpspec_gauss (maxkpno)
character*120 kpnospec
logical wrt_scr, smooth, doas
logical first_sun /.true./
save

if (ntype .eq. 1) then
  if (first_sun) then
    first_sun = .false.
    if (wrt_scr) write (*, '(5x,a)') 'wavelength reference spectrum.'
    read (*, '(a)') kpnospec
    open (unit = 25, file = kpnospec, status='old')
    i = 1
10  read (25, *, end = 20) kppos (i), kpspec (i)
      kpspec (i) = kpspec (i) !/ 1.d14
      i = i + 1
      go to 10
20  close (unit = 25)
    nsun = i - 1
  end if ! on first_sun
else if (ntype .eq. 3) then
! The following kluge keeps me from extra re-writing.
  nsun = npoints
  do i = 1, nsun
    sunpos (i) = pos (i)
    sunspec (i) = database (1, i)
    do j = 1, 10
      refspec (i, j) = database (j + 1, i)
    end do
  end do
! jmax used below for reference spectra calculations. It is the number of
! reference spectra.
  jmax = (npars - 12) / 4
  if (smooth) then
    do i = 1, nsun
      suntemp (i) = sunspec (i)
      do j = 1, 10
        reftemp (i, j) = refspec (i, j)
      end do
    end do
    do i = 3, npoints - 2
      sunspec (i) = 0.375d0 * suntemp (i) + 0.25d0 * (suntemp (i + 1) + &
        suntemp (i - 1)) + 0.0625d0 * (suntemp (i + 2) + suntemp (i - 2))
      do j = 1, 10
        refspec (i, j) = 0.375d0 * reftemp (i, j) + 0.25d0 * &
        (reftemp (i + 1, j) + reftemp (i - 1, j)) + 0.0625d0 * &
        (reftemp (i + 2, j) + reftemp (i - 2, j))
      end do
    end do
  end if ! on smoothing.
  return ! ntype = 3, database loaded.
end if ! on ntype.

! Calculate the spectrum: First do the shift and squeeze. Shift by var (11),
! squeeze by 1 + var (12); do in absolute sense, to make it easy to back-convert
! radiance data.
do i = 1, nsun
  if (ntype .eq. 1 .or. ntype .eq. 2) then
    kppos_ss (i) = kppos (i) * (1.d0 + var (12)) + var (11)
  else
    sunpos_ss (i) = sunpos (i) * (1.d0 + var (12)) + var (11)
  end if
end do

! Broadening and re-sampling of solar spectrum.
if (ntype .eq. 1 .or. ntype .eq. 2) then
! Case for wavelength fitting of irradiance and radiance.
! Broaden the solar reference by the hw1e value.
  hw1e = var (10)
  call gauss (kppos_ss, kpspec, kpspec_gauss, nsun, hw1e)
! Re-sample the solar reference spectrum to the radiance grid.
  call spline (kppos_ss, kpspec_gauss, nsun, d2sun)
  call splint (kppos_ss, kpspec_gauss, d2sun, nsun, npoints, pos, sunspec_ss)
else
! Re-sample the solar reference spectrum to the radiance grid.
  call spline (sunpos_ss, sunspec, nsun, d2sun)
  call splint (sunpos_ss, sunspec, d2sun, nsun, npoints, pos, sunspec_ss)
end if

! Add up the contributions, with solar intensity as var (9), trace species
! beginning at var (13), to include possible linear and Beer's law forms. Do
! these as linear-Beer's-linear. In order to do DOAS I need to be careful to
! include just linear contributions, which I already high-pass filtered.

! DOAS here - the spectrum to be fitted needs to be re-defined.
if (doas .and. ntype .eq. 4) then
  do i = 1, npoints
!   For DOAS, var (9) should == 1., and not be varied.
    fit (i) = var (9) * dlog (sunspec_ss (i))
!   Ring adjustment.
    fit (i) = fit (i) + var (13) * (refspec (i, 1) / sunspec_ss (i))
    do j = 2, jmax
      fit (i) = fit (i) + var (4 * j + 9) * refspec (i, j)
    end do
  end do
else
! Doing BOAS.
  do i = 1, npoints
    fit (i) = var (9) * sunspec_ss (i)
!   Initial add-on contributions.
    do j = 1, jmax
      fit (i) = fit (i) + var (4 * j + 9) * refspec (i, j)
    end do
!   Beer's law contributions.
    do j = 1, jmax
      fit (i) = fit (i) * dexp (- var (4 * j + 10) * refspec (i, j))
    end do
!   Final add-on contributions.
    do j = 1, jmax
      fit (i) = fit (i) + var (4 * j + 11) * refspec (i, j)
    end do
  end do
end if

! Add the scaling.
do i = 1, npoints
  del = pos (i) - avg
  fit (i) = fit (i) * (var (5) + del * var (6) + (del**2) * var (7) + &
  (del**3) * var (8))
end do

! Add baseline parameters.
do i = 1, npoints
  del = pos (i) - avg
  fit (i) = fit (i) + var (1) + del * var (2) + (del**2) * var (3) + &
  (del**3) * var (4)
end do

return
end subroutine spectrum
!
subroutine mrqmin (x, y, sig, ndata, a, ma, lista, mfit, covar, alpha, chisq, &
  funcs, alamda, delchi, smooth, avg, diff, wrt_scr, ntype, database, doas)

! Levenberg-Marquardt minimization adapted from Numerical Recipes.

implicit real*8 (a - h, o - z)
external funcs
parameter (mmax = 64)
parameter (maxpts = 2000)
dimension x (maxpts), y (maxpts), sig (maxpts), a (mmax), lista (mmax), &
  covar (mmax, mmax), alpha (mmax, mmax), atry (mmax), beta (mmax), da (mmax), &
  diff (mmax), index (mmax), database (11, maxpts)
logical wrt_scr, smooth, doas, first_call /.true./
real*8 identity (mmax, mmax), inverse (mmax, mmax)
save

if (first_call) then
  first_call = .false.
  identity = 0.d0
  do j = 1, mmax
    identity (j, j) = 1.d0
  end do
end if
if (alamda .lt. 0.) then
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
  call mrqcof (x, y, sig, ndata, a, ma, lista, mfit, alpha, beta, chisq, &
    smooth, avg, diff, wrt_scr, ntype, database, doas)
  ochisq = chisq
  do j = 1, ma
    atry (j) = a (j)
  end do
end if
do j = 1, mfit
  do k = 1, mfit
    covar (j, k) = alpha (j, k)
  end do
  covar (j, j) = alpha (j, j) * (1. + alamda)
  da (j) = beta (j)
end do
call ludcmp (covar, mfit, mmax, index, even)
call lubksb (covar, mfit, mmax, index, da)

if (alamda .eq. 0.) then
! Calculate inverse of curvature matrix to obtain covariance matrix.
  inverse = identity
  do j = 1, mfit
!   call lubksb (covar, mfit, mmax, index, inverse (1, j))
    call lubksb (covar, mfit, mmax, index, inverse (:, j))
!   Luke Valin correction 31mar2014; answers are exactly the same
  end do
  covar = inverse
  return
end if
do j = 1, mfit
  atry (lista (j)) = a (lista (j)) + da (j)
end do
call mrqcof (x, y, sig, ndata, atry, ma, lista, mfit, covar, da, chisq, &
  smooth, avg, diff, wrt_scr, ntype, database, doas)
if (chisq .lt. ochisq) then
  alamda = 0.1 * alamda
  ochisq = chisq
  do j = 1, mfit
    do k = 1, mfit
      alpha (j, k) = covar (j, k)
    end do
    beta (j) = da (j)
    a (lista (j)) = atry (lista (j))
  end do
else
  alamda = 10. * alamda
  chisq = ochisq
end if

return
end subroutine mrqmin
!
subroutine mrqcof (x, y, sig, ndata, a, ma, lista, mfit, alpha, beta, chisq, &
  smooth, avg, diff, wrt_scr, ntype, database, doas)

! Calculates chi-squared gradient and curvature matrices for the Levenberg-
! Marquardt minimization. From Numerical Recipes.

implicit real * 8 (a - h, o - z)
parameter (mmax = 64)
parameter (maxpts = 2000)
dimension y (maxpts), x (maxpts), sig (maxpts), alpha (mmax, mmax), &
  beta (mmax), dyda (mmax, maxpts), lista (mmax), a (mmax), ymod (maxpts), &
  diff (mmax)
logical wrt_scr, smooth, doas
dimension database (11, maxpts)
save
do j = 1, mfit
  do k = 1, j
    alpha (j, k) = 0.
  end do
  beta (j) = 0.
end do
chisq = 0.
print*, ntype
call funcs (x, ndata, a, ymod, dyda, ma, lista, mfit, smooth, avg, diff, &
  wrt_scr, ntype, database, doas)
print*, 'After funcs'
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
!
subroutine funcs (pos, npoints, vars, ymod, dyda, npars, lista, nvaried, &
  smooth, avg, diff, wrt_scr, ntype, database, doas)

! Calculates the spectrum and its derivatives for mrqcof.

implicit real*8 (a - h, o - z)
parameter (mmax = 64)
parameter (maxpts = 7000)
dimension lista (mmax)
dimension pos (maxpts), ymod (maxpts), dyplus (maxpts)
dimension diff (mmax), dyda (mmax, maxpts), vars (mmax)
logical wrt_scr, smooth, doas
dimension database (11, maxpts)
save

! Calculate the spectrum.
call spectrum (npoints, npars, smooth, avg, pos, ymod, vars, wrt_scr, ntype, &
database, doas)

! Calculate the derivatives by finite differences. dyplus is y (var + finite
! difference in variable)
do i = 1, nvaried
  var0 = vars (lista (i))
  vars (lista (i)) = vars (lista (i)) + diff (lista (i))
  call spectrum (npoints, npars, smooth, avg, pos, dyplus, vars, wrt_scr, &
    ntype, database, doas)
  do j = 1, npoints
    dyda (lista (i), j) = (dyplus (j) - ymod (j)) / diff (lista (i))
  end do
  vars (lista (i)) = var0
end do

return
end
!
SUBROUTINE write_input (npoints, npars, n_solar_pars, nvaried, nvar_sun, &
  iterate_rad, iterate_sun, weight_rad, weight_sun, if_varied, if_var_sun, &
  initial, init_sun, nfirstfit, update_pars, provar, autodiff, ll_rad, lu_rad, &
  ll_sun, lu_sun, phase, wrt_scr, smooth, renorm, szamax, szamin, &
  cldmax, latmax, latmin, doas)

IMPLICIT NONE

INTEGER*4, PARAMETER :: mmax = 64
INTEGER*4, PARAMETER :: maxpts = 7000

LOGICAL, INTENT(IN) :: iterate_sun, iterate_rad, weight_sun, weight_rad, autodiff, renorm, &
  update_pars, wrt_scr, smooth, doas
LOGICAL, DIMENSION(mmax), INTENT(IN) :: if_varied, if_var_sun
REAL*8, DIMENSION(mmax), INTENT(IN) :: initial, init_sun
REAL*8, INTENT(IN) :: latmax, latmin, szamax, szamin, phase, provar
INTEGER*4, INTENT(IN) :: cldmax, ll_sun, lu_sun, ll_rad, lu_rad, nfirstfit, nvar_sun, &
     nvaried, npoints, npars, n_solar_pars

INTEGER*4 :: i
save

! Summarizes the inputs for later reference.
if (doas) then
  write (22, '(1x, t2, a)') 'DOAS fitting'
else
  write (22, '(1x, t2, a)') 'BOAS fitting'
end if
write (22, '(1x, t2, a, t33, i3, t36, a)') 'maximum number of parameters =', &
  mmax, '.'
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
  write (22, *) 'with 0.001 times initial variable.'
end if
write (22, '(1x, t2, a, t51, 1pe9.3)') &
  'convergence for relative parameter change set at ', provar
write (22, *) 'times the differentiation parameter for each variable.'
write (22, *) ' renorm = ', renorm
write (22, *) ' update_pars = ', update_pars
write (22, *) ' wrt_scr = ', wrt_scr
write (22, *) ' smooth = ', smooth
write (22, '(a, f7.2)') '  szamax =', szamax
write (22, '(a, f7.2)') '  szamin =', szamin
write (22, '(a, I3)') '  cldmax =', cldmax
write (22, '(a, f7.2)') '  latmax =', latmax
write (22, '(a, f7.2)') '  latmin =', latmin
write (22, '(a, i6)') '  nfirst =', nfirstfit
write (22, '(a, f8.4)') '  phase =', phase

write (22, *)
write (22,'(1x, a)') 'Solar Fitting'
write (22,'(1x, a, /)') '                              baseline parameters'
write (22,'(1x, a)') '    parameter     initial value  varied'
write (22,'(1x,t1,a,t20,1pe14.7,t37,l)') ' baseline offset', init_sun (1), &
  if_var_sun (1)
write (22,'(1x,t1,a,t20,1pe14.7,t37,l)') ' baseline tilt', init_sun (2), &
  if_var_sun (2)
write (22,'(1x,t1,a,t20,1pe14.7,t37,l)') ' baseline quad.', init_sun (3), &
  if_var_sun (3)
write (22,'(1x,t1,a,t20,1pe14.7,t37,l, /)') ' baseline cubic', init_sun (4), &
  if_var_sun (4)
write (22,'(1x, a, /)') '                              scaling parameters'
write (22,'(1x, a)') '    parameter     initial value  varied ?'
write (22,'(1x,t1,a,t20,1pe14.7,t37,l)') ' scaling factor', init_sun (5), &
  if_var_sun (5)
write (22,'(1x,t1,a,t20,1pe14.7,t37,l)') ' scaling tilt', init_sun (6), &
  if_var_sun (6)
write (22,'(1x,t1,a,t20,1pe14.7,t37,l)') ' scaling quad.', init_sun (7), &
  if_var_sun (7)
write (22,'(1x,t1,a,t20,1pe14.7,t37,l, /)') ' scaling cubic', init_sun (8), &
  if_var_sun (8)

write (22,'(1x, a, /)') '                             spectral parameters'
write (22,'(1x, a)') '    parameter     initial value  varied ?'
do i = 9, n_solar_pars, 4
  write (22,'(1x,t1,a,t20,1pe14.7,t37,l)') ' parameter', init_sun (i), &
    if_var_sun (i)
  write (22,'(1x,t1,a,t20,1pe14.7,t37,l)') ' parameter', init_sun (i + 1), &
    if_var_sun (i + 1)
  write (22,'(1x,t1,a,t20,1pe14.7,t37,l)') ' parameter', init_sun (i + 2), &
    if_var_sun (i + 2)
  write (22,'(1x,t1,a,t20,1pe14.7,t37,l)') ' parameter', init_sun (i + 3), &
    if_var_sun (i + 3)
end do

write (22, *)
write (22,'(1x, a, /)') 'Radiance Fitting'
write (22,'(1x, a, /)') '                              baseline parameters'
write (22,'(1x, a)') '    parameter     initial value  varied'
write (22,'(1x,t1,a,t20,1pe14.7,t37,l)') ' baseline offset', initial (1), &
  if_varied (1)
write (22,'(1x,t1,a,t20,1pe14.7,t37,l)') ' baseline tilt', initial (2), &
  if_varied (2)
write (22,'(1x,t1,a,t20,1pe14.7,t37,l)') ' baseline quad.', initial (3), &
  if_varied (3)
write (22,'(1x,t1,a,t20,1pe14.7,t37,l, /)') ' baseline cubic', initial (4), &
  if_varied (4)
write (22,'(1x, a, /)') '                              scaling parameters'
write (22,'(1x, a)') '    parameter     initial value  varied ?'
write (22,'(1x,t1,a,t20,1pe14.7,t37,l)') ' scaling factor', initial (5), &
  if_varied (5)
write (22,'(1x,t1,a,t20,1pe14.7,t37,l)') ' scaling tilt', initial (6), &
  if_varied (6)
write (22,'(1x,t1,a,t20,1pe14.7,t37,l)') ' scaling quad.', initial (7), &
  if_varied (7)
write (22,'(1x,t1,a,t20,1pe14.7,t37,l)') ' scaling cubic', initial (8), &
  if_varied (8)

write (22, *)
write (22,'(1x, a, /)') '                             spectral parameters'
write (22,'(1x, a)') '    parameter     initial value  varied ?'
do i = 9, npars, 4
  write (22,'(1x,t1,a,t20,1pe14.7,t37,l)') ' parameter', initial (i), &
    if_varied (i)
  write (22,'(1x,t1,a,t20,1pe14.7,t37,l)') ' parameter', initial (i + 1), &
    if_varied (i + 1)
  write (22,'(1x,t1,a,t20,1pe14.7,t37,l)') ' parameter', initial (i + 2), &
    if_varied (i + 2)
  write (22,'(1x,t1,a,t20,1pe14.7,t37,l,/)') ' parameter', initial (i + 3), &
    if_varied (i + 3)
end do
write (22, *)

return
end subroutine write_input
!
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
!
subroutine gauss (pos, spec, specmod, npoints, hw1e)

! Convolves input spectrum with Gaussian slit function of specified hw1e
! (half-width at 1/e intensity). For evenly-spaced spectra.

implicit real*8 (a - h, o - z)
parameter (maxkpno = 7000)
dimension pos (maxkpno), spec (maxkpno)
dimension specmod (maxkpno), slit (maxkpno)
save

if (hw1e .eq. 0) then
  write (*, *) ' no gaussian convolution applied.'
  specmod = spec
  return
end if

emult = - 1.d0 / (hw1e**2)
delpos = pos (2) - pos (1)

! Apply slit function convolution.
! Calculate slit function values out to 0.001 times x0 value, normalize so that
! sum = 1.
slitsum = 1.d0
slit0 = 1.d0
do i = 1, 1000
  slit (i) = exp (emult * (delpos * i)**2)
  slitsum = slitsum + 2.d0 * slit (i)
  if (slit (i) / slit0 .le. 0.001) then
    nslit = i
!   write (*, *) ' nslit = ', nslit
    go to 40
  end if
end do
40 continue
slit0 = slit0 / slitsum
  do i = 1, nslit
    slit (i) = slit (i) / slitsum
  end do

! Convolve spectrum, reflect at endpoints.
do i = 1, npoints
  specmod (i) = slit0 * spec (i)
  do j = 1, nslit
    nlo = i - j
    nhi = i + j
    if (nlo .le. 0) nlo = - nlo + 1
    if (nhi .gt. npoints) nhi = 2 * npoints - nhi
    specmod (i) = specmod (i) + slit (j) * (spec (nlo) + spec (nhi))
  end do
end do

! This legacy code corrects for nlo = 0. Not currently needed, but might be with
! another compiler.
! do i = 1, npoints
!   specmod (i) = slit0 * spec (i)
!   do j = 1, nslit
!     nlo = i - j
!     nhi = i + j
!     if (nlo .lt. 0) nlo = - nlo
!     if (nhi .gt. npoints) nhi = 2 * npoints - nhi
!     if (nlo .ne. 0) then
!       specmod (i) = specmod (i) + slit (j) * (spec (nlo) + spec (nhi))
!     else
!       specmod (i) = specmod (i) + slit (j) * spec (nhi)
!     end if
!   end do
! end do

return
end subroutine gauss
!
subroutine gauss_uneven (pos, spec, specmod, npoints, hw1e)

! Convolves input spectrum with Gaussian slit function of specified hw1e
! (half-width at 1/e intensity). For unevenly-spaced spectra. Does not treat
! spectral end-points carefully: we usually have at least several nanometers of
! margin so this is not an issue.

implicit real*8 (a - h, o - z)
parameter (maxkpno = 7000)
dimension pos (maxkpno), spec (maxkpno), specmod (maxkpno), temp (maxkpno)
specmod = 0.d0
temp = 0.d0

if (hw1e .eq. 0) then
  write (*, *) ' no gaussian convolution applied.'
  specmod = spec
  return
end if

emult = - 1.d0 / (hw1e**2)

do i = 1, npoints
  sum = 0.d0
  do j = 1, npoints
    slit = dexp (emult * (pos (j) - pos (i))**2)
    temp (j) = spec (i) * slit
    sum = sum + slit
  end do
  specmod = specmod + temp / sum
end do

return
end subroutine gauss_uneven
!
subroutine spline(x, y, n, y2)

! Cubic spline derivative calculation from Numerical Recipes. Modified to always
! use "natural" boundary conditions (2nd derivatives = 0 at boundaries).

implicit real*8 (a - h, o - z)
parameter (maxkpno = 7000)
dimension x (n), y (n), y2 (n), u (maxkpno)
save

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
!
subroutine splint (xa, ya, y2a, n, np, x, y)

! Cubic spline interpolation from Numerical Recipes, using results of
! subroutine spline.

implicit real*8 (a - h, o - z)
parameter (maxpts = 2000)
dimension x (maxpts), y (maxpts)
dimension xa (n), ya (n), y2a (n)
save

do i = 1, np
  klo = 1
  khi = n
1 if (khi - klo .gt. 1) then
    k = (khi + klo) / 2
    if (xa (k) .gt. x (i)) then
      khi = k
    else
      klo = k
    end if
    go to 1
  end if
  h = xa (khi) - xa (klo)
  if (h .eq. 0.) print*, 'bad xa input.'
  a = (xa (khi) - x (i)) / h
  b = (x (i) - xa (klo)) / h
  y (i) = a * ya (klo) + b * ya (khi) + ((a**3 - a) * y2a (klo) + (b**3 - b) * &
    y2a (khi)) * (h**2) / 6.
end do

return
end subroutine splint
!
subroutine undersample (wav, nwav, underspec, hw1e, fraction, wrt_scr)

! Convolves input spectrum with Gaussian slit function of specified HW1e (half-
! width at 1/e of maximum intensity), and samples at a particular input phase to
! give the undersampling spectrum, that is, the correction for not Nyquist
! sampling the spectra (see Applied Optics 44, 1296-1304, 2005 for more detail).
! This version calculates both phases of the undersampling spectrum, phase1 -
! i.e., underspec (1, i) - being the more common in radiance spectra. "KPNO"
! data means the high spectral resolution solar reference spectrum (now, most
! usually, a portion of the SAO2010 Solar Irradiance Reference Spectrum,
! available at http://www.cfa.harvard.edu/atmosphere/).

implicit real*8 (a - h, o - z)
parameter (maxpts = 7000)
parameter (maxkpno = 7000)
character*120 overspec
logical wrt_scr
dimension pos (maxkpno), spec (maxkpno), specmod (maxkpno), d2spec (maxkpno)
dimension slit (maxpts), wav (maxpts), wav2 (maxpts), over (maxpts), &
  under (maxpts), resample (maxpts), d2res (maxpts), underspec (2, maxpts)
save

if (wrt_scr) write (*, '(a)') 'enter standard solar spectrum file.'
read (*, '(a)') overspec
open (unit = 27, file = overspec, status = 'old')
! read in the KPNO data
i = 1
20 read (27, *, end = 30) pos (i), spec (i)
  spec (i) = spec (i) / 1.d14
  i = i + 1
  go to 20
30 npoints = i - 1
delpos = pos (2) - pos (1)
close (unit = 27)

! Apply slit function convolution. Calculate slit function values out to 0.001
! times x0 value, normalize so that sum = 1.
emult = - 1.d0 / (hw1e**2)

slitsum = 1.d0
slit0 = 1.d0
do i = 1, maxpts
  slit (i) = exp (emult * (delpos * i)**2)
  slitsum = slitsum + 2.d0 * slit (i)
  if (slit (i) / slit0 .le. 0.001) then
    nslit = i
    go to 40
  end if
end do
40 continue
slit0 = slit0 / slitsum
do i = 1, nslit
  slit (i) = slit (i) / slitsum
end do

! Convolve spectrum, reflect at endpoints.
do i = 1, npoints
  specmod (i) = slit0 * spec (i)
  do j = 1, nslit
    nlo = i - j
    nhi = i + j
    if (nlo .lt. 0) nlo = - nlo
    if (nhi .gt. npoints) nhi = 2 * npoints - nhi
    specmod (i) = specmod (i) + slit (j) * (spec (nlo) + spec (nhi))
  end do
end do

! Legacy code.
! do i = 1, npoints
!   specmod (i) = slit0 * spec (i)
!   do j = 1, nslit
!     nlo = i - j
!     nhi = i + j
!     if (nlo .lt. 0) nlo = - nlo
!     if (nhi .gt. npoints) nhi = 2 * npoints - nhi
!     if (nlo .ne. 0) then
!       specmod (i) = specmod (i) + slit (j) * (spec (nlo) + spec (nhi))
!     else
!       specmod (i) = specmod (i) + slit (j) * spec (nhi)
!     end if
!   end do
! end do

! Phase1 calculation: Calculate spline derivatives for KPNO data.
call spline (pos, specmod, npoints, d2spec)

! Calculate solar spectrum at radiance positions
call splint (pos, specmod, d2spec, npoints, nwav, wav, resample)

! Calculate spline derivatives for resampled data.
call spline (wav, resample, nwav, d2res)

! Calculate solar spectrum at radiance + phase positions, original and
! resampled.
wav2 (1) = wav (1)
do i = 2, nwav
  wav2 (i) = (1.d0 - fraction) * wav (i - 1) + fraction * wav (i)
end do
call splint (pos, specmod, d2spec, npoints, nwav, wav2, over)
call splint (wav, resample, d2res, nwav, nwav, wav2, under)
over (1) = 0.
do i = 1, nwav
  underspec (1, i) = over (i) - under (i)
  resample (i) = over (i)
end do

! Phase2 calculation: Calculate spline derivatives for KPNO data.
call spline (pos, specmod, npoints, d2spec)

! Calculate spline derivatives for resampled data.
call spline (wav2, resample, nwav, d2res)

! Calculate solar spectrum at radiance positions, original and resampled.
call splint (pos, specmod, d2spec, npoints, nwav, wav, over)
call splint (wav2, resample, d2res, nwav, nwav, wav, under)
do i = 1, nwav
  underspec (2, i) = over (i) - under (i)
end do

return
end subroutine undersample
!
subroutine dataspline (wav, nwav, database, wrt_scr)

! Opens reference spectrum files.
! Convolves reference spectra with supplied slit width if requested.
! Samples reference spectra to standard grid determined from wavelength
! calibration of a radiance spectrum.
! Loads spectra into the database file.

implicit real*8 (a - h, o - z)
parameter (maxpts = 7000)
parameter (maxkpno = 7000)

! Input spectra.
dimension xring (maxkpno), yring (maxkpno), deriv2ring (maxkpno)
dimension xgas1 (maxkpno), ygas1 (maxkpno), deriv2gas1 (maxkpno)
dimension xgas2 (maxkpno), ygas2 (maxkpno), deriv2gas2 (maxkpno)
dimension xgas3 (maxkpno), ygas3 (maxkpno), deriv2gas3 (maxkpno)
dimension xgas4 (maxkpno), ygas4 (maxkpno), deriv2gas4 (maxkpno)
dimension xgas5 (maxkpno), ygas5 (maxkpno), deriv2gas5 (maxkpno)
dimension xgas6 (maxkpno), ygas6 (maxkpno), deriv2gas6 (maxkpno)
dimension xcom (maxkpno), ycom (maxkpno), deriv2com (maxkpno)
dimension xtemp (maxkpno), ytemp (maxkpno)

! Spectrum arrays.
dimension wav (maxpts), database (11, maxpts)
dimension temp1 (maxpts), temp2 (maxpts), temp3 (maxpts), temp4 (maxpts), &
  temp5 (maxpts), temp6 (maxpts), temp7 (maxpts), temp8 (maxpts)
character*120 :: spec1, spec2, spec3, spec4, spec5, spec6, spec7, spec8

! Convolution switches.
logical :: if_conv_ring, if_conv_gas1, if_conv_gas2, if_conv_gas3, if_conv_gas4, &
  if_conv_gas5, if_conv_gas6, if_conv_com

! Verbose mode.
logical wrt_scr

! Edlen 1966 standard air parameters. Gas6 is normally O2-O2, which used to
! require air-vacuum correction. Parameters saved here just in case.
parameter (c1 = 1.0000834213d0)
parameter (c2 = 2.406030d-2)
parameter (c3 = 130.d0)
parameter (c4 = 1.5997d-4)
parameter (c5 = 38.9d0)

save

! Open reference spectrum files.
if (wrt_scr) write (*, '(5x, a)') 'enter input Ring spectrum'
read (*, '(a)') spec1
read (*, *) ring_mult, if_conv_ring, hw1e_ring
if (wrt_scr) write (*, '(5x, a)') 'enter gas1 spectrum'
read (*, '(a)') spec2
read (*, *) gas1_mult, if_conv_gas1, hw1e_gas1
if (wrt_scr) write (*, '(5x, a)') 'enter gas2 spectrum'
read (*, '(a)') spec3
read (*, *) gas2_mult, if_conv_gas2, hw1e_gas2
if (wrt_scr) write (*, '(5x, a)') 'enter gas3 spectrum'
read (*, '(a)') spec4
read (*, *) gas3_mult, if_conv_gas3, hw1e_gas3
if (wrt_scr) write (*, '(5x, a)') 'enter gas4 spectrum'
read (*, '(a)') spec5
read (*, *) gas4_mult, if_conv_gas4, hw1e_gas4
if (wrt_scr) write (*, '(5x, a)') 'enter gas5 spectrum'
read (*, '(a)') spec6
read (*, *) gas5_mult, if_conv_gas5, hw1e_gas5
if (wrt_scr) write (*, '(5x, a)') 'enter gas6 spectrum'
read (*, '(a)') spec7
read (*, *) gas6_mult, if_conv_gas6, hw1e_gas6
if (wrt_scr) write (*, '(5x, a)') 'enter common mode spectrum'
read (*, '(a)') spec8
read (*, *) com_mult, if_conv_com, hw1e_com

open (unit = 11, file = spec1, status = 'old')
open (unit = 12, file = spec2, status = 'old')
open (unit = 13, file = spec3, status = 'old')
open (unit = 14, file = spec4, status = 'old')
open (unit = 15, file = spec5, status = 'old')
open (unit = 16, file = spec6, status = 'old')
open (unit = 17, file = spec7, status = 'old')
open (unit = 18, file = spec8, status = 'old')

! Convolve and spline Ring spectrum.
i = 1
20 read (11, *, end = 30) xring (i), yring (i)
  yring (i) = yring (i) * ring_mult
  i = i + 1
  go to 20
30 npoints = i - 1
if (if_conv_ring) then
  call gauss_uneven (xring, yring, ytemp, npoints, hw1e_ring)
  yring = ytemp
end if
if (xring (npoints) .lt. xring (1)) then
  do i = 1, npoints
    xtemp (i) = xring (i)
    ytemp (i) = yring (i)
  end do
  do i = 1, npoints
    xring (i) = xtemp (npoints + 1 - i)
    yring (i) = ytemp (npoints + 1 - i)
  end do
end if
call spline (xring, yring, npoints, deriv2ring)
call splint (xring, yring, deriv2ring, npoints, nwav, wav, temp1)

! Convolve and spline gas1 spectrum.
i = 1
40 read (12, *, end = 50) xgas1 (i), ygas1 (i)
  ygas1 (i) = ygas1 (i) * gas1_mult
  i = i + 1
  go to 40
50 npoints = i - 1
if (if_conv_gas1) then
  call gauss_uneven (xgas1, ygas1, ytemp, npoints, hw1e_gas1)
  ygas1 = ytemp
end if
if (xgas1 (npoints) .lt. xgas1 (1)) then
  do i = 1, npoints
    xtemp (i) = xgas1 (i)
    ytemp (i) = ygas1 (i)
  end do
  do i = 1, npoints
    xgas1 (i) = xtemp (npoints + 1 - i)
    ygas1 (i) = ytemp (npoints + 1 - i)
  end do
end if
call spline (xgas1, ygas1, npoints, deriv2gas1)
call splint (xgas1, ygas1, deriv2gas1, npoints, nwav, wav, temp2)

! Convolve and spline gas2 spectrum.
i = 1
60 read (13, *, end = 70) xgas2 (i), ygas2 (i)
  ygas2 (i) = ygas2 (i) * gas2_mult
  i = i + 1
  go to 60
70 npoints = i - 1
if (if_conv_gas2) then
  call gauss_uneven (xgas2, ygas2, ytemp, npoints, hw1e_gas2)
  ygas2 = ytemp
end if
if (xgas2 (npoints) .lt. xgas2 (1)) then
  do i = 1, npoints
    xtemp (i) = xgas2 (i)
    ytemp (i) = ygas2 (i)
  end do
  do i = 1, npoints
    xgas2 (i) = xtemp (npoints + 1 - i)
    ygas2 (i) = ytemp (npoints + 1 - i)
  end do
end if
call spline (xgas2, ygas2, npoints, deriv2gas2)
call splint (xgas2, ygas2, deriv2gas2, npoints, nwav, wav, temp3)

! Convolve and spline gas3 spectrum.
i = 1
80 read (14, *, end = 90) xgas3 (i), ygas3 (i)
  ygas3 (i) = ygas3 (i) * gas3_mult
  i = i + 1
  go to 80
90 npoints = i - 1
if (if_conv_gas3) then
  call gauss_uneven (xgas3, ygas3, ytemp, npoints, hw1e_gas3)
  ygas3 = ytemp
end if
if (xgas3 (npoints) .lt. xgas3 (1)) then
  do i = 1, npoints
    xtemp (i) = xgas3 (i)
    ytemp (i) = ygas3 (i)
  end do
  do i = 1, npoints
    xgas3 (i) = xtemp (npoints + 1 - i)
    ygas3 (i) = ytemp (npoints + 1 - i)
  end do
end if
call spline (xgas3, ygas3, npoints, deriv2gas3)
call splint (xgas3, ygas3, deriv2gas3, npoints, nwav, wav, temp4)

! Convolve and spline gas4 spectrum.
i = 1
100 read (15, *, end = 110) xgas4 (i), ygas4 (i)
  ygas4 (i) = ygas4 (i) * gas4_mult
  i = i + 1
  go to 100
110 npoints = i - 1
if (if_conv_gas4) then
  call gauss_uneven (xgas4, ygas4, ytemp, npoints, hw1e_gas4)
  ygas4 = ytemp
end if
if (xgas4 (npoints) .lt. xgas4 (1)) then
  do i = 1, npoints
    xtemp (i) = xgas4 (i)
    ytemp (i) = ygas4 (i)
  end do
  do i = 1, npoints
    xgas4 (i) = xtemp (npoints + 1 - i)
    ygas4 (i) = ytemp (npoints + 1 - i)
  end do
end if
call spline (xgas4, ygas4, npoints, deriv2gas4)
call splint (xgas4, ygas4, deriv2gas4, npoints, nwav, wav, temp5)

! Convolve and spline gas5 spectrum.
i = 1
120 read (16, *, end = 130) xgas5 (i), ygas5 (i)
  ygas5 (i) = ygas5 (i) * gas5_mult
  i = i + 1
  go to 120
130 npoints = i - 1
if (if_conv_gas5) then
  call gauss_uneven (xgas5, ygas5, ytemp, npoints, hw1e_gas5)
  ygas5 = ytemp
end if
if (xgas5 (npoints) .lt. xgas5 (1)) then
  do i = 1, npoints
    xtemp (i) = xgas5 (i)
    ytemp (i) = ygas5 (i)
  end do
  do i = 1, npoints
    xgas5 (i) = xtemp (npoints + 1 - i)
    ygas5 (i) = ytemp (npoints + 1 - i)
  end do
end if
call spline (xgas5, ygas5, npoints, deriv2gas5)
call splint (xgas5, ygas5, deriv2gas5, npoints, nwav, wav, temp6)

! Convolve and spline gas6 spectrum. Previously multiplied by 1.d40; another
! 1.d5 here.
i = 1
140 read (17, *, end = 150) xgas6 (i), ygas6 (i)
  ygas6 (i) = ygas6 (i) * gas6_mult
!! correct to vacuum wavelength
!  sigsq = 1.d6 / (xgas6 (i) * xgas6 (i))
!  rindex = c1 + c2 /(c3 - sigsq) + c4 / (c5 - sigsq)
!  wvac = xgas6 (i) * rindex
!! iterate once
!  sigsq = 1.d6 / (wvac * wvac)
!  rindex = c1 + c2 /(c3 - sigsq) + c4 / (c5 - sigsq)
!  xgas6 (i) = xgas6 (i) * rindex
  i = i + 1
  go to 140
150 npoints = i - 1
if (if_conv_gas6) then
  call gauss_uneven (xgas6, ygas6, ytemp, npoints, hw1e_gas6)
  ygas6 = ytemp
end if
if (xgas6 (npoints) .lt. xgas6 (1)) then
  do i = 1, npoints
    xtemp (i) = xgas6 (i)
    ytemp (i) = ygas6 (i)
  end do
  do i = 1, npoints
    xgas6 (i) = xtemp (npoints + 1 - i)
    ygas6 (i) = ytemp (npoints + 1 - i)
  end do
end if
call spline (xgas6, ygas6, npoints, deriv2gas6)
call splint (xgas6, ygas6, deriv2gas6, npoints, nwav, wav, temp7)

! Convolve and spline common mode spectrum if it is far off in wavelength grid.
! If it is not far off, use the original grid.
i = 1
160 read (18, *, end = 170) xcom (i), ycom (i)
  ycom (i) = ycom (i) * com_mult
  i = i + 1
  go to 160
170 npoints = i - 1
if (if_conv_com) then
  call gauss_uneven (xcom, ycom, ytemp, npoints, hw1e_com)
  ycom = ytemp
end if
if (abs (xcom (1) - wav (1)) .le. 0.02 .and. npoints .eq. nwav) then
  do i = 1, npoints
    temp8 (i) = ycom (i)
  end do
else
  if (xcom (npoints) .lt. xcom (1)) then
    do i = 1, npoints
      xtemp (i) = xcom (i)
      ytemp (i) = ycom (i)
    end do
    do i = 1, npoints
      xcom (i) = xtemp (npoints + 1 - i)
      ycom (i) = ytemp (npoints + 1 - i)
    end do
  end if
  call spline (xcom, ycom, npoints, deriv2com)
  call splint (xcom, ycom, deriv2com, npoints, nwav, wav, temp8)
end if

! Load results into the database file. Note ordering of spectra:
! 1 = Irradiance
! 2 = Ring
! 3 = gas1
! 4 = gas2
! 5 = gas3
! 6 = gas4
! 7 = gas5
! 8 = gas6
! 9 = Undersampling, phase 1
! 10 = Undersampling, phase 2
! 11 = Common mode
do i = 1, nwav
  database (2, i) = temp1 (i)
  database (3, i) = temp2 (i)
  database (4, i) = temp3 (i)
  database (5, i) = temp4 (i)
  database (6, i) = temp5 (i)
  database (7, i) = temp6 (i)
  database (8, i) = temp7 (i)
  database (11, i) = temp8 (i)
end do

close (unit = 11)
close (unit = 12)
close (unit = 13)
close (unit = 14)
close (unit = 15)
close (unit = 16)
close (unit = 17)
close (unit = 18)

return
end subroutine dataspline
!
subroutine subtract_cubic (pos, npoints, database, ll_rad, lu_rad)

! "High-pass filtering" for reference spectra, as required for DOAS.

implicit real*8 (a - h, o - z)
parameter (maxpts = 7000)
dimension pos (maxpts), database (11, maxpts), temp (maxpts), ptemp (maxpts)
dimension par (4), covar (4, 4), lista (4)
save

! Remember order of parameters (to satisfy lfit).
do i = 1, 4
  lista (i) = i
end do

! Find limits for polynomial fitting, with ~1 nm overlap.
do i = 1, npoints
  if (pos (i) .ge. pos (ll_rad) - 1.d0) then
    nlower = i
!   write (*, *) 'subtract_cubic lower:', i, pos (i), pos (ll_rad) - 1.
    go to 10
  end if
end do
10 continue
do i = nlower + 1, npoints
  if (pos (i) .gt. pos (lu_rad) + 1.d0) then
    nupper = i - 1
!   write (*, *) 'subtract_cubic upper:', i - 1, pos (i - 1), pos (lu_rad) + 1.
    go to 20
  end if
end do
20 continue
nfitted = nupper - nlower + 1

! Find average position over fitted region.
avg = 0.d0
do i = 1, nfitted
  avg = avg + pos (i + nlower - 1)
end do
avg = avg / nfitted

! Load temporary position file. Re-define positions in order to fit about mean
! position.
do i = 1, nfitted
  ptemp (i) = pos (i + nlower - 1) - avg
end do

! Load and fit database spectra nos. 2-11.
do i = 2, 11
  do j = 1, nfitted
    temp (j) = database (i, j + nlower - 1)
  end do
  call lfit (ptemp, temp, nfitted, par, lista, covar, chisq)
! Re-load database with high-pass filtered data, over whole spectral region.
  do j = 1, npoints
    x = pos (j) - avg
    database (i, j) = database (i, j) - (par (1) + par (2) * x + par (3) * &
    x**2 + par (4) * x**3)
  end do
end do

return
end subroutine subtract_cubic
!
subroutine subtract_cubic_meas (pos, npoints, spec, ll_rad, lu_rad)

! "High-pass filtering" for atmospheric spectra, as required for DOAS.

implicit real*8 (a - h, o - z)
parameter (maxpts = 7000)
dimension pos (maxpts), spec (maxpts), temp (maxpts), ptemp (maxpts)
dimension par (4), covar (4, 4), lista (4)
save

! Remember order of parameters (to satisfy lfit).
do i = 1, 4
  lista (i) = i
end do

! Find limits for polynomial fitting, with ~1 nm overlap.
do i = 1, npoints
  if (pos (i) .ge. pos (ll_rad) - 1.d0) then
    nlower = i
!   write (*, *) ' subtract_cubic_meas lower:', i, pos (i), pos (ll_rad) - 1.
    go to 10
  end if
end do
10 continue
do i = nlower + 1, npoints
  if (pos (i) .gt. pos (lu_rad) + 1.d0) then
    nupper = i - 1
!   write (*, *) ' subtract_cubic_meas upper:', i - 1, pos (i - 1), &
!     pos (lu_rad) + 1.
    go to 20
  end if
end do
20 continue
nfitted = nupper - nlower + 1

! Find average position over fitted region.
avg = 0.d0
do i = 1, nfitted
  avg = avg + pos (i + nlower - 1)
end do
avg = avg / nfitted

! Load temporary position file. Re-define positions in order to fit about mean
! position.
do i = 1, nfitted
  ptemp (i) = pos (i + nlower - 1) - avg
end do

! Load and fit spectrum
do i = 1, nfitted
  temp (i) = spec (i + nlower - 1)
end do
call lfit (ptemp, temp, nfitted, par, lista, covar, chisq)
! Re-load spec with high-pass filtered data, over whole spectral region.
do i = 1, npoints
  x = pos (i) - avg
  spec (i) = spec (i) - (par (1) + par (2) * x + par (3) * x**2 + par (4) * &
  x**3)
end do

return
end subroutine subtract_cubic_meas
!
subroutine lfit (x, y, ndata, a, lista, covar, chisq)

! Linear least-squares fitting, from Numerical Recipes.

implicit real*8 (a - h, o - z)
parameter (maxpts = 7000)
dimension x (maxpts), y (maxpts), a (4), lista (4), covar (4, 4), beta (4), &
  afunc (4), index (4)
save

! Dimension for the exact number of parameters:
ncvm = 4
mfit = 4
ma = 4

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
    print*, 'improper set in lista'
  end if
end do
if (kk .ne. (ma + 1)) print*, 'improper set in lista'
do j = 1, mfit
  do k = 1, mfit
    covar (j, k) = 0.
  end do
  beta (j) = 0.
end do
do i = 1,ndata
  afunc (1) = 1.
  afunc (2) = x (i)
  afunc (3) = x (i)**2
  afunc (4) = x (i)**3
  ym = y (i)
  if(mfit .lt. ma) then
    do j = mfit + 1, ma
      ym = ym - a (lista (j)) * afunc (lista (j))
    end do
  end if
  do j = 1, mfit
    wt = afunc (lista (j))
    do k = 1, j
      covar (j, k) = covar (j, k) + wt * afunc (lista (k))
    end do
    beta(j) = beta(j) + ym * wt
  end do
end do
if (mfit .gt. 1) then
  do j = 2, mfit
    do k = 1, j - 1
      covar (k, j) = covar (j, k)
    end do
  end do
end if
call ludcmp (covar, mfit, ncvm, index, even)
call lubksb (covar, mfit, ncvm, index, beta)
do j = 1, mfit
  a (lista (j)) = beta (j)
end do
chisq = 0.
do i = 1, ndata
  afunc (1) = 1.
  afunc (2) = x (i)
  afunc (3) = x (i)**2
  afunc (4) = x (i)**3
  sum = 0.
  do j = 1, ma
    sum = sum + a(j) * afunc (j)
  end do
  chisq = chisq + (y (i) - sum)**2
end do

return
end subroutine lfit
!
subroutine ludcmp (a, n, np, indx, d)

! LU decomposition, from Numerical Recipes.

implicit real*8 (a - h, o - z)
parameter (nmax = 100, tiny = 1.0d-20)
dimension a (np, np), indx (n), vv (nmax)
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
!
subroutine lubksb (a, n, np, indx, b)

! LU forward- and back-substitution, from Numerical Recipes.

implicit real*8 (a - h, o - z)
dimension a (np, np), indx (n), b (n)

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
