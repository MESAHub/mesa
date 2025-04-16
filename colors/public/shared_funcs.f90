
MODULE shared_func
  USE const_def, ONLY: dp
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: loadsed, dilute_flux

CONTAINS






     SUBROUTINE loadsed(directory, index, wavelengths, flux)
       CHARACTER(LEN=*), INTENT(IN) :: directory
       INTEGER, INTENT(IN) :: index
       REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: wavelengths, flux

       CHARACTER(LEN=512) :: line
       INTEGER :: unit, n_rows, status, i
       REAL(DP) :: temp_wavelength, temp_flux

       ! Open the file
       unit = 20
       OPEN(unit, FILE=TRIM(directory), STATUS='OLD', ACTION='READ', IOSTAT=status)
       IF (status /= 0) THEN
         PRINT *, "Error: Could not open file ", TRIM(directory)
         STOP
       END IF

       ! Skip header lines
       DO
         READ(unit, '(A)', IOSTAT=status) line
         IF (status /= 0) THEN
           PRINT *, "Error: Could not read the file", TRIM(directory)
           STOP
         END IF
         IF (line(1:1) /= "#") EXIT
       END DO

       ! Count rows in the file
       n_rows = 0
       DO
         READ(unit, '(A)', IOSTAT=status) line
         IF (status /= 0) EXIT
         n_rows = n_rows + 1
       END DO

       ! Allocate arrays
       ALLOCATE(wavelengths(n_rows))
       ALLOCATE(flux(n_rows))

       ! Rewind to the first non-comment line
       REWIND(unit)
       DO
         READ(unit, '(A)', IOSTAT=status) line
         IF (status /= 0) THEN
           PRINT *, "Error: Could not rewind file", TRIM(directory)
           STOP
         END IF
         IF (line(1:1) /= "#") EXIT
       END DO

       ! Read and parse data
       i = 0
       DO
         READ(unit, *, IOSTAT=status) temp_wavelength, temp_flux
         IF (status /= 0) EXIT
         i = i + 1
         ! Convert f_lambda to f_nu
         wavelengths(i) = temp_wavelength
         flux(i) = temp_flux
       END DO

       CLOSE(unit)

     END SUBROUTINE loadsed



  !---------------------------------------------------------------------------
  ! Apply dilution factor to convert surface flux to observed flux
  !---------------------------------------------------------------------------
  SUBROUTINE dilute_flux(surface_flux, R, d, calibrated_flux)
    REAL(dp), INTENT(IN)  :: surface_flux(:)
    REAL(dp), INTENT(IN)  :: R, d  ! R = stellar radius, d = distance (both in the same units, e.g., cm)
    REAL(dp), INTENT(OUT) :: calibrated_flux(:)

    ! Check that the output array has the same size as the input
    IF (SIZE(calibrated_flux) /= SIZE(surface_flux)) THEN
      PRINT *, "Error in dilute_flux: Output array must have the same size as input array."
      STOP 1
    END IF

    ! Apply the dilution factor (R/d)^2 to each element
    calibrated_flux = surface_flux * ((R / d)**2)
  END SUBROUTINE dilute_flux


END MODULE shared_func