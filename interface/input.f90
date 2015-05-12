

!=======================================================================
!
!   Interface input call
!
!
!
!=======================================================================


 SUBROUTINE input_parameter(lon, lat, z, p, X , deep, moist, pertt, zcoords,dt)

    IMPLICIT NONE

!-----------------------------------------------------------------------
!     input parameters asked to the user
!-----------------------------------------------------------------------
    INTEGER, INTENT(OUT)  :: &
    deep,       & ! Deep (1) or Shallow (0) test case
    moist,      & ! Moist (1) or Dry (0) test case
    pertt         ! Perturbation type

    REAL(8), INTENT(OUT)  :: &
    lon,        & ! Longitude (radians)
    lat,        & ! Latitude (radians)
    X,          & ! Earth scaling parameter
    dt,         & ! size of physics time step in second
    z,          & ! Altitude (m) 
    p             ! Pressure (Pa)



    INTEGER, INTENT(OUT) :: zcoords     ! 1 if z coordinates are specified
                                        ! 0 if p coordinates are specified

print *,"Deep"
read  *,deep
print *,"moist"
read  *,moist
print *,"pertt"
read  *,pertt
print *,"X"
read  *,X
print *,"type 1 if z is specified, 0 if p is specified"
read  *,zcoords
print *,"longitude"
read  *,lon
!print *,"latitude"
read  *,lat
!print *,"size of physics time step in second"
read  *,dt


if (zcoords .eq. 1) then
    print *,"z"
    !read  *,z
else
print *,"p"
read  *,p

end if

END SUBROUTINE input_parameter






