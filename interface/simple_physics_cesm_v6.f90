Subroutine tphyssim (dtime, state, tend, surf_state )
!----------------------------------------------------------------------- 
! 
! Purpose: Simple Physics Package
!
! Author: K. A. Reed
! 
! Description: Includes large-scale precipitation, a call to convection and
!              turbulence (PBL and surface). The processes are time
!              split in that order
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid            , only: plev,plat,plevp
   use ppgrid
   use phys_grid         , only: get_lat_all_p, get_rlat_all_p
   use physics_types     , only: physics_state, physics_tend, physics_ptend, physics_dme_adjust
   use geopotential      , only: geopotential_t
   use cam_history,        only: outfld
   use physconst,          only: gravit, cappa, rair, cpair, latvap, rh2o, epsilo, karman, zvir, omega, rearth, pi
   use camsrfexch,         only: cam_out_t
   use cam_diagnostics,    only: diag_phys_writeout
   use perf_mod

   implicit none

!
! Input arguments
!
   real(r8), intent(in) :: dtime                   ! model timestep
!
! Output arguments
!
   type(physics_state), intent(inout) :: state
   type(physics_tend ), intent(inout) :: tend
   type(cam_out_t), intent(inout) :: surf_state
!
!---------------------------Local workspace-----------------------------
!
   type(physics_ptend)   :: ptend       ! indivdual parameterization tendencies

!===============================================================================
! Test Case Flags  -  NEEDS TO BE SET FOR DCMIP 2016
!===============================================================================
  
   integer :: dcmip_test_option = 2     ! = 1 for Moist Baroclinic Wave Test
                                        ! = 2 for Tropical Cyclone Test     
   logical :: RJ2012_precip = .TRUE.    ! Turn on (TRUE) or off (FALSE) Reed and Jablonowski (2012)
                                        ! large-scale condensation scheme
                                        ! Should be set to FALSE when running
                                        ! with Kessler physics 
   logical :: TC_PBL_mod = .FALSE.      ! Turn on (TRUE) or off (FALSE) George Bryan's 
                                        ! tropical cyclone PBL modification
                                        ! Should be set to FALSE wehn running
                                        ! Reed and Jablonowski (2012) configuration
!===============================================================================

! Integers for loops

   integer :: lchnk                     ! Chunk identifier
   integer :: ncol                      ! Number of atmospheric columns
   integer  i,k                         ! Longitude, level indices
    

! Simple Physics Specific Constants 

!+++++++
   real(r8) Tsurf(pcols)                ! Sea Surface Temperature
   real(r8) lat(pcols)                  ! Latitude for SST calculation
!+++++++

   real(r8) SST_TC                      ! Sea Surface Temperature for tropical cyclone test
   real(r8) T0                          ! Control temp for calculation of qsat
   real(r8) e0                          ! Saturation vapor pressure at T0 for calculation of qsat
   real(r8) rhow                        ! Density of Liquid Water
   real(r8) p0                          ! Constant for calculation of potential temperature
   real(r8) Cd0                         ! Constant for calculating Cd from Smith and Vogl 2008
   real(r8) Cd1                         ! Constant for calculating Cd from Smith and Vogl 2008
   real(r8) Cm                          ! Constant for calculating Cd from Smith and Vogl 2008
   real(r8) v20                         ! Threshold wind speed for calculating Cd from Smith and Vogl 2008
   real(r8) C                           ! Drag coefficient for sensible heat and evaporation
   real(r8) T00                         ! Horizontal mean T at surface for moist baro test
   real(r8) u0                          ! Zonal wind constant for moist baro test
   real(r8) latw                        ! halfwidth for  for baro test
   real(r8) eta0                        ! Center of jets (hybrid) for baro test
   real(r8) etav                        ! Auxiliary variable for baro test
   real(r8) q0                          ! Maximum specific humidity for baro test
   real(r8) kappa                       ! von Karman constant

! Temporary variables for tendency calculations

   real(r8) tmp                         ! Temporary
   real(r8) qsat                        ! Saturation vapor pressure
   real(r8) qsats                       ! Saturation vapor pressure of SST

! Variables for Boundary Layer Calculation

   real(r8) wind(pcols)                 ! Magnitude of Wind
   real(r8) Cd(pcols)                   ! Drag coefficient for momentum
   real(r8) Km(pcols,pver+1)            ! Eddy diffusivity for boundary layer calculations 
   real(r8) Ke(pcols,pver+1)            ! Eddy diffusivity for boundary layer calculations
   real(r8) rho                         ! Density at lower/upper interface
   real(r8) za(pcols)                   ! Heights at midpoints of first model level
   real(r8) zi(pcols,pver+1)            ! Heights at model interfaces
   real(r8) dlnpint                     ! Used for calculation of heights
   real(r8) pbltop                      ! Top of boundary layer
   real(r8) zpbltop                     ! Top of boundary layer for George Bryan Modifcation
   real(r8) pblconst                    ! Constant for the calculation of the decay of diffusivity 
   real(r8) CA(pcols,pver)              ! Matrix Coefficents for PBL Scheme 
   real(r8) CC(pcols,pver)              ! Matrix Coefficents for PBL Scheme 
   real(r8) CE(pcols,pver+1)            ! Matrix Coefficents for PBL Scheme
   real(r8) CAm(pcols,pver)             ! Matrix Coefficents for PBL Scheme
   real(r8) CCm(pcols,pver)             ! Matrix Coefficents for PBL Scheme
   real(r8) CEm(pcols,pver+1)           ! Matrix Coefficents for PBL Scheme
   real(r8) CFu(pcols,pver+1)           ! Matrix Coefficents for PBL Scheme
   real(r8) CFv(pcols,pver+1)           ! Matrix Coefficents for PBL Scheme
   real(r8) CFt(pcols,pver+1)           ! Matrix Coefficents for PBL Scheme
   real(r8) CFq(pcols,pver+1)           ! Matrix Coefficents for PBL Scheme

! Variables for Convection Calculation

   real(r8) t_ref(pcols,pver)           ! Temperature reference profile
   real(r8) q_ref(pcols,pver)           ! Specific humidity reference profile
   real(r8) bmflag(pcols)               ! Flag for which type of convection
   real(r8) klzb(pcols)                 ! Model level of LZB
   real(r8) cape(pcols)                 ! Convectively available potential energy
   real(r8) cin(pcols)                  ! Convective inhibition
   real(r8) capeflag(pcols)             ! Flag that says why CAPE is zero
   real(r8) invtau_bm_q(pcols)          ! 1 / B-M relaxation timescale for q 
   real(r8) invtau_bm_t(pcols)          ! 1 / B-M relaxation timescale for T
   real(r8) tdel(pcols,pver)            ! Temperature difference
   real(r8) qdel(pcols,pver)            ! Specific humidity difference

! Precipitation Arrays

   real(r8) prectot(pcols)              ! Total precipitation (convective + large-scale)

! Variable for Dry Mass Adjustment 

   real(r8) qini(pcols,pver)            ! Initial specific humidity

!
!-----------------------------------------------------------------------
!

   lchnk = state%lchnk
   ncol  = state%ncol

!===============================================================================
!
! Local Constants for Simple Physics
!
!===============================================================================
      C        = 0.0011_r8      ! From Simth and Vogl 2008
      SST_TC   = 302.15_r8      ! Constant Value for SST
      T0       = 273.16_r8      ! control temp for calculation of qsat
      e0       = 610.78_r8      ! saturation vapor pressure at T0 for calculation of qsat
      rhow     = 1000.0_r8      ! Density of Liquid Water 
      Cd0      = 0.0007_r8      ! Constant for Cd calc. Simth and Vogl 2008
      Cd1      = 0.000065_r8    ! Constant for Cd calc. Simth and Vogl 2008
      Cm       = 0.002_r8       ! Constant for Cd calc. Simth and Vogl 2008
      v20      = 20.0_r8        ! Threshold wind speed for calculating Cd from Smith and Vogl 2008
      p0       = 100000.0_r8    ! Constant for potential temp calculation
      pbltop   = 85000._r8      ! Top of boundary layer in p
      zpbltop  = 1000._r8       ! Top of boundary layer in z
      pblconst = 10000._r8      ! Constant for the calculation of the decay of diffusivity
      T00      = 288.0_r8         ! Horizontal mean T at surface for moist baro test
      u0       = 35.0_r8          ! Zonal wind constant for moist baro test
      latw     = 2.0_r8*pi/9.0_r8 ! Halfwidth for  for baro test
      eta0     = 0.252_r8         ! Center of jets (hybrid) for baro test
      etav     = (1._r8-eta0)*0.5_r8*pi ! Auxiliary variable for baro test
      q0       = 0.021_r8         ! Maximum specific humidity for baro test
      kappa    = 0.4_r8         ! von Karman constant

!===============================================================================
!
! Definition of local arrays
!
!===============================================================================
!
! Calculate hydrostatic height
!
     do i=1,ncol
        dlnpint = (state%lnpint(i,pver+1) - state%lnpint(i,pver))
        za(i) = rair/gravit*state%t(i,pver)*(1._r8+zvir*state%q(i,pver,1))*0.5_r8*dlnpint
        zi(i,pver+1) = 0.0_r8
     end do
!
! Set Intial Total Precipitation to Zero
!
     do i=1,ncol
        prectot(i) = 0._r8
     end do

!
! Set Initial Specific Humidity
!
     qini(:ncol,:pver) = state%q(:ncol,:pver,1)
!
! Set Sea Surface Temperature (constant for tropical cyclone)
! Tsurf needs to be dependent on latitude for moist baroclinic wave test
! Tsurf needs to be constant for tropical cyclone test
!
     if (dcmip_test_option .eq. 1) then ! Moist Baroclinic Wave Test
        call get_rlat_all_p(lchnk, ncol, lat)
        do i=1,pcols
           Tsurf(i) = (T00 + pi*u0/rair * 1.5_r8 * sin(etav) * (cos(etav))**0.5_r8 *                 &
                     ((-2._r8*(sin(lat(i)))**6 * ((cos(lat(i)))**2 + 1._r8/3._r8) + 10._r8/63._r8)* &
                     u0 * (cos(etav))**1.5_r8  +                                                    &
                     (8._r8/5._r8*(cos(lat(i)))**3 * ((sin(lat(i)))**2 + 2._r8/3._r8) - pi/4._r8)*rearth*omega*0.5_r8 ))/ &
                     (1._r8+zvir*q0*exp(-(lat(i)/latw)**4))

        end do
     else if (dcmip_test_option .eq. 2) then ! Tropical Cyclone Test
        do i=1,pcols
           Tsurf(i) = SST_TC
        end do
     end if

!===============================================================================
!
! Large-Scale Condensation and Precipitation
!
! note: Scheme is set forth by the ECMWF (without cloud stage)
!
!===============================================================================

      if (RJ2012_precip) then
!
! Calculate Tendencies
!
      do k=1,pver
         do i=1,ncol
            qsat = epsilo*e0/state%pmid(i,k)*exp(-latvap/rh2o*((1._r8/state%t(i,k))-1._r8/T0))
            if (state%q(i,k,1) > qsat) then
               tmp  = 1._r8/dtime*(state%q(i,k,1)-qsat)/(1._r8+(latvap/cpair)*(epsilo*latvap*qsat/(rair*state%t(i,k)**2)))
               tend%dtdt(i,k) = tend%dtdt(i,k)+latvap/cpair*tmp
               tend%dqdt(i,k) = tend%dqdt(i,k)-tmp
               surf_state%precl(i) = surf_state%precl(i)+tmp*state%pdel(i,k)/(gravit*rhow)
            end if
         end do
      end do
!
! Update moisture and temperature fields from Larger-Scale Precipitation Scheme
!
      do k=1,pver
         do i=1,ncol
            state%t(i,k) =  state%t(i,k) + tend%dtdt(i,k)*dtime
            state%q(i,k,1) =  state%q(i,k,1) + tend%dqdt(i,k)*dtime
         end do
      end do

!===============================================================================
! Send variables to history file
!===============================================================================
      call diag_phys_writeout(state)
    
      end if
     
!===============================================================================
!
! Turbulence Scheme for Momentum, Heat and Moisture
!
! note: The Boundary layer is assumed to only be in the lower 
!       lowermost model levels, including a parameterized
!       surface flux.  So, here we parameterize the PBL for the 
!       lowermost model levels below 850 hPa.  We are using Simplified Ekman 
!       theory (constant Ke) but Ke is calculated each time step
!       and in each column.  First, the fields are updated with the
!       surface flux at the lowermost model level and then the implicit 
!       PBL scheme is implemented.
!
!===============================================================================
!
! Compute magnitude of the wind and drag coeffcients for turbulence scheme
!
     do i=1,ncol
        wind(i) = sqrt(state%u(i,pver)**2+state%v(i,pver)**2)
     end do
     do i=1,ncol
        if( wind(i) .lt. v20) then
           Cd(i) = Cd0+Cd1*wind(i) 
        else
           Cd(i) = Cm
        endif
     end do

     if (TC_PBL_mod) then !Bryan TC PBL Modification 
     do k=pver,1,-1
        do i=1,ncol
           dlnpint = (state%lnpint(i,k+1) - state%lnpint(i,k))
           zi(i,k) = zi(i,k+1)+rair/gravit*state%t(i,k)*(1._r8+zvir*state%q(i,k,1))*dlnpint
           if( zi(i,k) .le. zpbltop) then
              Km(i,k) = kappa*sqrt(Cd(i))*wind(i)*zi(i,k)*(1._r8-zi(i,k)/zpbltop)*(1._r8-zi(i,k)/zpbltop)
              Ke(i,k) = kappa*sqrt(C)*wind(i)*zi(i,k)*(1._r8-zi(i,k)/zpbltop)*(1._r8-zi(i,k)/zpbltop) 
           else
              Km(i,k) = 0.0_r8
              Ke(i,k) = 0.0_r8
           end if 
        end do
     end do     
     else ! Reed and Jablonowski (2012) Configuration
     do k=1,pver
        do i=1,ncol
           if( state%pint(i,k) .ge. pbltop) then
              Km(i,k) = Cd(i)*wind(i)*za(i) 
              Ke(i,k) = C*wind(i)*za(i)
           else
              Km(i,k) = Cd(i)*wind(i)*za(i)*exp(-(pbltop-state%pint(i,k))**2/(pblconst)**2)
              Ke(i,k) = C*wind(i)*za(i)*exp(-(pbltop-state%pint(i,k))**2/(pblconst)**2)
           end if
        end do
     end do
     end if
!
! Compute surface fluxes using an implicit approach
! note: this only occurs in lowermost model level
!
     do i=1,ncol
        qsats = epsilo*e0/state%ps(i)*exp(-latvap/rh2o*((1._r8/Tsurf(i))-1._r8/T0))
        tend%dudt(i,pver) = tend%dudt(i,pver) + (state%u(i,pver)/(1._r8+Cd(i)*wind(i)*dtime/za(i))-state%u(i,pver))/dtime
        tend%dvdt(i,pver) = tend%dvdt(i,pver) + (state%v(i,pver)/(1._r8+Cd(i)*wind(i)*dtime/za(i))-state%v(i,pver))/dtime
        state%u(i,pver) = state%u(i,pver)/(1._r8+Cd(i)*wind(i)*dtime/za(i))
        state%v(i,pver) = state%v(i,pver)/(1._r8+Cd(i)*wind(i)*dtime/za(i))
        tend%dtdt(i,pver) = tend%dtdt(i,pver) +((state%t(i,pver)+C*wind(i)*Tsurf(i)*dtime/za(i))/(1._r8+C*wind(i)*dtime/za(i))-state%t(i,pver))/dtime 
        state%t(i,pver) = (state%t(i,pver)+C*wind(i)*Tsurf(i)*dtime/za(i))/(1._r8+C*wind(i)*dtime/za(i))  
        tend%dqdt(i,pver) = tend%dqdt(i,pver) +((state%q(i,pver,1)+C*wind(i)*qsats*dtime/za(i))/(1._r8+C*wind(i)*dtime/za(i))-state%q(i,pver,1))/dtime
        state%q(i,pver,1) = (state%q(i,pver,1)+C*wind(i)*qsats*dtime/za(i))/(1._r8+C*wind(i)*dtime/za(i))
     end do


!
! Calculate Diagonal Variables for Implicit PBL Scheme
!
      do k=1,pver-1
         do i=1,ncol
            rho = (state%pint(i,k+1)/(rair*(state%t(i,k+1)+state%t(i,k))/2.0_r8))
            CAm(i,k) = state%rpdel(i,k)*dtime*gravit*gravit*Km(i,k+1)*rho*rho/(state%pmid(i,k+1)-state%pmid(i,k))    
            CCm(i,k+1) = state%rpdel(i,k+1)*dtime*gravit*gravit*Km(i,k+1)*rho*rho/(state%pmid(i,k+1)-state%pmid(i,k))
            CA(i,k) = state%rpdel(i,k)*dtime*gravit*gravit*Ke(i,k+1)*rho*rho/(state%pmid(i,k+1)-state%pmid(i,k))
            CC(i,k+1) = state%rpdel(i,k+1)*dtime*gravit*gravit*Ke(i,k+1)*rho*rho/(state%pmid(i,k+1)-state%pmid(i,k))
         end do
      end do
      do i=1,ncol
         CAm(i,pver) = 0._r8
         CCm(i,1) = 0._r8
         CEm(i,pver+1) = 0._r8
         CA(i,pver) = 0._r8
         CC(i,1) = 0._r8
         CE(i,pver+1) = 0._r8
         CFu(i,pver+1) = 0._r8
         CFv(i,pver+1) = 0._r8
         CFt(i,pver+1) = 0._r8
         CFq(i,pver+1) = 0._r8 
      end do
      do i=1,ncol
         do k=pver,1,-1
            CE(i,k) = CC(i,k)/(1._r8+CA(i,k)+CC(i,k)-CA(i,k)*CE(i,k+1)) 
            CEm(i,k) = CCm(i,k)/(1._r8+CAm(i,k)+CCm(i,k)-CAm(i,k)*CEm(i,k+1))
            CFu(i,k) = (state%u(i,k)+CAm(i,k)*CFu(i,k+1))/(1._r8+CAm(i,k)+CCm(i,k)-CAm(i,k)*CEm(i,k+1))
            CFv(i,k) = (state%v(i,k)+CAm(i,k)*CFv(i,k+1))/(1._r8+CAm(i,k)+CCm(i,k)-CAm(i,k)*CEm(i,k+1))
            CFt(i,k) = ((p0/state%pmid(i,k))**(rair/cpair)*state%t(i,k)+CA(i,k)*CFt(i,k+1))/(1._r8+CA(i,k)+CC(i,k)-CA(i,k)*CE(i,k+1)) 
            CFq(i,k) = (state%q(i,k,1)+CA(i,k)*CFq(i,k+1))/(1._r8+CA(i,k)+CC(i,k)-CA(i,k)*CE(i,k+1))
       end do
      end do
!
! Calculate the updated temperaure and specific humidity and wind tendencies
!
! First we need to calculate the tendencies at the top model level
!
      do i=1,ncol
            tend%dudt(i,1)  = tend%dudt(i,1)+(CFu(i,1)-state%u(i,1))/dtime
            tend%dvdt(i,1)  = tend%dvdt(i,1)+(CFv(i,1)-state%v(i,1))/dtime
            state%u(i,1)    = CFu(i,1)
            state%v(i,1)    = CFv(i,1)
            tend%dtdt(i,1)  = tend%dtdt(i,1)+(CFt(i,1)*(state%pmid(i,1)/p0)**(rair/cpair)-state%t(i,1))/dtime
            state%t(i,1)    = CFt(i,1)*(state%pmid(i,1)/p0)**(rair/cpair)
            tend%dqdt(i,1)  = tend%dqdt(i,1)+(CFq(i,1)-state%q(i,1,1))/dtime
            state%q(i,1,1)  = CFq(i,1)
      end do

      do i=1,ncol
         do k=2,pver
            tend%dudt(i,k)  = tend%dudt(i,k)+(CEm(i,k)*state%u(i,k-1)+CFu(i,k)-state%u(i,k))/dtime
            tend%dvdt(i,k)  = tend%dvdt(i,k)+(CEm(i,k)*state%v(i,k-1)+CFv(i,k)-state%v(i,k))/dtime
            state%u(i,k)    = CEm(i,k)*state%u(i,k-1)+CFu(i,k) 
            state%v(i,k)    = CEm(i,k)*state%v(i,k-1)+CFv(i,k)
            tend%dtdt(i,k)  = tend%dtdt(i,k)+((CE(i,k)*state%t(i,k-1)*(p0/state%pmid(i,k-1))**(rair/cpair)+CFt(i,k))*(state%pmid(i,k)/p0)**(rair/cpair)-state%t(i,k))/dtime 
            state%t(i,k)    = (CE(i,k)*state%t(i,k-1)*(p0/state%pmid(i,k-1))**(rair/cpair)+CFt(i,k))*(state%pmid(i,k)/p0)**(rair/cpair)
            tend%dqdt(i,k)  = tend%dqdt(i,k)+(CE(i,k)*state%q(i,k-1,1)+CFq(i,k)-state%q(i,k,1))/dtime
            state%q(i,k,1)  = CE(i,k)*state%q(i,k-1,1)+CFq(i,k)
         end do
      end do

!===============================================================================
!
! Compute Total Precipitation
!
!===============================================================================
      do i=1,ncol
         prectot(i) = surf_state%precl(i) + surf_state%precc(i)
      end do

!===============================================================================
!
! Dry Mass Adjustment 
!
!===============================================================================
 call physics_dme_adjust(state, tend, qini, dtime)

!! Archive idealized temperature tendency and precipitation
!
   call outfld('PRECL   ',surf_state%precl      ,pcols   ,lchnk       )
   call outfld('PRECC   ',surf_state%precc      ,pcols   ,lchnk       )
   call outfld('PRECT   ',prectot               ,pcols   ,lchnk       )

   return
end subroutine tphyssim

