program vertical_coordinate_generator

  ! This program generates a distribution of sigma levels
  ! that satisfies the following criteria::
  !   a) approx total number of levels (N)
  !   b) the lowest half-model level (sig1)
  !   c) number of levels desired in a stable pbl
  !   d) minimum number of stratospheric levels (above sigma = 0.1)
  ! Then, it converts the sigma levels to height levels and outputs both.
  !
  ! Mostly written by Joseph Olson (NOAA/GSL), but some of the first-draft
  ! code was written by Gemini.

  implicit none

  integer, parameter :: dp = kind(1.0d0)
  ! Physical constants or parameters
  real(dp), parameter :: rd   = 287.05    ! Gas constant for dry air [J/(kg K)]
  real(dp), parameter :: g    = 9.81      ! Acceleration due to gravity [m/s^2]
  real(dp), parameter :: psfc = 101325.0  ! 1013.25 mb in Pascals
  real(dp), parameter :: zsfc = 0.0       ! Sea-level height

  !specified grid specs
  integer :: nsig_cos    !approximate number of analytical (cosine) levels (suggest: 50-150)
  integer :: nsig_pbl    !exact number of levels desired in a stable boundary layer (suggest: 8-15)
  real(dp):: z_spbl      !approx height of stable pblg (m) (suggest: 200-400)
  integer :: nsig_p1     !number of sigma levels abive sigma = 0.1 (suggest: >15)
  real(dp):: sig1        !top pf lowest model layer (suggest: 0.998-0.997)
  real(dp):: alfa1       !alfa1 > 1 increases clustering near the surface; alfa1 < 1 increases clustering near the top
  real(dp):: ptop        !pressure (Pa) at model top, needed for sigma to z conversion
  real(dp):: ddelz_max   !maximum delta(delta-z) for poor vertical resolution and/or high ptop (suggest: 300 m)
  
  !derived variables
  real(dp):: dsig1       ! 1.0 - sig1
  integer :: nsig_top    !number of sigma levels used for high-resolution near top of model, derived from nsig_p05
  real(dp):: sbl_top     !sigma at top of stable pbl + delta (=0.005, for blending reasons)
  real(dp):: blend_top   !top of blending zone
  real(dp):: mdsig_high  !mean delta-sigma in the stratosphere
  real(dp):: mdsig_pbl   !mean delta-sigma in the stable pbl
  real(dp):: min_dsig
  real(dp):: min_dsig_trop
  real(dp):: min_dsig_pbl
  real(dp):: midsig      !sigma at the middle of the blending regions
  real(dp):: dsigma_int  !interpolated delta-sigma
  real(dp):: dsigma_int2 !interpolated delta-sigma

  !1d arrays
  real(dp),dimension(1:500):: sigma        !(final) output sigma levels
  real(dp),dimension(1:500):: sigma_cos    !analytical (cosine) distribution with nsig_cos levels
  real(dp),dimension(1:500):: sigma_pbl    !levels in the stable pbl
  real(dp),dimension(1:500):: sigma_top    !levels near the top of the model
  real(dp),dimension(1:500):: dsigma       !(final) output sigma layer depths
  real(dp),dimension(1:500):: dsigma_cos   !analytical (cosine) sigma layer depths
  real(dp),dimension(1:500):: dsigma_pbl   !sigma layer depths within the stable pbl
  real(dp),dimension(1:500):: dsigma_top   !sigma layer depths near top of model
  real(dp),dimension(1:500):: tv_mean      !mean virtual temperature profile (K)
  real(dp),dimension(1:500):: p            !pressure profile (Pa)
  real(dp),dimension(1:500):: z            !height (m)
  real(dp),dimension(1:500):: delz
  
  !miscellaneous
  real(dp):: pi,delta,x,dsig,wt,wt_int,rsig,ddelz,ddelz0,ddelz1
  integer :: k,k3,ks,ki,nlevs_top,k_high,nz,limit_ddelz,nddelz,kddelz
  logical,parameter::verbose=.false.       !extra output for debugging
  
 ! Define namelist group
  namelist /vcg_config/ sig1, nsig_cos, nsig_pbl, z_spbl, nsig_p1, alfa1, ptop, ddelz_max

  pi         = 4.0_dp * atan(1.0_dp)

  !set default sigma level specs
  sig1       = 0.9975_dp    !first sigma level above the surface (top of first model layer)
  nsig_cos   = 65           !number of layers in the analytical (cosine) sigma levels
  nsig_pbl   = 10           !number of layers within the stable pbl (0 to z_spbl)
  z_spbl     = 300._dp      !height of stable boundary layer (meters)
  nsig_p1    = 12           !minimum number of layers above sigma = 0.05
  alfa1      = 1.0_dp       !skewness factor
  ptop       = 200.         !pressure at model top (Pascals)
  ddelz_max  = 300.         !maximum allowable increase in delta-z (meters)
  
  ! Open and read namelist file
  open(unit=10, file='vcg.nml', status='old', action='read')
  read(10, nml=vcg_config)
  close(10)
  
  !initialize arrays to zero
  sigma      = 0._dp
  dsigma     = 0._dp
  sigma_cos  = 0._dp
  dsigma_cos = 0._dp
  sigma_pbl  = 0._dp
  dsigma_pbl = 0._dp
  sigma_top  = 0._dp
  dsigma_top = 0._dp

  !deriver variables based on input specs
  dsig1      = 1.0_dp - sig1
  sbl_top    = 1.0_dp - (z_spbl * 1e-4)
  blend_top  = sbl_top - 0.2_dp
  print*,"sbl_top=",sbl_top," blend_top=",blend_top
  
  !----------------------------------------------------------------
  print*," Generating analytical (cosine) sigma levels:"
  do k = 1, nsig_cos+1
     !sigma(k) = 0.5_dp * (1.0_dp + cos(pi * real(k-1,dp) / real(nsig_cos,dp)))
     !Sharper mid-level stretching:
     x = (real(k-1,dp) / real(nsig_cos,dp))**alfa1
     sigma_cos(k) = 0.5_dp * (1.0_dp + cos(pi * x))
  enddo

  if(verbose)print*," Computing the corresponding delta-sigmas..."
  do k = 1, nsig_cos
     dsigma_cos(k) = sigma_cos(k) - sigma_cos(k+1)
  enddo

  print*," Generating the near-surface sigma levels in the stable pbl:"
  if(verbose)print*," Calculate the mean delta-sigma within the stable pbl:"
  mdsig_pbl = (1.0_dp - sbl_top)/real(nsig_pbl,dp)
  if (dsig1 > mdsig_pbl) then
     print*,"Need to lower sig1 to:", 1.0_dp-mdsig_pbl
     print*," or decrease the resolution in the pbl."
     sig1  = 1.0_dp-mdsig_pbl
     dsig1 = mdsig_pbl
  endif
  
  sigma_pbl(1)  = 1.0_dp
  sigma_pbl(2)  = sig1
  dsigma_pbl(1) = sigma_pbl(1) - sigma_pbl(2)

  wt = 1._dp/real(nsig_pbl, dp)
  dsigma_pbl(2) = (1._dp-wt)*dsig1 + wt*(mdsig_pbl + (mdsig_pbl-dsig1))
  delta =  dsigma_pbl(2) - dsigma_pbl(1) 
  
  do k=3,500
     !specify linear increase of dsigma_pbl with height
     dsigma_pbl(k) = dsigma_pbl(k-1) + delta 
     sigma_pbl(k)=sigma_pbl(k-1) - dsigma_pbl(k-1)
     if (sigma_pbl(k) < 0._dp) then
        sigma_pbl(k) = 0._dp
        dsigma_pbl(k-1)=sigma_pbl(k-1)
        exit
     endif
  enddo

  print*," Now blending the cosine levels with the pbl levels:"
  sigma(1:nsig_pbl)    = sigma_pbl(1:nsig_pbl)
  dsigma(1:nsig_pbl-1) = dsigma_pbl(1:nsig_pbl-1)

  midsig      = 0.5_dp*(sbl_top + blend_top)
  do k=2,500
     !interpolate dsigma_cos to nearest estimated sigma 
     do ki = 2, nsig_cos
        if (sigma_cos(ki) < sigma(k-1)) then
           wt_int = (sigma_cos(ki-1)-sigma(k-1))/(sigma_cos(ki-1)-sigma_cos(ki))
           dsigma_int = (1._dp-wt_int)*dsigma_cos(ki-1) + wt_int*dsigma_cos(ki)
           if(verbose)print*,"------------------------------------------------------"
           if(verbose)print*,"ki=",ki," wt_int=",wt_int," dsigma_int=",dsigma_int
           if(verbose)print*,"dsigma_cos(ki-1)=",dsigma_cos(ki-1)," dsigma_cos(ki)=",dsigma_cos(ki)
           exit
        endif
     enddo
     min_dsig_trop = max(dsigma_int, 0.0004)
     !interpolate dsigma_pbl to nearest estimated sigma
     do ki = 2, 500
        if (sigma_pbl(ki) < sigma(k-1)) then
           wt_int = (sigma_pbl(ki-1)-sigma(k-1))/(sigma_pbl(ki-1)-sigma_pbl(ki))
           dsigma_int2 = (1._dp-wt_int)*dsigma_pbl(ki-1) + wt_int*dsigma_pbl(ki)
           if(verbose)print*,"------------------------------------------------------"
           if(verbose)print*,"ki=",ki," wt_int=",wt_int," dsigma_int2=",dsigma_int2
           if(verbose)print*,"dsigma_pbl(ki-1)=",dsigma_pbl(ki-1)," dsigma_pbl(ki)=",dsigma_pbl(ki)
           exit
        endif
     enddo
     min_dsig_pbl = max(dsigma_int2, 0.0004)
     wt          = 0.5_dp*TANH((sigma(k-1) - midsig)/(0.20_dp*(sbl_top - blend_top))) + 0.5_dp
     if(verbose)print '(A2,I4,A4,F7.5,A12,f7.5)',"k=",k," wt=",wt," sigma(k-1)=",sigma(k-1)
     if(verbose)print '(A15,F7.5,A15,f7.5)',"--min_dsig_pbl=",min_dsig_pbl," min_dsig_trop=",min_dsig_trop
     dsigma(k-1) = wt*dsigma_pbl(k-1) + (1._dp-wt)*min_dsig_trop
     dsigma(k-1) = wt*min_dsig_pbl + (1._dp-wt)*min_dsig_trop
     sigma(k)    = sigma(k-1) - dsigma(k-1)
     
     if (sigma(k) <= 0.0_dp) then
        sigma(k)    = 0.0_dp
        dsigma(k)   = 0.0_dp
        !even out spacing for the top two sigma layers
        rsig        = dsigma(k-3) / (dsigma(k-3) + dsigma(k-2))
        dsigma(k-2) = sigma(k-2) * rsig
        dsigma(k-1) = sigma(k-2) * (1._dp - rsig)
        sigma(k-1)  = sigma(k-2) - dsigma(k-2)
        k3 = k
        exit
     endif
  enddo
  
  print*," Check to make sure adequate levels exist above sigma=0.1"
  min_dsig      = 0.0003_dp
  nlevs_top     = 0
  do k = nsig_pbl+1,k3
     if (sigma(k) < 0.1) nlevs_top = nlevs_top + 1
  enddo
  if (nlevs_top >= nsig_p1) then
     print*," --the number of levels above sigma=0.1 is already adequate--no need to modify"
  else
     print*," --Must add some levels above sigma=0.1. Only had",nlevs_top
     print*," --Generating analytical (cosine) sigma levels for top of mmodel:"
     nsig_top = int(5.333_dp*(real(nsig_p1,dp) - 2.2_dp))
     do k = 1, nsig_top+1
        x = (real(k-1,dp) / real(nsig_top,dp))   !note alfa1 is removed for better control
        sigma_top(k) = 0.5_dp * (1.0_dp + cos(pi * x))
     enddo
     if(verbose)print*," Computing the corresponding delta-sigmas..."
     do k = 1, nsig_top
        dsigma_top(k) = sigma_top(k) - sigma_top(k+1)
     enddo

     midsig      = 0.5_dp*(0.15_dp + 0.05_dp)
     do k=2,500
        !interpolate dsigma_top to nearest estimated sigma
        do ki = 2, nsig_top
           if (sigma_top(ki) < sigma(k-1)) then
              wt_int     = (sigma_top(ki-1)-sigma(k-1))/(sigma_top(ki-1)-sigma_top(ki))
              dsigma_int = (1._dp-wt_int)*dsigma_top(ki-1) + wt_int*dsigma_top(ki)
              if(verbose)print*,"------------------------------------------------------"
              if(verbose)print*,"ki=",ki," wt_int=",wt_int," dsigma_int=",dsigma_int
              if(verbose)print*,"dsigma_top(ki-1)=",dsigma_top(ki-1)," dsigma_top(ki)=",dsigma_top(ki)
              exit
          endif
       enddo
       min_dsig_trop = max(dsigma_int, min_dsig)
       wt            = 0.5_dp*TANH((sigma(k-1) - midsig)/(0.25_dp*(0.15_dp - 0.05_dp))) + 0.5_dp
       dsigma(k-1)   = wt*dsigma(k-1) + (1._dp-wt)*min_dsig_trop
       sigma(k)      = sigma(k-1) - dsigma(k-1)

       if (sigma(k) <= 0.0_dp) then
          sigma(k)    = 0.0_dp
          dsigma(k-1) = sigma(k-1)
          dsigma(k)   = 0.0_dp
          k3 = k
          exit
       endif
    enddo
    nlevs_top     = 0
    do k = nsig_pbl+1,k3
       if (sigma(k) < 0.1) nlevs_top = nlevs_top + 1
    enddo
    print*," --Now we have",nlevs_top,"levels"
  endif

  !----------------------------------------------------------------------
  ! convert to height coordinates for MPAS
  !----------------------------------------------------------------------
  ! Define dummy mean virtual temperatures for the column [K]
  nz = max(nsig_cos+1,k3)
  do k = 1, nz
     tv_mean(k) = (724.0_dp - 40.0_dp*sigma(k) + 241.0_dp*((sigma(k)-0.1)**2 + 0.0004_dp))/3.0_dp
  enddo
       
  ! Call the conversion subroutine
  call calc_height_from_sigma(nz, sigma(1:nz), psfc, ptop, zsfc, tv_mean(1:nz), p(1:nz), z(1:nz))  

  !----------------------------------------------------------------------
  ! perform check on d(delz) to limit excessive growth if it surpasses 300 m
  !----------------------------------------------------------------------
  delz(1) = 0._dp
  limit_ddelz = 0
  nddelz = 0
  do k = 2, nz
     delz(k) = z(k) - z(k-1)
     ddelz   = delz(k) - delz(k-1)
     if ((ddelz > ddelz_max) .and. (limit_ddelz == 0) .and. (k .ne. nz)) then
        limit_ddelz = 1
        ddelz0 = delz(k-1) - delz(k-2)
        ddelz1 = delz(k)   - delz(k-1)
        kddelz = k
        print*,"Warning: large increase in delta-z is detected."
        print*,"--recommend increasing nsig_cos and/or reducing ptop."
        print'(I4,3F12.6)',k-1,z(k-1),delz(k-1),ddelz0
     endif
     !impose limit for all levels above the first levels that pierces the limit
     if (limit_ddelz == 1) then
        nddelz = nddelz + 1
        ddelz  = ddelz0 + (ddelz1 - ddelz0) * (1._dp - 1._dp/(1._dp + real(nddelz,dp)))
     endif
     delz(k) = delz(k-1) + ddelz
     z(k)    = z(k-1) + delz(k)
     if (limit_ddelz .gt. 0) then
        print'(I4,3F12.6)',k,z(k),delz(k),ddelz
     endif
  enddo
  
  !----------------------------------------------------------------------
  ! Output results
  !----------------------------------------------------------------------
  print *, '-------------------------------------------------------------'
  print *, ' actual number of sigma (interface) levels needed:',k3
  print *, '   k  sigma_cos    sigma     dsigma_cos     dsigma    height'
  print *, '-------------------------------------------------------------'
  do k = 1, nz
     print '(I4,4F12.8,F10.2)', k, sigma_cos(k),sigma(k),dsigma_cos(k),dsigma(k),z(k)
     !with commas
     !print '(I4,A1,F14.8,A1,F14.8,A1,F14.8,A1,F14.8)', k,",",sigma_cos(k),",",sigma(k),",",dsigma_cos(k),",",dsigma(k)
     !print '(F15.8)',sigma(k)
  end do
  print*
  print*,"sigma:"
  do k = 1, nz
     print '(F15.8)',sigma(k)
  end do
  print*
  print*,"sigma for WRF:"
  do k = 1, nz-1
     write(*,'(F11.8,A1)',advance='NO' ),sigma(k),","
  end do
  write(*,'(F11.8,A1)'),sigma(k),","
  print*
  print*,"delta-sigma:"
  do k = 1, nz
     print '(F15.8)',dsigma(k)
  end do
  print*
  print*,"height"
  do k = 1, nz
     print '(F15.6)',z(k)
  end do

  contains
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------

  subroutine calc_height_from_sigma(nz, sigma, psfc, ptop, zsfc, tv_mean, p, z)
    ! Arguments
        integer, parameter :: dp = kind(1.0d0)
        integer, intent(in) :: nz                               ! Number of vertical levels
        real(dp), intent(in), dimension(nz) :: sigma            ! Sigma levels (fraction, 0 to 1)
        real(dp), intent(in) :: psfc                            ! Surface pressure [Pa]
        real(dp), intent(in) :: ptop                            ! Model top pressure [Pa]
        real(dp), intent(in) :: zsfc                            ! Surface elevation / terrain height [m]
        real(dp), intent(in), dimension(nz) :: tv_mean          ! Mean virtual temp from surface to level [K]
        
        real(dp), intent(out), dimension(nz) :: p               ! Calculated pressure at sigma levels [Pa]
        real(dp), intent(out), dimension(nz) :: z               ! Calculated height at sigma levels [m]
        
        ! Local variables
        integer :: k

        ! Loop over each vertical level
        do k = 1, nz
            ! Step 1: Calculate the pressure at the current sigma level
            p(k) = sigma(k) * (psfc - ptop) + ptop
            
            ! Step 2: Calculate the geopotential height using the Hypsometric Equation
            ! Note: Avoid log(1) edge cases if a sigma level is exactly at the surface
            if (p(k) >= psfc) then
                z(k) = zsfc
            else
                z(k) = zsfc + ((rd * tv_mean(k)) / g) * log(psfc / p(k))
            end if
        end do

  end subroutine calc_height_from_sigma

end program vertical_coordinate_generator

