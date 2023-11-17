!-------------------------------------------------------------------------------------------
!Copyright (c) 2013-2016 by Wolfgang Kurtz .and. Guowei He (Forschungszentrum Juelich GmbH)
!
!This file is part of TerrSysMP-PDAF
!
!TerrSysMP-PDAF is free software: you can redistribute it .and./or modify
!it under the terms of the GNU Lesser General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.
!
!TerrSysMP-PDAF is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU LesserGeneral Public License for more details.
!
!You should have received a copy of the GNU Lesser General Public License
!along with TerrSysMP-PDAF.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------------------
!obs_op_pdaf.F90: TerrSysMP-PDAF implementation of routine
!                 'obs_op_pdaf' (PDAF online coupling)
!-------------------------------------------------------------------------------------------

!$Id: obs_op_pdaf.F90 1441 2013-10-04 10:33:42Z lnerger $
!BOP
!
! !ROUTINE: obs_op_pdaf --- Implementation of observation operator
!
! !INTERFACE:
SUBROUTINE obs_op_pdaf(step, dim_p, dim_obs_p, state_p, m_state_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/SEIK/EnKF/ETKF/ESTKF
!
! The routine is called during the analysis step.
! It has to perform the operation of the
! observation operator acting on a state vector.
! For domain decomposition, the action is on the
! PE-local sub-domain of the state .and. has to
! provide the observed sub-state for the PE-local
! domain.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
   USE mod_assimilation, &
        ONLY: obs_index_p, obs_p
#if (defined CLMFIVE)
   use clm_instMod, only : soilstate_inst
   use clm_varcon  , only : zsoi
   use enkf_clm_mod, only : clmupdate_swc, clmcrns_bd, clmcrns_nflux
   use mod_parallel_model, only: mype_world
#endif

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step               ! Currrent time step
  INTEGER, INTENT(in) :: dim_p              ! PE-local dimension of state
  INTEGER, INTENT(in) :: dim_obs_p          ! Dimension of observed state
  REAL, INTENT(in)    :: state_p(dim_p)     ! PE-local model state
  REAL, INTENT(out) :: m_state_p(dim_obs_p) ! PE-local observed state
  integer :: i, j, z, n
! !CALLING SEQUENCE:
! Called by: PDAF_seek_analysis   (as U_obs_op)
! Called by: PDAF_seik_analysis, PDAF_seik_analysis_newT
! Called by: PDAF_enkf_analysis_rlm, PDAF_enkf_analysis_rsm
!EOP

! *** local variables ***

! CRNS implementation based on Schrön et al. 2017
#if (defined CLMFIVE)
  REAL :: weights_r1(920), weights_r2(920), weights_r3(920)
  Real :: weights_layer(8)
  Integer :: nweights(8)
  Real  :: d86_r1, d86_r2, d86_r3
  REAL :: r1 = 1.0
  REAL :: r2 = 20.0
  REAL :: r3 = 85.0
  REAL :: bd, y 
  REAL :: sum_r1, sum_r2, sum_r3, totw, avesm
! based on W. James Shuttleworth and Rafael Rosolem - January/2012
  REAL :: vwclat = 0.0753  ! Volumetric "lattice" water content (m3/m3)
  REAL :: Nflux = 510.5173790200   ! High energy neutron flux, N (-)
  REAL :: alpha = 0.2392421548 ! Ratio of Fast Neutron Creaton Factor (soil to water), alpha (-)
  REAL :: L1 = 161.98621864    ! High energy soil attenuation length
  REAL :: L2 = 129.14558985    ! High energy water attenuation length
  REAL :: L3 = 107.82204562    ! fast neutron soil attenuation length
  REAL :: L4 = 3.1627190566    ! fast neutron water attenuation length
  REAL :: zdeg, zrad, ideg, costheta, dtheta ! angle variables
  REAL :: h2odens = 1000.0 !Density of water (g/cm3)
  REAL :: pi = 3.14159265359
  REAL :: totflux = 0.0     ! Total flux of above-ground fast neutrons
  REAL :: vwc = 0.0 ! temp. variable for volumetric water content
  REAL, dimension(:), allocatable :: wetsoidens  ! Density of wet soil layer (g/cm3)
  REAL, dimension(:), allocatable :: wetsoimass  ! Mass of wet soil layer (g)
  REAL, dimension(:), allocatable :: isoimass    ! Integrated dry soil mass above layer (g)
  REAL, dimension(:), allocatable :: iwatmass    ! Integrated water mass above layer (g)
  REAL, dimension(:), allocatable :: iwetsoimass ! Integrated wet soil mass above layer (g)
  REAL, dimension(:), allocatable :: hiflux      ! High energy neutron flux
  REAL, dimension(:), allocatable :: fastpot     ! Fast neutron source strength of layer
  REAL, dimension(:), allocatable :: h2oeffdens  ! "Effective" density of water in layer (g/cm3)
  REAL, dimension(:), allocatable :: h2oeffmass  ! "Effective" mass of water in layer (g)
  REAL, dimension(:), allocatable :: ih2oeffmass ! Integrated water mass above layer (g)
  REAL, dimension(:), allocatable :: idegrad     ! Integrated neutron degradation factor (-)
  REAL, dimension(:), allocatable :: fastflux    ! Contribution to above-ground neutron flux
  REAL, dimension(:), allocatable :: normfast    ! Normalized contribution to neutron flux (-) [weighting factors]
  REAL, dimension(:), allocatable :: inormfast   ! Cumulative fraction of neutrons (-)
  
  REAL :: zthick = 0.01 ! Thickness of the crns soil layers has to be consistent with nlayers (cm)
  Integer :: nlayers = 3200 ! total number of crns soil layers
  Integer :: angle, angledz, maxangle ! loop indices for an integration interval
#endif
! *********************************************
! *** Perform application of measurement    ***
! *** operator H on vector orclmcrns_bd matrix column ***
! *********************************************

#if (defined CLMFIVE)
  ! CLMUPDATE_SWC == 3 is for average SWC from Cosmic ray neutron sensor
  if(clmupdate_swc.eq.3) then
    DO i = 1, dim_obs_p
      ! CRNS implementation based on Schrön et al. 2017
      ! First calculate d86 for 3 different radius values
      ! Bulk density average for the 8 considered layers
      IF (clmcrns_bd > 0.0) THEN
        bd = clmcrns_bd
      ELSE
        bd = 0.0
        DO j = 1, 8
          bd = bd + soilstate_inst%bd_col(obs_index_p(i),j) ! bulk density
        END DO
        bd = bd / 8.0 * 0.001 ! average and convert from kg/m^3 to g/cm^3
      ENDIF

      ! CRNS observed value
      y = obs_p(i) ! CRNS observation
      ! Penetration depth calculations D86(bd, r, y)
      d86_r1 = (1/bd*(8.321+0.14249*(0.96655+exp(-0.01*r1))*(20+y)/(0.0429+y))) 
      d86_r2 = (1/bd*(8.321+0.14249*(0.96655+exp(-0.01*r2))*(20+y)/(0.0429+y)))
      d86_r3 = (1/bd*(8.321+0.14249*(0.96655+exp(-0.01*r3))*(20+y)/(0.0429+y)))
      ! Then calculate the weights for thin (1mm) slices of the layers for 85cm
      sum_r1 = 0.0
      sum_r2 = 0.0
      sum_r3 = 0.0
      DO j = 1, 920 ! depth in mm but in calculation used in cm:
        weights_r1(j) = exp(-2*(j/10)/d86_r1)
        weights_r2(j) = exp(-2*(j/10)/d86_r2)
        weights_r3(j) = exp(-2*(j/10)/d86_r3)

        sum_r1 = sum_r1 + weights_r1(j)
        sum_r2 = sum_r2 + weights_r2(j)
        sum_r3 = sum_r3 + weights_r3(j)
      END DO
      ! Normalize the weights:
      weights_r1(:) = weights_r1(:) / sum_r1
      weights_r2(:) = weights_r2(:) / sum_r2
      weights_r3(:) = weights_r3(:) / sum_r3
      ! assign average weights to each layer
      nweights(:) = 0
      weights_layer(:) = 0.0
      ! z index for different layers, manually here, could be done better with model layer depth
      DO j = 1, 920
        IF (j > 680) then 
          z = 8
        ELSEIF (j < 680 .and. j > 480) then 
          z = 7
        ELSEIF (j < 480 .and. j > 320) then 
          z = 6
        ELSEIF (j < 320 .and. j > 200) then 
          z = 5
        ELSEIF (j < 200 .and. j > 120) then 
          z = 4
        ELSEIF (j < 120 .and. j > 60) then 
          z = 3
        ELSEIF (j < 60 .and. j > 20) then 
          z = 2
        ELSEIF (j < 20) then 
          z = 1
        ENDIF
        weights_layer(z) = weights_layer(z) + weights_r1(j) + weights_r2(j) + weights_r3(j)
        nweights(z) = nweights(z) + 1 
      END DO
      ! Normalize the weights
      totw = 0.0
      DO j = 1, 8
        weights_layer(j) = weights_layer(j) / (nweights(j) * 3)
        totw = totw + weights_layer(j)
      END DO
      weights_layer(:) = weights_layer(:) / totw
      ! Finally use the weights to calculate the weighted average of the state variable
      avesm = 0.0
      DO j = 1, 8
        avesm = avesm + weights_layer(j) * state_p(obs_index_p(i) + (j-1)) 
        ! THIS MIGHT NOT WORK FOR MULTI GRIDCELLS
        ! IT assumes that obs_index_p(i) for obs i is the index of 
        ! the first layer of the gridcell where obs i is
      END DO
      ! Assign new average as the state variable
      m_state_p(i) = avesm
    ! end loop over observations
    END DO
 ! CLMUPDATE_SWC == 4 is for assimilating neutron counts directly.
  else if (clmupdate_swc.eq.4) then
    ! COSMIC Operator for CRNS based on W. James Shuttleworth and Rafael Rosolem - January/2012
    ! Init temporary arrays
    allocate(wetsoidens(nlayers),wetsoimass(nlayers),&
             iwetsoimass(nlayers),hiflux(nlayers),fastpot(nlayers),&
             h2oeffdens(nlayers),h2oeffmass(nlayers),ih2oeffmass(nlayers),&
             idegrad(nlayers),fastflux(nlayers),normfast(nlayers),&
             inormfast(nlayers),isoimass(nlayers),iwatmass(nlayers))
    do i = 1,nlayers
      wetsoidens(i)  = 0.0
      wetsoimass(i)  = 0.0
      iwetsoimass(i) = 0.0
      isoimass(i)    = 0.0
      iwatmass(i)    = 0.0
      hiflux(i)      = 0.0
      fastpot(i)     = 0.0
      h2oeffdens(i)  = 0.0
      h2oeffmass(i)  = 0.0
      ih2oeffmass(i) = 0.0
      idegrad(i)     = 0.0
      fastflux(i)    = 0.0
      normfast(i)    = 0.0
      inormfast(i)   = 0.0  
    enddo
    totflux = 0.0
    ! Do the calculations
    DO i = 1, dim_obs_p
      IF (clmcrns_bd > 0.0) THEN
        bd = clmcrns_bd
      ELSE
        bd = 0.0
        DO j = 1, 8
          bd = bd + soilstate_inst%bd_col(obs_index_p(i),j) ! bulk density
        END DO
        bd = bd / 8.0 * 0.001 ! average and convert from kg/m^3 to g/cm^3
      ENDIF

      ! Angle distribution parameters
      ideg = 0.5 ! angle interval 
      angledz = nint(ideg*10.0) ! angle loop
      maxangle = 900 - angledz ! to create integers with no remainder
      dtheta = ideg*(pi/180.0)

      if (real(angledz) /= ideg*10.0) then
        write(*,*) 'ERROR ideg*10.0 must result in an integer - it results in ', ideg*10
      endif

      DO j = 1, nlayers
      ! Volumetric water content for the layer:
        z = 0
        IF (j < 3200 .and. j > 2000) then
          z = 5
        ELSEIF (j < 2000 .and. j > 1200) then
          z = 4
        ELSEIF (j < 1200 .and. j > 600) then
          z = 3
        ELSEIF (j < 600 .and. j > 200) then
          z = 2
        ELSEIF (j < 200) then
          z = 1
        ENDIF

        vwc = state_p(obs_index_p(i) + (z-1))
        ! THIS MIGHT NOT WORK FOR MULTI GRIDCELLS
        ! IT assumes that obs_index_p(i) for obs i is the index of 
        ! the first layer of the gridcell where obs i is

        ! High energy neutron downward flux
        ! The integration is now performed at the node of each layer (i.e., center of the layer)
        h2oeffdens(j) = ((vwc+vwclat)*h2odens)/1000.0
        IF (j > 1) THEN
          ! Assuming an area of 1 cm2
          isoimass(j) = isoimass(j-1) + bd*(0.5*zthick)*1.0 &
                                      + bd*(0.5*zthick)*1.0
           
          ! Assuming an area of 1 cm2
          iwatmass(j) = iwatmass(j-1) + h2oeffdens(j-1)*(0.5*zthick)*1.0 &
                                      + h2oeffdens(j)*(0.5*zthick)*1.0
        ELSE
          ! Assuming an area of 1 cm2
          isoimass(j) = bd*(0.5*zthick)*1.0 
          iwatmass(j) = h2oeffdens(j)*(0.5*zthick)*1.0
        ENDIF

        IF (clmcrns_nflux > 0.0) THEN
          hiflux(j) = clmcrns_nflux*exp(-(isoimass(j)/L1 + iwatmass(j)/L2))
        ELSE
          hiflux(j)  = Nflux*exp(-(isoimass(j)/L1 + iwatmass(j)/L2) )
        ENDIF

        fastpot(j) = zthick*hiflux(j)*(alpha*bd + h2oeffdens(j))
    
        ! This second loop needs to be done for the distribution of angles for fast neutron release
        ! the intent is to loop from 0 to 89.5 by 0.5 degrees - or similar.
        ! Because Fortran loop indices are integers, we have to divide the indices by 10 - you get the idea.  
    
        do angle=0,maxangle,angledz
          zdeg     = real(angle)/10.0   ! 0.0  0.5  1.0  1.5 ...
          zrad     = (zdeg*pi)/180.0
          costheta = cos(zrad)
    
          ! Angle-dependent low energy (fast) neutron upward flux
          fastflux(j) = fastflux(j) + fastpot(j)*exp(-(isoimass(j)/L3 + iwatmass(j)/L4)/costheta)*dtheta
        enddo
        ! After contribution from all directions are taken into account,
        ! need to multiply fastflux by 2/PI
        fastflux(j) = (2.0/pi)*fastflux(j)
    
        ! Low energy (fast) neutron upward flux
        totflux = totflux + fastflux(j)
        ENDDO ! End of totflux calculations
        ! Assign totflux for the observation to the transformed state vector
        m_state_p(i) = totflux
    END DO ! End for loop over observations 

    deallocate(wetsoidens, wetsoimass, iwetsoimass, hiflux,&
               fastpot, h2oeffdens, h2oeffmass, ih2oeffmass, idegrad, fastflux,&
               normfast, inormfast, isoimass, iwatmass)
  else
    DO i = 1, dim_obs_p
      m_state_p(i) = state_p(obs_index_p(i))
    END DO
  endif
#else
  DO i = 1, dim_obs_p
    m_state_p(i) = state_p(obs_index_p(i))
  END DO
#endif
END SUBROUTINE obs_op_pdaf
