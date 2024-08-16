!   Copyright 2016 University of Auckland.

!   This file is part of Waiwera.

!   Waiwera is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.

!   Waiwera is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU Lesser General Public License for more details.

!   You should have received a copy of the GNU Lesser General Public License
!   along with Waiwera.  If not, see <http://www.gnu.org/licenses/>.

module IFC67_module
  !!  IFC-67 industrial thermodynamic formulation, as described by:
  !!
  !! **"A formulation of the thermodynamic properties of ordinary water substance"**, *International Formulation
  !! Committee*, D&uuml;sseldorf, Germany, 1967.

  ! These are higher-speed versions developed by Mike O'Sullivan for AUTOUGH2.

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use thermodynamics_module

  implicit none
  private

!------------------------------------------------------------------------
! Constants
!------------------------------------------------------------------------

  type(critical_point_type), parameter :: critical = critical_point_type( &
       647.3_dp, 374.15_dp, 22.12e6_dp, 322.0_dp)

!------------------------------------------------------------------------
! Saturation curve type
!------------------------------------------------------------------------

  type, public, extends(saturation_type) :: IFC67_saturation_type
     !! IFC-67 saturation curve calculations.
     private
     PetscReal :: &
          A1 = -7.691234564_dp,   A2 = -2.608023696e1_dp, &
          A3 = -1.681706546e2_dp, A4 =  6.423285504e1_dp, &
          A5 = -1.189646225e2_dp, A6 =  4.167117320_dp, &
          A7 =  2.097506760e1_dp, A8 =  1.0e9_dp, A9 = 6.0_dp
     contains
       private
       procedure, public :: temperature => saturation_temperature
       procedure, public :: pressure => saturation_pressure
  end type IFC67_saturation_type

  !------------------------------------------------------------------------
  ! Region 1 (liquid water) type
  !------------------------------------------------------------------------

  type, public, extends(region_type) :: IFC67_region1_type
     !! IFC-67 region 1 (pure water) type.
     private
     PetscReal :: &
          A1 = 6.824687741e3_dp,    A2 = -5.422063673e2_dp, &
          A3 = -2.096666205e4_dp,   A4 = 3.941286787e4_dp, &
          A5 = -13.466555478e4_dp,  A6 = 29.707143084e4_dp, &
          A7 = -4.375647096e5_dp,   A8 = 42.954208335e4_dp, &
          A9 = -27.067012452e4_dp,  A10 = 9.926972482e4_dp, &
          A11 = -16.138168904e3_dp, A12 = 7.982692717_dp, &
          A13 = -2.616571843e-2_dp, A14 = 1.522411790e-3_dp, &
          A15 = 2.284279054e-2_dp,  A16 = 2.421647003e2_dp, &
          A17 = 1.269716088e-10_dp, A18 = 2.074838328e-7_dp, &
          A19 = 2.174020350e-8_dp,  A20 = 1.105710498e-9_dp, &
          A21 = 1.293441934e1_dp,   A22 = 1.308119072e-5_dp, &
          A23 = 6.047626338e-14_dp
     PetscReal :: &
          SA1 = 8.438375405e-1_dp,  SA2 = 5.362162162e-4_dp, &
          SA3 = 1.72_dp,            SA4 = 7.342278489e-2_dp, &
          SA5 = 4.975858870e-2_dp,  SA6 = 6.537154300e-1_dp, &
          SA7 = 1.150e-6_dp,        SA8 = 1.51080e-5_dp, &
          SA9 = 1.41880e-1_dp,      SA10 = 7.002753165_dp, &
          SA11 = 2.995284926e-4_dp, SA12 = 2.040e-1_dp
     PetscReal :: max_temperature
   contains
     private
     procedure, public :: init => region1_init
     procedure, public :: destroy => region1_destroy
     procedure, public :: properties => region1_properties
     procedure, public :: viscosity => region1_viscosity
  end type IFC67_region1_type

!------------------------------------------------------------------------
! Region 2 (steam) type
!------------------------------------------------------------------------

  type, public, extends(region_type) :: IFC67_region2_type
     !! IFC-67 region 2 (steam) type.
     private
     PetscReal :: &
          B0 = 16.83599274_dp,      B01 = 28.56067796_dp, &
          B02 = -54.38923329_dp,    B03 = 0.4330662834_dp, &
          B04 = -0.6547711697_dp,   B05 = 8.565182058e-2_dp, &
          B11 = 6.670375918e-2_dp,  B12 = 1.388983801_dp, &
          B21 = 8.390104328e-2_dp,  B22 = 2.614670893e-2_dp, &
          B23 = -3.373439453e-2_dp, B31 = 4.520918904e-1_dp, &
          B32 = 1.069036614e-1_dp,  B41 = -5.975336707e-1_dp, &
          B42 = -8.847535804e-2_dp, B51 = 5.958051609e-1_dp, &
          B52 = -5.159303373e-1_dp, B53 = 2.075021122e-1_dp, &
          B61 = 1.190610271e-1_dp,  B62 = -9.867174132e-2_dp, &
          B71 = 1.683998803e-1_dp,  B72 = -5.809438001e-2_dp, &
          B81 = 6.552390126e-3_dp,  B82 = 5.710218649e-4_dp, &
          B90 = 1.936587558e2_dp,   B91 = -1.388522425e3_dp, &
          B92 = 4.126607219e3_dp,   B93 = -6.508211677e3_dp, &
          B94 = 5.745984054e3_dp,   B95 = -2.693088365e3_dp, &
          B96 = 5.235718623e2_dp
     PetscReal :: &
          SB = 7.633333333e-1_dp, SB61 = 4.006073948e-1_dp, &
          SB71 = 8.636081627e-2_dp, SB81 = -8.532322921e-1_dp, &
          SB82 = 3.460208861e-1_dp
   contains
     private
     procedure, public :: init => region2_init
     procedure, public :: destroy => region2_destroy
     procedure, public :: properties => region2_properties
     procedure, public :: viscosity => region2_viscosity
  end type IFC67_region2_type

!------------------------------------------------------------------------
! IFC-67 thermodynamics type
!------------------------------------------------------------------------

  type, extends(thermodynamics_type), public :: IFC67_type
     !! IFC-67 thermodynamics type.
     private
   contains
     private
     procedure, public :: init => IFC67_init
     procedure, public :: destroy => IFC67_destroy
     procedure, public :: phase_composition => IFC67_phase_composition
  end type IFC67_type

!------------------------------------------------------------------------

contains

!------------------------------------------------------------------------
! IFC67 class
!------------------------------------------------------------------------

  subroutine IFC67_init(self, json, logfile)
    !! Constructs IFC-67 thermodynamics object.

    use fson
    use fson_value_m, only: TYPE_OBJECT
    use fson_mpi_module
    use logfile_module

    class(IFC67_type), intent(in out) :: self
    type(fson_value), pointer, intent(in), optional :: json !! JSON input object
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    PetscInt :: i, thermo_type
    PetscBool :: defaults
    PetscBool, parameter :: default_extrapolate = PETSC_FALSE

    self%name = 'IFC-67'

    self%max_temperature = 800._dp
    self%max_pressure = 100.e6_dp
    self%critical = critical

    allocate(IFC67_saturation_type :: self%saturation)
    call self%saturation%init(self)

    self%num_regions = 2
    allocate(IFC67_region1_type :: self%water)
    allocate(IFC67_region2_type :: self%steam)
    allocate(self%region(self%num_regions))

    call self%region(1)%set(self%water)
    call self%region(2)%set(self%steam)

    defaults = PETSC_TRUE
    if (present(json)) then
       if (fson_has_mpi(json, "thermodynamics")) then
          thermo_type = fson_type_mpi(json, "thermodynamics")
          if (thermo_type == TYPE_OBJECT) then
             call fson_get_mpi(json, "thermodynamics.extrapolate", &
                  default_extrapolate, self%extrapolate, logfile)
             defaults = PETSC_FALSE
          end if
       end if
    end if

    if (defaults) then
       self%extrapolate = default_extrapolate
       if (present(logfile)) then
          call logfile%write(LOG_LEVEL_INFO, 'input', 'default', &
               logical_keys = ['thermodynamics.extrapolate'], &
               logical_values = [default_extrapolate])
       end if
    end if

    do i = 1, self%num_regions
       call self%region(i)%ptr%init(self)
    end do

  end subroutine IFC67_init

!------------------------------------------------------------------------

  subroutine IFC67_destroy(self)
    !! Destroys IFC-67 thermodynamics object.

    class(IFC67_type), intent(in out) :: self
    ! Locals:
    PetscInt :: i

    do i = 1, self%num_regions
       call self%region(i)%ptr%destroy()
    end do

    deallocate(self%region)
    deallocate(self%water, self%steam)
    deallocate(self%saturation)

  end subroutine IFC67_destroy

!------------------------------------------------------------------------

  PetscInt function IFC67_phase_composition(self, region, pressure, &
       temperature) result(phases)
    !! Returns phase composition integer for given region, pressure
    !! and temperature. Here the bits represent:
    !! 0: liquid
    !! 1: vapour

    class(IFC67_type), intent(in) :: self
    PetscInt, intent(in) :: region
    PetscReal, intent(in) :: pressure, temperature

    select case(region)
    case(1)   ! liquid water
       phases = int(b'01')
    case(2)   ! dry steam
       phases = int(b'10')
    case(4)   ! two-phase
       phases = int(b'11')
    case default
       phases = int(b'00')
    end select

  end function IFC67_phase_composition

!------------------------------------------------------------------------
! Region 1 (liquid water)
!------------------------------------------------------------------------

  subroutine region1_init(self, thermo)
    !! Initializes IFC-67 region 1 object.

    class(IFC67_region1_type), intent(in out) :: self
    class(thermodynamics_type), intent(in), target :: thermo
    ! Locals:
    PetscReal, parameter :: default_max_temperature = 350._dp
    PetscReal, parameter :: extrapolated_max_temperature = 360._dp

    self%name = 'water'
    self%thermo => thermo

    if (thermo%extrapolate) then
       self%max_temperature = extrapolated_max_temperature
    else
       self%max_temperature = default_max_temperature
    end if

  end subroutine region1_init

!------------------------------------------------------------------------

  subroutine region1_destroy(self)

    !! Destroys IFC-67 region 1 object.

    class(IFC67_region1_type), intent(in out) :: self

  end subroutine region1_destroy

!------------------------------------------------------------------------

  subroutine region1_properties(self, param, props, err)
    !! Calculates density and internal energy of liquid water as a function of
    !! pressure (Pa) and temperature (deg C).
    !!
    !! Returns err = 1 if called outside its operating range (t<=self%max_temperature, p<=100 MPa).

    class(IFC67_region1_type), intent(in out) :: self
    PetscReal, intent(in) :: param(:) !! Primary variables (pressure, temperature)
    PetscReal, intent(out):: props(:) !! (density, internal energy)
    PetscInt, intent(out) :: err !! error code
    ! Locals:
    PetscReal :: AA1, BB1, BB2, CC1, CC2, CC4, CC8, CC10
    PetscReal :: CZ, DD1, DD2, DD4, EE1, EE3, ENTR, H
    PetscReal :: PAR1, PAR2, PAR3, PAR4, PAR5
    PetscReal :: PNMR, PNMR2, PNMR3, PNMR4, PRT1, PRT2, PRT3, PRT4, PRT5
    PetscReal :: SNUM, TKR, TKR2, TKR3, TKR4, TKR5, TKR6, TKR7, TKR8, TKR9, TKR10, &
         TKR11, TKR18, TKR19, TKR20
    PetscReal :: V, VMKR, Y, YD, Z, ZP, D, U
    
    associate (p => param(1), t => param(2))

      if ((t <= self%max_temperature).and.(p <= self%thermo%max_pressure)) then

         TKR = (t + tc_k) / self%thermo%critical%temperature_k
         TKR2 = TKR * TKR
         TKR3 = TKR * TKR2
         TKR4 = TKR2 * TKR2
         TKR5 = TKR2 * TKR3
         TKR6 = TKR4 * TKR2
         TKR7 = TKR4 * TKR3
         TKR8 = TKR4 * TKR4
         TKR9 = TKR4 * TKR5
         TKR10 = TKR4 * TKR6
         TKR11 = TKR * TKR10
         TKR18 = TKR8 * TKR10
         TKR19 = TKR8 * TKR11
         TKR20 = TKR10 * TKR10
         PNMR = p / self%thermo%critical%pressure
         PNMR2 = PNMR * PNMR
         PNMR3 = PNMR * PNMR2
         PNMR4 = PNMR * PNMR3
         Y = 1.0_dp - self%SA1 * TKR2 - self%SA2 / TKR6
         ZP = self%SA3 * Y * Y - 2.0_dp * self%SA4 * TKR + &
              2.0_dp * self%SA5 * PNMR
         if (ZP >= 0._dp) then
            Z = Y + sqrt(ZP)
            CZ = Z ** (5.0_dp / 17.0_dp)
            PAR1 = self%A12 * self%SA5 / CZ
            CC1 = self%SA6 - TKR
            CC2 = CC1 * CC1
            CC4 = CC2 * CC2
            CC8 = CC4 * CC4
            CC10 = CC2 * CC8
            AA1 = self%SA7 + TKR19
            PAR2 = self%A13 + self%A14 * TKR + self%A15 * TKR2 + &
                 self%A16 * CC10 + self%A17 / AA1
            PAR3 = (self%A18 + 2._dp * self%A19 * PNMR + &
                 3._dp * self%A20 * PNMR2) / (self%SA8 + TKR11)
            DD1 = self%SA10 + PNMR
            DD2 = DD1 * DD1
            DD4 = DD2 * DD2
            PAR4 = self%A21 * TKR18 * (self%SA9 + TKR2) * &
                 (-3.0_dp / DD4 + self%SA11)
            PAR5 = 3.0_dp * self%A22 * (self%SA12 - TKR) * &
                 PNMR2 + 4.0_dp * self%A23 / TKR20 * PNMR3
            VMKR = PAR1 + PAR2 - PAR3 - PAR4 + PAR5
            V = VMKR * 3.17e-3_dp
            D = 1.0_dp / V
            YD = -2.0_dp * self%SA1 * TKR + 6.0_dp * self%SA2 / TKR7
            SNUM =  self%A10 + self%A11 * TKR
            SNUM = SNUM * TKR + self%A9
            SNUM = SNUM * TKR + self%A8
            SNUM = SNUM * TKR + self%A7
            SNUM = SNUM * TKR + self%A6
            SNUM = SNUM * TKR + self%A5
            SNUM = SNUM * TKR + self%A4
            SNUM = SNUM * TKR2 - self%A2
            PRT1 = self%A12 * (Z * (17.0_dp * (Z / 29.0_dp - Y / 12.0_dp) + &
                 5.0_dp * TKR * YD / 12.0_dp) + &
                 self%SA4 * TKR - (self%SA3 - 1.0_dp) * TKR * Y * YD) / CZ
            PRT2 = PNMR * (self%A13 - self%A15 * TKR2 + self%A16 * &
                 (9.0_dp * TKR + self%SA6) * CC8 * CC1 + &
                 self%A17 * (19.0_dp * TKR19 + AA1) / (AA1 * AA1))
            BB1 = self%SA8 + TKR11
            BB2 = BB1 * BB1
            PRT3 = (11.0_dp * TKR11 + BB1) / BB2 * (self%A18 * PNMR + &
                 self%A19 * PNMR2 + self%A20 * PNMR3)
            EE1 = self%SA10 + PNMR
            EE3 = EE1 * EE1 * EE1
            PRT4 = self%A21 * TKR18 * (17.0_dp * self%SA9 + 19.0_dp * TKR2) * &
                 (1.0_dp / EE3 + self%SA11 * PNMR)
            PRT5 = self%A22 * self%SA12 * PNMR3 + 21.0_dp * self%A23 / TKR20 * PNMR4
            ENTR = self%A1 * TKR - SNUM + PRT1 + PRT2 - PRT3 + PRT4 + PRT5
            H = ENTR * 70120.4_dp
            U = H - p * V

            props(1) = D ! density
            props(2) = U ! internal energy
            err = 0

         else
            err = 1
         end if
      else
         err = 1
      end if

    end associate

  end subroutine region1_properties

!------------------------------------------------------------------------

  subroutine region1_viscosity(self, temperature, pressure, density, viscosity)
    !! Calculates water viscosity. Density is a dummy argument and is not used.

    class(IFC67_region1_type), intent(in out) :: self
    PetscReal, intent(in) :: temperature  !! Fluid temperature (\(^\circ C\))
    PetscReal, intent(in) :: pressure     !! Fluid pressure (not used)
    PetscReal, intent(in) :: density      !! Fluid density (\(kg. m^{-3}\))
    PetscReal, intent(out) :: viscosity   !! Viscosity (\(kg.m^{-1}.s^{-1}\))
    ! Locals:
    PetscReal :: ex, phi, am, ps
    PetscInt :: err

    ex = 247.8_dp / (temperature + 133.15_dp)
    phi = 1.0467_dp * (temperature - 31.85_dp)
    call self%thermo%saturation%pressure(temperature, ps, err)
    am = 1.0_dp + phi * (pressure - ps) * 1.0e-11_dp
    viscosity = 1.0e-7_dp * am * 241.4_dp * 10.0_dp ** ex

  end subroutine region1_viscosity

!------------------------------------------------------------------------
! Region 2 (steam)
!------------------------------------------------------------------------

  subroutine region2_init(self, thermo)
    !! Initializes IFC-67 region 2 object.

    class(IFC67_region2_type), intent(in out) :: self
    class(thermodynamics_type), intent(in), target :: thermo

    self%name = 'steam'
    self%thermo => thermo

  end subroutine region2_init

!------------------------------------------------------------------------

  subroutine region2_destroy(self)
    !! Destroys IFC-67 region 2 object.

    class(IFC67_region2_type), intent(in out) :: self

  end subroutine region2_destroy

!------------------------------------------------------------------------

  subroutine region2_properties(self, param, props, err)
    !! Calculates density and internal energy of dry steam as a function of
    !! pressure (Pa) and temperature (deg C).
    !!
    !! Returns err = 1 if called outside its operating range (t<=800 deg C, p<=100 MPa).

    class(IFC67_region2_type), intent(in out) :: self
    PetscReal, intent(in) :: param(:) !! Primary variables (pressure, temperature)
    PetscReal, intent(out):: props(:)  !! (density, internal energy)
    PetscInt, intent(out) :: err !! error code
    ! Locals:
    PetscReal :: BETA, BETA2, BETA3, BETA4, BETA5, BETA6, BETA7, BETAL, CHI2
    PetscReal :: DBETAL, EPS2, H, OS1, OS2, OS5, OS6, OS7, R, R2, R4, R6, R10, RI1
    PetscReal :: SC, SD1, SD12, SD2, SD22, SD3, SD32, SN, SN6, SN7, SN8
    PetscReal ::  THETA, THETA2, THETA3, THETA4, V, D, U
    PetscReal :: X, X2, X3, X4, X5, X6, X8, X10, X11, X14, X18, X19, X24, X27

    associate (p => param(1), t => param(2))

      ! Check input:
      if ((t <= self%thermo%max_temperature).and.(p <= self%thermo%max_pressure)) then

         THETA = (T + tc_k) / self%thermo%critical%temperature_k
         BETA = P / self%thermo%critical%pressure
         RI1 = 4.260321148_dp
         X = exp(self%SB * (1.0_dp - THETA))
         X2 = X * X
         X3 = X2 * X
         X4 = X3 * X
         X5 = X4 * X
         X6 = X5 * X
         X8 = X6 * X2
         X10 = X6 * X4
         X11 = X10 * X
         X14 = X10 * X4
         X18 = X14 * X4
         X19 = X18 * X
         X24 = X18 * X6
         X27 = X24 * X3

         THETA2 = THETA * THETA
         THETA3 = THETA2 * THETA
         THETA4 = THETA3 * THETA

         BETA2 = BETA * BETA
         BETA3 = BETA2 * BETA
         BETA4 = BETA3 * BETA
         BETA5 = BETA4 * BETA
         BETA6 = BETA5 * BETA
         BETA7 = BETA6 * BETA

         BETAL = 15.74373327_dp - 34.17061978_dp * THETA + 19.31380707_dp * THETA2
         DBETAL = -34.17061978_dp + 38.62761414_dp * THETA
         R = BETA / BETAL
         R2 = R * R
         R4 = R2 * R2
         R6 = R4 * R2
         R10 = R6 * R4

         CHI2 = RI1 * THETA / BETA
         SC = (self%B11 * X10 + self%B12) * X3
         CHI2 = CHI2 - SC
         SC = self%B21 * X18 + self%B22 * X2 + self%B23 * X
         CHI2 = CHI2 - 2._dp * BETA * SC
         SC = (self%B31 * X8 + self%B32) * X10
         CHI2 = CHI2 - 3._dp * BETA2 * SC
         SC = (self%B41 * X11 + self%B42) * X14
         CHI2 = CHI2 - 4._dp * BETA3 * SC
         SC = (self%B51 * X8 + self%B52 * X4 + self%B53) * X24
         CHI2 = CHI2 - 5._dp * BETA4 * SC

         SD1 = 1.0_dp / BETA4 + self%SB61 * X14
         SD2 = 1.0_dp / BETA5 + self%SB71 * X19
         SD3 = 1.0_dp / BETA6 + (self%SB81 * X27 + self%SB82) * X27
         SD12 = SD1 * SD1
         SD22 = SD2 * SD2
         SD32 = SD3 * SD3

         SN = (self%B61 * X + self%B62) * X11
         CHI2 = CHI2 - SN / SD12 * 4._dp / BETA5
         SN = (self%B71 * X6 + self%B72) * X18
         CHI2 = CHI2 - SN / SD22 * 5._dp / BETA6
         SN = (self%B81 * X10 + self%B82) * X14
         CHI2 = CHI2 - SN / SD32 * 6._dp / BETA7
         SC = self%B96
         SC = SC * X + self%B95
         SC = SC * X + self%B94
         SC = SC * X + self%B93
         SC = SC * X + self%B92
         SC = SC * X + self%B91
         SC = SC * X + self%B90
         CHI2 = CHI2 + 11.0_dp * R10 * SC
         V = CHI2 * 0.00317_dp
         D = 1.0_dp / V

         OS1 = self%SB * THETA
         EPS2 = self%B0 * THETA - (-self%B01 + self%B03 * THETA2 + &
              2.0_dp * self%B04 * THETA3 + 3.0_dp * self%B05 * THETA4)
         SC = (self%B11 * (1.0_dp + 13.0_dp * OS1) * X10 + self%B12 * &
              (1.0_dp + 3.0_dp * OS1)) * X3
         EPS2 = EPS2 - BETA * SC
         SC = self%B21 * (1.0_dp + 18.0_dp * OS1) * X18 + self%B22 * &
              (1.0_dp + 2.0_dp * OS1) * X2 + self%B23 * (1.0_dp + OS1) * X
         EPS2 = EPS2 - BETA2 * SC
         SC = (self%B31 * (1.0_dp + 18.0_dp * OS1) * X8 + self%B32 * &
              (1.0_dp + 10.0_dp * OS1)) * X10
         EPS2 = EPS2 - BETA3 * SC
         SC = (self%B41 * (1.0_dp + 25.0_dp * OS1) * X11 + self%B42 * &
              (1.0_dp + 14.0_dp * OS1)) * X14
         EPS2 = EPS2 - BETA4 * SC
         SC = (self%B51 * (1.0_dp + 32.0_dp * OS1) * X8 + self%B52 * &
              (1.0_dp + 28.0_dp * OS1) * X4 + self%B53 * &
              (1.0_dp + 24.0_dp * OS1)) * X24
         EPS2 = EPS2 - BETA5 * SC

         SN6 = 14.0_dp * self%SB61 * X14
         SN7 = 19.0_dp * self%SB71 * X19
         SN8 = (54.0_dp * self%SB81 * X27 + 27.0_dp * self%SB82) * X27
         OS5 = 1.0_dp + 11.0_dp * OS1 - OS1 * SN6 / SD1
         SC = (self%B61 * X * (OS1 + OS5) + self%B62 * OS5) * (X11 / SD1)
         EPS2 = EPS2 - SC
         OS6 = 1.0_dp + 24.0_dp * OS1 - OS1 * SN7 / SD2
         SC = (self%B71 * X6 * OS6 + self%B72 * (OS6 - 6.0_dp * OS1)) * &
              (X18 / SD2)
         EPS2 = EPS2 - SC
         OS7 = 1.0_dp + 24.0_dp * OS1 - OS1 * SN8 / SD3
         SC = (self%B81 * X10 * OS7 + self%B82 * (OS7 - 10.0_dp *  OS1)) * &
              (X14 / SD3)
         EPS2 = EPS2 - SC
         OS2 = 1.0_dp + THETA * 10.0_dp * DBETAL / BETAL
         SC = (OS2 + 6.0_dp * OS1) * self%B96
         SC = SC * X + (OS2 + 5.0_dp * OS1) * self%B95
         SC = SC * X + (OS2 + 4.0_dp * OS1) * self%B94
         SC = SC * X + (OS2 + 3.0_dp * OS1) * self%B93
         SC = SC * X + (OS2 + 2.0_dp * OS1) * self%B92
         SC = SC * X + (OS2 + OS1) * self%B91
         SC = SC * X + OS2 * self%B90
         EPS2 = EPS2 + BETA * R10 * SC
         H = EPS2 * 70120.4_dp
         U = H - P * V

         props(1) = D ! density
         props(2) = U ! internal energy
         err = 0

      else
         err = 1
      end if

    end associate

  end subroutine region2_properties

!------------------------------------------------------------------------

  subroutine region2_viscosity(self, temperature, pressure, density, viscosity)
    !! Calculates water viscosity. Pressure is a dummy argument and is not used.

    class(IFC67_region2_type), intent(in out) :: self
    PetscReal, intent(in) :: temperature !! Fluid temperature (\(^\circ C\))
    PetscReal, intent(in) :: pressure  !! Fluid pressure (not used)
    PetscReal, intent(in) :: density   !! Fluid density (\(kg. m^{-3}\))
    PetscReal, intent(out) :: viscosity !! Viscosity (\(kg.m^{-1}.s^{-1}\))
    ! Locals:
    PetscReal :: v1

    v1 = 0.407_dp * temperature + 80.4_dp
    if (temperature <= 350.0_dp) then
       viscosity = 1.0e-7_dp * (v1 - density * &
       (1858.0_dp - 5.9_dp * temperature) * 1.0e-3_dp)
    else
        viscosity = 1.0e-7_dp * (v1 + density * &
        (0.353_dp + density * (676.5e-6_dp + density * 102.1e-9_dp)))
    end if

  end subroutine region2_viscosity

!------------------------------------------------------------------------
  ! Saturation curve
!------------------------------------------------------------------------

subroutine saturation_pressure(self, t, p, err)
  !! Calculates saturation pressure as a function of temperature.
  !! Returns err = 1 if called outside its operating range (1 <= t <= critical temperature).

  class(IFC67_saturation_type), intent(in) :: self
  PetscReal, intent(in) :: t  !! Fluid temperature (\(^\circ C\))
  PetscReal, intent(out):: p  !! Fluid pressure (\(kg. m. s^{-1}\))
  PetscInt, intent(out) :: err  !! Error code
  ! Locals:
  PetscReal :: PC, SC, TC, X1, X2

  if ((t >= 1._dp).and.(t <= self%thermo%critical%temperature)) then
     TC = (t + tc_k) / self%thermo%critical%temperature_k
     X1 = 1._dp - TC
     X2 = X1 * X1
     SC = self%A5 * X1 + self%A4
     SC = SC * X1 + self%A3
     SC = SC * X1 + self%A2
     SC = SC * X1 + self%A1
     SC = SC * X1
     PC = exp(SC / (TC * (1._dp + self%A6 * X1 + self%A7 * X2)) - X1 / (self%A8 * X2 + self%A9))
     p = PC * self%thermo%critical%pressure
     err = 0
  else
     err = 1
  end if

end subroutine saturation_pressure

!------------------------------------------------------------------------

subroutine saturation_temperature(self, p, t, err)
  !! Calculates saturation temperature (deg C) as a function of pressure.
  !! Returns err = 1 if called outside its operating range (611.213 Pa <= p <= critical pressure).

  use utils_module, only: newton1d

  class(IFC67_saturation_type), intent(in) :: self
  PetscReal, intent(in) :: p  !! Fluid pressure (\(kg. m. s^{-1}\))
  PetscReal, intent(out):: t  !! Fluid temperature (\(^\circ C\))
  PetscInt, intent(out) :: err !! Error code
  ! Locals:
  PetscInt, parameter :: maxit = 200
  PetscReal, parameter :: ftol = 1.e-10_dp, xtol = 1.e-10_dp
  PetscReal, parameter :: inc = 1.e-8_dp

  if ((p >= 0.0061e5_dp) .and. (p <= self%thermo%critical%pressure)) then

     ! Initial estimate:
     t = max(4606.0_dp / (24.02_dp - dlog(p)) - tc_k, 5._dp)

     call newton1d(f, t, ftol * p, xtol, maxit, inc, err)

  else
     err = 1
  end if

contains

  PetscReal function f(x, err)
    PetscReal, intent(in) :: x
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscReal :: ps

    call self%pressure(x, ps, err)
    f = p - ps

  end function f

end subroutine saturation_temperature

!------------------------------------------------------------------------

end module IFC67_module
