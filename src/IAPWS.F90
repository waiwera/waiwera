module IAPWS_module

  !! IAPWS-97 industrial thermodynamic formulation, as described by:
  !!
  !! Wagner, W., Cooper, J.R., Dittman, A., Kijima, J., Kretzschmar, H.-J., Kruse, A., Mares, R., 
  !! Oguchi, K., Sato, H., Stocker, I., Sifner, O., Takaishi, Y., Tanishita, I., Trubenbach, J. and
  !! Willkommen, Th., 1997.  **The IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of
  !! Water and Steam.**  *Trans. ASME* 150(122), 150-182.
  !!
  !! The viscosity function is described by:
  !! IAPWS, 2008.  **Release on the IAPWS Formulation 2008 for the Viscosity of Ordinary Water Substance.**

  use kinds_module
  use powertable_module
  use thermodynamics_module

  implicit none
  private

#include <petsc/finclude/petscsys.h>

!------------------------------------------------------------------------
! Saturation curve type
!------------------------------------------------------------------------

  type, public, extends(saturation_type) :: IAPWS_saturation_type
     !! IAPWS-97 saturation curve calculations.
     private
     PetscReal :: pstar = 1.0e6_dp
     PetscReal :: n(10) = [ &
           0.11670521452767e4_dp, -0.72421316703206e6_dp, -0.17073846940092e2_dp,  &
           0.12020824702470e5_dp, -0.32325550322333e7_dp,  0.14915108613530e2_dp,  &
          -0.48232657361591e4_dp,  0.40511340542057e6_dp, -0.23855557567849_dp, &
           0.65017534844798e3_dp]
     contains
       private
       procedure, public :: temperature => saturation_temperature
       procedure, public :: pressure => saturation_pressure
  end type IAPWS_saturation_type

!------------------------------------------------------------------------
! Viscosity type
!------------------------------------------------------------------------

  type :: IAPWS_viscosity_type
     !! IAPWS viscosity type.
     private
     PetscReal :: mustar = 1.0e-6_dp
     PetscReal :: h0(0:3) = [ &
          1.67752_dp, 2.20462_dp, 0.6366564_dp, -0.241605_dp]
     PetscReal :: h1(21) = [ &
           5.20094e-1_dp,  8.50895e-2_dp, -1.08374_dp,    -2.89555e-1_dp,  2.22531e-1_dp,  &
           9.99115e-1_dp,  1.88797_dp,     1.26613_dp,     1.20573e-1_dp, -2.81378e-1_dp,  &
          -9.06851e-1_dp, -7.72479e-1_dp, -4.89837e-1_dp, -2.57040e-1_dp,  1.61913e-1_dp, &
           2.57399e-1_dp, -3.25372e-2_dp,  6.98452e-2_dp,  8.72102e-3_dp, -4.35673e-3_dp,  &
          -5.93264e-4_dp]
     PetscInt :: I(21) = [0,1,2,3,0,1,2,3,5,0,1,2,3,4,0,1,0,3,4,3,5]
     PetscInt :: J(21) = [0,0,0,0,1,1,1,1,1,2,2,2,2,2,3,3,4,4,5,6,6]
     PetscInt :: K(4)  = [0,1,2,3]
     type(powertable_type) :: pi, pj, pk
   contains
     private
     procedure, public :: init => viscosity_init
     procedure, public :: destroy => viscosity_destroy
  end type IAPWS_viscosity_type

!------------------------------------------------------------------------
! IAPWS region type
!------------------------------------------------------------------------

  type, extends(region_type) :: IAPWS_region_type
     !! IAPWS-97 region type- just implements viscosity method, which is common
     !! to all regions
     private
     type(IAPWS_viscosity_type), public :: visc
   contains
     private
     procedure, public :: init => region_init
     procedure, public :: destroy => region_destroy
     procedure, public :: properties => region_properties
     procedure, public :: viscosity => region_viscosity
  end type IAPWS_region_type

!------------------------------------------------------------------------
  ! Region 1 (liquid water) type
!------------------------------------------------------------------------

  type, public, extends(IAPWS_region_type) :: IAPWS_region1_type
     !! IAPWS-97 region 1 (liquid water) type.
     private
     PetscReal :: pstar = 16.53e6_dp, tstar = 1386._dp
     PetscReal :: n(34) = [ &                
           0.14632971213167_dp,     -0.84548187169114_dp,    -0.37563603672040e1_dp,   &
           0.33855169168385e1_dp,   -0.95791963387872_dp,     0.15772038513228_dp,     &
          -0.16616417199501e-1_dp,   0.81214629983568e-3_dp,  0.28319080123804e-3_dp,  &
          -0.60706301565874e-3_dp,  -0.18990068218419e-1_dp, -0.32529748770505e-1_dp,  &
          -0.21841717175414e-1_dp,  -0.52838357969930e-4_dp, -0.47184321073267e-3_dp,  &
          -0.30001780793026e-3_dp,   0.47661393906987e-4_dp, -0.44141845330846e-5_dp,  &
          -0.72694996297594e-15_dp, -0.31679644845054e-4_dp, -0.28270797985312e-5_dp,  &
          -0.85205128120103e-9_dp,  -0.22425281908000e-5_dp, -0.65171222895601e-6_dp,  &
          -0.14341729937924e-12_dp, -0.40516996860117e-6_dp, -0.12734301741641e-8_dp,  &
          -0.17424871230634e-9_dp,  -0.68762131295531e-18_dp, 0.14478307828521e-19_dp, &
           0.26335781662795e-22_dp, -0.11947622640071e-22_dp, 0.18228094581404e-23_dp, &
          -0.93537087292458e-25_dp]
     PetscInt :: I(34) = [ &           
           0,  0,  0,  0,  0,  0, 0, 0, 1, 1, 1, 1, 1, 1, &
           2,  2,  2,  2,  2,  3, 3, 3, 4, 4, 4, 5, 8, 8, &
          21, 23, 29, 30, 31, 32]
     PetscInt :: J(34) = [ &            
           -2, -1,   0,   1,   2,   3,   4,   5, -9, -7, -1,  0,  1, &
            3, -3,   0,   1,   3,  17,  -4,   0,  6, -5, -2, 10, -8,  &
          -11, -6, -29, -31, -38, -39, -40, -41]
     PetscReal :: nI(34), nJ(34)
     PetscInt :: I_1(34), J_1(34)
     type(powertable_type) :: pi, pj
   contains
     private
     procedure, public :: init => region1_init
     procedure, public :: destroy => region1_destroy
     procedure, public :: properties => region1_properties
  end type IAPWS_region1_type

!------------------------------------------------------------------------
  ! Region 2 (steam) type
!------------------------------------------------------------------------

  type, public, extends(IAPWS_region_type) :: IAPWS_region2_type
     !! IAPWS-97 region 2 (steam) type.
     private
     PetscReal :: pstar = 1.0e6_dp, tstar = 540.0_dp
     PetscReal :: n0(9) = [ &
          -0.96927686500217e1_dp,   0.10086655968018e2_dp, -0.56087911283020e-2_dp, &
           0.71452738081455e-1_dp, -0.40710498223928_dp,    0.14240819171444e1_dp,  &
          -0.43839511319450e1_dp,  -0.28408632460772_dp,    0.21268463753307e-1_dp]
     PetscReal :: n(43) = [ &
          -0.17731742473213e-2_dp,  -0.17834862292358e-1_dp,  -0.45996013696365e-1_dp,  -0.57581259083432e-1_dp,  &
          -0.50325278727930e-1_dp,  -0.33032641670203e-4_dp,  -0.18948987516315e-3_dp,  -0.39392777243355e-2_dp,  &
          -0.43797295650573e-1_dp,  -0.26674547914087e-4_dp,   0.20481737692309e-7_dp,   0.43870667284435e-6_dp,  &
          -0.32277677238570e-4_dp,  -0.15033924542148e-2_dp,  -0.40668253562649e-1_dp,  -0.78847309559367e-9_dp,  &
           0.12790717852285e-7_dp,   0.48225372718507e-6_dp,   0.22922076337661e-5_dp,  -0.16714766451061e-10_dp, &
          -0.21171472321355e-2_dp,  -0.23895741934104e2_dp,   -0.59059564324270e-17_dp, -0.12621808899101e-5_dp,  &
          -0.38946842435739e-1_dp,   0.11256211360459e-10_dp, -0.82311340897998e1_dp,    0.19809712802088e-7_dp,  &
           0.10406965210174e-18_dp, -0.10234747095929e-12_dp, -0.10018179379511e-8_dp,  -0.80882908646985e-10_dp, &
           0.10693031879409_dp,     -0.33662250574171_dp,      0.89185845355421e-24_dp,  0.30629316876232e-12_dp, &
          -0.42002467698208e-5_dp,  -0.59056029685639e-25_dp,  0.37826947613457e-5_dp,  -0.12768608934681e-14_dp, &
           0.73087610595061e-28_dp,  0.55414715350778e-16_dp, -0.94369707241210e-6_dp ]
     PetscInt :: J0(9) = [0, 1, -5, -4, -3, -2, -1, 2, 3]
     PetscInt :: I(43) = [ &
           1,  1,  1,  1,  1,  2,  2,  2,  2,  2,  3,  3,  3,  3,  3, &
           4,  4,  4,  5,  6,  6,  6,  7,  7,  7,  8,  8,  9, 10, 10, &
          10, 16, 16, 18, 20, 20, 20, 21, 22, 23, 24, 24, 24]
     PetscInt :: J(43) = [ &
           0,  1,  2,  3,  6,  1,  2,  4,  7, 36,  0,  1,  3, 6,&
          35,  1,  2,  3,  7,  3, 16, 35,  0, 11, 25,  8, 36,  &
          13,  4, 10, 14, 29, 50, 57, 20, 35, 48, 21, 53,      &
          39, 26, 40, 58]
     PetscReal :: n0J0(9), nI(43), nJ(43)
     PetscInt :: J0_1(9), I_1(43), J_1(43)
     type(powertable_type) :: pj0, pi, pj
   contains
     private
     procedure, public :: init => region2_init
     procedure, public :: destroy => region2_destroy
     procedure, public :: properties => region2_properties
  end type IAPWS_region2_type

!------------------------------------------------------------------------
  ! Region 3 (supercritical) type
!------------------------------------------------------------------------

  type, public, extends(IAPWS_region_type) :: IAPWS_region3_type
     !! IAPWS-97 region 3 (supercritical) type.
     private
     PetscReal :: dstar = dcritical, tstar = tcriticalk
     PetscReal :: n(40) = [ &
           0.10658070028513e1_dp, -0.15732845290239e2_dp,   0.20944396974307e2_dp,  -0.76867707878716e1_dp,  &
           0.26185947787954e1_dp, -0.28080781148620e1_dp,   0.12053369696517e1_dp,  -0.84566812812502e-2_dp, &
          -0.12654315477714e1_dp, -0.11524407806681e1_dp,   0.88521043984318_dp,    -0.64207765181607_dp,    &
           0.38493460186671_dp,   -0.85214708824206_dp,     0.48972281541877e1_dp,  -0.30502617256965e1_dp,  &
           0.39420536879154e-1_dp, 0.12558408424308_dp,    -0.27999329698710_dp,     0.13899799569460e1_dp,  &
          -0.20189915023570e1_dp, -0.82147637173963e-2_dp, -0.47596035734923_dp,     0.43984074473500e-1_dp, &
          -0.44476435428739_dp,    0.90572070719733_dp,     0.70522450087967_dp,     0.10770512626332_dp,    &
          -0.32913623258954_dp,   -0.50871062041158_dp,    -0.22175400873096e-1_dp,  0.94260751665092e-1_dp, &
           0.16436278447961_dp,   -0.13503372241348e-1_dp, -0.14834345352472e-1_dp,  0.57922953628084e-3_dp, &
           0.32308904703711e-2_dp, 0.80964802996215e-4_dp, -0.16557679795037e-3_dp, -0.44923899061815e-4_dp]
     PetscInt :: I(40) = [ &
          0,  0,  0,   0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, &
          3,  3,  3,   3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 8, 9, &
          9 , 10, 10, 11]
     PetscInt :: J(40) = [ &
           0, 0,  1, 2,  7, 10, 12, 23, 2,  6, 15, 17,  0, 2, 6,  7, 22,  &
          26, 0,  2, 4, 16, 26,  0,  2, 4, 26,  1,  3, 26, 0, 2, 26,  2, &
          26, 2, 26, 0,  1, 26]
     PetscReal :: nI(40), nJ(40)
     PetscInt :: I_1(40), J_1(40)
     type(powertable_type) :: pi, pj
   contains
     private
     procedure, public :: init => region3_init
     procedure, public :: destroy => region3_destroy
     procedure, public :: properties => region3_properties
  end type IAPWS_region3_type

!------------------------------------------------------------------------
  ! Region 2/3 boundary type
!------------------------------------------------------------------------

  type, public :: IAPWS_boundary23_type
     !! IAPWS-97 boundary between regions 2 and 3.
     private
     PetscReal :: pstar = 1.0e6_dp
     PetscReal :: n(5) = [ &
          0.34805185628969e3_dp, -0.11671859879975e1_dp, 0.10192970039326e-2_dp, &
          0.57254459862746e3_dp,  0.13918839778870e2_dp]
     contains
       private
       procedure, public :: temperature => boundary23_temperature
       procedure, public :: pressure => boundary23_pressure
    end type IAPWS_boundary23_type

!------------------------------------------------------------------------
  ! IAPWS thermodynamics type
!------------------------------------------------------------------------

  type, extends(thermodynamics_type), public :: IAPWS_type
     !! IAPWS thermodynamics type.
     private
     type(IAPWS_boundary23_type), public :: boundary23
   contains
     private
     procedure, public :: init => IAPWS_init
     procedure, public :: destroy => IAPWS_destroy
     procedure, public :: phase_composition => IAPWS_phase_composition
  end type IAPWS_type

!------------------------------------------------------------------------

contains

!------------------------------------------------------------------------
  ! IAPWS class
!------------------------------------------------------------------------

  subroutine IAPWS_init(self)
    !! Initializes IAPWS object.

    class(IAPWS_type), intent(in out) :: self
    ! Locals:
    PetscInt :: i

    self%name = "IAPWS-97"

    allocate(IAPWS_saturation_type :: self%saturation)

    self%num_regions = 3
    allocate(IAPWS_region1_type :: self%water)
    allocate(IAPWS_region2_type :: self%steam)
    allocate(IAPWS_region3_type :: self%supercritical)
    allocate(self%region(self%num_regions))

    call self%region(1)%set(self%water)
    call self%region(2)%set(self%steam)
    call self%region(3)%set(self%supercritical)

    do i = 1, self%num_regions
       call self%region(i)%ptr%init()
    end do

  end subroutine IAPWS_init

!------------------------------------------------------------------------

  subroutine IAPWS_destroy(self)
    !! Destroys IAPWS object.

    class(IAPWS_type), intent(in out) :: self
    ! Locals:
    PetscInt :: i
    
    do i = 1, self%num_regions
       call self%region(i)%ptr%destroy()
    end do
    deallocate(self%region)
    deallocate(self%water, self%steam, self%supercritical)
    deallocate(self%saturation)

  end subroutine IAPWS_destroy

!------------------------------------------------------------------------

  PetscInt function IAPWS_phase_composition(self, region, pressure, &
       temperature) result(phases)
    !! Returns phase composition integer for given region, pressure
    !! and temperature. Here the bits represent:
    !! 0: liquid
    !! 1: vapour
    !! 2: supercritical

    class(IAPWS_type), intent(in) :: self
    PetscInt, intent(in) :: region
    PetscReal, intent(in) :: pressure, temperature
    ! Locals:
    PetscReal :: saturation_pressure
    PetscInt :: ierr

    phases = b'000'

    if (region == 4) then
       phases = b'011'
    else
       if (temperature <= tcritical) then
          select case(region)
             case (1)
                phases = b'001'
             case (2)
                phases = b'010'
             case (3)
                call self%saturation%pressure(temperature, &
                     saturation_pressure, ierr)
                if (ierr == 0) then
                   if (pressure >= saturation_pressure) then
                      phases = b'001'
                   else
                      phases = b'010'
                   end if
                else
                   phases = -1 ! error in saturation pressure
                end if
             end select
       else
          if (pressure <= pcritical) then
             phases = b'010'
          else
             phases = b'100'
          end if
       end if
    end if

  end function IAPWS_phase_composition

!------------------------------------------------------------------------
  ! Abstract region type
!------------------------------------------------------------------------

  subroutine region_init(self)
    !! Initializes abstract IAPWS-97 region object.
    
    class(IAPWS_region_type), intent(in out) :: self

    call self%visc%init()

  end subroutine region_init

!------------------------------------------------------------------------

  subroutine region_destroy(self)
    !! Destroys abstract IAPWS-97 region object.

    class(IAPWS_region_type), intent(in out) :: self

    call self%visc%destroy()

  end subroutine region_destroy

!------------------------------------------------------------------------

  subroutine region_properties(self, param, props, err)
    !! Dummy properties routine for abstract IAPWS-97 region object- to be
    !! overridden for specific regions.

    class(IAPWS_region_type), intent(in out) :: self
    PetscReal, intent(in) :: param(:) !! Primary variables
    PetscReal, intent(out):: props(:) !! Properties
    PetscInt, intent(out) :: err !! error code

    continue
    
  end subroutine region_properties

!------------------------------------------------------------------------

  subroutine region_viscosity(self, temperature, pressure, density, viscosity)
    !! IAPWS viscosity routine.

    class(IAPWS_region_type), intent(in out) :: self
    PetscReal, intent(in) :: temperature !! Temperature
    PetscReal, intent(in) :: pressure    !! Pressure (not used)
    PetscReal, intent(in) :: density     !! Density
    PetscReal, intent(out) :: viscosity  !! Viscosity

    ! Locals:
    PetscReal:: del, tk, tau
    PetscReal:: mu0, mu1, s0, s1

    tk = temperature + tc_k
    tau = tk / tcriticalk
    del = density / dcritical

    call self%visc%pk%compute(1._dp / tau)
    call self%visc%pi%compute(self%visc%pk%power(1) - 1._dp)
    call self%visc%pj%compute(del - 1._dp)

    ! Viscosity in dilute-gas limit:
    s0 = dot_product(self%visc%h0, self%visc%pk%power(0:3))
    mu0 = 100._dp * dsqrt(tau) / s0

    ! Contribution due to finite density:
    s1 = sum(self%visc%pi%power(self%visc%I) * self%visc%h1 * self%visc%pj%power(self%visc%J))
    mu1 = exp(del*s1)

    viscosity = self%visc%mustar * mu0 * mu1

  end subroutine region_viscosity

!------------------------------------------------------------------------
  ! Region 1 (liquid water)
!------------------------------------------------------------------------

  subroutine region1_init(self)
    !! Initializes IAPWS region 1 object.

    class(IAPWS_region1_type), intent(in out) :: self

    call self%IAPWS_region_type%init()

    self%name = 'water'

    self%nI = self%n * self%I
    self%nJ = self%n * self%J
    self%I_1 = self%I - 1
    self%J_1 = self%J - 1

    ! Configure power tables:
    call self%pi%configure(self%I)
    call self%pi%configure(self%I_1)

    call self%pj%configure(self%J)
    call self%pj%configure(self%J_1)

  end subroutine region1_init

!------------------------------------------------------------------------

  subroutine region1_destroy(self)
    !! Destroys IAPWS region 1 object.

    class(IAPWS_region1_type), intent(in out) :: self

    call self%pi%destroy()
    call self%pj%destroy()

    call self%IAPWS_region_type%destroy()

  end subroutine region1_destroy

!------------------------------------------------------------------------

  subroutine region1_properties(self, param, props, err)
    !! Calculates density and internal energy of liquid water as a function of
    !! pressure (Pa) and temperature (deg C).
    !!
    !! Returns err = 1 if called outside its operating range (t<=350 deg C, p<=100 MPa).

    class(IAPWS_region1_type), intent(in out) :: self
    PetscReal, intent(in) :: param(:) !! Primary variables (pressure, temperature)
    PetscReal, intent(out):: props(:) !! (density, internal energy)
    PetscInt, intent(out) :: err !! error code
    ! Locals:
    PetscReal:: tk, rt, pi, tau, gampi, gamt
    
    associate (p => param(1), t => param(2))

      ! Check input:
      if ((t <= 350.0_dp).and.(p <= 100.e6_dp)) then
         !      
         tk = t + tc_k
         rt = rconst * tk
         pi = p / self%pstar
         tau = self%tstar / tk

         call self%pi%compute(7.1_dp - pi)
         call self%pj%compute(tau - 1.222_dp)

         gampi = -sum(self%nI * self%pi%power(self%I_1) * self%pj%power(self%J))
         gamt = sum(self%nJ * self%pi%power(self%I) * self%pj%power(self%J_1))

         props(1) = self%pstar / (rt * gampi)      ! density
         props(2) = rt * (tau * gamt - pi * gampi) ! internal energy
         err = 0

      else
         err = 1
      end if

    end associate

  end subroutine region1_properties

!------------------------------------------------------------------------
  ! Region 2 (steam)
!------------------------------------------------------------------------

  subroutine region2_init(self)
    !! Initializes IAPWS region 2 object.

    class(IAPWS_region2_type), intent(in out) :: self

    call self%IAPWS_region_type%init()

    self%name = 'steam'

    self%n0J0 = self%n0 * self%J0
    self%nI = self%n * self%I
    self%nJ = self%n * self%J
    self%J0_1 = self%J0 - 1
    self%I_1 = self%I - 1
    self%J_1 = self%J - 1

    ! Configure power tables:
    call self%pj0%configure(self%J0)
    call self%pj0%configure(self%J0_1)

    call self%pi%configure(self%I)
    call self%pi%configure(self%I_1)
    call self%pi%configure([-1]) ! also need pi%power(-1) to calculate gampi

    call self%pj%configure(self%J)
    call self%pj%configure(self%J_1)

  end subroutine region2_init

!------------------------------------------------------------------------

  subroutine region2_destroy(self)
    !! Destroys IAPWS region 2 object.

    class(IAPWS_region2_type), intent(in out) :: self

    call self%pj0%destroy()
    call self%pi%destroy()
    call self%pj%destroy()

    call self%IAPWS_region_type%destroy()

  end subroutine region2_destroy

!------------------------------------------------------------------------

  subroutine region2_properties(self, param, props, err)
    !! Calculates density and internal energy of dry steam as a function of
    !! pressure (Pa) and temperature (deg C).
    !!
    !! Returns err = 1 if called outside its operating range (t<=800 deg C, p<=100 MPa).

    class(IAPWS_region2_type), intent(in out) :: self
    PetscReal, intent(in) :: param(:) !! Primary variables (pressure, temperature)
    PetscReal, intent(out):: props(:)  !! (density, internal energy)
    PetscInt, intent(out) :: err  !! error code
    ! Locals:
    PetscReal:: tk, rt, pi, tau, gampir, gamt0, gamtr, gampi

    associate (p => param(1), t => param(2))

      ! Check input:
      if ((t <= 800.0_dp).and.(p <= 100.e6_dp)) then

         tk = t + tc_k
         rt = rconst * tk
         pi = p / self%pstar
         tau = self%tstar / tk

         call self%pj0%compute(tau)
         call self%pi%compute(pi)
         call self%pj%compute(tau - 0.5_dp)

         gamt0  = sum(self%n0J0 * self%pj0%power(self%J0_1))
         gampir = sum(self%nI * self%pi%power(self%I_1) * self%pj%power(self%J))
         gamtr  = sum(self%nJ * self%pi%power(self%I)   * self%pj%power(self%J_1))

         gampi = self%pi%power(-1) + gampir

         props(1) = self%pstar / (rt * gampi)                 ! density
         props(2) = rt * (tau * (gamt0 + gamtr) - pi * gampi) ! internal energy
         err = 0

      else
         err = 1
      end if

    end associate

  end subroutine region2_properties

!------------------------------------------------------------------------
  ! Region 3 (supercritical)
!------------------------------------------------------------------------

  subroutine region3_init(self)
    !! Initializes IAPWS region 3 object.

    class(IAPWS_region3_type), intent(in out) :: self

    call self%IAPWS_region_type%init()

    self%name = 'supercritical'

    self%nI = self%n * self%I
    self%nJ = self%n * self%J
    self%I_1 = self%I - 1
    self%J_1 = self%J - 1

    ! Configure power tables:
    call self%pi%configure(self%I)
    call self%pi%configure(self%I_1)

    call self%pj%configure(self%J)
    call self%pj%configure(self%J_1)

  end subroutine region3_init

!------------------------------------------------------------------------

  subroutine region3_destroy(self)
    !! Destroyes IAPWS region 3 object.

    class(IAPWS_region3_type), intent(in out) :: self

    call self%pi%destroy()
    call self%pj%destroy()

    call self%IAPWS_region_type%destroy()

  end subroutine region3_destroy

!------------------------------------------------------------------------

  subroutine region3_properties(self, param, props, err)
    !! Calculates pressure and internal energy of supercritical water/steam
    !! as a function of density and temperature (deg C).
    !!
    !! Returns err = 1 if resulting pressure is outside its operating range (p<=100 MPa).

    class(IAPWS_region3_type), intent(in out) :: self
    PetscReal, intent(in) :: param(:) !! Primary variables (density, temperature)
    PetscReal, intent(out):: props(:)  !! (pressure, internal energy)
    PetscInt, intent(out) :: err   !! Error code
    ! Locals:
    PetscReal:: tk, rt, delta, tau
    PetscReal:: phidelta, phitau

    associate (d => param(1), t => param(2))

      tk = t + tc_k
      rt = rconst * tk
      tau = self%tstar / tk
      delta = d / self%dstar

      call self%pi%compute(delta)
      call self%pj%compute(tau)

      phidelta = self%n(1) * self%pi%power(-1) + &
           sum(self%nI * self%pi%power(self%I_1) * self%pj%power(self%J))
      phitau = sum(self%nJ * self%pi%power(self%I)   * self%pj%power(self%J_1))

      props(1) = d * rt * delta * phidelta ! pressure
      props(2) = rt * tau * phitau         ! internal energy
      if (props(1) > 100.0e6_dp) then
         err = 1
      else
         err = 0
      end if

    end associate

  end subroutine region3_properties

!------------------------------------------------------------------------
  ! Viscosity
!------------------------------------------------------------------------

  subroutine viscosity_init(self)
    !! Initializes viscosity object.

    class(IAPWS_viscosity_type), intent(in out) :: self

    ! Configure power tables:
    call self%pi%configure(self%I)
    call self%pj%configure(self%J)
    call self%pk%configure(self%k)

  end subroutine viscosity_init

!------------------------------------------------------------------------

  subroutine viscosity_destroy(self)
    !! Destroys viscosity object.

    class(IAPWS_viscosity_type), intent(in out) :: self

    call self%pi%destroy()
    call self%pj%destroy()
    call self%pk%destroy()

  end subroutine viscosity_destroy

!------------------------------------------------------------------------
  ! Saturation curve
!------------------------------------------------------------------------

  subroutine saturation_pressure(self, t, p, err)
    !! Calculates saturation pressure as a function of temperature.
    !! Returns err = 1 if called outside its operating range (0 <= t <= critical temperature).

    class(IAPWS_saturation_type), intent(in) :: self
    PetscReal, intent(in) :: t  !! Fluid temperature (\(^\circ C\))
    PetscReal, intent(out):: p  !! Fluid pressure (\(kg. m. s^{-1}\))
    PetscInt, intent(out) :: err  !! Error code
    ! Locals:
    PetscReal:: tk
    PetscReal:: theta, theta2, a, b, c, x

    if ((t >= 0._dp).and.(t <= tcritical)) then
       tk = t + tc_k      
       theta = tk + self%n(9) / (tk - self%n(10))
       theta2 = theta * theta
       a = theta2 + self%n(1) * theta + self%n(2)
       b = self%n(3) * theta2 + self%n(4) * theta + self%n(5)
       c = self%n(6) * theta2 + self%n(7) * theta + self%n(8)
       x = 2._dp * c / (-b + dsqrt(b*b - 4._dp*a*c))
       x = x * x
       p = self%pstar * x * x
       err = 0
    else
       err = 1
    end if

  end subroutine saturation_pressure

!------------------------------------------------------------------------

  subroutine saturation_temperature(self, p, t, err)
    !! Calculates saturation temperature (deg C) as a function of pressure.
    !! Returns err = 1 if called outside its operating range (611.213 Pa <= p <= critical pressure).

    class(IAPWS_saturation_type), intent(in) :: self
    PetscReal, intent(in) :: p  !! Fluid pressure (\(kg. m. s^{-1}\))
    PetscReal, intent(out):: t  !! Fluid temperature (\(^\circ C\))
    PetscInt, intent(out) :: err !! Error code
    ! Locals:
    PetscReal:: beta, beta2, d, e, f, g, x

    if ((p >= 611.213_dp).and.(p <= pcritical)) then
       beta2 = dsqrt(p / self%pstar)
       beta = dsqrt(beta2)
       e = beta2 + self%n(3) * beta + self%n(6)
       f = self%n(1) * beta2 + self%n(4) * beta + self%n(7)
       g = self%n(2) * beta2 + self%n(5) * beta + self%n(8)
       d = 2.0_dp * g / (-f - dsqrt(f*f - 4._dp*e*g))
       x = self%n(10) + d
       t = 0.5_dp * (self%n(10) + d - dsqrt(x*x - 4._dp * (self%n(9) + self%n(10) * d))) - tc_k
       err = 0
    else
       err = 1
    end if

  end subroutine saturation_temperature

!------------------------------------------------------------------------
  ! Region 2/3 boundary
!------------------------------------------------------------------------

  subroutine boundary23_pressure(self, t, p)
    !! Calculates the pressure p (Pa) on the boundary between regions 2 and 3,
    !! given a temperature t (deg C).

    class(IAPWS_boundary23_type), intent(in) :: self
    PetscReal, intent(in) :: t  !! Fluid temperature (\(^\circ C\))
    PetscReal, intent(out):: p  !! Fluid pressure (\(kg. m. s^{-1}\))
    ! Local variable:      
    PetscReal:: tk

    tk = t + tc_k
    p  = self%pstar * (self%n(1) + tk * (self%n(2) + tk * self%n(3))) 

  end subroutine boundary23_pressure

!-----------------------------------------------------------------------

  subroutine boundary23_temperature(self, p, t)
    !! Calculates the temperature t (deg C) on the boundary between regions 2 and 3,
    !! given a pressure p (Pa).

    class(IAPWS_boundary23_type), intent(in) :: self
    PetscReal, intent(in) :: p  !! Fluid pressure (\(kg. m. s^{-1}\))
    PetscReal, intent(out):: t  !! Fluid temperature (\(^\circ C\))

    t = self%n(4) + dsqrt((p/self%pstar - self%n(5)) / self%n(3)) - tc_k 

  end subroutine boundary23_temperature

!------------------------------------------------------------------------

end module IAPWS_module
