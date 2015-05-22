module powertable_module
  !! Module for efficiently computing a table of powers of a real number.
  !! The algorithm uses repeated squaring and a lookup table, and avoids calculating
  !! any powers that are not needed.

  use kinds_module

  implicit none
  private

#include <petsc-finclude/petscdef.h>

  type, private :: product_pointer
     private
     PetscReal, pointer :: fac1, fac2
     PetscReal, pointer :: prod
   contains
     procedure, public :: set => product_pointer_set
  end type product_pointer

  type, public :: powertable
     !! Type for efficiently calculating a table of powers of a
     !! double precision number.
     private
     PetscReal, allocatable, public :: power(:) !! Computed powers
     PetscInt, allocatable :: product(:,:)
     type(product_pointer), allocatable :: powerlist(:)
     PetscInt, public :: lower = 0, upper = 0 !! Lower and upper bounds of powers
     PetscInt, allocatable :: required(:)
     PetscInt :: powerlist_size
   contains
     private
     procedure, public  :: configure => powertable_configure
     procedure, public  :: compute => powertable_compute 
     procedure, public  :: destroy => powertable_destroy
     procedure :: product_configured => powertable_product_configured
     procedure :: configure_product => powertable_configure_product
     procedure :: set_powerlist => powertable_set_powerlist
  end type powertable

contains

!------------------------------------------------------------------------

  subroutine product_pointer_set(self, lower, upper, arr, i1, i2, iprod)
    !! Sets up product pointer into an array. The pointers fac1 and fac2 point
    !! to the factors in the array, and prod points to the resulting product.

    class(product_pointer), intent(in out) :: self
    PetscInt, intent(in) :: lower, upper
    PetscReal, intent(in), target :: arr(lower:upper)
    PetscInt, intent(in) :: i1, i2, iprod

    self%fac1 => arr(i1)
    self%fac2 => arr(i2)
    self%prod => arr(iprod)
    
  end subroutine product_pointer_set

!------------------------------------------------------------------------

  PetscBool function powertable_product_configured(self, i)
    !! Returns true if product combination has already been computed
    !! for the specified power.

    class(powertable), intent(in out) :: self
    PetscInt, intent(in) :: i

    powertable_product_configured = &
         (abs(i) <= 1) .or. &
         ((self%product(1,i) /= 0) .and. (self%product(2,i) /= 0))

  end function powertable_product_configured

!------------------------------------------------------------------------

  recursive subroutine powertable_configure_product(self, i)
    !! Configures product combination for the specifed power.

    class(powertable), intent(in out) :: self
    PetscInt, intent(in) :: i
    ! Locals:
    PetscInt :: s, j, i2, c

    if ((abs(i) > 1) .and. (.not.(self%product_configured(i)))) then
       ! First look for combinations of existing powers:
       s = sign(1, i)
       i2 = nint(i/2.)
       do c = i2, s, -s
          j = i-c
          if ((self%product_configured(c) .and. (self%product_configured(j)))) then
             self%product(1,i) = c
             self%product(2,i) = j
             exit
          end if
       end do
       if (.not.(self%product_configured(i))) then
          ! Use repeated squaring:
          if (mod(i,2) == 0) then
             call self%configure_product(i2)
             self%product(1:2,i) = i2
             self%required(i2) = 2
          else
             j = i-s
             call self%configure_product(j)
             self%product(1,i) = s
             self%product(2,i) = j
             self%required(j) = 2
          end if
       end if
    end if

  end subroutine powertable_configure_product

!------------------------------------------------------------------------

  subroutine powertable_configure(self, powers)
    !! Calculates the most efficient way of computing the required powers,
    !! adding any extra powers needed for intermediate steps.

    class(powertable), intent(in out) :: self
    PetscInt, intent(in), dimension(:) :: powers
    ! Locals:
    PetscInt :: min_power, max_power, i, old_lower, old_upper
    PetscInt :: old_required(self%lower:self%upper), s, u
    PetscBool :: enlarge

    old_lower = 0; old_upper =0

    ! Always store zeroth and first powers: 
    min_power = min(minval(powers), 0)
    max_power = max(maxval(powers), 1)

    old_required = 0

    enlarge = allocated(self%power)

    if (enlarge) then ! enlarge table:
       old_lower = self%lower
       old_upper = self%upper
       old_required = self%required
       self%lower = min(self%lower, min_power)
       self%upper = max(self%upper, max_power)
       deallocate(self%power, self%product, self%required)
    else
       self%lower = min_power
       self%upper = max_power
    end if

    allocate(self%power(self%lower : self%upper), &
         self%product(2,self%lower : self%upper), &
         self%required(self%lower : self%upper))

    self%power = 0._dp
    self%power(0) = 1._dp
    self%product = 0
    self%required = 0

    ! Copy old required table (only powers actually specified,
    ! not secondary values with required = 2)
    if (enlarge) then
       do i = old_lower, old_upper
          if (old_required(i) == 1) then
             self%required(i) = old_required(i)
          end if
       end do
    end if

    ! Set new required values:
    do i = lbound(powers,1), ubound(powers,1)
       if (abs(powers(i)) > 1) then
          self%required(powers(i)) = 1
       end if
    end do

    ! Configure products:
    do s = 1, -1, -2
       if (s > 0) then
          u = self%upper
       else
          u = self%lower
       end if
       do i = s*2, u, s
          if (self%required(i) > 0) then
             call self%configure_product(i)
          end if
       end do
    end do

    call self%set_powerlist()

  end subroutine powertable_configure

!------------------------------------------------------------------------

  subroutine powertable_set_powerlist(self)
    !! Generates list of which powers are required, in the order in which they
    !! should be computed. (Powers 0,1 are not included, as they are always
    !! required, as is -1 if negative powers are required.)

    class(powertable), intent(in out) :: self
    ! Locals:
    PetscInt :: p, n, s, u

    if (allocated(self%powerlist)) then
       deallocate(self%powerlist)
    end if
    self%powerlist_size = count(self%required > 0)
    allocate(self%powerlist(self%powerlist_size))
    
    n = 0
    do s = 1, -1, -2
       if (s > 0) then
          u = self%upper
       else
          u = self%lower
       end if
       do p = s*2, u, s
          if (self%required(p) > 0) then
             n = n + 1
             call self%powerlist(n)%set(self%lower, self%upper, self%power,&
                  self%product(1,p), self%product(2,p), p)
          end if
       end do
    end do

  end subroutine powertable_set_powerlist

!------------------------------------------------------------------------

  subroutine powertable_destroy(self)
    !! Destroys a powertable object.

    class(powertable), intent(in out) :: self

    deallocate(self%power, self%product, self%powerlist, &
         self%required)

  end subroutine powertable_destroy

!------------------------------------------------------------------------

  subroutine powertable_compute(self, val)
    !! Computes the table of powers for the specified value.

    class(powertable), intent(in out) :: self
    PetscReal, intent(in) :: val !! value to compute powers of
    ! Locals:
    PetscInt :: n

    self%power(1) = val
    if (self%lower < 0) then
       self%power(-1) = 1._dp / val
    end if

    do n = 1, self%powerlist_size
       self%powerlist(n)%prod = self%powerlist(n)%fac1 * self%powerlist(n)%fac2
    end do

  end subroutine powertable_compute

!------------------------------------------------------------------------

end module powertable_module
