module kinds_module
  !! Numeric kinds module. This contains kind definitions for commonly used
  !! numeric data types.
  
  implicit none
  private
  
  ! double precision:
  integer, parameter, public :: dp = kind(0.d0) !! double precision kind

  ! Quiet NaN constant:
  real(dp), parameter, public :: qnan_dp = transfer(Z'FFF8000000000000', 1._dp)

end module kinds_module
