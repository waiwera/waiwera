module kinds_module
  !! Numeric kinds module. This contains kind definitions for commonly used
  !! numeric data types.
  
  implicit none
  private
  
  ! double precision:
  integer, parameter, public :: dp = kind(0.d0) !! double precision kind

end module kinds_module
