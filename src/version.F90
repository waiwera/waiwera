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

module version_module
  !! Module for software version information.

#include <petsc/finclude/petscsys.h>

  use petscsys

  implicit none
  private
  
  PetscInt, parameter :: max_version_string_length = 64
  character(len = max_version_string_length), public :: &
       waiwera_version = "1.1.0"

end module version_module
