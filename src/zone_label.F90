!   Copyright 2017 University of Auckland.

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

module zone_label_module
  !! Module for a zone labelling.

#include <petsc/finclude/petscsys.h>

  use petscsys

  implicit none
  private

  PetscInt, parameter, public :: max_zone_name_length = 80

  public :: zone_label_name

contains

!------------------------------------------------------------------------

  function zone_label_name(zone_name) result(name)
    !! Returns name of zone label associated with a zone name.

    character(*), intent(in) :: zone_name
    character(:), allocatable :: name

    name = 'zone_' // trim(zone_name)

  end function zone_label_name

!------------------------------------------------------------------------

end module zone_label_module
