!> \file mod_types_basic.f90
!> \brief Module of elfe3d containing type definitions
!> \details Provides fundamental derived types used across elfe3D for geometric
!> \details coordinates, element locations, and array index management.
!!
!> original module in emilia inversion software (Kalscheuer 2008, 2010)
!!
!> written by Thomas Kalscheuer, 09/11/2009
!> Copyright (C) Thomas Kalscheuer 2010
!>
!>  This file is part of elfe3D.
!> 
!>  Licensed under the Apache License, Version 2.0 (the "License"); 
!>  you may not use this file except in compliance with the License.  
!>  You may obtain a copy of the License at
!> 
!>      https://www.apache.org/licenses/LICENSE-2.0
!> 
!>  Unless required by applicable law or agreed to in writing, software
!>  distributed under the License is distributed on an "AS IS" BASIS, 
!>  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or 
!>  implied.  See the
!>  License for the specific language governing permissions and 
!>  limitations under the License. 

module mod_types_basic

  use mod_constant

  implicit none

  !> \brief Two-dimensional point coordinates
  type :: point_2D
     real(kind=dp) :: x, y
  end type point_2D

  !> \brief Three-dimensional point coordinates
  type :: point_3D
     real(kind=dp) :: x, y, z
  end type point_3D

  !> \brief Complex-valued three-dimensional point coordinates
  type :: cpoint_3D
     complex(kind=dp) :: x, y, z
  end type cpoint_3D


  !> \brief Two-dimensional local coordinates used in surface or element definitions
  type :: loc_2D
     real(kind=dp) :: y, z
  end type loc_2D

  !> \brief Three-dimensional local coordinates used in surface or element definitions
  type :: loc_3D
     real(kind=dp) :: x, y, z
  end type loc_3D


  !> \brief Two-dimensional integer index for array addressing
  type :: index_2D
     integer :: iy, iz
  end type index_2D

  !> \brief Three-dimensional integer index for array addressing
  type :: index_3D
     integer :: ix, iy, iz
  end type index_3D

end module mod_types_basic
