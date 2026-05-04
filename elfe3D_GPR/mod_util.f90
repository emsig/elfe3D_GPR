!> \file mod_util.f90
!> \brief Module of elfe3d containing useful utility routines
!> \details Contains I/O helpers, text processing, allocation checking, and
!> \details string-handling support functions used throughout the application.
!>
!> original module in emilia inversion software (Kalscheuer 2008, 2010)
!>
!> written by Thomas Kalscheuer, last change: 16/11/2009
!>
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

module mod_util

  use mod_types_basic

  implicit none

contains

  !-----------------------------------------------------------------
  !> \brief Write a message to the specified unit and optionally to screen
  !> \param[in] iu Fortran unit number where the message is written
  !> \param[in] msg Text message to output
  !-----------------------------------------------------------------
  subroutine Write_Message(iu, msg)

    ! input variables
    integer, intent(in) :: iu
    character(len=*), intent(in) :: msg

    ! local variables
    ! for inquire statement
    logical :: unit_open
    character(len=9) :: unit_action

    ! output message to unit iu
    inquire (unit=iu, opened=unit_open, action=unit_action)
    if ( (unit_open .eqv. .true.) .and. &
         (scan('WRITE',unit_action) .gt. 0) ) then
       write (unit=iu,fmt='(a)') trim(msg)
    end if

    ! output message to terminal
    if (LOGFILE_SCREEN .eqv. .true.) &
         write (unit=*,fmt='(a)') trim(msg)

  end subroutine Write_Message


  !-----------------------------------------------------------------
  !> \brief Write an error message to the specified unit and terminate execution
  !> \param[in] iu Fortran unit number where the error message is written
  !> \param[in] error_msg Error text to output before program termination
  !-----------------------------------------------------------------
  subroutine Write_Error_Message(iu, error_msg)

    ! input variables
    integer, intent(in) :: iu
    character(len=*), intent(in) :: error_msg

    ! local variables
    ! for inquire statement
    logical :: unit_open
    character(len=9) :: unit_action

    ! output error message to unit iu, ...
    inquire (unit=iu, opened=unit_open, action=unit_action)
    if ( (unit_open .eqv. .true.) .and. &
         (scan('WRITE',unit_action) .gt. 0) ) then
       write (unit=iu,fmt='(a,/,a,/,a)') trim(error_msg) , &
            "Correct input files and restart!!!", &
            "Quitting..."
    end if

    ! output error message to terminal, and ...
    if (LOGFILE_SCREEN .eqv. .true.) then
       write (unit=*,fmt='(a,/,a,/,a)') trim(error_msg) , &
            "Correct input files and restart!!!", &
            "Quitting..."
    end if

    ! stop program
    stop 1

  end subroutine Write_Error_Message


  !-----------------------------------------------------------------
  !> \brief Validate a required input variable and report problems
  !> \param[in] iu Fortran unit number to write validation messages
  !> \param[in] error_msg Name or description of the invalid input variable
  !-----------------------------------------------------------------
  subroutine Check_Input(iu, error_msg)

    ! input variables
    integer, intent(in) :: iu
    character(len=*), intent(in) :: error_msg

    ! local variables
    ! for inquire statement
    logical :: unit_open
    character(len=9) :: unit_action

    ! output error message to unit iu, ...
    inquire (unit=iu, opened=unit_open, action=unit_action)
    if ( (unit_open .eqv. .true.) .and. &
         (scan('WRITE',unit_action) .gt. 0) ) then
       write (unit=iu,fmt='(2a)') &
            "!!!!!! Incorrect input in variable: ", trim(error_msg)
       write (unit=iu,fmt='(a)') &
            "!!!!!! Correct input file and restart!"
    end if

    ! output error message to terminal, and ...
    if (LOGFILE_SCREEN .eqv. .true.) then
       write (unit=*,fmt='(2a)') &
            "!!!!!! Incorrect input in variable: ", trim(error_msg)
       write (unit=*,fmt='(a)') &
            "!!!!!! Correct input file and restart!"
    end if

    ! stop program
    stop

  end subroutine Check_Input


  !-----------------------------------------------------------------
  !> \brief Report an unrecognized input item and ask for correction
  !> \param[in] iu Fortran unit number to write validation messages
  !> \param[in] error_msg Name or description of the unknown input
  !-----------------------------------------------------------------
  subroutine Check_Unknown_Input(iu, error_msg)

    ! input variables
    integer, intent(in) :: iu
    character(len=*), intent(in) :: error_msg

    ! local variables
    ! for inquire statement
    logical :: unit_open
    character(len=9) :: unit_action

    ! output error message to unit iu, ...
    inquire (unit=iu, opened=unit_open, action=unit_action)
    if ( (unit_open .eqv. .true.) .and. &
         (scan('WRITE',unit_action) .gt. 0) ) then
       write (unit=iu,fmt='(2a)') &
            "!!!!!! Unknown input variable: ", trim(error_msg)
       write (unit=iu,fmt='(a)') &
            "!!!!!! Correct input file and restart!"
    end if

    ! output error message to terminal, and ...
    if (LOGFILE_SCREEN .eqv. .true.) then
       write (unit=*,fmt='(2a)') &
            "!!!!!! Unknown input variable: ", trim(error_msg)
       write (unit=*,fmt='(a)') &
            "!!!!!! Correct input file and restart!"
    end if

    ! stop program
    stop

  end subroutine Check_Unknown_Input


  !-----------------------------------------------------------------
  !> \brief Check allocation status and report failures
  !> \param[in] iu Fortran unit used for error messages
  !> \param[in] allo_stat Allocation return code from memory allocation
  !> \param[in] ctext Context string describing the allocation call
  !-----------------------------------------------------------------
  subroutine allocheck(iu, allo_stat, ctext)

    ! input variables
    integer, intent(in) :: iu, allo_stat
    character(len=*), intent(in) :: ctext

    ! local variables
    ! for inquire statement
    logical :: unit_open
    character(len=9) :: unit_action

    ! if allocation failed, ...
    if (allo_stat .ne. 0) then
       ! output error message to unit iu, ...
       inquire (unit=iu, opened=unit_open, action=unit_action)
       if ( (unit_open .eqv. .true.) .and. &
            (scan('WRITE',unit_action) .gt. 0) ) then
          write (iu,9000)
          write (iu,'(a,/,a,i4,/,a)') ctext, &
               "Allocation status: ", allo_stat, &
               "Quitting..."
          write (iu,9000)
       end if

       ! output error message to terminal, and ...
       if (LOGFILE_SCREEN .eqv. .true.) then
          write (*,9000)
          write (*,'(a,/,a,i4,/,a)') ctext, &
               "Allocation status: ", allo_stat, &
               "Quitting..."
          write (*,9000)
       end if

       ! stop program
       stop
    end if

9000 format('------------------------------------------------------')

  end subroutine allocheck


  !----------------------------------------------------------------------
  !> \brief Convert a character string to lower case
  !> \param[in,out] string Input/output string modified in place
  subroutine lower(string)

    ! input variables
    character(len = *) :: string

    ! local variables
    integer :: i, ismall, iBIG

    ! initialize
    ismall = iachar('a')
    iBIG = iachar('A')

    do i = 1,len(string)
       select case (string(i:i))
       case ('A':'Z')
          string(i:i) = achar( iachar(string(i:i)) + ismall - iBIG )
       end select
    end do

  end subroutine lower


  !----------------------------------------------------------------------
  !> \brief Convert a character string to upper case
  !> \param[in,out] string Input/output string modified in place
  subroutine upper(string)


    ! input variables
    character(len = *) :: string

    ! local variables
    integer :: i, ismall, iBIG

    ! initialize
    ismall = iachar('a')
    iBIG = iachar('A')

    do i = 1,len(string)
       select case (string(i:i))
       case ('a':'z')
          string(i:i) = achar( iachar(string(i:i)) + iBIG - ismall )
       end select
    end do

  end subroutine upper


  !----------------------------------------------------------------------
  !> \brief Find a substring within a string ignoring case
  !> \param[in] str1 Input search string
  !> \param[in] str2 Substring to locate within str1
  !> \return Position of the first occurrence of str2 in str1, or 0 if not found
  integer function findstr(str1,str2)

    ! input data
    character(len = *), intent(in) :: str1, str2
    ! returns the position of str2 in str1.  Ignores CASE.
    ! returns 0 IF str2 not found in str1

    ! local data
    integer :: i, j, capdif
    logical :: same

    capdif= ichar('a')-ichar('A')

    outer: do i= 1, len(str1)-len(str2)+1
       inner: do j = 1,len(str2)

          same= str1(i+j-1:i+j-1) .eq. str2(j:j)        .or. &
               
               'A'.le.str2(j:j) .and. str2(j:j).le.'Z' .and. &
               ichar(str1(i+j-1:i+j-1)) .eq. ichar(str2(j:j))+capdif .or. &
               
               'a'.le.str2(j:j) .and. str2(j:j).le.'z' .and. &
               ichar(str1(i+j-1:i+j-1)) .eq. ichar(str2(j:j)) - capdif

          if( .not. same ) cycle outer
       end do inner
       findstr = i
       return
    end do outer

    findstr = 0

  end function findstr


  !----------------------------------------------------------------------
  !> \brief Count words in a string
  !> \details Words are separated by spaces, tabs, or commas.
  !> \param[in] string Input string to parse
  !> \return Number of words found in the string
  integer function nword(string)

    ! input data
    character(len = *), intent(in) :: string

    ! local data
    integer :: i, ntmp
    logical :: wasblk

    wasblk = .true.
    ntmp = 0
    do i = 1,len(string)
       if ( string(i:i) == ' ' .or. string(i:i) == ',' ) then

          ! current CHARACTER is blank
          wasblk = .true.
       else
          if (wasblk) ntmp = ntmp + 1
          wasblk = .false.
       end if
    end do

    nword = ntmp

  end function nword


  !----------------------------------------------------------------------
  !> \brief Find first non-blank character index of the specified word
  !> \details Words are separated by spaces, tabs, or commas.
  !> \param[in] string Input string containing words
  !> \param[in] iwrd Word number to locate (1-based)
  !> \return Index of the first character of the iwrd-th word, or len(string) if not found
  integer function begwrd(string,iwrd)

    ! input data
    integer, intent(in) :: iwrd
    character(len = *), intent(in) :: string

    ! local data
    integer :: i, nword
    logical :: wasblk

    ! initialise
    wasblk = .true.
    nword = 0
    begwrd = -1

    do i = 1,len(string)
       if ( string(i:i) == ' ' .or. string(i:i) == ',' ) then
          ! current CHARACTER is blank
          wasblk = .true.
       else
          if (wasblk) nword = nword + 1
          wasblk = .false.
          if (nword == iwrd) then
             begwrd = i
             return
          end if
       end if
    end do

    begwrd = len(string)

  end function begwrd


  !----------------------------------------------------------------------
  !> \brief Find last non-blank character index of the specified word
  !> \details Words are separated by spaces, tabs, or commas.
  !> \param[in] string Input string containing words
  !> \param[in] iwrd Word number to locate (1-based)
  !> \return Index of the last character of the iwrd-th word, or len(string) if not found
  integer function endwrd(string,iwrd)

    ! input data
    integer, intent(in) :: iwrd
    character(len = *), intent(in) :: string

    ! local data
    integer :: i, nword
    logical :: wasblk
    intrinsic len

    ! initialise
    wasblk = .true.
    nword = 0
    endwrd = -1

    do i = 1,len(string)
       if ( string(i:i) == ' ' .or. string(i:i) == ',' ) then
          ! current CHARACTER is blank
          wasblk = .true.
          if (nword == iwrd) return
       else
          if (wasblk) nword = nword + 1
          wasblk = .false.
          if (nword == iwrd) endwrd = i
       end if
    end do
    endwrd = len(string)

  end function endwrd


  !----------------------------------------------------------------------
  !> \brief Sort an integer array using bubble sort
  !! \param[in,out] a Integer array to sort in place
  subroutine bubblesort_int(a)
    ! input/output variables
    integer, dimension(:), intent(inout) :: a

    ! local variables
    integer :: i, j, n
    logical :: sorted
    integer :: temp

    n = size(a)
    sorted = .false.
    j = 0

    do while (.not. sorted .and. j < n-1)
       sorted = .true.
       j = j+1
       do i = 1,n-j
          if (a(i) > a(i+1)) then
             temp = a(i)
             a(i) = a(i+1)
             a(i+1) = temp
             sorted = .false.
          end if
       end do
    end do
  end subroutine bubblesort_int


!----------------------------------------------------------------------
  !> \brief Sort a real array using bubble sort
  !! \param[in,out] a Real array to sort in place
  subroutine bubblesort_real(a)
    ! input/output variables
    real(kind=dp), dimension(:), intent(inout) :: a

    ! local variables
    integer :: i, j, n
    logical :: sorted
    real(kind=dp) :: temp

    n = size(a)
    sorted = .false.
    j = 0

    do while (.not. sorted .and. j < n-1)
       sorted = .true.
       j = j+1
       do i = 1,n-j
          if (a(i) > a(i+1)) then
             temp = a(i)
             a(i) = a(i+1)
             a(i+1) = temp
             sorted = .false.
          end if
       end do
    end do
  end subroutine bubblesort_real

  !----------------------------------------------------------------------
  ! bubble sort for array with real values
  subroutine bubblesort_real_2d(a)
    ! input/output variables
    type(loc_2D), dimension(:), intent(inout) :: a

    ! local variables
    integer :: i, j, n
    logical :: sorted
    type(loc_2D) :: temp

    n = size(a)
    sorted = .false.
    j = 0

    ! sort only after y-positions
    do while (.not. sorted .and. j < n-1)
       sorted = .true.
       j = j+1
       do i = 1,n-j
          if (a(i)%y > a(i+1)%y) then
             temp = a(i)
             a(i) = a(i+1)
             a(i+1) = temp
             sorted = .false.
          end if
       end do
    end do
    
!!$    ! for elements with equal y-positions sort after z-position
!!$    sorted = .false.
!!$    j = 0
!!$    do while (.not. sorted .and. j < n-1)
!!$       sorted = .true.
!!$       j = j+1
!!$       do i = 1,n-j
!!$          if ((a(i)%y .eq. a(i+1)%y) .and.) then
!!$             temp = a(i)
!!$             a(i) = a(i+1)
!!$             a(i+1) = temp
!!$             sorted = .false.
!!$          end if
!!$       end do
!!$    end do

  end subroutine bubblesort_real_2d

  !----------------------------------------------------------------------
  ! subroutine unique_real stores elements of input array a in ascending
  ! order and removes duplicate entries. The result is saved in an array
  ! a_u.
  subroutine unique_real(a, N_u, a_u)
    ! input variables
    real(kind=dp), dimension(:), intent(in) :: a

    ! output variables
    ! number of unique elements
    integer, intent(out) :: N_u
    ! array with unique matrix elements
    real(kind=dp), dimension(:), allocatable, intent(out) :: a_u

    ! local variables
    ! number of elements in a
    integer :: N_a
    ! temporary a
    real(kind=dp), dimension(:), allocatable :: a_tmp
    ! allocation status
    integer :: allo_stat
    ! counters
    integer :: i, j

    ! initialise
    N_a = size(a)

    ! allocate temporary array to store sorted results in
    allocate (a_tmp(N_a), stat = allo_stat)
    call allocheck(log_unit, allo_stat, &
         "unique_real: error allocating a_tmp")
    a_tmp = a

    ! sort into ascending order
    call bubblesort_real(a_tmp)

    ! initialise counters ...
    ! ... over unique elements
    i = 1
    ! ... to next non-duplicate element
    j = 2
    do while (i <= N_a .and. j <= N_a)
       ! find next entry that differs from i-th entry
       do while (j <= N_a)
          if (abs(a_tmp(j)-a_tmp(i)) .lt. 1E-06) then
             j = j+1
          else
             exit
          end if
       end do

       if (j <= N_a) then
          ! move found entry into (i+1)-th position
          a_tmp(i+1) = a_tmp(j)
          ! update i
          i = i+1
          ! update j
          j = j+1
       end if
    end do

    ! store array with unique elements for output
    N_u = i
    allocate(a_u(N_u), stat = allo_stat)
    call allocheck(log_unit, allo_stat, &
         "unique_real: error allocating a_u")
    a_u = a_tmp(1:N_u)

    if (allocated(a_tmp)) deallocate(a_tmp)

  end subroutine unique_real

  
  !----------------------------------------------------------------------
  ! subroutine unique_real stores elements of input array a in ascending
  ! order and removes duplicate entries. The result is saved in an array
  ! a_u.
  subroutine unique_real_2D(a, N_u, a_u)
    ! input variables
    type(loc_2D), dimension(:), intent(in) :: a

    ! output variables
    ! number of unique elements
    integer, intent(out) :: N_u
    ! array with unique matrix elements
    type(loc_2D), dimension(:), allocatable, intent(out) :: a_u

    ! local variables
    ! number of elements in a
    integer :: N_a
    ! temporary a
    type(loc_2D), dimension(:), allocatable :: a_tmp
    ! allocation status
    integer :: allo_stat
    ! counters
    integer :: i, j

    ! initialise
    N_a = size(a)

    ! allocate temporary array to store sorted results in
    allocate (a_tmp(N_a), stat = allo_stat)
    call allocheck(log_unit, allo_stat, "unique_real: error allocating a_tmp")
    a_tmp = a

    ! sort into ascending order
    call bubblesort_real_2D(a_tmp)

    ! initialise counters ...
    ! ... over unique elements
    i = 1
    ! ... to next non-duplicate element
    j = 2
    do while (i <= N_a .and. j <= N_a)
       ! find next entry that differs from i-th entry
       do while (j <= N_a)
          if (abs(a_tmp(j)%y-a_tmp(i)%y) .lt. 1E-06) then
             j = j+1
          else
             exit
          end if
       end do

       if (j <= N_a) then
          ! move found entry into (i+1)-th position
          a_tmp(i+1) = a_tmp(j)
          ! update i
          i = i+1
          ! update j
          j = j+1
       end if
    end do

    ! store array with unique elements for output
    N_u = i
    allocate(a_u(N_u), stat = allo_stat)
    call allocheck(log_unit, allo_stat, "unique_real: error allocating a_u")
    a_u = a_tmp(1:N_u)

    if (allocated(a_tmp)) deallocate(a_tmp)

  end subroutine unique_real_2D

  !---------------------------------------------------------------------
  !> TETRAHEDRON GEOMETRY
  !> new in version 1.1.0
  !> subroutine for calculating the centroids of all mesh tetrahedra with  
  !> nodes n1, n2 and n3
  !---------------------------------------------------------------------
  subroutine tetrahedra_centroids(n1, n2, n3, n4, M, elem_centr)
    ! INPUT
    ! x,y,z coordinates of element nodes
    real(kind=dp), dimension(:,:), intent(in) :: n1, n2, n3, n4
    ! number of elements in mesh
    integer :: M

    ! OUTPUT
    ! x,y,z coordinates of element centroids
    real (kind=dp), allocatable, dimension(:,:), intent(out) :: &
                                                              elem_centr

    ! LOCAL VARIABLES
    ! allocation status
    integer :: allo_stat

    !-------------------------------------------------------------------
    ! allocate output array
    allocate (elem_centr(M,3), stat = allo_stat)
    call allocheck(log_unit, allo_stat, &
                    "error allocating element centroid array")

    ! initialise output array
    elem_centr = 0.0_dp
    
    ! calculate centroids
    elem_centr(:,1) = (n1(:,1)+n2(:,1)+n3(:,1)+n4(:,1))/4.0_dp
    elem_centr(:,2) = (n1(:,2)+n2(:,2)+n3(:,2)+n4(:,2))/4.0_dp
    elem_centr(:,3) = (n1(:,3)+n2(:,3)+n3(:,3)+n4(:,3))/4.0_dp
   !--------------------------------------------------------------------
   end subroutine tetrahedra_centroids


  !----------------------------------------------------------------------

end module mod_util
