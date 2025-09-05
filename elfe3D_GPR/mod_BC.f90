!> @brief
!> Module of elfe3d containing subroutines for boundary
!>  conditions-related calculations
!!
!>  written by Paula Rulff and Chaitanya Dinesh Singh, 24/06/2019
!!
!> So far only of Dirichlet-type boundary condition is implemented.
!> Also includes perfect electric conductor options, and 
!> the PML complex stretching functions.
!!
!>  Copyright (C) Paula Rulff and Chaitanya Dinesh Singh
!>
!>  This file is part of elfe3D_GPR
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

module BC

  use mod_util
  use define_model

  implicit none

contains

  !---------------------------------------------------------------------
  !> @brief
  !> subroutine for detecting surface edges
  !---------------------------------------------------------------------
  subroutine detect_surface_edges (M, el2ed, el2neigh, &
                                       x_start, x_end, &
                                       y_start, y_end, &
                                       z_start, z_end, &
                                       s_edges, num_s_edges)

    ! INPUT
    integer, intent(in) :: M
    integer, dimension(:,:), intent(in) :: el2ed, el2neigh
    real(kind=dp), dimension(:,:) :: x_start, x_end, &
                                     y_start, y_end, &
                                     z_start, z_end

    ! OUTPUT
    integer, allocatable, dimension(:), intent(out) :: s_edges
    integer, intent(out) :: num_s_edges

    ! LOCAL variables
    ! model dimension coordinates
    real(kind=dp) :: x_min, x_max, y_min, y_max, z_min, z_max
    integer :: l,i,j, ineigh
    integer :: allo_stat, no_ele
    integer, dimension(1000000) :: dummy_s_edges 
    ! required to increase the size of dummy_s_edges from 40k to 1m for finest mesh

    !------------------------------------------------------------------- 
    call define_model_size(x_min, x_max, y_min, y_max, z_min, z_max)

    dummy_s_edges = 0

    ! set value for not-existing neighbour-element at boundary
    ! (as in .neigh file)
    no_ele = -1

    ! Detect surface edge numbers for given mesh:

    l = 1
      do i = 1,M
        do ineigh = 1,4 ! loop over element faces
          ! if one face of the current element is a surface face
          if (el2neigh(i,ineigh) .eq. no_ele) then
              do j = 1,6    
                  if (( &

                      ! edges of constant y
                      abs(y_start(i,j) - y_min) <= eps_dp .and. &
                          x_start(i,j) .ge. x_min .and. &
                          x_start(i,j) .le. x_max .and. &
                          z_start(i,j) .ge. z_min .and. &
                          z_start(i,j) .le. z_max .and. &
                      abs(y_end(i,j) - y_min) <= eps_dp .and. &
                          x_end(i,j) .ge. x_min .and. &
                          x_end(i,j) .le. x_max .and. &
                          z_end(i,j) .ge. z_min .and. &
                          z_end(i,j) .le. z_max .or. &
                      
                      abs(y_start(i,j) - y_max) <= eps_dp .and. &
                          x_start(i,j) .ge. x_min .and. &
                          x_start(i,j) .le. x_max .and. &
                          z_start(i,j) .ge. z_min .and. &
                          z_start(i,j) .le. z_max .and. &
                      abs(y_end(i,j) - y_max) <= eps_dp .and. &
                          x_end(i,j) .ge. x_min .and. &
                          x_end(i,j) .le. x_max .and. &
                          z_end(i,j) .ge. z_min .and. &
                          z_end(i,j) .le. z_max .or. &

                      ! edges of constant x
                      abs(x_start(i,j) - x_min) <= eps_dp .and. &
                          y_start(i,j) .ge. y_min .and. &
                          y_start(i,j) .le. y_max .and. &
                          z_start(i,j) .ge. z_min .and. &
                          z_start(i,j) .le. z_max .and. &
                      abs(x_end(i,j) - x_min) <= eps_dp .and. &
                          y_end(i,j) .ge. y_min .and. &
                          y_end(i,j) .le. y_max .and. &
                          z_end(i,j) .ge. z_min .and. &
                          z_end(i,j) .le. z_max .or. &
                      
                      abs(x_start(i,j) - x_max) <= eps_dp .and. &
                          y_start(i,j) .ge. y_min .and. &
                          y_start(i,j) .le. y_max .and. &
                          z_start(i,j) .ge. z_min .and. &
                          z_start(i,j) .le. z_max .and. &
                      abs(x_end(i,j) - x_max) <= eps_dp .and. &
                          y_end(i,j) .ge. y_min .and. &
                          y_end(i,j) .le. y_max .and. &
                          z_end(i,j) .ge. z_min .and. &
                          z_end(i,j) .le. z_max .or. &
                      
                      ! edges of constant z
                      abs(z_start(i,j) - z_min) <= eps_dp .and. &
                          x_start(i,j) .ge. x_min .and. &
                          x_start(i,j) .le. x_max .and. &
                          y_start(i,j) .ge. y_min .and. &
                          y_start(i,j) .le. y_max .and. &
                      abs(z_end(i,j) - z_min) <= eps_dp .and. &
                          x_end(i,j) .ge. x_min .and. &
                          x_end(i,j) .le. x_max .and. &
                          y_end(i,j) .ge. y_min .and. &
                          y_end(i,j) .le. y_max .or. &
                      
                      abs(z_start(i,j) - z_max) <= eps_dp .and. &
                          x_start(i,j) .ge. x_min .and. &
                          x_start(i,j) .le. x_max .and. &
                          y_start(i,j) .ge. y_min .and. &
                          y_start(i,j) .le. y_max .and. &
                      abs(z_end(i,j) - z_max) <= eps_dp .and. &
                          x_end(i,j) .ge. x_min .and. &
                          x_end(i,j) .le. x_max .and. &
                          y_end(i,j) .ge. y_min .and. &
                          y_end(i,j) .le. y_max ) .and. &   
                      

                      ! write in array if it is not yey there
                      ((any(dummy_s_edges .eq. el2ed(i,j))) .eqv. &
                        .false.)) then
                          dummy_s_edges(l) = el2ed(i,j)
                          l = l+1
                  end if ! if edge is a surface edge

            end do ! edge loop
          end if ! if element has a surface face
        end do ! face loop
      end do ! element loop


    ! call Write_Message (log_unit, 'Number of surface edges:')
    ! print *, l-1

    num_s_edges = l-1

    ! Allocate array of length num_s_edges with surface edges for 
    ! storing detected surface edges
    allocate (s_edges(num_s_edges), stat = allo_stat) 
    call allocheck(log_unit, allo_stat, "error allocating array")

    s_edges = dummy_s_edges(1:num_s_edges)

  !---------------------------------------------------------------------
  end subroutine detect_surface_edges


  !---------------------------------------------------------------------
  !> @brief
  !> subroutine for applying Dirichlet boundary condition - fast.
  !> new in elfe3D_GPR, @CS
  !---------------------------------------------------------------------
  subroutine apply_DBC_fast(num_s_edges, s_edges, NNZ, Agcoo, Agcol, Agrow)

    ! INPUT
    integer, intent(in) :: num_s_edges
    integer, dimension(:), intent(in) :: s_edges
  
    ! IN & OUTPUT
    integer, intent(inout) :: NNZ
    complex(kind=dp), dimension(:), intent(inout) :: Agcoo
    integer, dimension(:), intent(inout) :: Agcol, Agrow

    ! LOCAL variables
    integer :: N_nodes, i_max, j_max
    integer :: t, i, j, l
    complex(kind=dp), allocatable, dimension(:) :: Agcoo_dummy
    integer, allocatable, dimension(:) ::  jAgcoo_dummy, iAgcoo_dummy
    integer :: allo_stat
    logical, allocatable :: is_boundary(:) ! stores boundary mask over all entries of system matrix
    integer, allocatable :: diag_idx(:)    ! stores COO index of the first diagonal entry corresponding to boundary edges

    ! N_nodes is the maximum edge number that possibly has a boundary condition
    N_nodes = maxval(s_edges)
    i_max = maxval(Agrow)
    j_max = maxval(Agcol)
    N_nodes = max(N_nodes, max(i_max,j_max)) ! just being safe

    ! Build boundary mask and reset diag_idx
    allocate(is_boundary(N_nodes))
    is_boundary = .false.
    do l = 1, num_s_edges
      if (s_edges(l)>=1 .and. s_edges(l)<=N_nodes) then ! just being safe
        is_boundary(s_edges(l)) = .true.
      end if
    end do

    allocate(diag_idx(N_nodes))
    diag_idx = 0    ! zero mean haven’t seen a diagonal yet

    ! Mark zeros for off diagonals and duplicates,
    ! and pick first diagonal where found to be boundary DOF
    do t = 1, NNZ
      i = Agrow(t)
      j = Agcol(t)

      if (i == j .and. is_boundary(i)) then
        if (diag_idx(i) == 0) then
          diag_idx(i) = t    ! remember first diag
        else
          Agcoo(t) = cmplx(0.0_dp, 0.0_dp, kind=dp)   ! duplicate diagonal
        end if

      else if (is_boundary(i) .or. is_boundary(j)) then
        Agcoo(t) = cmplx(0.0_dp, 0.0_dp, kind=dp)     ! off‐diagonal
      end if
    end do

    ! Now enforce diagonal=1 on each boundary DOF
    do l = 1, num_s_edges
      i = s_edges(l)
      if (diag_idx(i) > 0) then
        Agcoo(diag_idx(i)) = cmplx(1.0_dp, 0.0_dp, kind=dp)
      end if
    end do

     ! Delete zero entries in Agcoo

     allocate (Agcoo_dummy(NNZ), jAgcoo_dummy(NNZ), iAgcoo_dummy(NNZ), &
               stat = allo_stat)
     call allocheck(log_unit, allo_stat, "error allocating arrays")
     Agcoo_dummy = (0.0_dp,0.0_dp) 
     jAgcoo_dummy = 0 
     iAgcoo_dummy = 0 

     ! Use dummy arrays to create Agcoo, Agrow, Agcol 
     ! without zero entries:
     t = 1
     do l = 1,NNZ
      if (abs(Agcoo(l) - cmplx(0.0_dp, 0.0_dp, kind=dp)) &
         .gt. eps_dp) then

        Agcoo_dummy(t) = Agcoo(l)
        jAgcoo_dummy(t) = Agcol(l)
        iAgcoo_dummy(t) = Agrow(l)

        t = t+1

      end if

     end do

     NNZ = t-1
     Agcoo = (0.0_dp,0.0_dp)
     Agrow = 0
     Agcol = 0 

     Agcoo(1:NNZ) = Agcoo_dummy(1:NNZ)
     Agrow(1:NNZ) = iAgcoo_dummy(1:NNZ)
     Agcol(1:NNZ) = jAgcoo_dummy(1:NNZ)

     deallocate (Agcoo_dummy, jAgcoo_dummy, iAgcoo_dummy)
    !-------------------------------------------------------------------
  end subroutine apply_DBC_fast


  !---------------------------------------------------------------------
  !> @brief
  !> subroutine for applying Dirichlet boundary condition
  !---------------------------------------------------------------------
  subroutine apply_BC(num_s_edges, s_edges, NNZ, Agcoo, Agcol, Agrow)

    ! INPUT
    integer, intent(in) :: num_s_edges
    integer, dimension(:), intent(in) :: s_edges
  
    ! IN & OUTPUT
    integer, intent(inout) :: NNZ
    complex(kind=dp), dimension(:), intent(inout) :: Agcoo
    integer, dimension(:), intent(inout) :: Agcol, Agrow
    
    ! LOCAL variables
    complex(kind=dp), allocatable, dimension(:) :: Agcoo_dummy
    integer, allocatable, dimension(:) ::  jAgcoo_dummy, iAgcoo_dummy
    ! model dimension coordinates
    integer :: l,t, num_diag_entry
    integer :: allo_stat

    !-------------------------------------------------------------------
    ! call Write_Message (log_unit, 'Implementing Dirichlet BC')
     ! l = index of boundary (surface) edges
     ! t = index of edges in Agcoo

      do l = 1,num_s_edges

        num_diag_entry = 0

        do t = 1,NNZ
           ! detect diagonal entries
           if (Agcol(t) .eq. s_edges(l) .and. &
               Agrow(t) .eq. s_edges(l)) then
              if (num_diag_entry .eq. 0) then
                ! set diagonal element to (1.0,0.0)
                Agcoo(t) = cmplx(1.0_dp, 0.0_dp, kind=dp)
                !counter = counter + 1
                num_diag_entry = 1
              else if (num_diag_entry .gt. 0) then
                ! set duplicate diagonal element to (0.0,0.0)
                Agcoo(t) = cmplx(0.0_dp, 0.0_dp, kind=dp)
              end if

           else if (Agrow(t) .eq. s_edges(l) .and. &
                    Agcol(t) .ne. s_edges(l)) then
              ! set row to zero/delete the entries
              Agcoo(t) = cmplx(0.0_dp, 0.0_dp, kind=dp)

           else if (Agrow(t) .ne. s_edges(l) .and. &
                    Agcol(t) .eq. s_edges(l)) then
              ! set column to zero/delete the entries
              Agcoo(t) = cmplx(0.0_dp, 0.0_dp, kind=dp)
           end if

        end do

     end do

     ! Delete zero entries in Agcoo

     allocate (Agcoo_dummy(NNZ), jAgcoo_dummy(NNZ), iAgcoo_dummy(NNZ), &
               stat = allo_stat)
     call allocheck(log_unit, allo_stat, "error allocating arrays")
     Agcoo_dummy = (0.0_dp,0.0_dp) 
     jAgcoo_dummy = 0 
     iAgcoo_dummy = 0 

     ! Use dummy arrays to create Agcoo, Agrow, Agcol 
     ! without zero entries:
     t = 1
     do l = 1,NNZ
      if (abs(Agcoo(l) - cmplx(0.0_dp, 0.0_dp, kind=dp)) &
         .gt. eps_dp) then

        Agcoo_dummy(t) = Agcoo(l)
        jAgcoo_dummy(t) = Agcol(l)
        iAgcoo_dummy(t) = Agrow(l)

        t = t+1

      end if

     end do

     NNZ = t-1
     Agcoo = (0.0_dp,0.0_dp)
     Agrow = 0
     Agcol = 0 

     Agcoo(1:NNZ) = Agcoo_dummy(1:NNZ)
     Agrow(1:NNZ) = iAgcoo_dummy(1:NNZ)
     Agcol(1:NNZ) = jAgcoo_dummy(1:NNZ)

     deallocate (Agcoo_dummy, jAgcoo_dummy, iAgcoo_dummy)
    !-------------------------------------------------------------------
  end subroutine apply_BC
  

  !---------------------------------------------------------------------
  !> @brief
  !> subroutine for detecting perfect electric conductor (PEC) edges
  !> so far only exactly vertical boreholes, z negative downwards
  !--------------------------------------------------------------------
  subroutine detect_PEC_edges (E, nd, ed2nd, num_PEC, &
                               x_start, y_start, z_start, &
                               x_end, y_end, z_end, &
                               PEC_edges, num_PEC_edges)

    ! INPUT
    integer, intent(in) :: E, num_PEC
    real(kind=dp), dimension(:,:), intent(in) :: nd 
    integer, dimension(:,:), intent(in) :: ed2nd
    real(kind=dp),dimension(:) :: x_start, x_end, &
                                  y_start, y_end, &
                                  z_start, z_end

    ! OUTPUT
    integer, allocatable, dimension(:), intent(out) :: PEC_edges
    integer, intent(out) :: num_PEC_edges

    ! LOCAL variables
    integer :: l,i,p
    integer :: allo_stat
    integer, dimension(E) :: dummy_PEC_edges
    !-------------------------------------------------------------------
    dummy_PEC_edges = 0
    l = 0

    num_PEC_loop: do p = 1,num_PEC
      global_edge_loop: do i = 1,E  

        if(abs(nd(ed2nd(i,1),2) - y_start(p)) <= eps_dp .and. & ! y
           abs(nd(ed2nd(i,1),1) - x_start(p)) <= eps_dp .and. & ! x
           ! node 1 on PEC edge
           ! if z negative downwards
           nd(ed2nd(i,1),3) .le. z_start(p) .and. &
           nd(ed2nd(i,1),3) .ge. z_end(p) .and. &! z

           abs(nd(ed2nd(i,2),2) - y_start(p)) <= eps_dp .and. & ! y
           abs(nd(ed2nd(i,2),1) - x_start(p)) <= eps_dp .and. & ! x
           ! node 2 on PEC edge
           ! if z negative downwards
           nd(ed2nd(i,2),3) .le. z_start(p) .and. &
           nd(ed2nd(i,2),3) .ge. z_end(p)) then !z

              ! write edge number in dummy array
              l = l+1
              dummy_PEC_edges(l) = i
              
        end if
      end do global_edge_loop 
    end do num_PEC_loop


    call Write_Message (log_unit, 'Number of PEC edges:')
    print *, l
    call Write_Message (log_unit, 'PEC coordinates:')
    print*, x_start, x_end, y_start, y_end, z_start, z_end

    num_PEC_edges = l

    ! Allocate array of length num_s_edges with surface edges for 
    ! storing detected PEC edges
    allocate (PEC_edges(num_PEC_edges), stat = allo_stat) 
    call allocheck(log_unit, allo_stat, &
                  "error allocating array PEC_edges")

    PEC_edges = dummy_PEC_edges(1:num_PEC_edges)
    !print*, PEC_edges
  !---------------------------------------------------------------------
  end subroutine detect_PEC_edges


  !---------------------------------------------------------------------
  !> @brief
  !> subroutine for applying PEC boundary conditions (E = 0 on PEC edges)
  !---------------------------------------------------------------------
  subroutine apply_PEC_BC(num_PEC_edges, PEC_edges, NNZ, &
                          Agcoo, Agcol, Agrow)

    ! INPUT
    integer, intent(in) :: num_PEC_edges
    integer, dimension(:), intent(in) :: PEC_edges
  
    ! IN & OUTPUT
    integer, intent(inout) :: NNZ
    complex(kind=dp), dimension(:), intent(inout) :: Agcoo
    integer, dimension(:), intent(inout) :: Agcol, Agrow
    
    ! LOCAL variables
    complex(kind=dp), allocatable, dimension(:) :: Agcoo_dummy
    integer, allocatable, dimension(:) ::  jAgcoo_dummy, iAgcoo_dummy
    integer :: l,t, num_diag_entry 
    integer :: allo_stat

    !-------------------------------------------------------------------
    call Write_Message (log_unit, 'Implementing PECs in system matrix')
    call Write_Message (log_unit, 'Number of PEC edges:')
    print *, num_PEC_edges

     ! l = index of PEC edges
     ! t = index of edges in Agcoo

      do l = 1,num_PEC_edges

        num_diag_entry = 0

        do t = 1,NNZ
           ! detect diagonal entries
           if (Agcol(t) .eq. PEC_edges(l) .and. &
               Agrow(t) .eq. PEC_edges(l)) then
              if (num_diag_entry .eq. 0) then
                ! set diagonal element to (1.0,0.0)
                Agcoo(t) = cmplx(1.0_dp, 0.0_dp, kind=dp)
                !counter = counter + 1
                num_diag_entry = 1
              else if (num_diag_entry .gt. 0) then
                ! set duplicate diagonal element to (0.0,0.0)
                Agcoo(t) = cmplx(0.0_dp, 0.0_dp, kind=dp)
              end if

           else if (Agrow(t) .eq. PEC_edges(l) .and. &
                    Agcol(t) .ne. PEC_edges(l)) then
              ! set row to zero/delete the entries
              Agcoo(t) = cmplx(0.0_dp, 0.0_dp, kind=dp)

           else if (Agrow(t) .ne. PEC_edges(l) .and. &
                    Agcol(t) .eq. PEC_edges(l)) then
              ! set column to zero/delete the entries
              Agcoo(t) = cmplx(0.0_dp, 0.0_dp, kind=dp)
           end if

        end do

     end do



     ! Delete zero entries in Agcoo
     allocate (Agcoo_dummy(NNZ), jAgcoo_dummy(NNZ), iAgcoo_dummy(NNZ), &
               stat = allo_stat)
     call allocheck(log_unit, allo_stat, "error allocating arrays")
     Agcoo_dummy = (0.0_dp,0.0_dp) 
     jAgcoo_dummy = 0 
     iAgcoo_dummy = 0 


     ! Use dummy arrays to create Agcoo, Agrow, Agcol 
     ! without zero entries:
     t = 1
     do l = 1,NNZ
      if (abs(Agcoo(l) - cmplx(0.0_dp, 0.0_dp, kind=dp)) &
         .gt. eps_dp) then

        Agcoo_dummy(t) = Agcoo(l)
        jAgcoo_dummy(t) = Agcol(l)
        iAgcoo_dummy(t) = Agrow(l)

        t = t+1

      end if

     end do

     NNZ = t-1
     Agcoo = (0.0_dp,0.0_dp)
     Agrow = 0
     Agcol = 0 

     Agcoo(1:NNZ) = Agcoo_dummy(1:NNZ)
     Agrow(1:NNZ) = iAgcoo_dummy(1:NNZ)
     Agcol(1:NNZ) = jAgcoo_dummy(1:NNZ)

     deallocate (Agcoo_dummy, jAgcoo_dummy, iAgcoo_dummy)
    !-------------------------------------------------------------------
  end subroutine apply_PEC_BC

  !---------------------------------------------------------------------
  !> @brief
  !> new in elfe3D_GPR, @CS:
  !> subroutine for detecting PML elements and where they are located
  !----------------------------------------------------------------------
  subroutine detect_PML_element_type(M, nd, el2nd, PML_thickness, num_PML_elements, &
    PML_elements, PML_type, num_non_PML_elements, non_PML_elements)
    ! INPUT
    integer, intent(in) :: M
    integer, dimension(:,:), intent(in) :: el2nd
    real(kind=dp), dimension(:,:), intent(in) :: nd
    real(kind=dp), intent(in) :: PML_thickness
    
    ! OUTPUT
    integer, intent(out) :: num_PML_elements
    integer, allocatable, dimension(:), intent(out) :: PML_elements
    integer, allocatable, dimension(:,:), intent(out) :: PML_type
    integer, intent(out) :: num_non_PML_elements
    integer, allocatable, dimension(:), intent(out) :: non_PML_elements

    ! LOCAL variables
    integer :: l, i
    logical :: found
    integer :: allo_stat
    real(kind=dp) :: x_min, x_max, y_min, y_max, z_min, z_max
    integer, dimension(M) :: dummy_PML_elements ! element numbers
    integer, dimension(M,3) :: dummy_PML_type ! +1 or -1 for each direction
    !-------------------------------------------------------------------

    call Write_Message (log_unit, 'Detecting PML elements')
    call define_model_size(x_min, x_max, y_min, y_max, z_min, z_max)

    dummy_PML_elements = 0
    dummy_PML_type = 0
    i = 0

    global_PML_classification_loop: do l = 1,M

      found = .false.

      ! Step 1: Classify which region of the PML the element is in
      if (nd(el2nd(l,1),1) .ge. x_max - PML_thickness .and. &
        nd(el2nd(l,2),1) .ge. x_max - PML_thickness .and. &
        nd(el2nd(l,3),1) .ge. x_max - PML_thickness .and. &
        nd(el2nd(l,4),1) .ge. x_max - PML_thickness) then

          if (found .eqv. .false.) then
            i = i+1
            found = .true.
          end if
          dummy_PML_elements(i) = l
          dummy_PML_type(i,1) = 1 ! in the right

      else if (nd(el2nd(l,1),1) .le. x_min + PML_thickness .and. &
           nd(el2nd(l,2),1) .le. x_min + PML_thickness .and. &
           nd(el2nd(l,3),1) .le. x_min + PML_thickness .and. &
           nd(el2nd(l,4),1) .le. x_min + PML_thickness) then
          
          if (found .eqv. .false.) then
            i = i+1
            found = .true.
          end if
          dummy_PML_elements(i) = l
          dummy_PML_type(i,1) = -1 ! in the left
          found = .true.
      end if ! if element is in x-extents of the PML region
      
      if (nd(el2nd(l,1),2) .ge. y_max - PML_thickness .and. &
        nd(el2nd(l,2),2) .ge. y_max - PML_thickness .and. &
        nd(el2nd(l,3),2) .ge. y_max - PML_thickness .and. &
        nd(el2nd(l,4),2) .ge. y_max - PML_thickness) then
          
          if (found .eqv. .false.) then
            i = i+1
            found = .true.
          end if
          dummy_PML_elements(i) = l
          dummy_PML_type(i,2) = 1 ! in the front

      else if (nd(el2nd(l,1),2) .le. y_min + PML_thickness .and. &
           nd(el2nd(l,2),2) .le. y_min + PML_thickness .and. &
           nd(el2nd(l,3),2) .le. y_min + PML_thickness .and. &
           nd(el2nd(l,4),2) .le. y_min + PML_thickness) then

          if (found .eqv. .false.) then
            i = i+1
            found = .true.
          end if
          dummy_PML_elements(i) = l
          dummy_PML_type(i,2) = -1 ! in the back
      end if ! if element is in y-extents of the PML region
      
      if (nd(el2nd(l,1),3) .ge. z_max - PML_thickness .and. &
        nd(el2nd(l,2),3) .ge. z_max - PML_thickness .and. &
        nd(el2nd(l,3),3) .ge. z_max - PML_thickness .and. &
        nd(el2nd(l,4),3) .ge. z_max - PML_thickness) then

          if (found .eqv. .false.) then
            i = i+1
            found = .true.
          end if
          dummy_PML_elements(i) = l
          dummy_PML_type(i,3) = 1 ! in the top

      else if (nd(el2nd(l,1),3) .le. z_min + PML_thickness .and. &
           nd(el2nd(l,2),3) .le. z_min + PML_thickness .and. &
           nd(el2nd(l,3),3) .le. z_min + PML_thickness .and. &
           nd(el2nd(l,4),3) .le. z_min + PML_thickness) then

          if (found .eqv. .false.) then
            i = i+1
            found = .true.
          end if
          dummy_PML_elements(i) = l
          dummy_PML_type(i,3) = -1 ! in the bottom
      end if ! if element is in the z-extents of PML region
      
      end do global_PML_classification_loop

    ! Step 2: Remove duplicates in dummy_PML_elements
    ! and create the final PML_elements array
    ! and the PML_type array

    num_PML_elements = i
    call Write_Message (log_unit, 'Number of PML elements:')
    print *, num_PML_elements

    ! Allocate array of length num_PML_elements with PML elements for
    ! storing detected PML elements
    allocate (PML_elements(num_PML_elements), stat = allo_stat)
    call allocheck(log_unit, allo_stat, &
                  "error allocating array PML_elements")
    allocate (PML_type(num_PML_elements,3), stat = allo_stat)
    call allocheck(log_unit, allo_stat, &
                  "error allocating array PML_type")
    
    PML_elements = dummy_PML_elements(1:num_PML_elements)
    PML_type = dummy_PML_type(1:num_PML_elements,:)

    ! Step 3: Create non-PML elements array

    num_non_PML_elements = M - num_PML_elements
    call Write_Message (log_unit, 'Number of non-PML elements:')
    print *, num_non_PML_elements

    allocate (non_PML_elements(num_non_PML_elements), stat = allo_stat)
    call allocheck(log_unit, allo_stat, &
                  "error allocating array non_PML_elements")
                  
    i = 0
    non_PML_elements = 0
    non_PML_elements_loop: do l = 1,M
      if (any(PML_elements .eq. l) .eqv. .false.) then
        i = i+1
        non_PML_elements(i) = l
      end if
    end do non_PML_elements_loop
    call Write_Message (log_unit, 'To Check, number of non-PML elements:')
    print *, i

    ! Print if any of the PML_type is 0
    do l = 1,num_PML_elements
      if (PML_type(l,1) == 0 .and. PML_type(l,2) == 0 .and. &
          PML_type(l,3) == 0) then
        call Write_Message (log_unit, 'Could not classify PML element:')
        print *, PML_elements(l)
      end if
    end do

  end subroutine detect_PML_element_type


  !-----------------------------------------------------------------------------
  !> @brief
  !> new in elfe3D_GPR, @CS
  !> subroutine for evaluating the Exact PML stretching functions based on the 
  !> anisotropic formulation of the PML.
  !> Simple decay used.
  !-----------------------------------------------------------------------------
  subroutine evaluate_anisotropic_PML_simple(num_PML_elements, elem_centr, PML_thickness, &
    k_bg, PML_elements, PML_type, LAMBDA, LAMBDA_inv)
    ! INPUT
    integer, intent(in) :: num_PML_elements
    integer, dimension(:), intent(in) :: PML_elements
    integer, dimension(:,:), intent(in) :: PML_type
    complex(kind=dp), dimension(:), intent(in) :: k_bg
    real(kind=dp), dimension(:,:), intent(in) :: elem_centr
    real(kind=dp), intent(in) :: PML_thickness

    ! OUTPUT
    complex(kind=dp), allocatable, dimension(:,:), intent(out) :: LAMBDA, LAMBDA_inv
    
    ! LOCAL variables
    integer :: i, PML_i
    integer :: allo_stat
    real(kind=dp) :: x_min, x_max, y_min, y_max, z_min, z_max
    complex(kind=dp) :: k_bg_m
    complex(kind=dp) :: s_x, s_y, s_z, s_temp !local stretching functions
    real(kind=dp), dimension(3) :: elem_centr_PML_i
    !---------------------------------------------------------------
    
    call define_model_size(x_min, x_max, y_min, y_max, z_min, z_max)

    allocate (LAMBDA(num_PML_elements,3), LAMBDA_inv(num_PML_elements,3), stat = allo_stat)

    call allocheck(log_unit, allo_stat, &
                  "error allocating array PML_stretching")

    ! initialize the stretching factors, with 1.0, or no stretching
    LAMBDA = cmplx(1.0_dp, 0.0_dp, kind=dp)
    LAMBDA_inv = cmplx(1.0_dp, 0.0_dp, kind=dp)
    k_bg_m = k_bg(minloc(abs(k_bg), dim=1))

    PML_global_loop: do i = 1,num_PML_elements ! loop over PML elements

      ! get the element number and the background wavenumber
      PML_i = PML_elements(i)

      ! initialise the stretching factors with 1.0, or no stretching
      s_x = cmplx(1.0_dp, 0.0_dp, kind=dp)
      s_y = cmplx(1.0_dp, 0.0_dp, kind=dp)
      s_z = cmplx(1.0_dp, 0.0_dp, kind=dp)
      s_temp = cmplx(0.0_dp, 0.0_dp, kind=dp)
      elem_centr_PML_i = elem_centr(PML_i,:)

      ! for the right and left sides of the PML
      if (PML_type(i,1) .eq. 1) then
        s_temp = cmplx(0.0_dp, (1.0_dp/(PML_thickness - abs(elem_centr_PML_i(1) - (x_max - PML_thickness)))), kind=dp) / k_bg_m
        s_x = 1.0_dp - s_temp
      else if (PML_type(i,1) .eq. -1) then
        s_temp = cmplx(0.0_dp, (1.0_dp/(PML_thickness - abs(elem_centr_PML_i(1) - (x_min + PML_thickness)))), kind=dp) / k_bg_m
        s_x = 1.0_dp - s_temp
      end if
      
      ! for the front and back sides of the PML
      if (PML_type(i,2) .eq. 1) then
          s_temp = cmplx(0.0_dp, (1.0_dp/(PML_thickness - abs(elem_centr_PML_i(2) - (y_max - PML_thickness)))), kind=dp) / k_bg_m
          s_y = 1.0_dp - s_temp
      else if (PML_type(i,2) .eq. -1) then
          s_temp = cmplx(0.0_dp, (1.0_dp/(PML_thickness - abs(elem_centr_PML_i(2) - (y_min + PML_thickness)))), kind=dp) / k_bg_m
          s_y = 1.0_dp - s_temp
      end if

      ! for the top and bottom sides of the PML
      if (PML_type(i,3) .eq. 1) then
          s_temp = cmplx(0.0_dp, (1.0_dp/(PML_thickness - abs(elem_centr_PML_i(3) - (z_max - PML_thickness)))), kind=dp) / k_bg_m
          s_z = 1.0_dp - s_temp
      else if (PML_type(i,3) .eq. -1) then
          s_temp = cmplx(0.0_dp, (1.0_dp/(PML_thickness - abs(elem_centr_PML_i(3) - (z_min + PML_thickness)))), kind=dp) / k_bg_m
          s_z = 1.0_dp - s_temp
      end if

      LAMBDA(i,1) = s_y * s_z / s_x
      LAMBDA(i,2) = s_x * s_z / s_y
      LAMBDA(i,3) = s_x * s_y / s_z
      LAMBDA_inv(i,1) = s_x / (s_y * s_z)
      LAMBDA_inv(i,2) = s_y / (s_x * s_z)
      LAMBDA_inv(i,3) = s_z / (s_x * s_y)

    end do PML_global_loop

  end subroutine evaluate_anisotropic_PML_simple


  !-----------------------------------------------------------------------------
  !> @brief
  !> new in elfe3D_GPR, @CS
  !> subroutine for evaluating the Exact PML stretching functions based on the 
  !> anisotropic formulation of the PML.
  !> Logarithmic decay used.
  !-----------------------------------------------------------------------------
  subroutine evaluate_anisotropic_PML_log(num_PML_elements, elem_centr, PML_thickness, &
    k_bg, PML_elements, PML_type, LAMBDA, LAMBDA_inv)
    
    ! INPUT
    integer, intent(in) :: num_PML_elements
    integer, dimension(:), intent(in) :: PML_elements
    integer, dimension(:,:), intent(in) :: PML_type
    complex(kind=dp), dimension(:), intent(in) :: k_bg
    real(kind=dp), dimension(:,:), intent(in) :: elem_centr
    real(kind=dp), intent(in) :: PML_thickness

    ! OUTPUT
    complex(kind=dp), allocatable, dimension(:,:), intent(out) :: LAMBDA, LAMBDA_inv
    
    ! LOCAL variables
    integer :: i, PML_i
    integer :: allo_stat
    real(kind=dp) :: x_min, x_max, y_min, y_max, z_min, z_max
    complex(kind=dp) :: s_x, s_y, s_z, s_temp !local stretching functions
    complex(kind=dp) :: k_bg_m
    real(kind=dp), dimension(3) :: elem_centr_PML_i
    !---------------------------------------------------------------
    
    call define_model_size(x_min, x_max, y_min, y_max, z_min, z_max)

    allocate (LAMBDA(num_PML_elements,3), LAMBDA_inv(num_PML_elements,3), stat = allo_stat)

    call allocheck(log_unit, allo_stat, &
                  "error allocating array PML_stretching")

    ! initialize the stretching factors, with 1.0, or no stretching
    LAMBDA = cmplx(1.0_dp, 0.0_dp, kind=dp)
    LAMBDA_inv = cmplx(1.0_dp, 0.0_dp, kind=dp)

    PML_global_loop: do i = 1,num_PML_elements ! loop over PML elements

      ! get the element number
      PML_i = PML_elements(i)
      k_bg_m = k_bg(PML_i)

      ! initialise the stretching factors with 1.0, or no stretching
      s_x = cmplx(1.0_dp, 0.0_dp, kind=dp)
      s_y = cmplx(1.0_dp, 0.0_dp, kind=dp)
      s_z = cmplx(1.0_dp, 0.0_dp, kind=dp)
      s_temp = cmplx(0.0_dp, 0.0_dp, kind=dp)
      elem_centr_PML_i = elem_centr(PML_i,:)

      ! for the right and left sides of the PML
      if (PML_type(i,1) .eq. 1) then
        s_temp = cmplx(0.0_dp, log(1.0_dp - ((abs(elem_centr_PML_i(1) - (x_max - PML_thickness)))/PML_thickness)), kind=dp) / k_bg_m
        s_x = 1.0_dp - s_temp
      else if (PML_type(i,1) .eq. -1) then
        s_temp = cmplx(0.0_dp, log(1.0_dp - ((abs(elem_centr_PML_i(1) - (x_min + PML_thickness)))/PML_thickness)), kind=dp) / k_bg_m
        s_x = 1.0_dp - s_temp
      end if
      
      ! for the front and back sides of the PML
      if (PML_type(i,2) .eq. 1) then
          s_temp = cmplx(0.0_dp, log(1.0_dp - ((abs(elem_centr_PML_i(2) - (y_max - PML_thickness)))/PML_thickness)), kind=dp) / k_bg_m
          s_y = 1.0_dp - s_temp
      else if (PML_type(i,2) .eq. -1) then
          s_temp = cmplx(0.0_dp, log(1.0_dp - ((abs(elem_centr_PML_i(2) - (y_min + PML_thickness)))/PML_thickness)), kind=dp) / k_bg_m
          s_y = 1.0_dp - s_temp
      end if

      ! for the top and bottom sides of the PML
      if (PML_type(i,3) .eq. 1) then
          s_temp = cmplx(0.0_dp, log(1.0_dp - ((abs(elem_centr_PML_i(3) - (z_max - PML_thickness)))/PML_thickness)), kind=dp) / k_bg_m
          s_z = 1.0_dp - s_temp
      else if (PML_type(i,3) .eq. -1) then
          s_temp = cmplx(0.0_dp, log(1.0_dp - ((abs(elem_centr_PML_i(3) - (z_min + PML_thickness)))/PML_thickness)), kind=dp) / k_bg_m
          s_z = 1.0_dp - s_temp
      end if

      LAMBDA(i,1) = s_y * s_z / s_x
      LAMBDA(i,2) = s_x * s_z / s_y
      LAMBDA(i,3) = s_x * s_y / s_z
      LAMBDA_inv(i,1) = s_x / (s_y * s_z)
      LAMBDA_inv(i,2) = s_y / (s_x * s_z)
      LAMBDA_inv(i,3) = s_z / (s_x * s_y)

    end do PML_global_loop

  end subroutine evaluate_anisotropic_PML_log


  !-----------------------------------------------------------------------------
  !> @brief
  !> new in elfe3D_GPR, @CS
  !> subroutine for evaluating the Exact PML stretching functions based on the 
  !> anisotropic formulation of the PML.
  !> Exponential decay used.
  !-----------------------------------------------------------------------------
  subroutine evaluate_anisotropic_PML_exp(num_PML_elements, elem_centr, PML_thickness, &
    k_bg, PML_elements, PML_type, LAMBDA, LAMBDA_inv)
    ! INPUT
    integer, intent(in) :: num_PML_elements
    integer, dimension(:), intent(in) :: PML_elements
    integer, dimension(:,:), intent(in) :: PML_type
    complex(kind=dp), dimension(:), intent(in) :: k_bg
    real(kind=dp), dimension(:,:), intent(in) :: elem_centr
    real(kind=dp), intent(in) :: PML_thickness

    ! OUTPUT
    complex(kind=dp), allocatable, dimension(:,:), intent(out) :: LAMBDA, LAMBDA_inv
    
    ! LOCAL variables
    integer :: i, PML_i
    integer :: allo_stat
    real(kind=dp) :: x_min, x_max, y_min, y_max, z_min, z_max
    complex(kind=dp) :: s_x, s_y, s_z, s_temp !local stretching functions
    complex(kind=dp) :: k_bg_m
    real(kind=dp), dimension(3) :: elem_centr_PML_i
    !---------------------------------------------------------------
    
    call define_model_size(x_min, x_max, y_min, y_max, z_min, z_max)

    allocate (LAMBDA(num_PML_elements,3), LAMBDA_inv(num_PML_elements,3), stat = allo_stat)

    call allocheck(log_unit, allo_stat, &
                  "error allocating array PML_stretching")

    ! initialize the stretching factors, with 1.0, or no stretching
    LAMBDA = cmplx(1.0_dp, 0.0_dp, kind=dp)
    LAMBDA_inv = cmplx(1.0_dp, 0.0_dp, kind=dp)

    PML_global_loop: do i = 1,num_PML_elements ! loop over PML elements

      ! get the element number
      PML_i = PML_elements(i)
      k_bg_m = k_bg(PML_i)

      ! initialise the stretching factors with 1.0, or no stretching
      s_x = cmplx(1.0_dp, 0.0_dp, kind=dp)
      s_y = cmplx(1.0_dp, 0.0_dp, kind=dp)
      s_z = cmplx(1.0_dp, 0.0_dp, kind=dp)
      s_temp = cmplx(0.0_dp, 0.0_dp, kind=dp)
      elem_centr_PML_i = elem_centr(PML_i,:)

      ! for the right and left sides of the PML
      if (PML_type(i,1) .eq. 1) then
        s_temp = cmplx(0.0_dp, (exp(5 * abs(elem_centr_PML_i(1) - (x_max - PML_thickness))/PML_thickness)), kind=dp) / k_bg_m
        s_x = 1.0_dp - s_temp
      else if (PML_type(i,1) .eq. -1) then
        s_temp = cmplx(0.0_dp, (exp(5 * abs(elem_centr_PML_i(1) - (x_min + PML_thickness))/PML_thickness)), kind=dp) / k_bg_m
        s_x = 1.0_dp - s_temp
      end if
      
      ! for the front and back sides of the PML
      if (PML_type(i,2) .eq. 1) then
          s_temp = cmplx(0.0_dp, (exp(5 * abs(elem_centr_PML_i(2) - (y_max - PML_thickness))/PML_thickness)), kind=dp) / k_bg_m
          s_y = 1.0_dp - s_temp
      else if (PML_type(i,2) .eq. -1) then
          s_temp = cmplx(0.0_dp, (exp(5 * abs(elem_centr_PML_i(2) - (y_min + PML_thickness))/PML_thickness)), kind=dp) / k_bg_m
          s_y = 1.0_dp - s_temp
      end if

      ! for the top and bottom sides of the PML
      if (PML_type(i,3) .eq. 1) then
          s_temp = cmplx(0.0_dp, (exp(5 * abs(elem_centr_PML_i(3) - (z_max - PML_thickness))/PML_thickness)), kind=dp) / k_bg_m
          s_z = 1.0_dp - s_temp
      else if (PML_type(i,3) .eq. -1) then
          s_temp = cmplx(0.0_dp, (exp(5 * abs(elem_centr_PML_i(3) - (z_min + PML_thickness))/PML_thickness)), kind=dp) / k_bg_m
          s_z = 1.0_dp - s_temp
      end if

      LAMBDA(i,1) = s_y * s_z / s_x
      LAMBDA(i,2) = s_x * s_z / s_y
      LAMBDA(i,3) = s_x * s_y / s_z
      LAMBDA_inv(i,1) = s_x / (s_y * s_z)
      LAMBDA_inv(i,2) = s_y / (s_x * s_z)
      LAMBDA_inv(i,3) = s_z / (s_x * s_y)

    end do PML_global_loop

  end subroutine evaluate_anisotropic_PML_exp


  !-----------------------------------------------------------------------------
  !> @brief
  !> new in elfe3D_GPR, @CS
  !> subroutine for evaluating the Exact PML stretching functions based on the 
  !> anisotropic formulation of the PML.
  !> Polynomial decay used.
  !-----------------------------------------------------------------------------
  subroutine evaluate_anisotropic_PML_poly(num_PML_elements, elem_centr, PML_thickness, &
    k_bg, PML_elements, PML_type, LAMBDA, LAMBDA_inv)
    ! INPUT
    integer, intent(in) :: num_PML_elements
    integer, dimension(:), intent(in) :: PML_elements
    integer, dimension(:,:), intent(in) :: PML_type
    complex(kind=dp), dimension(:), intent(in) :: k_bg
    real(kind=dp), dimension(:,:), intent(in) :: elem_centr
    real(kind=dp), intent(in) :: PML_thickness

    ! OUTPUT
    complex(kind=dp), allocatable, dimension(:,:), intent(out) :: LAMBDA, LAMBDA_inv
    
    ! LOCAL variables
    integer :: i, PML_i
    integer :: allo_stat
    integer :: m ! order of polynomial
    real(kind=dp) :: x_min, x_max, y_min, y_max, z_min, z_max
    complex(kind=dp) :: s_x, s_y, s_z, s_temp !local stretching functions
    complex(kind=dp) :: k_bg_m
    real(kind=dp), dimension(3) :: elem_centr_PML_i
    !---------------------------------------------------------------
    
    call define_model_size(x_min, x_max, y_min, y_max, z_min, z_max)

    allocate (LAMBDA(num_PML_elements,3), LAMBDA_inv(num_PML_elements,3), stat = allo_stat)

    call allocheck(log_unit, allo_stat, &
                  "error allocating array PML_stretching")

    ! initialize the stretching factors, with 1.0, or no stretching
    LAMBDA = cmplx(1.0_dp, 0.0_dp, kind=dp)
    LAMBDA_inv = cmplx(1.0_dp, 0.0_dp, kind=dp)
    m = 4 ! order of polynomial

    PML_global_loop: do i = 1,num_PML_elements ! loop over PML elements

      ! get the element number
      PML_i = PML_elements(i)
      k_bg_m = k_bg(PML_i)

      ! initialise the stretching factors with 1.0, or no stretching
      s_x = cmplx(1.0_dp, 0.0_dp, kind=dp)
      s_y = cmplx(1.0_dp, 0.0_dp, kind=dp)
      s_z = cmplx(1.0_dp, 0.0_dp, kind=dp)
      s_temp = cmplx(0.0_dp, 0.0_dp, kind=dp)
      elem_centr_PML_i = elem_centr(PML_i,:)

      ! for the right and left sides of the PML
      if (PML_type(i,1) .eq. 1) then
        s_temp = cmplx(0.0_dp, (sqrt(c_0)*(abs(elem_centr_PML_i(1) - (x_max - PML_thickness))/PML_thickness)**m), kind=dp) / k_bg_m
        s_x = 1.0_dp - s_temp
      else if (PML_type(i,1) .eq. -1) then
        s_temp = cmplx(0.0_dp, (sqrt(c_0)*(abs(elem_centr_PML_i(1) - (x_min + PML_thickness))/PML_thickness)**m), kind=dp) / k_bg_m
        s_x = 1.0_dp - s_temp
      end if
      
      ! for the front and back sides of the PML
      if (PML_type(i,2) .eq. 1) then
          s_temp = cmplx(0.0_dp, (sqrt(c_0)*(abs(elem_centr_PML_i(2) - (y_max - PML_thickness))/PML_thickness)**m), kind=dp) / k_bg_m
          s_y = 1.0_dp - s_temp
      else if (PML_type(i,2) .eq. -1) then
          s_temp = cmplx(0.0_dp, (sqrt(c_0)*(abs(elem_centr_PML_i(2) - (y_min + PML_thickness))/PML_thickness)**m), kind=dp) / k_bg_m
          s_y = 1.0_dp - s_temp
      end if

      ! for the top and bottom sides of the PML
      if (PML_type(i,3) .eq. 1) then
          s_temp = cmplx(0.0_dp, (sqrt(c_0)*(abs(elem_centr_PML_i(3) - (z_max - PML_thickness))/PML_thickness)**m), kind=dp) / k_bg_m
          s_z = 1.0_dp - s_temp
      else if (PML_type(i,3) .eq. -1) then
          s_temp = cmplx(0.0_dp, (sqrt(c_0)*(abs(elem_centr_PML_i(3) - (z_min + PML_thickness))/PML_thickness)**m), kind=dp) / k_bg_m
          s_z = 1.0_dp - s_temp
      end if

      LAMBDA(i,1) = s_y * s_z / s_x
      LAMBDA(i,2) = s_x * s_z / s_y
      LAMBDA(i,3) = s_x * s_y / s_z
      LAMBDA_inv(i,1) = s_x / (s_y * s_z)
      LAMBDA_inv(i,2) = s_y / (s_x * s_z)
      LAMBDA_inv(i,3) = s_z / (s_x * s_y)

    end do PML_global_loop

  end subroutine evaluate_anisotropic_PML_poly


end module BC
