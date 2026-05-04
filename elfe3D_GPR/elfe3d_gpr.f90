!> \file elfe3d_gpr.f90
!> \brief Main program of elfe3D_GPR
!!
!> \details Main executable for 3D ground-penetrating radar (GPR) forward modeling.
!> \details Reads TetGen mesh input files, assembles the finite-element system,
!> \details applies boundary conditions, optionally performs refinement, solves
!> \details the linear system with MUMPS, and computes field components at receivers.
!!
!> \author Paula Rulff and Chaitanya Dinesh Singh
!> \since 2018
!!
!> Developed by Paula Rulff and Chaitanya Dinesh Singh, started 16/07/2018
!> Last change: July 2025
!!
!> Program to run 3D forward modelling for GPR setups
!> using tetrahedral meshes and first-order finite-element approximation
!!
!> The program is designed 
!> to read a mesh from ASCII node, mesh and element files 
!> created with TetGen,
!!
!> to calculate linear-edge based interpolation functions
!!
!> to build global system matrix and source term,
!!
!> to apply boundary conditions,
!!
!> to optionally run refinement based on a dual problem,
!!
!> to solve the system of equations with the direct solver MUMPS
!!
!> to calculate electric and magnetic fields components 
!> at given positions.
!!
!! Specify modelling input parameters in ./in/elfe3D_input.txt
!!
!> Compile with: the provided Makefile, tested with gfortran compiler.
!!
!> Provide mesh files in tetgen format.    
!!
!>
!>
!>  Copyright (C) Paula Rulff and Chaitanya Dinesh Singh 2024
!>  This file is part of elfe3D_GPR.
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

program elfe3d_gpr

  ! Include modules
  use mod_util
  use vector_products
  use define_model
  use read_mesh
  use calculate_matrices
  use model_parameters
  use interp_functions
  use calculate_local_left
  use calculate_global_source
  use BC
  use sparse_matrix_operations
  use solvers
  use calculate_tf
  use error_estimates
  use tetgen_operations
  use omp_lib


  implicit none 

  !---------------------------------------------------------------------
  ! Variable definitions
  !---------------------------------------------------------------------
  ! Counter variables
  integer :: i,j,l,p_l
  integer :: numfreq
  integer :: refStep


  ! Mesh filenames
  character(len = 255) :: NodeFile, EdgeFile, ElementFile, NeighFile, &
                          NewVolumeFile
  character(len = 50) :: MeshFileName

  ! Output filenames
  character(len = 500) :: EFile, HFile

  ! Output field arrays
  complex(kind=dp), allocatable, dimension(:,:,:) :: EFields
  complex(kind=dp), allocatable, dimension(:,:,:) :: HFields

  ! Allocation status check
  integer :: allo_stat

  ! Error check
  integer :: ierr

  ! Number of frequencies
  integer :: Nfreq 

  ! Frequencies
  real(kind=dp), allocatable, dimension(:) :: freq 

  ! Angular frequency 
  real(kind=dp) :: w

  ! Number of nodes
  integer :: N

  ! Number of edges
  integer :: E

  ! Number of elements
  integer :: M

  ! Declaration of array for node coordinates N times 3 (x,y,z)
  real(kind=dp), allocatable, dimension(:,:) :: nd

  ! Declaration of array for edge nodes E times 2
  integer, allocatable, dimension(:,:) :: ed2nd

  ! Declaration of array for element nodes M times 4
  integer, allocatable, dimension(:,:) :: el2nd

  ! Declaration of array for element edge lengths M times 6
  real(kind=dp), allocatable, dimension(:,:) :: el2edl

  ! Declaration of array for edge lengths
  real(kind=dp), allocatable, dimension(:) :: edl

  ! Declaration of array for element edges M times 6
  integer, allocatable, dimension(:,:) :: el2ed

  ! Declaration of array for element neighbours M times 4
  integer, allocatable, dimension(:,:) :: el2neigh

  ! Array of local edge signs M times 6
  real(kind=dp), allocatable, dimension(:,:) :: ed_sign
  
  ! Declaration of array for node boundary marker N times 1
  integer, allocatable, dimension(:) :: nodemarker
  
  ! Declaration of array for edge boundary marker E times 1
  integer, allocatable, dimension(:) :: edgemarker
  
  ! Declaration of array for element attributes M times 1
  integer, allocatable, dimension(:) :: eleattr

  !> \note New in elfe3D_GPR: element centroids, array of size M times 3
  real (kind=dp), allocatable, dimension(:,:) :: elem_centr

  ! Declaration of coordinate matrices and linear interpolation
  ! function arrays
  real(kind=dp), allocatable, dimension(:,:) :: x, y, z, a, b, c, d, &
       a_start, a_end, b_start, b_end, c_start, c_end, d_start, d_end
  
  ! Source start and end coordinates
  real(kind=dp) :: sx_start, sy_start, sz_start, sx_end, sy_end, sz_end

  ! Declaration of element volumes
  real(kind=dp), allocatable, dimension(:) :: ve

  ! Region attributes, model parameters
  real(kind=dp), allocatable, dimension(:) :: rho, mu, eps, mu_r

  !> \note New in elfe3D_GPR: High-frequency specific effective permittivity
  complex(kind=dp), allocatable, dimension(:) :: epsilon_r_eff

  !> \note New in elfe3D_GPR: High-frequency specific wavenumber
  real(kind=dp) :: k_0 ! free space wave number
  complex(kind=dp), allocatable, dimension(:) :: k_bg ! complete mesh wavenumbers

  ! Declaration of local stiffness and mass matrix K, MM
  !> \note New in elfe3D_GPR: made complex for complex material parameters
  complex(kind=dp), allocatable, dimension(:,:) :: K, MM

  ! Declaration of local left-handside matrix AA
  complex(kind=dp), allocatable, dimension(:,:) :: AA

  ! Declaration of current strength vector Is and
  ! interpolation functions on start and end nodes of all elements edges
  real(kind=dp), allocatable, dimension(:,:) :: Is, &
       x_start, x_end, y_start, y_end, z_start, z_end

  ! Declaration of global source vectors of primal (Bg) and dual (Cg) 
  ! problems
  ! and RHS matrix (E x 2) containing Bg and Cg for solver routine
  complex(kind=dp), allocatable, dimension(:) :: Bg, Cg
  complex(kind=dp), allocatable, dimension(:,:) :: RHS

  ! Declaration of vector for source calculation
  real(kind=dp), allocatable, dimension(:) :: source
  ! Total source length
  real(kind=dp) :: source_length

  ! Source midpoint (x,y,z)
  real(kind=dp), dimension(3) :: midp_source

  ! Surface edges vector and PEC edges vector
  integer, allocatable, dimension(:) :: s_edges, PEC_edges

  !> \note New in elfe3D_GPR: PML parameters
  integer :: PML_present
  real(kind=dp) :: PML_thickness
  integer :: PML_decay_type ! 1: Simple Reciprocal, 2: Logarithmic, 3: Exponential, 4: Polynomial

  !> \note New in elfe3D_GPR: PML elements and type
  integer :: num_PML_elements, num_non_PML_elements
  integer, allocatable, dimension(:) :: PML_elements      ! element numbers containing PML
  integer, allocatable, dimension(:) :: non_PML_elements  ! element numbers not containing PML
  integer, allocatable, dimension(:,:) :: PML_type        ! which of the three components get stretched at which element

  !> \note New in elfe3D_GPR: PML stretching functions
  !> For all elements, all components of the tensor
  complex(kind=dp), allocatable, dimension(:,:) :: LAMBDA, LAMBDA_inv

  !> \note New in elfe3D_GPR: complex stretched factors for local matrices
  !> Per element, all components of the tensor
  complex(kind=dp), allocatable, dimension(:) :: LAMBDA_l, LAMBDA_inv_l
  ! check per element to see if an element is in the PML or not
  integer :: stretched_element

  ! Matrix arrays
  ! Declaration of compressed system matrix Agcsr (values)
  complex(kind=dp), allocatable, dimension(:) :: Agcsr
  ! Declaration of compressed system matrix Agcsr (column & pointer)
  integer, allocatable, dimension(:) :: jAgcsr,iAgcsr

  ! Declaration of compressed system matrix Agcsr in 
  ! upper triangle format (values)
  complex(kind=dp), allocatable, dimension(:) :: Agcsr_o, Agcsr_oo
  ! Declaration of compressed system matrix Agcsr in 
  ! upper triangle format (column and pointer)
  integer, allocatable, dimension(:) :: jAgcsr_o,iAgcsr_o, jAgcsr_oo
  ! for MUMPS in COO format
  integer, allocatable, dimension(:) :: iAgcoo_o
  ! Variables for assembling global Matrix in COO format
  complex(kind=dp), allocatable, dimension(:) :: Agcoo
  integer, allocatable, dimension(:) :: Agrow, Agcol
  ! test for matrix plot
  ! Number of nonzero elements, counter
  integer :: NNZ, NNZ_c, NNZ_new, icount

  ! Declaration of solution vector S (primal problem) and 
  ! W (dual problem)
  ! and Solution matrix (E x 2) containing S and Wg for solver routines
  complex (kind=dp), allocatable, dimension(:) :: S,Wg
  complex (kind=dp), allocatable, dimension(:,:) :: Solution

  ! Variables to calculate field components at receiver
  ! Elements containing receiver
  integer, allocatable, dimension(:) :: rec1_el 
  ! Vector containing Ex, Ey and Ez component at receiver site
  complex(kind=dp), dimension(3) :: E_rec1 
  ! Vector containing Hx,Hy and Hz component at receiver site
  complex(kind=dp), dimension(3) :: H_rec1 
  !> \note New in elfe3D_GPR: array containing domain E field components at element centroids
  complex(kind=dp), allocatable, dimension(:,:) :: domain_Efields
  !> \note New in elfe3D_GPR: array containing domain H field components at element centroids
  complex(kind=dp), allocatable, dimension(:,:) :: domain_Hfields
  
  ! Receiver coordinates
  real(kind=dp), allocatable, dimension(:) :: u1, v1, w1 
  ! Number of receivers
  integer :: num_rec

  ! Source parameters
  ! Source type and direction
  integer :: CSTYPE, direction
  ! Source moment
  real(kind=dp) :: p
  ! Variables for influence sources (RHS of dual problem)
  ! edges of the element containing receiver
  integer, dimension(6) :: rec1_ed 
  complex(kind=dp), dimension(6) :: infl_source

  ! Number of source edges and number of surface edges
  integer:: num_source_edges, num_s_edges

  ! PEC 
  ! (perfect electric conductor for borehole casing present (PEC = 1))
  integer :: PEC, num_PEC, num_PEC_edges
  real(kind=dp) , allocatable, dimension(:,:) :: PEC_start, PEC_end


  ! Variables for the refinement loop:
  ! max number of refinement steps, max number of unknowns
  integer :: maxRefSteps, maxUnknowns
  ! accuracy tolerance for refinement, elements to be refined based
  ! on accuracy above factor beta
  real(kind=dp) :: accuracyTol, betaRef
  ! termination criterion of refinement loop
  integer :: terminationCrit
  ! write out vtk files (1/0)
  integer :: vtk
  
  ! Variables for error estimator calculation:
  ! Elemental residuals primal and dual problem
  real(kind=dp), allocatable, dimension(:) :: res_s, res_w
  ! Elemental face jumps primal and dual problem
  real(kind=dp), allocatable, dimension(:) :: fjJ_s, fjJ_w, fjH_s, fjH_w 
  ! Elemental error estimator
  real(kind=dp), allocatable, dimension(:) :: errorEst
  ! Elemental error erstimator for primal and dual problem
  real(kind=dp), allocatable, dimension(:) :: eta_s, eta_w
  ! Sum of all elemental error estimators, average relative error 
  ! estimate, average error estimate at receiver sites
  real(kind=dp) :: SumerrorEst, RelErr, RecErr
  ! residuals (1), residuals and face jumps J (2), 
  ! residuals and face jumps J and B (3)
  integer :: errorEst_method

  ! Variables for refinement:
  ! Number of unknowns of previous refinement step
  integer :: E_before = 1 
  !integer :: mi
  integer :: refStrategy
  character(len = 50) :: StringStep, StringEnding

  !> \note New in elfe3D_GPR: input information - domain field components to be saved in VTK file
  integer :: fields_vtk

  ! time testing
  real(kind=dp) :: start, finish, seconds,seconds_solve

  ! solver type
  integer :: solver
  ! Variables for MUMPS and pseudo MPI
  logical :: flag

  !---------------------------------------------------------------------
  call cpu_time(start)  ! CPU time measurement start
  seconds = omp_get_wtime ( ) ! Wall time measurments start
  !---------------------------------------------------------------------
  call Write_Message (log_unit, '***************************************')
  call Write_Message (log_unit, 'GPR forward modelling with elfe3D_GPR')
  call Write_Message (log_unit, '***************************************')
  !---------------------------------------------------------------------
  ! Reading solver type
  call define_solver(solver)
  call Write_Message (log_unit, &
          'Your forward problem is solved with')
  select case (solver)
    case(1)
      print*, 'PARDISO'
    case(2)
      print*, 'MUMPS'
  end select
  
  ! Reading refinement information
  call Write_Message (log_unit, '*************************************')
  call Write_Message (log_unit, 'Reading refinement information')
  call define_refinement(maxRefSteps, maxUnknowns, accuracyTol, &
                            betaRef, vtk, errorEst_method, refStrategy)
  
  call Write_Message (log_unit, 'Your mesh will be refined')
  print *,'for  ',maxRefSteps,'refinement steps'
  print *,'or until',maxUnknowns,'unknowns'
  print *,'with error estim. method',errorEst_method
  print *,'with refinement strategy',refStrategy


  ! Initialise refinement step
  refStep = 0
  
  ! termination criterion of loop (0 = continue, 1 = terminate)
  terminationCrit = 0

  ! Open Output files
  !> \note New in elfe3D_GPR: fields_vtk controls VTK output generation
  call define_output(EFile, HFile, num_rec, fields_vtk)
  open (unit = (50+1), file = trim(EFile)//".txt") ! electric fields
  open (unit = (50+2), file = trim(HFile)//".txt") ! magnetic fields
  ! electric fields ordered along receiver line
  open (unit = (50+3), file = trim(EFile)//"_receiver_line.txt") 
  ! magnetic fields ordered along receiver line
  open (unit = (50+4), file = trim(HFile)//"_receiver_line.txt") 

  !> \note New in elfe3D_GPR
  call Write_Message (log_unit, &
     'Your output files will be generated in: /out')
  if (fields_vtk == 1) then
    call Write_Message (log_unit, &
     'Field components in the domain will be in *.vtk file in: /in')
  end if

  if (vtk == 1) then
    ! Refinement-info 
    open(unit=9999, file = 'out/elfe3D_refinement_Info.txt') 
  end if
  
  !---------------------------------------------------------------------
  ! Loop for refinement
  !---------------------------------------------------------------------
  refinement_loop: do while (terminationCrit == 0)

  ! increment current refinement step:
  refStep = refStep + 1

  !---------------------------------------------------------------------
  ! I: read mesh files:
  !---------------------------------------------------------------------
  ! Reading model definitions
  call Write_Message (log_unit, '*************************************')
  call Write_Message (log_unit, 'Reading model definitions')
  call define_mesh(NodeFile, EdgeFile, ElementFile, NeighFile, &
                                                  refStep, MeshFileName)

  ! Reading coordinates for all nodes from file nodes.txt 
  ! output in nd (N times 3)
  call read_nodes(NodeFile, N, nd, nodemarker)

  ! Reading edges with corresponding node numbers from file edges.txt
  call read_edges(EdgeFile, E, ed2nd, edgemarker)

  ! Reading elements with corresponding node numbers from file elem.txt
  call read_elements(ElementFile, M, el2nd, eleattr)

  ! Reading element-neighbours with corresponding element numbers 
  ! from file neigh.txt
  call read_neigh(NeighFile, M, el2neigh)
  
  call Write_Message (log_unit, 'Your input meshfile is called')
  print *, MeshFileName
  call Write_Message (log_unit, 'Your mesh has')
  print *,N,'Nodes'
  print *,M,'Elements'
  print *,E,'Edges (dof)'

  !---------------------------------------------------------------------
  ! II: assemble connectivity arrays
  !---------------------------------------------------------------------
  call Write_Message (log_unit, '*************************************')
  call Write_Message (log_unit, 'Assembling connectivity arrays')
  ! Creating element2edge connectivity array
  call calc_connect_matrix(el2nd, ed2nd, E, M, el2ed)

  ! Calculate arrays for coordinates and _start & _end 
  ! coordinates for edges
  ! and set up coordinate_start and _end array
  call calc_coord_matrices(M, nd, el2nd, x, y, z, &
       x_start, x_end, y_start, y_end, z_start, z_end)
  ! Create edge sign array
  call calc_edge_signs(M, el2nd, ed_sign)
  
  ! --------------------------------------------------------------------
  ! III: calculate linear interpolation functions
  !---------------------------------------------------------------------
  ! Calculating coefficients for linear interpolation functions beween 
  ! nodes, put them into arrays and
  ! set up start and end (node 1 and 2) arrays for a,b,c,d  
  call Write_Message (log_unit, &
    'Calculating interpolation coefficients')
  call calc_abcd (M, x, y, z, a, b, c, d, &
       a_start, a_end, b_start, b_end, c_start, c_end, d_start, d_end)

  !---------------------------------------------------------------------
  ! IV: obtain model properties
  !---------------------------------------------------------------------
  ! Calculate volume of each element
  call Write_Message (log_unit, 'Obtaining model properties')
  call calc_vol (M, x, y, z, ve)

  ! Calculate length of each edge for each element (M times 6 matrix)
  call calc_edge_length (M, E, nd, ed2nd, x, y, z, el2edl,edl)
  ! Coordinate matrices x,y,z are not needed anymore at this point

  ! Get region attributes from eleattr and assign to model parameters
  call read_model_param(eleattr, M, rho, mu, eps)

  !> \note New in elfe3D_GPR: high-frequency specific effective permittivity
  allocate (epsilon_r_eff(M), stat = allo_stat)
  call allocheck(log_unit, allo_stat, &
                 "error allocating array epsilon_r_eff")

  !> \note New in elfe3D_GPR: high-frequency specific wave number
  allocate (k_bg(M), stat = allo_stat)
  call allocheck(log_unit, allo_stat, &
                 "error allocating array k_bg")

  call Write_Message (log_unit, '*************************************')
  ! Read source definitions
  call define_source(sx_start, sy_start, sz_start, &
                     sx_end, sy_end, sz_end, &
                     CSTYPE, direction, p, midp_source)
  call Write_Message (log_unit, 'Source parameters:')
  print *, 'CSTYPE'                  , CSTYPE
  print *, 'source-start coordinates:'
  print *, sx_start, sy_start, sz_start
  print *, 'source-end coordinates  '
  print *, sx_end, sy_end, sz_end
  print *, 'source moment            ', p

  ! Allocating global vector source: indicating location (edges), 
  ! type and direction of source
  allocate (source(E), stat = allo_stat)
  call allocheck(log_unit, allo_stat, "Error allocating array source")

  source = (0.0_dp)

  select case (CSTYPE)
    case (0) ! HED_x
      call Write_Message (log_unit, 'HED in x-direction')
      call HED_x (E, direction, nd, ed2nd, &
                  sx_start, sy_start, sz_start, sx_end, source)

    case (1) ! HED_y
      call Write_Message (log_unit, 'HED in y-direction')
      call HED_y (E, direction, nd, ed2nd, &
                  sx_start, sy_start, sz_start, sy_end, source)

    case (2) ! loop source
      call Write_Message (log_unit, 'Loop source')
      call loop_source (E, direction, midp_source, nd, ed2nd, &
                        sx_start, sy_start, sz_start, sx_end, source)

    case (3) ! arbitrary HED_x
      call Write_Message (log_unit, 'Arbitrary HED in x-direction')
      call arbitrary_HED_x (E, direction, nodemarker, nd, ed2nd, source)

    case (4) ! arbitrary HED_y
      call Write_Message (log_unit, 'Arbitrary HED in y-direction')
      call arbitrary_HED_y (E, direction, nodemarker, nd, ed2nd, source)

    case (5) ! straight source segment in any direction
      call Write_Message (log_unit, 'Straight source segment')
      call straight_source_segment (E, direction, nd, ed2nd, &
           sx_start, sy_start, sz_start, sx_end, sy_end, sz_end, source)

    case (6) ! segmented line source
      call Write_Message (log_unit, 'Segmented line source')
      call segmented_source(E, direction, nd, ed2nd, CSTYPE, source)

    case (7) ! segmented loop source
      call Write_Message (log_unit, 'Segmented loop source')
      call segmented_source(E, direction, nd, ed2nd, CSTYPE, source)

  end select


  ! Detect number of source edges
  num_source_edges = nint(sum(abs(source)))
  call Write_Message (log_unit, 'Number of source edges:')
  print *,num_source_edges

  call Write_Message (log_unit, 'Source edges and nodes:')
  l = 0
  do i = 1,E
    if (source(i) .ne. 0)  then
     l = l+1
     print *,i
     print *, ed2nd(i,:)
    end if
  end do

  call Write_Message (log_unit, '*************************************')
  ! Detect surface edges for imposing BC later on on these edges
  call Write_Message (log_unit, 'Detecting surface edges')
  call detect_surface_edges(M, el2ed, el2neigh, x_start, x_end, &
                            y_start, y_end, z_start, z_end, &
                            s_edges, num_s_edges)

  ! Check if PEC borehole is present
  call Write_Message (log_unit, 'Check if PEC present')
  call define_PEC (PEC, num_PEC, PEC_start, PEC_end)
  if (PEC == 1) then
    call detect_PEC_edges(E, nd, ed2nd, num_PEC, &
                          PEC_start(:,1), PEC_start(:,2), &
                          PEC_start(:,3), &
                          PEC_end(:,1), PEC_end(:,2), PEC_end(:,3), &
                          PEC_edges, num_PEC_edges)
  end if

  !> \note New in elfe3D_GPR: check if PML is present and detect whether it lies on faces, edges, or vertices
  call Write_Message (log_unit, 'Check if PML present')
  call define_PML (PML_present, PML_thickness, PML_decay_type)
  
  if (PML_present == 1) then
    call detect_PML_element_type(M, nd, el2nd, PML_thickness, num_PML_elements, &
    PML_elements, PML_type, num_non_PML_elements, non_PML_elements)
  else 
    call Write_Error_Message(log_unit, 'Cannot proceed GPR problem without the PML defined.')
  end if

  deallocate(x,y,z)
  
  ! Reading frequencies from mod_define_model
  call define_freq(Nfreq, freq)

  call Write_Message (log_unit, '*************************************')
  ! Reading receiver coordinates from mod_define_model
  call Write_Message (log_unit, 'Checking elements for all receivers')
  call define_rec(M, num_rec, a, b, c, d, Ve, u1, v1, w1, rec1_el)

  ! Allocating global matrix Agcoo in coordinate format
  allocate (Agcoo(M*36), Agrow(M*36), Agcol(M*36), stat = allo_stat)
  call allocheck(log_unit, allo_stat, "Error allocating array Agcoo") 

  ! Allocate output arrays for electric an magnetic fields
  allocate (EFields(Nfreq,num_rec,3), HFields(Nfreq,num_rec,3), &
            stat = allo_stat)
  call allocheck(log_unit, allo_stat, &
                 "Error allocating arrays EFields and HFields") 
  ! initialise
  Efields = ZEROW
  Hfields = ZEROW
  epsilon_r_eff = ZEROW
  k_bg = ZEROW

  ! --------------------------------------------------------------------
  ! Loop over frequencies:
  ! --------------------------------------------------------------------
  frequency_loop: do numfreq = 1,Nfreq

    ! defining all required frequency dependent electromagnetic parameters
    !> \note New in elfe3D_GPR: define all frequency-dependent electromagnetic parameters except w, which was previously defined
     w = 2.0_dp*pi*freq(numfreq)
     k_0 = w*sqrt(mu_0*epsilon_0)
     epsilon_r_eff = cmplx(eps/epsilon_0, -(1.0_dp/(w*rho*epsilon_0)), kind=dp)
     mu_r = mu/mu_0
     k_bg = k_0 * sqrt(mu_r*epsilon_r_eff)

     call Write_Message (log_unit, '**********************************')
     print *, 'Loop for frequency',freq(numfreq),'[Hz]'
     call Write_Message (log_unit, '**********************************')

     !> \note New in elfe3D_GPR: evaluate the stretching effect in the PML at the current frequency
     if (PML_present == 1) then

      call Write_Message (log_unit, 'Evaluating PML complex stretching factors')

      call tetrahedra_centroids(nd(el2nd(:,1),:), nd(el2nd(:,2),:), &
                                    nd(el2nd(:,3),:), nd(el2nd(:,4),:), &
                                    M, elem_centr)
      
     ! Assign PML stretching based on different decay types
     select case (PML_decay_type)
      case (1)
        print *, 'Simple Reciprocal PML decay type'
        call evaluate_anisotropic_PML_simple(num_PML_elements, elem_centr, PML_thickness, &
          k_bg, PML_elements, PML_type, LAMBDA, LAMBDA_inv)
      case (2)
        print *, 'Logarithmic PML decay type'
        call evaluate_anisotropic_PML_log(num_PML_elements, elem_centr, PML_thickness, &
          k_bg, PML_elements, PML_type, LAMBDA, LAMBDA_inv)
      case (3)
        print *, 'Exponential PML decay type'
        call evaluate_anisotropic_PML_exp(num_PML_elements, elem_centr, PML_thickness, &
          k_bg, PML_elements, PML_type, LAMBDA, LAMBDA_inv)
      case (4)
        print *, 'Polynomial PML decay type'
        call evaluate_anisotropic_PML_poly(num_PML_elements, elem_centr, PML_thickness, &
          k_bg, PML_elements, PML_type, LAMBDA, LAMBDA_inv)
      case default
        call Write_Error_Message(log_unit, 'PML decay type not recognised. Please check your input parameters.')
      end select
      
     end if

     ! -----------------------------------------------------------------
     ! V: assemble global system matrix A in sparse COO format
     ! -----------------------------------------------------------------
     ! Local matrix allocation:
     allocate (AA(6,6), K(6,6), MM(6,6), stat = allo_stat)
     call allocheck(log_unit, allo_stat, "error allocating array AA, K, MM")

     AA = ZEROW
     K = ZEROW
     MM = ZEROW

     ! Initialising global matrix Agcoo in coordinate format to zero
     Agrow = 0
     Agcol = 0
     Agcoo = ZEROW
     NNZ = 1

     !> \note New in elfe3D_GPR: allocate local stretching factors for PML computations
     if (PML_present == 1) then
       allocate(LAMBDA_l(3), LAMBDA_inv_l(3), stat = allo_stat)
       call allocheck(log_unit, allo_stat, "error allocating local stretching factors")
     end if

     call Write_Message (log_unit, 'Setting up system matrix')

     ! Loop for calculating all local AA matrices and assembling 
     ! them to global matrix in COO format
     matrix_element_loop: do l = 1,M

        stretched_element = 0
        !> \note New in version 1.2.0: PML element segmentation and handling
        if (PML_present == 1) then


          !> \note New in elfe3D_GPR: initialise local stretching factors to
          !> default values corresponding to no stretching
          LAMBDA_l = ZEROW
          LAMBDA_inv_l = ZEROW

          ! evaluating local matrices for PML elements
          do p_l = 1,num_PML_elements
            if (l == PML_elements(p_l)) then

              stretched_element = 1
              LAMBDA_l = LAMBDA(p_l,1:3)
              LAMBDA_inv_l = LAMBDA_inv(p_l,1:3)
            
              !> \note New in elfe3D_GPR: PML-specific local stiffness matrix
              call calc_stiffness_matrix_anisotropic_PML (l, el2edl, Ve, ed_sign, &
                b_start, b_end, c_start, c_end, d_start, d_end, &
                mu_r, LAMBDA_inv_l, K)

              !> \note New in elfe3D_GPR: PML-specific local mass matrix
              call calc_mass_matrix_anisotropic_PML (l, el2edl, Ve, ed_sign, b, c, d, k_0, epsilon_r_eff, LAMBDA_l, MM)
              exit
            end if
          end do

          if (stretched_element == 0) then
            ! Calculating stiffness matrix Kij for a single element 
            !> \note Modified in elfe3D_GPR: evaluates the integral together with model parameter scaling
            call calc_stiffness_matrix (l, el2edl, Ve, ed_sign, &
              b_start, b_end, c_start, c_end, d_start, d_end, &
              mu_r, K)

            ! Calculating mass matrix MMij for a single element
            !> \note Modified in elfe3D_GPR: evaluates the integral together with model parameter scaling
            call calc_mass_matrix (l, el2edl, Ve, ed_sign, b, c, d, k_0, epsilon_r_eff, MM)
          end if

        else
          call Write_Error_Message(log_unit, &
            'PML is not defined. Cannot proceed GPR problem without the PML defined.')
        end if
          
        ! check if matrix contains NaN elements, 
        ! based on NAN is not equal to itself
        do concurrent (i = 1:6, j = 1:6)
         if (K(i,j) .ne. K(i,j)) print*, 'NaN elements in K'
         if (MM(i,j) .ne. MM(i,j)) print*, 'NaN elements in MM'
        end do
      

        ! Calculating local lefthandside matrix AAij = Kij-MMij
        ! for a single element 
        !> \note Simplified in elfe3D_GPR: model parameters are already included in the matrices
        AA = K - MM

        ! check if matrix contains NaN elements, 
        ! based on NAN is not equal to itself
        do concurrent (i = 1:6, j = 1:6)
         if (AA(i,j) .ne. AA(i,j)) print*, 'NaN elements in AA'
        end do

        ! Loop through local matrix for assembly in coordinate format, 
        ! duplicates are taken care of later on:
        do j = 1,6
           do i = 1,6
              ! if entry of AA is nonzero, 
              ! assemble it in global matrix Agcoo
              if (l .ge. 1 .and. l .le. M .and. &
                  abs(AA(i,j) - cmplx(0.0_dp, 0.0_dp, kind=dp)) &
                  .gt. eps_dp) then

                      Agcoo(NNZ) = AA(i,j)
                      Agrow(NNZ) = el2ed(l,i)
                      Agcol(NNZ) = el2ed(l,j)

                      NNZ = NNZ+1

              end if
           end do
        end do

     end do matrix_element_loop

     ! deallocate local matrices
     deallocate(K, MM, AA, LAMBDA, LAMBDA_inv, LAMBDA_l, LAMBDA_inv_l) 

     ! Number of nonzeros in Agcoo
     NNZ = NNZ-1

     do l = 1,NNZ
      ! check if matrix contains NaN elements, 
      ! based on NAN is not equal to itself
      if(Agcoo(l) .ne. Agcoo(l)) then
        call Write_Message (log_unit, &
                           'Matrix Agcoo contains NaN elements!!!') 
      end if
     end do

     !call Write_Message (log_unit, &
     !              'Non-zero elements in system matrix in COO format:')
     !print *, NNZ

     !------------------------------------------------------------------
     ! VI: calculate global RHS for primal BVP 
     !------------------------------------------------------------------
     call Write_Message (log_unit, 'Calculating global source vector')

     allocate (Bg(E), stat = allo_stat)
     call allocheck(log_unit, allo_stat , "error allocating array Bg")
     ! initialise
     Bg = (0.0_dp,0.0_dp)

     ! improved source term
     source_length = 0.0_dp
     source_length = sum(abs(source)*edl)
     call Write_Message (log_unit, 'Source length:')
     print*, source_length

     Bg = cmplx(0.0_dp, (edl*source*k_0*Z_0*p)/(source_length),kind=dp) * (-1.0_dp)

     !------------------------------------------------------------------
     ! VIa: calculate global RHS for adjoint BVP
     !------------------------------------------------------------------
     call Write_Message (log_unit, 'Calculating RHS of dual problem')

     allocate (Cg(E), stat = allo_stat)
     call allocheck(log_unit, allo_stat , "error allocating array Cg")
     ! initialise
     Cg = (0.0_dp,0.0_dp)

     ! Loop over all elements containing receivers
     do l = 1, num_rec
      call calculate_elemental_influence_source (u1(l), v1(l), w1(l), &
                                                 rec1_el(l), el2ed, &
                                                 a_start, a_end, &
                                                 b_start, b_end, &
                                                 c_start, c_end, &
                                                 d_start, d_end, &
                                                 el2edl, ed_sign, &
                                                 Ve, k_0, p, &
                                                 midp_source, rec1_ed, &
                                                 infl_source)

      ! Assemble influence source term
      do i = 1,6
         Cg(rec1_ed(i)) = infl_source(i) * (-1.0_dp)
      end do
     end do

     !------------------------------------------------------------------
     ! Allocate RHS (E x 2) so that PARDISO can solve for 2 RHS vectors
     allocate (RHS(E,2), stat = allo_stat)
     call allocheck(log_unit, allo_stat , "error allocating array RHS")
     ! initialise
     RHS = (0.0_dp,0.0_dp)
     ! Put Bg (RHS of primal problem) & CG (RHS of dual problem) in RHS
     RHS(:,1) = Bg
     RHS(:,2) = Cg

     ! -----------------------------------------------------------------
     ! VII: apply boundary conditions
     ! -----------------------------------------------------------------
     ! apply boundary condition on identified boundary edges s_edges
     call Write_Message (log_unit, 'Applying Dirichlet boundary condition')
     call apply_DBC_fast(num_s_edges, s_edges, NNZ, Agcoo, Agcol, Agrow)

     call Write_Message (log_unit, &
            'Number of non-zero elements in Agcoo after applying BC:')
     print *, NNZ

     if (PEC==1) then
        ! apply boundary condition E = 0 on identified PEC/borehole edges 
        call Write_Message (log_unit, 'Apply BC for PEC')
        call apply_PEC_BC(num_PEC_edges, PEC_edges, NNZ, &
                          Agcoo, Agcol, Agrow)

        call Write_Message (log_unit, &
          'Number of non-zero elements in Agcoo after applying PEC BC:')
        print *, NNZ
     end if

     ! -----------------------------------------------------------------
     ! IX: conversion of system matrix to CSR format
     ! -----------------------------------------------------------------
     ! Convert system matrix into CSR format using adapted SPARSEKIT 
     ! functions
     ! Set size of Agcsr matrix to size of nonzero elements of Ag

     ! conversion from COO into CSR format
     call compcoicsr (E,NNZ,1,Agcoo,Agrow,Agcol)

     allocate (Agcsr(NNZ), jAgcsr(NNZ), iAgcsr(E+1), stat = allo_stat)
     call allocheck(log_unit, allo_stat, &
                    "error allocating array Agcsr etc.")
     jAgcsr = Agrow(1:NNZ)
     iAgcsr = Agcol(1:E+1)
     Agcsr = Agcoo(1:NNZ)

     ! Sorting
     call csort (E, .true., Agcsr, jAgcsr, iAgcsr)

     ! Summing duplicates
     call csumdup (jAgcsr, iAgcsr, Agcsr)

     NNZ_c = count(Agcsr .ne. (0.0_dp, 0.0_dp))
     !call Write_Message (log_unit, &
     !'Number of non-zero elements in Agcsr in csr, summed dublicates:')
     !print *, NNZ_c

     ! Conversion into upper triangular CSR format
     ! using modified SPARSEKIT function GETU for complex matrices
     allocate (Agcsr_o(NNZ), jAgcsr_o(NNZ), iAgcsr_o(E+1), &
               stat = allo_stat)
     call allocheck(log_unit, allo_stat, &
                                    "error allocating array Agcsr etc.")

     jAgcsr_o = 0
     iAgcsr_o = 0 
     Agcsr_o = (0.0_dp,0.0_dp) 

     call compgetu (E,NNZ,Agcsr,jAgcsr,iAgcsr,Agcsr_o,jAgcsr_o,iAgcsr_o)

     NNZ_c = count(Agcsr_o .ne. (0.0_dp, 0.0_dp))
     !call Write_Message (log_unit, &
     ! 'Non-zero elements in system matrix in upper triangular format:')
     ! print *, NNZ_c

     ! Agcsr,jAgcsr,iAgcsr is not needed anymore, 
     ! is in upper trangle format in Agcsr_o etc
     deallocate(Agcsr,jAgcsr,iAgcsr)

     ! Assign Agcsr_o to smaller array by eliminating zero elements
     allocate (Agcsr_oo(NNZ_c), jAgcsr_oo(NNZ_c), stat = allo_stat)
     call allocheck(log_unit, allo_stat, &
                    "error allocating array Agcsr_oo etc.")

     jAgcsr_oo = 0
     Agcsr_oo = (0.0_dp,0.0_dp) 

     icount = 1
     do l = 1,NNZ
       ! check if matrix contains NaN elements, 
       ! based on NAN is not equal to itself
       if (Agcsr_o(l) .ne. Agcsr_o(l)) then
        call Write_Message (log_unit, &
                              'Matrix Agcsr_o contains NaN elements!!!')
       ! otherwise check, if matrix element is zero
       else if (abs(Agcsr_o(l) - cmplx(0.0_dp, 0.0_dp, kind=dp)) .gt. &
                eps_dp) then

        Agcsr_oo(icount) = Agcsr_o(l)
        jAgcsr_oo(icount) = jAgcsr_o(l)

        icount = icount+1

      end if

     end do

     if(icount-1 .ne. NNZ_c) then 
      call Write_Message (log_unit, &
                    'Matrix Agcsr_oo has an incorrect number of NNZ!!!')
     end if

     deallocate(Agcsr_o,jAgcsr_o)


     ! -----------------------------------------------------------------
     ! X: apply direct solver to solve primal and adjoint BVP
     ! -----------------------------------------------------------------
     call Write_Message (log_unit, '**********************************')
     ! Allocate Solution (E x 2) so that the solver can put solution 
     ! of primal and dual problem in array
     allocate (Solution(E,2), stat = allo_stat)
     call allocheck(log_unit, allo_stat , &
                    "error allocating array Solution")
     ! initialise
     Solution = (0.0_dp,0.0_dp)

     select case (solver)

      case (1) ! PARDISO

       call Write_Message (log_unit, &
                        'Solving primal and dual problems with PARDISO')

       seconds_solve = omp_get_wtime ( ) ! timing solve

       call PARDISO_solving (E, Agcsr_oo, jAgcsr_oo, iAgcsr_o, RHS, &
                             Solution)

       seconds_solve = omp_get_wtime ( ) - seconds_solve;

      case (2) ! MUMPS

       call Write_Message (log_unit, &
                          'Solving primal and dual problems with MUMPS')

       ! back conversion into COO format 
       ! -> change column index array iAgcsr_o to IAgcoo_o
       ! allocate space for column indices iAgcoo_o
       allocate (iAgcoo_o(NNZ_c), stat = allo_stat)
       call allocheck(log_unit, allo_stat, &
                      "error allocating array Agcsr etc.")
       iAgcoo_o = 0

       call compcsrcoo (E,1,NNZ_c,Agcsr_oo, jAgcsr_oo, iAgcsr_o, &
                        nnz_new,Agcsr_oo,iAgcoo_o,jAgcsr_oo,ierr)
       if (ierr == 1) then
           call Write_Error_Message(log_unit, &
            "Wrong number of nonzeros in compcsrcoo.")
       end if

       ! update number of nonzeros?
       if (NNZ_c .ne. NNZ_new) then
         print*, 'Number of nonzero changed from', NNZ_c, 'to', NNZ_new
       end if


       seconds_solve = omp_get_wtime ( ) ! timing solve

       select case(maxRefSteps)
       case(0)
         !solve only primal problem
         call MUMPS_solving (Agcsr_oo, jAgcsr_oo, iAgcoo_o, &
                                             RHS(:,1), Solution(:,1), 1)
       case default
         ! solve primal and dual problem
         print *, 'Solving primal and dual problem since multiple refinement steps'
         call MUMPS_solving_multiple_RHS (Agcsr_oo, jAgcsr_oo, &
                                             iAgcoo_o, RHS, Solution, 2)
       end select

       seconds_solve = omp_get_wtime ( ) - seconds_solve;


      case default
        call Write_Error_Message(log_unit, &
            "Wrong solver type.")
     end select

     ! Put the primal and dual solution from the Solution matrix into 
     ! two vectors S and Wg for later 
     allocate (S(E), Wg(E), stat = allo_stat)
     call allocheck(log_unit, allo_stat, &
                    "error allocating solution vector S and vector Wg")
     ! initialise solution vectors for primal (S) and dual (Wg) problem
     S = Solution(:,1)
     Wg = Solution(:,2)

     ! Deallocate the RHS matrix and the Solution matrix
     if (allocated(RHS)) deallocate(RHS)
     if (allocated(Solution)) deallocate(Solution)
     call Write_Message (log_unit, '**********************************')

     ! -----------------------------------------------------------------
     ! XII: error estimation
     ! -----------------------------------------------------------------
     select case(maxRefSteps)
     case(0)
       call Write_Message (log_unit, 'No error estimation!')
     case default
       call Write_Message (log_unit, 'Calculating error estimates')

        !$OMP PARALLEL SECTIONS

        !$OMP SECTION
          ! Calculate elemental residuals
          print *, "Calculating elemental residuals."
          call calculate_elemental_residuals(num_non_PML_elements, non_PML_elements, k_0, el2ed, &
                                             S, Wg, Cg, Bg, epsilon_r_eff, Ve, &
                                             el2edl, ed_sign, &
                                             a_start, a_end, b_start, &
                                             b_end, c_start, c_end, &
                                             d_start, d_end, &
                                             el2nd, nd, res_s, res_w)
    

        !$OMP SECTION
          ! Calculate face-jumps
          print *, "Evaluating face jumps."
          call calculate_face_jumps (M, w, el2ed, el2nd, el2neigh, &
                                     el2edl, ed_sign, nd, &
                                     a_start, a_end, b_start, b_end, &
                                     c_start, c_end, d_start, d_end, &
                                     epsilon_r_eff, mu, Ve, S, Wg, &
                                     fjJ_s, fjJ_w, fjH_s, fjH_w)
    
        !$OMP END PARALLEL SECTIONS  


       ! Calculate element error-estimates
       print *, "computing elemental error estimates."
       call compute_elemental_error_estimates (errorEst_method, M, &
                                               res_s, res_w, &
                                               fjJ_s, fjJ_w, &
                                               fjH_s, fjH_w, &
                                               errorEst, SumerrorEst, &
                                               eta_s, eta_w)

       ! Average relative error estimate:
       RelErr = SumerrorEst/M
       print *, "Relative Error: "
       print *, RelErr

       ! Relative error estimate at receivers:
       RecErr = 0.0_dp
       do l = 1, num_rec
        RecErr = RecErr + errorEst(rec1_el(l))
       end do
       RecErr = RecErr/num_rec
     end select

     !------------------------------------------------------------------
     ! Output refinement files if vtk = 1
     !------------------------------------------------------------------

    if (vtk == 1) then
       ! write element error estimates in .vtk files for paraview
       if (refStep > 0) then
        call write_errors_vtk (M, refStep, StringStep, StringEnding, &
                        MeshFileName, errorEst, eta_s, eta_w, &
                        Ve, fjJ_s, fjJ_w, fjH_s, fjH_w)
       end if

       ! Write Refinement-Info in file
       write(9999,*) refStep, E, RelErr, RecErr, &
                     sum(res_s)/M, sum(res_w)/M, &
                     sum(fjJ_s)/M, sum(fjJ_w)/M, &
                     sum(fjH_s)/M, sum(fjH_w)/M, &
                     sum(eta_s)/M, sum(eta_w)/M

     end if

     !------------------------------------------------------------------
     ! If either refinement steps, number of dof or error tolerance 
     ! levels are exceeded, stop refinement WHILE loop
     ! otherwise generate a new refined mesh and do the forward 
     ! modelling again with refined mesh
     !------------------------------------------------------------------
      ! termination of the loop
      if (E .ge. maxUnknowns) then
        terminationCrit = 1
        if (maxRefSteps > 0) then
         call Write_Message (log_unit, '------------------------------')
         call Write_Message (log_unit, 'Terminating refinement due to Exceeding maxUnknowns')
         print *, 'No of unknowns:', E
         print *, 'Refinement step:', refStep
         print *, 'at the current frequency of:', freq(numfreq)
        end if
      end if
     
      if (refStep > maxRefSteps) then
        terminationCrit = 1
        if (maxRefSteps > 0) then
         call Write_Message (log_unit, '-------------------------------')
         call Write_Message (log_unit, 'Terminating refinement due to Exceeding maxRefSteps')
         print *, 'No of unknowns:', E
         print *, 'Refinement step:', refStep
         print *, 'at the current frequency of:', freq(numfreq)
        end if
      end if
     
      if (RelErr .le. accuracyTol) then
       terminationCrit = 1
       if (maxRefSteps > 0) then
        call Write_Message (log_unit, '--------------------------------')
        call Write_Message (log_unit, 'Terminating refinement due to Relative Error being less than the required tolerance')
        print *, 'No of unknowns:', E
        print *, 'Refinement step:', refStep
        print *, 'at the current frequency of:', freq(numfreq)
       end if
      endif

      if (E == E_before) then
       terminationCrit = 1
       if (maxRefSteps > 0) then
        call Write_Message (log_unit, '-------------------------------')
        call Write_Message (log_unit, 'Terminating refinement due to previous step elements being same as now.')
        print *, 'No of unknowns:', E
        print *, 'refinement step:', refStep
        print *, 'at the current frequency of:', freq(numfreq)
       end if
      endif

     ! Save number of unknowns of current refinement step in E_before
     E_before = E
     
     ! IF the refinement is not stopped now, generate new mesh
     if (terminationCrit == 0) then
        call calculate_elemental_volume_constraints (M, refStep, Ve, &
                                                     betaRef, errorEst,&
                                                     MeshFileName, &
                                                     NewVolumeFile)

        call generate_new_mesh(NodeFile, ElementFile, refStep, &
                               refStrategy, maxRefSteps)

 
     ! -----------------------------------------------------------------
     ! XII: Calculate E- and H-Field response at receiver locations 
     ! ----------------------------------------------------------------- 
     ! AT FINAL REFINEMENT STEP
     else if (terminationCrit == 1) then
     
       receiver_loop: do l = 1, num_rec

        call calculate_fields (u1(l), v1(l), w1(l), rec1_el(l), &
                               el2ed, S, &
                               a_start, a_end, b_start, b_end, &
                               c_start, c_end, d_start, d_end,&
                               el2edl, ed_sign, Ve, w, mu, E_rec1, &
                               H_rec1)

        ! Write results in array
        Efields(numfreq,l,1) = E_rec1(1)
        Efields(numfreq,l,2) = E_rec1(2)
        Efields(numfreq,l,3) = E_rec1(3)

      
        ! Write results in array
        Hfields(numfreq,l,1) = H_rec1(1)
        Hfields(numfreq,l,2) = H_rec1(2)
        Hfields(numfreq,l,3) = H_rec1(3)

       end do receiver_loop

       !> \note New in elfe3D_GPR: obtain fields for the first frequency at all element centroids if VTK output is enabled
       if (fields_vtk == 1 .and. numfreq == 1) then
          call Write_Message(log_unit, '***********************************************************')
          call Write_Message(log_unit, 'Writing Volumetric Field Distribution in Entire Domain Now')

        ! calculate centroids of all mesh elements
        if (.not. allocated(elem_centr)) then
          call tetrahedra_centroids(nd(el2nd(:,1),:), nd(el2nd(:,2),:), &
                                    nd(el2nd(:,3),:), nd(el2nd(:,4),:), &
                                    M, elem_centr)
        end if
         ! allocate arrays for domain fields
         allocate (domain_Efields(M,3), domain_Hfields(M,3), &
                                                       stat = allo_stat)
         call allocheck(log_unit, allo_stat, &
                                       "error allocating domain fields")

         ! initialise
         domain_Efields = ZEROW
         domain_Hfields = ZEROW
         ! calculate fields at all element centriods
         fields_element_loop: do l = 1,M

          call calculate_fields (elem_centr(l,1), elem_centr(l,2), &
                                 elem_centr(l,3), l, &
                                 el2ed, S, &
                                 a_start, a_end, b_start, b_end, &
                                 c_start, c_end, d_start, d_end,&
                                 el2edl, ed_sign, Ve, w, mu, E_rec1, &
                                 H_rec1)
          ! assign element fields to domain field array
          domain_Efields(l,:) = E_rec1(:)
          domain_Hfields(l,:) = H_rec1(:)

         end do fields_element_loop

         ! write field components into .vtk file
         call write_vtk_fields (M, refStep, MeshFileName, &
                                domain_Efields, domain_Hfields)
         ! deallocate domain field arrays
         deallocate(domain_Efields, domain_Hfields)
         ! deallocate element centriod array
         if (allocated(elem_centr)) deallocate(elem_centr)

       end if
     
     end if !(end IF last refinement step)

     ! -----------------------------------------------------------------
     ! Cleaning in frequency loop variables
     ! -----------------------------------------------------------------
     deallocate(Bg)
     deallocate(S)
     deallocate(iAgcsr_o)
     if (allocated(Agcsr_o)) deallocate(Agcsr_o,jAgcsr_o)
     if (allocated(Agcsr_oo)) deallocate(Agcsr_oo,jAgcsr_oo)
     if (allocated(Cg)) deallocate(Cg)
     if (allocated(Wg)) deallocate(Wg)
     if (allocated(iAgcoo_o)) deallocate(iAgcoo_o)
     if (allocated(elem_centr)) deallocate(elem_centr)

  ! --------------------------------------------------------------------
  ! end loop over frequencies
  ! --------------------------------------------------------------------
  end do frequency_loop


  ! Write output files 
  ! receiver_loop
  do l = 1, num_rec
    ! frequency_loop
    do numfreq = 1, Nfreq
     ! E-fields
     write ((50+1),'(10(es15.8,2x))') u1(l), v1(l), w1(l), freq(numfreq), Efields(numfreq,l,:)
     ! H-fields
     write ((50+2),'(10(es15.8,2x))') u1(l), v1(l), w1(l), freq(numfreq), Hfields(numfreq,l,:)
    end do
  end do
  ! Write output files along receiver line
  ! frequency_loop
  do numfreq = 1, Nfreq
    ! receiver_loop
    do l = 1, num_rec
     ! E-fields
     write ((50+3),'(10(es15.8,2x))') u1(l), v1(l), w1(l), freq(numfreq), Efields(numfreq,l,:)
     ! H-fields
     write ((50+4),'(10(es15.8,2x))') u1(l), v1(l), w1(l), freq(numfreq), Hfields(numfreq,l,:)
    end do
  end do


     
  ! Clean allocated refinement loop variables
  if (allocated(res_w)) deallocate(res_w, res_s)
  if (allocated(fjJ_s)) deallocate(fjJ_s, fjJ_w, fjH_s, fjH_w)
  if (allocated(errorEst)) deallocate(errorEst)
  if (allocated(eta_s)) deallocate(eta_s)
  if (allocated(eta_w)) deallocate(eta_w)
     

  ! --------------------------------------------------------------------
  ! Cleaning
  ! --------------------------------------------------------------------
  deallocate(nd)
  deallocate(ed2nd)
  deallocate(el2nd, el2neigh)
  deallocate(el2edl,el2ed,ve)
  if(allocated(edl)) deallocate(edl)
  if (allocated(ed_sign)) deallocate(ed_sign)
  deallocate(a,b,c,d)
  deallocate(a_start,a_end,b_start,b_end,c_start,c_end,d_start,d_end)
  deallocate(rho,mu,eps)
  deallocate(epsilon_r_eff, k_bg)
  deallocate(PML_elements, PML_type)
  if (allocated(Is)) deallocate(Is)
  deallocate(x_start,x_end,y_start,y_end,z_start,z_end)
  deallocate(s_edges, source)
  if (allocated(PEC_edges)) deallocate(PEC_start, PEC_end, PEC_edges)
  deallocate(u1, v1, w1, rec1_el)
  deallocate(Agcoo, Agrow, Agcol)
  deallocate(freq)
  deallocate(nodemarker,edgemarker,eleattr)
  if (allocated(EFields)) deallocate(EFields, HFields)

  call Write_Message (log_unit, 'Allocated variables were deallocated')

  ! --------------------------------------------------------------------
  ! end while-loop for refinement
  ! --------------------------------------------------------------------
  end do refinement_loop


  !---------------------------------------------------------------------
  ! For MUMPS solving, MPI finalize
  ! Check if mpi is initialised
  select case (solver)
  case(2)
    flag = .false.
    call mpi_initialized(flag, ierr)
    if (flag .eqv. .true.) then
      ! Initialise mpi
      call MPI_FINALIZE(ierr)
    end if
  end select

  !---------------------------------------------------------------------
  ! Closing output files
   close (unit = (50+1)) ! electric fields
   close (unit = (50+2)) ! magnetic fields
   close (unit = (50+3)) ! electric fields receiver line
   close (unit = (50+4)) ! magnetic fields receiver line

  if (vtk == 1) then
    close(unit=9999) ! Refinement-info File
  end if

  !---------------------------------------------------------------------
  call Write_Message (log_unit, '*************************************')
  call Write_Message (log_unit, '*************************************')
  call Write_Message (log_unit, 'Modelling with elfe3D completed')
  call Write_Message (log_unit, 'Output files are located in: /out')
  call Write_Message (log_unit, '*************************************')
  call Write_Message (log_unit, '*************************************')

  call cpu_time(finish)  ! CPU time measurement stop
  print '("CPU time to run elfe3D = ",f8.3," seconds.")',finish-start
  seconds = omp_get_wtime ( ) - seconds;
  print '("Real time to run elfe3D = ",f9.3," seconds.")',seconds

end program elfe3d_gpr
