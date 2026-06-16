!> \file mod_tetgen_operations.f90
!> \brief Module of elfe3d for mesh operations involving TetGen
!> \details Provides routines for adaptive mesh refinement, volume constraint generation,
!> \details and external TetGen invocation to regenerate the mesh during the refinement cycle.
!!
!> written by Laura Maria Buntin & Paula Rulff, 19/11/2019
!!
!> Copyright (C) Laura Maria Buntin & Paula Rulff 2020
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

module tetgen_operations
  
  use mod_util

  implicit none

contains
   !--------------------------------------------------------------------
   !> \brief Calculate volume constraints for mesh refinement
   !> \param[in] M Number of elements in the mesh
   !> \param[in] refStep Current refinement step index
   !> \param[in] Ve Element volume array
   !> \param[in] betaRef Refinement threshold multiplier
   !> \param[in] errorEst Element error estimator array
   !> \param[in] StringName Base mesh file name used to derive the volume file name
   !> \param[in,out] NewVolumeFile Output volume constraint file name
   !--------------------------------------------------------------------
   subroutine calculate_elemental_volume_constraints (M, refStep, Ve, &
                                                      betaRef, &
                                                      errorEst, &
                                                      StringName, &
                                                      NewVolumeFile)
  
     ! INPUT
     ! number of elements
     integer, intent(in) :: M
     ! current refinementstep
     integer, intent(in) :: refStep
     ! element volume
     real(kind=dp), dimension(:), intent(in) :: Ve
     ! factor beta for refinement
     real(kind=dp), intent(in) :: betaRef
     ! elemental error estimator
     real(kind=dp), dimension(:), intent(in) :: errorEst
     ! array containing receiver-element numbers
     !integer, dimension(:), intent(in) :: rec1_el
     character(len = 50), intent(in) :: StringName
     
     ! OUTPUT
     ! FileName of file containing volume constraints
     character(len = 255), intent(inout) :: NewVolumeFile

     ! LOCAL variables
     integer :: mi
     ! new volume of elements
     real(kind=dp) :: newVol
     ! refinement threshold
     real(kind=dp) :: refThres
     ! counter: how many elements refined:
     integer :: noEleRef
     ! List, which elements are refined
     character(len = 50) :: StringStep, StringEnding
     ! file operation variables
     integer :: out_unit = 1
     integer :: opening
     
     
    !-------------------------------------------------------------------
    ! calculate volume constraints
    call Write_Message (log_unit, &
                    'Calculate volume constraints for the refined mesh')
    
    ! VolFileName 
    write(StringStep , *) refStep
    StringEnding = '.vol'
    NewVolumeFile = trim(adjustl(StringName))// &
                    trim(adjustl(StringStep))// &
                    trim(adjustl(StringEnding))
    
    ! open Volume file
    open(out_unit, file = trim(NewVolumeFile), status='new', &
                   action = 'write', iostat = opening)
    
    ! was opening successful?
    if (opening /= 0) then
        call Write_Error_Message(log_unit, &
            'calculate_elemental_volume_constraints: file ' // &
            trim(NewVolumeFile) // ' could not be opened')
    else
        
        ! determine the desired threshold for refinement
        ! (beta * maximum error estimator)
        !!! PR: changed to log10 here 08/21
        refThres = betaRef * log10(maxval(errorEst)) 
        
        ! counter of refined elements
        noEleRef = 0
        
        ! Write #elements in first line of .vol file
        write(out_unit,*) M

        ! Loop all elements and write #element new Vol in each line of
        ! .vol file
        do mi = 1,M
            
            ! if elemental error estimator is above threshold: refine
            !!! PR: changed to log10 here 08/21
            if (log10(errorEst(mi) )>= refThres) then 
                ! Do not refine rec_el
                !if (any(rec1_el == mi)) then
                !    newVol = 0 ! zero = no constraint
                !    write(out_unit,*) mi, newVol
                !else
                    newVol = Ve(mi) * 0.5
                    noEleRef = noEleRef + 1
                    write(out_unit,*) mi, newVol
                !end if
            else
                newVol = 0 ! zero = no constraint
                write(out_unit,*) mi, newVol
            end if
        
        end do
        
        ! print how many elements will be refined
        print*, noEleRef, 'of', M, 'elements will be refined'
        print*, maxval(errorEst), 'max error est'
        print*, refThres, 'refinement threshold'
        
        
        ! close file
        close(out_unit)
        
    end if
    !-------------------------------------------------------------------
    end subroutine calculate_elemental_volume_constraints

   !--------------------------------------------------------------------
   !> \brief Generate a new refined mesh by invoking TetGen
   !> \param[in] NodeFile Base node file name for the current mesh
   !> \param[in] ElementFile Base element file name for the current mesh
   !> \param[in] refStep Current refinement step index
   !> \param[in] refStrategy Refinement strategy selector
   !> \param[in] maxRefSteps Maximum number of refinement steps configured
   !--------------------------------------------------------------------
   subroutine generate_new_mesh (NodeFile, ElementFile, refStep, &
                                 refStrategy, maxRefSteps)
  
     ! INPUT
     ! Current mesh (.nd and .elem) and .vol file to create new mesh
     character(len = 255), intent(in) :: NodeFile, ElementFile
     integer, intent(in) :: refStep, refStrategy, maxRefSteps

     ! LOCAL variables
     integer :: i
     character(len = 255) :: TetGen, TetGenCommand
     character(len = 1) :: space

    !-------------------------------------------------------------------
    call Write_Message (log_unit, 'Refine the mesh calling TetGen')

    call execute_command_line ('pwd', exitstat=i)

    select case (refStrategy)
        case (0) ! constant, low quality factor

            TetGen = 'tetgen -rq1.6kAaen'

        case (1) ! maxRefSteps-1 with low quality , last high quality

            if (refStep < maxRefSteps) then
               TetGen = 'tetgen -rq1.6kAaen'
             else if (refStep == maxRefSteps) then
               TetGen = 'tetgen -rq1.2kAaen'
            end if


        case (2) ! increasing quality factor

            if (refStep == 1) then
               TetGen = 'tetgen -rq1.55kAaen'
             else if (refStep == 2) then
               TetGen = 'tetgen -rq1.5kAaen'
             else if (refStep == 3) then
               TetGen = 'tetgen -rq1.45kAaen'
             else if (refStep == 4) then
               TetGen = 'tetgen -rq1.4kAaen'
             else if (refStep == 5) then
               TetGen = 'tetgen -rq1.35kAaen'
             else if (refStep == 6) then
               TetGen = 'tetgen -rq1.3kAaen'
             else if (refStep == 7) then
               TetGen = 'tetgen -rq1.25kAaen'
             else if (refStep == 8) then
               TetGen = 'tetgen -rq1.2kAaen'
             else if (refStep .ge. 9) then
               TetGen = 'tetgen -rq1.15kAaen'
            end if

        ! increasing quality factor on mesh with detailled
        ! subsurface anomaly (-T and -d  option added)
        case (3) 

            if (refStep == 1) then
               TetGen = 'tetgen -rq1.6kAaenT0.00005dV'
             else if (refStep == 2) then
               TetGen = 'tetgen -rq1.55kAaenT0.000005dV'
             else if (refStep == 3) then
               TetGen = 'tetgen -rq1.5kAaenT0.00005dV'
             else if (refStep == 4) then
               TetGen = 'tetgen -rq1.45kAaenT0.000005dV'
             else if (refStep == 5) then
               TetGen = 'tetgen -rq1.4kAaenT0.000005dV'
             else if (refStep == 6) then
               TetGen = 'tetgen -rq1.35kAaenT0.00005dV'
             else if (refStep == 7) then
               TetGen = 'tetgen -rq1.3kAaenT0.00005dV'
             else if (refStep == 8) then
               TetGen = 'tetgen -rq1.25kAaenT0.00005dV'
             else if (refStep .ge. 9) then
               TetGen = 'tetgen -rq1.2kAaenT0.00005dV'
            end if

        case(4) ! for meshes from gmesh software

            if (refStep < maxRefSteps) then
               TetGen = 'tetgen -rq1.6keAafiVn'

             else if (refStep == maxRefSteps) then
               TetGen = 'tetgen -rq1.2keAafiVn'
            end if


    end select

    ! call TetGen in /in folder to refine previous mesh
    space = ' '
    TetGenCommand =  trim(adjustl(TetGen))//space// &
                     trim(adjustl(NodeFile))//space// &
                     trim(adjustl(ElementFile))
    print*, TetGenCommand
    call execute_command_line (TetGenCommand, exitstat=i)
    print *, "Exit status of TetGen was ", i

    !-------------------------------------------------------------------
    end subroutine generate_new_mesh
    
  end module tetgen_operations
