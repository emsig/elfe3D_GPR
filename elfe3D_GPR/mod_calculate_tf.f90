!> \file mod_calculate_tf.f90
!> \brief Module of elfe3D containing subroutines to calculate field components
!> \details Provides routines to compute electric and magnetic fields at receiver
!> \details locations from finite-element solution vectors and to write those
!> \details fields to VTK output for visualization.
!> \author Paula Rulff
!!
!> written by Paula Rulff, 24/06/2019
!!
!> Copyright (C) Paula Rulff 2025
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
module calculate_tf

  use mod_util
  use define_model

  implicit none

contains

  !---------------------------------------------------------------------
  !> \brief Subroutine for calculating electric and magnetic fields at a certain receiver location
  !> \param[in] u1 Cartesian x-coordinate of current receiver
  !> \param[in] v1 Cartesian y-coordinate of current receiver
  !> \param[in] w1 Cartesian z-coordinate of current receiver
  !> \param[in] rec1_el Current receiver element index
  !> \param[in] el2ed Element-to-edge mapping array
  !> \param[in] S Finite-element solution vector
  !> \param[in] a_start Start coefficients for shape functions
  !> \param[in] a_end End coefficients for shape functions
  !> \param[in] b_start Start coefficients for shape functions
  !> \param[in] b_end End coefficients for shape functions
  !> \param[in] c_start Start coefficients for shape functions
  !> \param[in] c_end End coefficients for shape functions
  !> \param[in] d_start Start coefficients for shape functions
  !> \param[in] d_end End coefficients for shape functions
  !> \param[in] el2edl Edge length values for local edges
  !> \param[in] ed_sign Edge orientation sign matrix
  !> \param[in] Ve Element volumes
  !> \param[in] w Angular frequency
  !> \param[in] mu Magnetic permeability per element
  !> \param[out] E_rec1 Electric field vector at the receiver
  !> \param[out] H_rec1 Magnetic field vector at the receiver
  !---------------------------------------------------------------------
  subroutine calculate_fields (u1, v1, w1, rec1_el, el2ed, S, &
                               a_start, a_end, b_start, b_end, &
                               c_start, c_end, d_start, d_end, &
                               el2edl, ed_sign, Ve, w, mu, &
                               E_rec1, H_rec1)
  
    real(kind=dp), intent(in) :: u1,v1,w1 
    integer, intent(in) :: rec1_el ! current receiver element
    integer, dimension(:,:), intent(in) :: el2ed
    complex (kind=dp), dimension(:), intent(in) :: S ! solution vector
    real(kind=dp), dimension(:,:), intent(in) :: a_start, a_end, &
                                                 b_start, b_end, &
                                                 c_start, c_end, &
                                                 d_start, d_end
    real(kind=dp), dimension(:,:), intent(in) :: el2edl
    real(kind=dp), dimension(:,:), intent(in) :: ed_sign
    real(kind=dp), dimension(:), intent(in) :: Ve
    real(kind=dp), intent(in) :: w
    real(kind=dp), dimension(:), intent(in) :: mu

    complex(kind=dp), dimension(3), intent(out) :: E_rec1, H_rec1
    

    ! LOCAL variables
    !> \brief Loop index
    integer :: l
    !> \brief Edges of the element containing receiver
    integer, dimension(6) :: rec1_ed 
    !> \brief Gradient of Lstart
    real(kind=dp), dimension(3) :: grad_Lstart ! grad Lstart
    !> \brief Gradient of Lend
    real(kind=dp), dimension(3) :: grad_Lend ! grad Lend
    !> \brief Nedelec basis functions of edges containing receiver (vector)
    real(kind=dp), dimension(3) :: N_rec1 
    !> \brief Factor -1/(iwmu)
    complex(kind=dp) :: factor   ! factor -1/(iwmu)
    !> \brief Curl of Nedelec basis function for one edge
    real(kind=dp), dimension(3) :: curl_N !
    !-------------------------------------------------------------------
    
     ! Find edges of receiver earth element
     rec1_ed = el2ed(rec1_el,:)

     ! Calculate electric fields:
     ! call Write_Message (log_unit, &
     !                    'Edge numbers of receiver earth element are:')
     ! do i = 1,6
     !    write(*,*) rec1_ed(i)
     ! end do

     E_rec1 = (/(0.0_dp,0.0_dp),(0.0_dp,0.0_dp),(0.0_dp,0.0_dp)/)
     N_rec1 = (/(0.0_dp),(0.0_dp),(0.0_dp)/)

     do l = 1,6

        ! Calculate grad Lstart and grad Lend vectors
        grad_Lstart = (/ b_start(rec1_el,l), &
                         c_start(rec1_el,l), &
                         d_start(rec1_el,l) /)

        grad_Lend = (/ b_end(rec1_el,l), &
                       c_end(rec1_el,l), &
                       d_end(rec1_el,l) /)

        ! Calculate Nedelec basis function for edge l
        N_rec1 = (((1.0_dp/(6.0_dp*Ve(rec1_el)))**2.0_dp)* &
                  (a_start(rec1_el,l) + b_start(rec1_el,l)*u1 &
                                      + c_start(rec1_el,l)*v1 &
                                      + d_start(rec1_el,l)*w1) &
                  *grad_Lend &
             
             - &
             
             ((1.0_dp/(6.0_dp*Ve(rec1_el)))**2.0_dp)* &
              (a_end(rec1_el,l) + b_end(rec1_el,l)*u1 &
                                + c_end(rec1_el,l)*v1 &
                                + d_end(rec1_el,l)*w1)*grad_Lstart) &
             
             * el2edl(rec1_el,l) &

             * ed_sign(rec1_el,l)


        ! Calculate E_rec1 contatining Ex, Ey and Ez 
        ! as the sum over all N times S
        ! of all edges of the element contatining the receiver
        E_rec1 =  E_rec1 + cmplx(N_rec1,kind=dp)*S(rec1_ed(l))

     end do

     ! Calculate magnetic fields:
     H_rec1 = (/(0.0_dp,0.0_dp),(0.0_dp,0.0_dp),(0.0_dp,0.0_dp)/)
     curl_N = (/(0.0_dp),(0.0_dp),(0.0_dp)/)
     factor = cmplx(0.0,(1/(w*mu(rec1_el))), kind=dp)

     do l = 1,6

        ! Calculate curl_N for every edge of receiver element
        curl_N = ((2.0_dp*el2edl(rec1_el,l))/ &
                  (6.0_dp*Ve(rec1_el))**2.0_dp)* &
                 (/(c_start(rec1_el,l)*d_end(rec1_el,l) &
                  - d_start(rec1_el,l)*c_end(rec1_el,l)), &
                   (d_start(rec1_el,l)*b_end(rec1_el,l) &
                  - b_start(rec1_el,l)*d_end(rec1_el,l)),&
                   (b_start(rec1_el,l)*c_end(rec1_el,l) &
                  - c_start(rec1_el,l)*b_end(rec1_el,l))/)&

                  * ed_sign(rec1_el,l)

        ! Calculate H_rec1 = -1/(iwmu)sum(Ej curl (Nj) and 
        ! add for all 6 edges
        H_rec1 =  H_rec1 + cmplx(S(rec1_ed(l))*curl_N, kind=dp)

     end do

     H_rec1 =  H_rec1*factor

   !--------------------------------------------------------------------
   end subroutine calculate_fields

   !---------------------------------------------------------------------
   !> \brief Subroutine for writing domain field components to VTK files
   !> \details Converts computed electric and magnetic field arrays into VTK
   !> \details output for visualization with ParaView. The output file name is
   !> \details generated from the mesh file name and the refinement step.
   !> \param[in] M Number of mesh points to write
   !> \param[in] refStep Refinement step index used for output file name
   !> \param[in] MeshFileName Base VTK output file name prefix
   !> \param[in] domain_Efields Complex electric field components at mesh points
   !> \param[in] domain_Hfields Complex magnetic field components at mesh points
   !---------------------------------------------------------------------
   subroutine write_vtk_fields (M, refStep, MeshFileName, &
                                domain_Efields, domain_Hfields)
   
     integer, intent(in) :: M, refStep
     complex(kind=dp), dimension(:,:), intent(in) :: domain_Efields, &
                                                    domain_Hfields


     ! LOCAL variables
     !> \brief Loop indices
     integer :: i, mi
     !> \brief File opening status
     integer :: opening
     !> \brief File length
     integer :: length
     !> \brief I/O status
     integer :: io
     !> \brief VTK file name
     character(len = 255) :: vtkFile
     !> \brief String for step
     character(len = 50) :: StringStep
     !> \brief String for ending
     character(len = 50) :: StringEnding
     !> \brief Mesh file name (local)
     character(len = 50) :: MeshFileName

    !--------------------------------------------------------------------
     ! initialise length to zero
     length = 0
     ! define name of vtk file
     write(StringStep , *) refStep
     StringEnding = ".vtk"
     vtkFile = trim(adjustl(MeshFileName))// &
               trim(adjustl(StringStep))//trim(adjustl(StringEnding))
     print*, 'Domain fields are written into ', vtkFile

     ! detect file length:
     ! open file
     open (999, file = trim(vtkFile), status='old', iostat = opening)
     ! was opening successful?
     if (opening /= 0) then
       call Write_Error_Message(log_unit, &
               ' file ' // trim(vtkFile) // ' could not be opened')
     else
         length = 0 !40829
         ! loop all lines
         do
             read(999, *, iostat = io)
             ! if end of file, exit loop
             if (io/=0) exit
             ! increment line counter
             length = length + 1
         end do
         close(unit=999)
     end if

     ! Write field components in file
     ! open file
     open (999, file = trim(vtkFile), status='old', iostat = opening)
     ! was opening successful?
     if (opening /= 0) then
        call Write_Error_Message(log_unit, &
               ' file ' // trim(vtkFile) // ' could not be opened')
     else

        do i = 1,length
            read(999,*)
        end do

        ! write new lines at the end of vtk file

        !-----------------------------Ex--------------------------------
        write(999,*)'SCALARS Real_Ex double 1'
        write(999,*)'LOOKUP_TABLE default'

        ! write real Ex component
        do mi = 1,M
            write(999,*) real(domain_Efields(mi,1))
        end do

        write(999,*)' '
        write(999,*)'SCALARS Imag_Ex double 1'
        write(999,*)'LOOKUP_TABLE default'

        ! write Imag Ex component
        do mi = 1,M
            write(999,*) aimag(domain_Efields(mi,1))
        end do

        write(999,*)' '
        write(999,*)'SCALARS Amplitude_Ex double 1'
        write(999,*)'LOOKUP_TABLE default'

        ! write Amplitude Ex
        do mi = 1,M
            write(999,*) abs(domain_Efields(mi,1))
        end do

        write(999,*)' '
        write(999,*)'SCALARS Phase_Ex double 1'
        write(999,*)'LOOKUP_TABLE default'

        ! write Phase Ex in degree
        do mi = 1,M
            write(999,*) datan2d(aimag(domain_Efields(mi,1)), &
                                  real(domain_Efields(mi,1)))
        end do


        !-----------------------------Ey--------------------------------
        write(999,*)' '
        write(999,*)'SCALARS Real_Ey double 1'
        write(999,*)'LOOKUP_TABLE default'

        ! write Real Ey
        do mi = 1,M
            write(999,*) real(domain_Efields(mi,2))
        end do

        write(999,*)' '
        write(999,*)'SCALARS Imag_Ey double 1'
        write(999,*)'LOOKUP_TABLE default'

        ! write Imag Ey
        do mi = 1,M
            write(999,*) aimag(domain_Efields(mi,2))
        end do

        write(999,*)' '
        write(999,*)'SCALARS Amplitude_Ey  double 1'
        write(999,*)'LOOKUP_TABLE default'

        ! write Amplitude Ey
        do mi = 1,M
            write(999,*) abs(domain_Efields(mi,2))
        end do

        write(999,*)' '
        write(999,*)'SCALARS Phase_Ey double 1'
        write(999,*)'LOOKUP_TABLE default'

        ! write Phase Ey in degrees
        do mi = 1,M
            write(999,*) atan2(aimag(domain_Efields(mi,2)), &
                                  real(domain_Efields(mi,2))) * 180.0_dp / pi ! rescale to degrees
        end do

        !-----------------------------Ez--------------------------------
        write(999,*)' '
        write(999,*)'SCALARS Real_Ez double 1'
        write(999,*)'LOOKUP_TABLE default'
        ! write Real Ez
        do mi = 1,M
            write(999,*) real(domain_Efields(mi,3))
        end do

        write(999,*)' '
        write(999,*)'SCALARS Imag_Ez double 1'
        write(999,*)'LOOKUP_TABLE default'

        ! write Imag Ez
        do mi = 1,M
            write(999,*) aimag(domain_Efields(mi,3))
        end do

        write(999,*)' '
        write(999,*)'SCALARS Amplitude_Ez double 1'
        write(999,*)'LOOKUP_TABLE default'

        ! write Amplitude Ez
        do mi = 1,M
            write(999,*) abs(domain_Efields(mi,3))
        end do

        write(999,*)' '
        write(999,*)'SCALARS Phase_Ez double 1'
        write(999,*)'LOOKUP_TABLE default'

        ! write Phase Ez in degrees
        do mi = 1,M
            write(999,*) atan2(aimag(domain_Efields(mi,3)), &
                                  real(domain_Efields(mi,3))) * 180.0_dp / pi ! rescale to degrees
        end do

        !-----------------------------Hx--------------------------------
        write(999,*)'SCALARS Real_Hx double 1'
        write(999,*)'LOOKUP_TABLE default'

        ! write real Hx component
        do mi = 1,M
            write(999,*) real(domain_Hfields(mi,1))
        end do

        write(999,*)' '
        write(999,*)'SCALARS Imag_Hx double 1'
        write(999,*)'LOOKUP_TABLE default'

        ! write Imag Hx component
        do mi = 1,M
            write(999,*) aimag(domain_Hfields(mi,1))
        end do

        write(999,*)' '
        write(999,*)'SCALARS Amplitude_Hx double 1'
        write(999,*)'LOOKUP_TABLE default'

        ! write Amplitude Ex
        do mi = 1,M
            write(999,*) abs(domain_Hfields(mi,1))
        end do

        write(999,*)' '
        write(999,*)'SCALARS Phase_Hx double 1'
        write(999,*)'LOOKUP_TABLE default'

        ! write Phase Hx in degree
        do mi = 1,M
            write(999,*) atan2(aimag(domain_Hfields(mi,1)), &
                                  real(domain_Hfields(mi,1))) * 180.0_dp / pi ! rescale to degrees
        end do


        !-----------------------------Ey--------------------------------
        write(999,*)' '
        write(999,*)'SCALARS Real_Hy double 1'
        write(999,*)'LOOKUP_TABLE default'

        ! write Real Hy
        do mi = 1,M
            write(999,*) real(domain_Hfields(mi,2))
        end do

        write(999,*)' '
        write(999,*)'SCALARS Imag_Hy double 1'
        write(999,*)'LOOKUP_TABLE default'

        ! write Imag Hy
        do mi = 1,M
            write(999,*) aimag(domain_Hfields(mi,2))
        end do

        write(999,*)' '
        write(999,*)'SCALARS Amplitude_Hy  double 1'
        write(999,*)'LOOKUP_TABLE default'

        ! write Amplitude Hy
        do mi = 1,M
            write(999,*) abs(domain_Hfields(mi,2))
        end do

        write(999,*)' '
        write(999,*)'SCALARS Phase_Hy double 1'
        write(999,*)'LOOKUP_TABLE default'

        ! write Phase Hy in degrees
        do mi = 1,M
            write(999,*) atan2(aimag(domain_Hfields(mi,2)), &
                                  real(domain_Hfields(mi,2))) * 180.0_dp / pi ! rescale to degrees
        end do

        !-----------------------------Hz--------------------------------
        write(999,*)' '
        write(999,*)'SCALARS Real_Hz double 1'
        write(999,*)'LOOKUP_TABLE default'
        ! write Real Hz
        do mi = 1,M
            write(999,*) real(domain_Hfields(mi,3))
        end do

        write(999,*)' '
        write(999,*)'SCALARS Imag_Hz double 1'
        write(999,*)'LOOKUP_TABLE default'

        ! write Imag Hz
        do mi = 1,M
            write(999,*) aimag(domain_Hfields(mi,3))
        end do

        write(999,*)' '
        write(999,*)'SCALARS Amplitude_Hz double 1'
        write(999,*)'LOOKUP_TABLE default'

        ! write Amplitude Hz
        do mi = 1,M
            write(999,*) abs(domain_Hfields(mi,3))
        end do

        write(999,*)' '
        write(999,*)'SCALARS Phase_Hz double 1'
        write(999,*)'LOOKUP_TABLE default'

        ! write Phase Hz in degrees
        do mi = 1,M
            write(999,*) atan2(aimag(domain_Hfields(mi,3)), &
                                  real(domain_Hfields(mi,3))) * 180.0_dp / pi ! rescale to degrees
        end do


        close(unit=999)
     end if
   !--------------------------------------------------------------------
   end subroutine write_vtk_fields

  end module calculate_tf
