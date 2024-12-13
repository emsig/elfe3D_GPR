!-------------------------------------------------------------------------------------------
! Program to write elemental model parameters (rho, mu_r, epsilon_r) into input .vtk file
!
! written by Paula Rulff, 23/04/2020
!
! compile with: gfortran -o write_model_para_in_vtk write_model_para_in_vtk.f90
! use with: ./write_model_para_in_vtk [name of file with region parameters] [name of .vtk file]
! e.g. ./write_model_para_in_vtk CSEM_regionpara.txt CSEM_input_model.1.vtk
!
! You have to define region attributes, when running tetgen before.
! The input [file with region parameters] must have the following structure:
	! # eleattr
	! 3
	! # eleattr rho mu_r epsilon_r
	! 1 100000000.0 1.0 0.0
	! 2 100 1.0 0.0
	! 3 10 1.0 0.0

!-------------------------------------------------------------------------------------------
program write_model_para_in_vtk
	implicit none

	! variable declaration

	! command line input
	character(len=100) :: ModParaFileName, vtkFileName

	! other variables
	character(len=10) :: string
	integer :: in_unit = 1
    integer :: opening, length, io
    integer :: i,j
    integer :: allo_stat
    integer :: skipcol

    integer :: num_attr, num_ele ! number of different element attributes, number of elements
    integer, allocatable, dimension(:) :: attr, eleattr ! existing attributes and element attributes
    real(8), allocatable, dimension(:) :: rho, mu_r, epsilon_r ! model parameters

	!-------------------------------------------------------------------------------------------
	! read command line input
	! make sure the right number of inputs have been provided
	if(COMMAND_ARGUMENT_COUNT().NE.2)THEN
 		 write(*,*)'ERROR, TWO COMMAND-LINE ARGUMENTS REQUIRED, STOPPING'
 	 stop
	endif
	! read in the two file names
	CALL GET_COMMAND_ARGUMENT(1,ModParaFileName)   
	CALL GET_COMMAND_ARGUMENT(2,vtkFileName)

	print *,'Your input filenames are:'
	print *, ModParaFileName
	print *, vtkFileName



	! open/read regionpara file and assign to arrays
	! open the file
    open (in_unit, file = trim(ModParaFileName), status='old', action = 'read', iostat = opening)

    ! was opening successful?
    if (opening /= 0) then
         write(*,*) 'file ' // trim(ModParaFileName) // ' could not be opened'
    else
    	! skip # line
    	read (in_unit, *)
    	! read number of element attributes from the second line
    	read (in_unit, *) num_attr
    	! skip # line
    	read (in_unit, *)
    	! allocate attr, rho, mu_r, epsilon_r arrays
    	allocate (attr(num_attr), rho(num_attr), mu_r(num_attr), epsilon_r(num_attr), stat = allo_stat)  
    	! read all following lines
    	do i = 1, num_attr
    		read (in_unit, *) attr(i), rho(i), mu_r(i), epsilon_r(i)	
    	end do

    	print *, 'attribute', attr
    	print *, 'rho', rho
    	print *, 'mu_r', mu_r
    	print *, 'epsilon_r', epsilon_r

    	! close modelpara file
    	close (unit = in_unit)



		! open/read .vtk file
		open (in_unit, file = trim(vtkFileName), status='old', iostat = opening)
		! was opening successful?
	    if (opening /= 0) then
	         write(*,*) 'file ' // trim(vtkFileName) // ' could not be opened'
	    else
	    	! detect file length:
	    	length = 0
	          	! loop all lines
	          	do
	            	read(in_unit, *, iostat = io)
	            	! if end of file, exit loop
	            	if (io/=0) exit
	             	! increment line counter
	            	length = length + 1
	          	end do

	        ! read number of elements from .vtk file
	        rewind(in_unit)
	        do while(string/='CELL_DATA')
	   			read(in_unit, *) string
			end do
			backspace(in_unit)
			read (in_unit, *) string, num_ele
			print *, 'number of elements', num_ele

			! skip two lines
	    	read (in_unit, *)
	    	read (in_unit, *)

	        ! allocate eleattr
	        allocate (eleattr(num_ele), stat = allo_stat)  

	        ! read elemental attrcibutes and store them in eleattr
	        do i = 1, num_ele
	    		read (in_unit, *) eleattr(i)
	    	end do

	    	! write header lines in the end of vtk file
	        write(in_unit,*)'SCALARS rho double 1'
	        write(in_unit,*)'LOOKUP_TABLE default'

	        ! write list of elemental rho in .vtk file
	        do i = 1, num_ele
	        	do j = 1, num_attr
	        		if (eleattr(i) .eq. attr(j)) then
	    				write (in_unit, *) rho(j)
	    			endif
	    		end do
	    	end do

	        ! write header lines in the end of vtk file
	        write(in_unit,*)'SCALARS mu_r double 1'
	        write(in_unit,*)'LOOKUP_TABLE default'

	        ! write list of elemental mu_r in .vtk file
	        do i = 1, num_ele
	        	do j = 1, num_attr
	        		if (eleattr(i) .eq. attr(j)) then
	    				write (in_unit, *) mu_r(j)
	    			endif
	    		end do
	    	end do
	        ! write header lines in the end of vtk file
	        write(in_unit,*)'SCALARS epsilon_r double 1'
	        write(in_unit,*)'LOOKUP_TABLE default'

	        ! write list of elemental epsilon_r in .vtk file
	        do i = 1, num_ele
	        	do j = 1, num_attr
	        		if (eleattr(i) .eq. attr(j)) then
	    				write (in_unit, *) epsilon_r(j)
	    			endif
	    		end do
	    	end do
	    	! close vtk file
	    	close (unit = in_unit)
	    	print *,'elemental model parameters written into .vtk file'
	    end if

	end if

	! deallocate allocated arrays
	if(allocated(attr)) deallocate(attr)
	if(allocated(rho)) deallocate(rho)
	if(allocated(mu_r)) deallocate(mu_r)
	if(allocated(epsilon_r)) deallocate(epsilon_r)
	if(allocated(eleattr)) deallocate(eleattr)

	

end program write_model_para_in_vtk