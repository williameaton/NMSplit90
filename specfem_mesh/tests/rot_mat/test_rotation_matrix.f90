
program test_W_matrix
    use params, only:  Wmat, rho_spl, nprocs, Vcen
    use allocation_module, only: allocate_if_unallocated, deallocate_if_allocated
    use mineos_model, only: mineos, mineos_ptr
    use modes, only: Mode, get_mode
    use specfem_mesh, only: create_SetMesh, SetMesh
    use piecewise_interpolation, only: InterpPiecewise, create_PieceInterp
    implicit none
    include "constants.h"

integer :: iproc, region, knot_lower, knot_upper, n1, n2, l1, l2, imode, cmdlength , num_args, i
character :: t1, t2
type(SetMesh)         :: sm
type(InterpPiecewise) :: interp
character(len=60)     :: out_name

character(len=:), allocatable :: arg


integer, dimension(31), parameter :: modeNs = (/0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 6, 8, 8, 9, 11, 11, 13, 13, 18, 18 /)
integer, dimension(31), parameter :: modeLs = (/2, 3, 4, 5, 6, 7, 8, 2, 3, 4, 5, 6, 7, 8, 9, 3, 4, 5, 6, 1, 2, 3, 1, 5, 3,  4,  5,  1,  2,  3,  4 /)

real(kind=8), allocatable :: gpsi(:,:,:,:,:)
real(kind=8) :: ggpsi(3,3)

! Read mineos model 
call mineos%process_mineos_model(.false.)
mineos_ptr => mineos



num_args = command_argument_count()
call get_command_argument(1, length=cmdlength)
allocate(character(len=cmdlength) :: arg)
call get_command_argument(1, arg)
if(cmdlength.eq.1)then 
    read(arg,'(i1)') iproc
elseif(cmdlength.eq.2)then
    read(arg,'(i2)') iproc
else 
    write(*,*)"ERROR IN LENGTH: HERE!!"
endif 

write(*,*)'iproc is ', iproc




do imode = 1, 1

    n1 = modeNs(imode)
    n2 = modeNs(imode)

    t1 = 'S'
    t2 = 'S'

    l1 = modeLs(imode)
    l2 = modeLs(imode)

    region = 0

    ! The matrix should be 2l + 1 from -m to m 
    allocate(Wmat(2*l1+1, 2*l2+1))
    allocate(Vcen(2*l1+1, 2*l2+1))
    Wmat = SPLINE_iZERO
    Vcen = SPLINE_iZERO


    !do iproc = 0, nprocs-1
    !do iproc = 1, 35
        ! Things that need to be done for each processor
        sm = create_SetMesh(iproc, region)
        call sm%read_proc_coordinates()
        call sm%load_ibool()
        call sm%setup_gll()
        call sm%compute_jacobian(.true.)
        call sm%compute_wglljac(.false.)

        call sm%compute_rtp_from_xyz(.true.)

        call sm%setup_global_coordinate_arrays(.true.)
        call sm%compute_rtp_from_xyz(.true.)
        call sm%get_unique_radii(.true.)
        call sm%compute_rotation_matrix()

        ! W matrix 
        !call compute_W_matrix(sm, sm%interp, n1, t1, l1, & 
        !                                     n2, t2, l2, & 
        !                                     .true.)

        !Vcen matrix
        call compute_Vcentrifugal(sm, sm%interp, n1, t1, l1, & 
                                                 n2, t2, l2, & 
                                                 .true.)

        !allocate(gpsi(3, sm%ngllx, sm%nglly, sm%ngllz, sm%nspec))
        !call compute_grad_centrifugal(sm, gpsi, ggpsi, myrank)
        !call compute_Vcen_matrix(sm, sm%interp, gpsi, n1, t1, l1, n2, t2, l2)


        !deallocate(gpsi)
        call sm%cleanup()

    !enddo 


    ! Constants so multiply after
    !Wmat = Wmat * OMEGA * iONE

    !write(out_name, '(a,i1,a,i1,a,i1,a,i1,a)')'rot_mat/Wmat_', n1, t1, l1, '_', n2, t2, l2, '_proc'//trim(arg)//'.txt'
    !call save_W_matrix(l1, l2, trim(out_name))

    write(out_name, '(a,i1,a,i1,a,i1,a,i1,a)')'rot_mat/VcenProcs/Vcen_', n1, t1, l1, '_', n2, t2, l2,'_proc'//trim(arg)//'.txt'

    call save_Vcen_matrix(l1, l2, trim(out_name))

    deallocate(Wmat)
    deallocate(Vcen)

    write(*,*)"Finished mode: ", imode

enddo 


end program test_W_matrix