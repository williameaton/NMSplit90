
program CCwhole
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
type(Mode)             :: mode_1 

character(len=60)     :: out_name
character(len=2)      :: nstr1 ,lstr1, nstr2 ,lstr2
character(len=1)      :: Tval
integer :: nstart, nstop, lfinish, lval, nval , lstart
logical :: toroidal 
character(len=:), allocatable :: arg

call mineos%process_mineos_model(.false.)
mineos_ptr => mineos


Tval = 'S'


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
    stop 
endif 

write(*,*)'iproc is ', iproc

region = 0 


n1 = 2
t1 = "S"
l1 = 4

n2 = 3
t2 = "T"
l2 = 2



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




allocate(Wmat(2*l1+1, 2*l2+1))

Wmat = SPLINE_iZERO


! W matrix 
call alt_W80_compute_W_matrix(sm, sm%interp, n1, t1, l1, & 
                                     n2, t2, l2, & 
                                     .true.)


! Constants so multiply after
Wmat = Wmat * OMEGA 

call buffer_int(nstr1, n1)
call buffer_int(lstr1, l1)
call buffer_int(nstr2, n2)
call buffer_int(lstr2, l2)

write(out_name, '(a)')'rot_mat/CrossCoupled/CCwhole/CC_Wmat_'//trim(nstr1)// trim(t1)// trim(lstr1)// "_"//trim(nstr2)// trim(t2)// trim(lstr2)// '_proc'//trim(arg)//'.txt'
call save_W_matrix(l1, l2, trim(out_name))


! write(out_name, '(a,i1,a,i1,a,i1,a,i1,a)')'rot_mat/VcenProcs/Vcen_', n1, t1, l1, '_', n2, t2, l2,'_proc'//trim(arg)//'.txt'
! call save_Vcen_matrix(l1, l2, trim(out_name))

deallocate(Wmat)
! deallocate(Vcen)



end program CCwhole