module specfem_mesh 
    use piecewise_interpolation, only: InterpPiecewise, create_PieceInterp
    use allocation_module, only: allocate_if_unallocated, deallocate_if_allocated
    use modes, only: Mode
    implicit none 
    include "constants.h"


    !interface map_local_global

    !end interface map_local_global

    !interface map_complex_vector

    !end interface map_complex_vector


    type :: SetMesh 
        ! Object holding mesh details for a specific set of the specfem
        ! mesh 
    
        ! Set id
        integer   :: iset
        integer   :: region 
        character :: set_str 

        ! GLL values
        integer :: nglob
        integer :: nspec
        integer :: ngllx
        integer :: nglly
        integer :: ngllz
        real(kind=CUSTOM_REAL), allocatable :: xi(:), wgll(:), dgll(:,:)
        real(kind=CUSTOM_REAL), allocatable :: wglljac(:,:,:,:)

        ! GLL coordinates 
        real(kind=CUSTOM_REAL), allocatable :: xstore(:,:,:,:),    & 
                                               ystore(:,:,:,:),    & 
                                               zstore(:,:,:,:),    &
                                               rstore(:,:,:,:),    & 
                                               thetastore(:,:,:,:),& 
                                               phistore(:,:,:,:)
        
        ! Global coordinates
        real(kind=CUSTOM_REAL), allocatable   :: x_glob(:), & 
                                                 y_glob(:), &
                                                 z_glob(:)

        ! Jacobians
        real(kind=CUSTOM_REAL), allocatable    :: jac(:,:,:,:,:,:)
        real(kind=CUSTOM_REAL), allocatable    :: detjac(:,:,:,:)
        real(kind=CUSTOM_REAL), allocatable    :: jacinv(:,:,:,:,:,:)


        ! ibool
        integer, allocatable                :: ibool(:,:,:,:)
        ! Unique mesh radii 
        ! Number of unique mesh radii 
        real(kind=CUSTOM_REAL), allocatable :: unique_r(:)
        integer                             :: n_unique_rad

        ! Map from GLL to its unique radius
        integer, allocatable                :: rad_id(:, :, :, :) 

        ! Interpolation map for the mesh 
        type(InterpPiecewise) :: interp
        real(kind=CUSTOM_REAL),    allocatable :: Rmat(:,:,:)

        ! Field variables: 
        complex(kind=SPLINE_REAL), allocatable :: disp1(:,:, :, :, :)
        complex(kind=SPLINE_REAL), allocatable :: disp2(:,:, :, :, :)
        complex(kind=SPLINE_REAL), allocatable :: gradphi_1(:,:, :, :, :)
        complex(kind=SPLINE_REAL), allocatable :: gradphi_2(:,:, :, :, :)
        complex(kind=SPLINE_REAL), allocatable :: strain1(:, :, :, :, :)
        complex(kind=SPLINE_REAL), allocatable :: strain2(:, :, :, :, :)
        complex(kind=SPLINE_REAL), allocatable :: gradS_1(:, :, :, :, :, :)
        complex(kind=SPLINE_REAL), allocatable :: gradS_2(:, :, :, :, :, :)
        complex(kind=SPLINE_REAL), allocatable :: devE_1(:, :, :, :, :, :)
        complex(kind=SPLINE_REAL), allocatable :: devE_2(:, :, :, :, :, :)
        complex(kind=SPLINE_REAL), allocatable :: gSR_1( :, :, :, :, :)
        complex(kind=SPLINE_REAL), allocatable :: gSR_2( :, :, :, :, :)

        ! Magnitude of g 
        real(kind=CUSTOM_REAL), allocatable :: gmag_at_r(:)


        contains 
             procedure :: get_unique_radii 
        
            !procedure :: map => map_local_global
            !procedure :: map_complex => map_complex_vector
            procedure map_complex_vector_4
            procedure map_complex_vector_8
            procedure map_local_global_double_precision
            procedure map_local_global_real_4
            procedure map_local_global_complex_4
            procedure map_local_global_complex_8
            
            procedure :: setup_mesh_sem_details
            procedure :: check_ibool_is_defined
            procedure :: load_ibool
            procedure :: read_proc_coordinates
            procedure :: read_integer_proc_variable
            procedure :: setup_global_coordinate_arrays
            procedure :: compute_rtp_from_xyz
            procedure :: compute_rotation_matrix
            procedure :: rotate_complex_vector_rtp_to_xyz
            procedure :: rotate_complex_matrix_rtp_to_xyz
            procedure :: rotate_complex_sym_matrix_rtp_to_xyz
            procedure :: compute_jacobian
            procedure :: recalc_jacobian_gll3D
            procedure :: save_get_mesh_radii_results
            procedure :: load_get_mesh_radii_results  
            procedure :: load_rho_spline
            procedure :: compute_mode_gradphi
            procedure :: compute_mode_grad_sr
            procedure :: compute_mode_displacement
            procedure :: compute_mode_strain
            procedure :: compute_mode_gradS
            procedure :: compute_strain_deviator
            procedure :: compute_strain_deviator_congj
            procedure :: save_mode_disp_binary
            procedure :: load_mode_disp_binary
            procedure :: save_mode_strain_binary
            procedure :: load_mode_strain_binary
            procedure :: load_jacobian
            procedure :: save_jacobian
            procedure :: setup_gll
            procedure :: compute_wglljac
            procedure :: save_wglljac
            procedure :: load_wglljac
            procedure :: save_global_xyz
            procedure :: load_global_xyz
            procedure :: save_elem_rtp
            procedure :: load_elem_rtp
            procedure :: compute_background_g

            procedure :: cleanup
    end type SetMesh

    contains


        subroutine check_ibool_is_defined(self)
            ! Checks that ibool is allocated and not just zero
            ! TODO: call ibool loading if not allocated? 
            implicit none 
            class(SetMesh) :: self
        
            if(.not. allocated(self%ibool))then 
                write(*,*)'ERROR: ibool is not allocated but is about to be used.'
                stop
            else 
                if (self%nglob.ne.maxval(self%ibool))then
                    write(*,*)'ERROR: nglob is not equal to the maximum ibool value'
                    write(*,*)'nglob    : ', self%nglob
                    write(*,*)'max ibool: ', maxval(self%ibool)
                    stop
                endif
            endif 
        end subroutine check_ibool_is_defined
        
        
        
        subroutine load_ibool(self)
            implicit none 
            class(SetMesh) :: self

            call allocate_if_unallocated(self%ngllx, self%nglly, self%ngllz, self%nspec,self%ibool)
            call self%read_integer_proc_variable(self%ibool, 'ibool')
            call self%check_ibool_is_defined()
        end subroutine
        


        
        subroutine read_proc_coordinates(self)
            use params, only: datadir, verbose, IIN
            implicit none 
            class(SetMesh) :: self

            ! Local variables
            integer :: ier
            character(len=550) :: binname
            
            double precision, allocatable :: xstore_dp(:,:,:,:)
            double precision, allocatable :: ystore_dp(:,:,:,:)
            double precision, allocatable :: zstore_dp(:,:,:,:)
        
            ! File name prefix: 
            write(binname,'(a,i0.6,a,i1,a)')trim(datadir)//'/proc',self%iset,'_'//'reg',self%region,'_'

            if(verbose.gt.1)then
                write(*,'(/,a,/)')'• Reading files from '//trim(datadir)
                write(*,'(a,i1)')'  -- region      : ', self%region
                write(*,'(a,i0.6,/)')'  -- processor id: ', self%iset
            endif
        
            ! Read processor info to get ngll and nspec
            open(unit=IIN,file=trim(binname)//'info.bin', &
            status='unknown',form='unformatted',action='read',iostat=ier)
            if (ier.ne.0)then 
                write(*,'(a, i0.6)')'Couldnt read info file for proc ', self%iset
                stop
            endif 
            read(IIN)self%nglob
            read(IIN)self%nspec
            read(IIN)self%ngllx
            read(IIN)self%nglly
            read(IIN)self%ngllz
            close(IIN)
                
            if(verbose.ge.3)then
                write(*,'(a)')'  -- Info: '
                write(*,*)'    --> nglob: ', self%nglob
                write(*,*)'    --> nspec: ', self%nspec
                write(*,'(a,i1)')'     --> ngllx: ', self%ngllx
                write(*,'(a,i1)')'     --> nglly: ', self%nglly
                write(*,'(a,i1)')'     --> ngllz: ', self%ngllz
            endif 
        
            ! Allocate mesh arrays:        
            call allocate_if_unallocated(self%ngllx, self%nglly, self%ngllz, self%nspec, xstore_dp)
            call allocate_if_unallocated(self%ngllx, self%nglly, self%ngllz, self%nspec, self%xstore)
            call allocate_if_unallocated(self%ngllx, self%nglly, self%ngllz, self%nspec, ystore_dp)
            call allocate_if_unallocated(self%ngllx, self%nglly, self%ngllz, self%nspec, self%ystore)
            call allocate_if_unallocated(self%ngllx, self%nglly, self%ngllz, self%nspec, zstore_dp)
            call allocate_if_unallocated(self%ngllx, self%nglly, self%ngllz, self%nspec, self%zstore)
        
            ! Open the x coordinate and load: 
            open(unit=IIN,file=trim(binname)//'xstore.bin', &
            status='unknown',form='unformatted',action='read',iostat=ier)
            if (ier.ne.0)then 
                write(*,'(a, i0.6)')'Couldnt read xstore file for proc ', self%iset
                stop
            endif 
            read(IIN)xstore_dp
            close(IIN)

            ! Open the y coordinate and load: 
            open(unit=IIN,file=trim(binname)//'ystore.bin', &
            status='unknown',form='unformatted',action='read',iostat=ier)
            if (ier.ne.0)then 
                write(*,'(a, i0.6)')'Couldnt read ystore file for proc ', self%iset
                stop
            endif 
            read(IIN)ystore_dp
            close(IIN)
            
            ! Open the z coordinate and load: 
            open(unit=IIN,file=trim(binname)//'zstore.bin', &
            status='unknown',form='unformatted',action='read',iostat=ier)
            if (ier.ne.0)then 
                write(*,'(a, i0.6)')'Couldnt read zstore file for proc ', self%iset
                stop
            endif 
            read(IIN)zstore_dp
            close(IIN)

            ! Cast from DP to CUSTOM_REAL (could also be DP)
            self%xstore(:,:,:,:) = real(xstore_dp(:,:,:,:), kind=CUSTOM_REAL)
            self%ystore(:,:,:,:) = real(ystore_dp(:,:,:,:), kind=CUSTOM_REAL)
            self%zstore(:,:,:,:) = real(zstore_dp(:,:,:,:), kind=CUSTOM_REAL)

            if(verbose.ge.3)then
                write(*,'(/,a)')'  -- X coordinates:'
                write(*,*)'     --> min. value: ', minval(self%xstore)
                write(*,*)'     --> max. value: ', maxval(self%xstore)

                write(*,'(/,a)')'  -- Y coordinates:'
                write(*,*)'     --> min. value: ', minval(self%ystore)
                write(*,*)'     --> max. value: ', maxval(self%ystore)
        
                write(*,'(/,a)')'  -- Z coordinates:'
                write(*,*)'     --> min. value: ', minval(self%zstore)
                write(*,*)'     --> max. value: ', maxval(self%zstore)
            endif 
            
            deallocate(xstore_dp)
            deallocate(ystore_dp)
            deallocate(zstore_dp)
        
            return 
        end subroutine read_proc_coordinates
        
        
        
        subroutine read_integer_proc_variable(self, variable, varname)
            use params, only: datadir, verbose, IIN
        
            implicit none 
            include "precision.h"
            class(SetMesh) :: self

            ! IO variables: 
            character(len=*) :: varname 
            integer          :: variable(self%ngllx, self%nglly, self%ngllz, self%nspec)
        
            ! Local variables
            integer :: ier
            character(len=250) :: binname
        
            ! File name prefix: 
            write(binname,'(a,i0.6,a,i1,a)')trim(datadir)//'/proc',self%iset,'_'//'reg',self%region,'_'//trim(varname)//'.bin'
            if(verbose.ge.3)then
                write(*,'(/,/,a)')'• Reading variable called '//trim(varname)
                write(*,'(a,i1)')'  -- data type : integer'
                write(*,'(a)')'  -- file name : '//trim(binname)
            endif 
        
            ! Open the variable file and load: 
            open(unit=IIN,file=trim(binname), &
            status='unknown',form='unformatted',action='read',iostat=ier)
            if (ier.ne.0)then 
                write(*,'(a, i0.6)')'Couldnt read "'//trim(varname)//'" file for proc ', self%iset
                stop
            endif 
            read(IIN)variable
        
            if(verbose.ge.2)then
                write(*,'(a)')'  -- '//trim(varname)//' :'
                write(*,*)'     --> min. value: ', minval(variable)
                write(*,*)'     --> max. value: ', maxval(variable)
            endif 
        
            close(IIN)
            
            return 
        end subroutine read_integer_proc_variable
        

        subroutine get_unique_radii(self, save)
            ! Determines a list of the unique radii of the GLL point
            ! Points each GLL point to its radial value in the list 
            use params, only: verbose, datadir
            implicit none 
            include "constants.h"
    
            class(SetMesh) :: self
            logical        :: save
        
            ! Local variables: 
            logical :: local_store_unq_radii

            real(kind=CUSTOM_REAL) :: rr(self%ngllx, self%nglly, self%ngllz,self%nspec) , tol, rgll
            real(kind=CUSTOM_REAL), allocatable :: unique_r_tmp(:)
            integer :: i, j, k, ispec, size_r, unique_id, ur
            logical :: match
        
            ! Compute radii
            call allocate_if_unallocated(self%ngllx, self%nglly, self%ngllz, self%nspec, self%rstore)
            self%rstore  = (self%xstore**TWO + self%ystore**TWO + self%zstore**TWO)**HALF
        
            ! Maximum number of unique radii is this
            size_r =  size(rr)
            allocate(unique_r_tmp(size_r))
            unique_r_tmp = -ONE
        
            ! Set tolerance for 'same radii' 
            tol = 1.0e-4
        
            ! Initial value
            unique_id = 1
        
            ! Allocate the array that will store the unique radius IDs: 
            call allocate_if_unallocated(self%ngllx, self%nglly, & 
                                         self%ngllz, self%nspec, & 
                                         self%rad_id)
        
            do ispec = 1, self%nspec
                do i = 1, self%ngllx
                    do j = 1, self%nglly
                        do k = 1, self%ngllz
                            ! Radius of this gll 
                            rgll = self%rstore(i, j, k, ispec)
                            match = .false. 
        
                            ! Loop over current unique radii: 
                            do ur = 1, unique_id 
                                if (abs(unique_r_tmp(ur) - rgll) < tol)then 
                                    ! already a 'unique' radius
                                    match = .true.
                                    self%rad_id(i,j,k,ispec) = ur
                                endif 
                            enddo ! ur 
        
                            ! No matches found
                            if (.not. match) then 
                                ! a new 'unique' radius 
                                self%rad_id(i,j,k,ispec)     = unique_id
                                unique_r_tmp(unique_id) = rgll
                                unique_id               = unique_id + 1
                            endif
                        enddo 
                    enddo 
                enddo 
            enddo 
        
        
            ! Allocate new unique radius array of the correct size
            ! It might be that the very last node tested in unique and so 
            ! unique_id is actually 1 too large -- test if last radius is -1 and if so cut it 
            if (unique_r_tmp(unique_id).eq.-ONE)then 
                self%n_unique_rad = unique_id - 1
            else
                self%n_unique_rad = unique_id
            endif
        
        
            ! Copy over the unique values to the smaller array and deallocate long.
            allocate(self%unique_r(self%n_unique_rad)) 
            self%unique_r(1:self%n_unique_rad) = unique_r_tmp(1:self%n_unique_rad)
            deallocate(unique_r_tmp)
        
            ! Print some details
            if(verbose.ge.2) then 
                write(*,'(/,a)')'• Determine unique mesh radii'
                write(*,'(a,es8.1)')  '  -- tolerance value       :', tol
                write(*,'(a,f8.2,a)') '  -- matches radii within  :', tol*SCALE_R, ' metres'
                write(*,'(a,i8,a,i8)')'  -- number of unique radii:', self%n_unique_rad, ' out of ', size_r
            endif 
        
            if(minval(self%rad_id).le.0 .or. maxval(self%rad_id).gt.self%n_unique_rad)then 
                write(*,*)'Error with rad_id:' 
                write(*,'(a,i8)')  '  -- min id       :', minval(self%rad_id)
                write(*,'(a,i8)')  '  -- max id       :', maxval(self%rad_id)
                write(*,'(a,i8)')  '  -- n_unique_rad :',self%n_unique_rad
                stop
            endif        
            
        
            self%interp = create_PieceInterp(self%n_unique_rad)
            self%interp%radial = self%unique_r
            call self%interp%setup()
            call self%interp%create_interpolation_radial_map()

            if(save) call self%save_get_mesh_radii_results()
        
        end subroutine get_unique_radii
    

        subroutine setup_global_coordinate_arrays(self, save)
            implicit none 
            logical :: save 
            class(SetMesh) :: self 

            call allocate_if_unallocated(self%nglob, self%x_glob)
            call allocate_if_unallocated(self%nglob, self%y_glob)
            call allocate_if_unallocated(self%nglob, self%z_glob)
            call self%map_local_global_double_precision(self%xstore, self%x_glob, 0)
            call self%map_local_global_double_precision(self%ystore, self%y_glob, 0)
            call self%map_local_global_double_precision(self%zstore, self%z_glob, 0)
        
            if(save)call self%save_global_xyz()
        end subroutine setup_global_coordinate_arrays
    


        subroutine compute_rtp_from_xyz(self, save)
            ! Converts the stored xyz to r theta phi
            implicit none
            include "constants.h" 
            logical :: save
            integer ispec, i, j, k
            class(SetMesh) :: self 

            call allocate_if_unallocated(self%ngllx,self%nglly,self%ngllz,self%nspec, self%rstore)
            call allocate_if_unallocated(self%ngllx,self%nglly,self%ngllz,self%nspec, self%thetastore)
            call allocate_if_unallocated(self%ngllx,self%nglly,self%ngllz,self%nspec, self%phistore)


            do ispec = 1, self%nspec
                do i = 1, self%ngllx
                    do j = 1, self%nglly
                        do k = 1, self%ngllz
                            self%rstore(i,j,k,ispec) = (self%xstore(i,j,k,ispec)**TWO + &
                                                self%ystore(i,j,k,ispec)**TWO + & 
                                                self%zstore(i,j,k,ispec)**TWO)**HALF
                                                
                            ! 0 <= phi <= 2pi
                            self%phistore(i,j,k,ispec) = atan2(self%ystore(i,j,k,ispec),&
                                                        self%xstore(i,j,k,ispec))   
                            if(self%phistore(i,j,k,ispec) .lt. ZERO)& 
                                self%phistore(i,j,k,ispec)  = TWO_PI + self%phistore(i,j,k,ispec) 

                            ! 0 <= theta <= pi                        
                            self%thetastore(i,j,k,ispec) = PI_OVER_TWO -  atan2(self%zstore(i,j,k,ispec),&
                                                                        (self%xstore(i,j,k,ispec)**TWO +&
                                                                        self%ystore(i,j,k,ispec)**TWO)**HALF) 
                        enddo 
                    enddo 
                enddo
            enddo 

            if(save)call self%save_elem_rtp()

        end subroutine compute_rtp_from_xyz



        subroutine compute_rotation_matrix(self)
            ! Computes rotation matrix R that converts a vector in r, theta, phi
            ! to cartesian coordinates
            ! R holds the normalised unit vectors r, theta, phi in cartesian coords
            ! as columns e.g. 
            !     ( r_x  θ_x  φ_x )
            ! R = ( r_y  θ_y  φ_y ) 
            !     ( r_z  θ_z  φ_z )

            use params, only : verbose, safety_checks
            use math, only: sinp, cosp, sqrtp, acosp, atan2p
            implicit none 
            include "constants.h"
            class(SetMesh) :: self 

            real(kind=CUSTOM_REAL) :: x, y, z, r, theta, phi, ct, cp, st, sp, norm
            integer :: iglob, p

            !TODO:  check allocations of x, y, z glob
            if(verbose.ge.2) write(*,'(/,a)')'• Computing rotation matrix'

            if(safety_checks)then 
                if(.not.allocated(self%x_glob))then 
                    write(*,*)'ERROR: X_glob not allocated. Stop.'
                    stop
                endif 
                if(.not.allocated(self%y_glob))then 
                    write(*,*)'ERROR: Y_glob not allocated. Stop.'
                    stop
                endif 
                if(.not.allocated(self%z_glob))then 
                    write(*,*)'ERROR: Z_glob not allocated. Stop.'
                    stop
                endif 
            endif 
            
            call allocate_if_unallocated(3, 3, self%nglob, self%Rmat)

            do iglob = 1, self%nglob
                ! Get x, y, z coordinates: 
                x = real(self%x_glob(iglob), kind=CUSTOM_REAL)
                y = real(self%y_glob(iglob), kind=CUSTOM_REAL)
                z = real(self%z_glob(iglob), kind=CUSTOM_REAL)

                r = sqrtp(x ** 2 + y ** 2 + z ** 2)

                if (r.eq.zero)then 
                    ! For now we will wont rotate it if the central GLL point
                    self%Rmat(:, :, iglob) = zero
                    do p = 1, 3
                        self%Rmat(p, p, iglob) = one
                    enddo 
                else
                    ! Note that theta and phi here are not necessarily the same
                    ! as in DT98 i.e. theta here will be between -pi/2 and pi/2
                    theta = acosp(z/r)
                    phi   = atan2p(y,x)
                    if (phi .lt. zero) phi = phi + TWO_PI

                    ct = real(cosp(theta), kind=CUSTOM_REAL)
                    cp = real(cosp(phi),   kind=CUSTOM_REAL)
                    st = real(sinp(theta), kind=CUSTOM_REAL)
                    sp = real(sinp(phi),   kind=CUSTOM_REAL)

                    ! Radial vector
                    !norm =((st*cp)**TWO + (st*sp)**TWO + ct**TWO)**half
                    self%Rmat(1, 1, iglob) = st*cp !/ norm
                    self%Rmat(2, 1, iglob) = st*sp !/ norm
                    self%Rmat(3, 1, iglob) = ct    !/ norm

                    ! Theta vector
                    !norm =((ct*cp)**TWO + (ct*sp)**TWO + st**TWO)**half
                    self%Rmat(1, 2, iglob) = ct*cp !/ norm
                    self%Rmat(2, 2, iglob) = ct*sp !/ norm
                    self%Rmat(3, 2, iglob) = -st   !/ norm

                    ! Phi vector
                    !norm =(sp*sp + cp*cp)**half
                    self%Rmat(1, 3, iglob) = -sp !/ norm
                    self%Rmat(2, 3, iglob) =  cp !/ norm
                    self%Rmat(3, 3, iglob) = zero
                endif
            enddo 
        end subroutine compute_rotation_matrix



        subroutine rotate_complex_vector_rtp_to_xyz(self, vector)
            implicit none  
            class(SetMesh) :: self 
            complex(kind=SPLINE_REAL) :: vector(3, self%ngllx, self%nglly, self%ngllz, self%nspec)

            integer :: i,j,k,ispec

            do ispec = 1, self%nspec
                do i = 1, self%ngllx
                    do j = 1, self%nglly
                        do k = 1, self%ngllz
                            vector(:,i,j,k,ispec) = matmul(real(self%Rmat(:,:,self%ibool(i,j,k,ispec)), kind=SPLINE_REAL), & 
                                                            vector(:,i,j,k,ispec))
                        enddo
                    enddo
                enddo 
            enddo
        end subroutine rotate_complex_vector_rtp_to_xyz



        subroutine rotate_complex_matrix_rtp_to_xyz(self, matrix)
            implicit none  
            class(SetMesh) :: self 
            complex(kind=SPLINE_REAL) :: matrix(3,3,self%ngllx, self%nglly, self%ngllz, self%nspec)
            complex(kind=SPLINE_REAL) :: RR(3,3)

            integer :: i,j,k,ispec

            do ispec = 1, self%nspec
                do i = 1, self%ngllx
                    do j = 1, self%nglly
                        do k = 1, self%ngllz
                            ! F' = R F R^t
                            RR = real(self%Rmat(:,:,self%ibool(i,j,k,ispec)), kind=SPLINE_REAL)
                            matrix(:,:,i,j,k,ispec) = matmul(matmul(RR,  matrix(:,:,i,j,k,ispec)), transpose(RR))
                        enddo
                    enddo
                enddo 
            enddo
        end subroutine rotate_complex_matrix_rtp_to_xyz



        subroutine rotate_complex_sym_matrix_rtp_to_xyz(self, symmat)
            ! Assumes the matrix is symmetric 3 x 3 in voigt ordering for matrix M
            ! 1  =  M_rr     4  =  M_tp
            ! 2  =  M_tt     5  =  M_rp
            ! 3  =  M_pp     6  =  M_rt
            implicit none 
            class(SetMesh) :: self

            complex(kind=SPLINE_REAL) :: symmat(6, self%ngllx, self%nglly, self%ngllz, self%nspec), M(3,3)
            real(kind=SPLINE_REAL) ::  R(3,3)
            integer :: i, j, k, ispec

            do ispec = 1, self%nspec
                do i = 1, self%ngllx
                    do j = 1, self%nglly
                        do k = 1, self%ngllz

                            R = real(self%Rmat(:,:, self%ibool(i,j,k,ispec)), kind=SPLINE_REAL)

                            ! Build back to 3 x 3
                            M(1,1) = symmat(1, i, j, k, ispec)
                            M(2,2) = symmat(2, i, j, k, ispec)
                            M(3,3) = symmat(3, i, j, k, ispec)

                            M(2,3) = symmat(4, i, j, k, ispec)
                            M(3,2) = symmat(4, i, j, k, ispec)

                            M(1,3) = symmat(5, i, j, k, ispec)
                            M(3,1) = symmat(5, i, j, k, ispec)

                            M(1,2) = symmat(6, i, j, k, ispec)
                            M(2,1) = symmat(6, i, j, k, ispec)

                            M = matmul(matmul(R, M), transpose(R))

                            ! Store back in symmat
                            symmat(1, i, j, k, ispec) = M(1,1) ! xx
                            symmat(2, i, j, k, ispec) = M(2,2) ! yy
                            symmat(3, i, j, k, ispec) = M(3,3) ! zz
                            symmat(4, i, j, k, ispec) = M(2,3) ! yz
                            symmat(5, i, j, k, ispec) = M(1,3) ! xz
                            symmat(6, i, j, k, ispec) = M(1,2) ! xy
                        enddo
                    enddo
                enddo 
            enddo
        end subroutine rotate_complex_sym_matrix_rtp_to_xyz



        subroutine compute_jacobian(self, save)
            use params, only: verbose
            use math, only: mat_inv
            implicit none 
            include "constants.h"
            class(SetMesh) :: self

            ! Save: 
            logical :: save

            ! Local variables: 
            integer :: i, j, s, t, n, p, ispec
            real(kind=CUSTOM_REAL) :: val, jl(3,3), tmp(3,3)


    
            if(verbose.ge.2)then
                write(*,'(/,a)')'• Setting up jacobian'
            endif 
    
            call allocate_if_unallocated(3, 3, self%ngllx, self%nglly, self%ngllz, self%nspec, self%jac)
            call allocate_if_unallocated(3, 3, self%ngllx, self%nglly, self%ngllz, self%nspec, self%jacinv)
            call allocate_if_unallocated(self%ngllx, self%nglly, self%ngllz, self%nspec, self%detjac)
    
            self%detjac = zero 
            self%jac    = zero 
    
            do ispec = 1, self%nspec
                do s = 1, self%ngllx
                    do t = 1, self%nglly
                        do n = 1, self%ngllz
                            ! Jij = del x_i/del Xi_j
    
                            ! Loop over x,y,z
                            do i = 1,3                      
                                
                                ! Loop over xi, eta, zeta
                                do j = 1, 3
                                    val = zero
    
                                    do p = 1, self%ngllx
                                        if (j.eq.1) then     ! xi 
                                            if     (i .eq. 1) then 
                                                val = val + self%xstore(p, t, n, ispec) * self%dgll(p, s)
                                            elseif (i .eq. 2) then 
                                                val = val + self%ystore(p, t, n, ispec) * self%dgll(p, s)
                                            else
                                                val = val + self%zstore(p, t, n, ispec) * self%dgll(p, s)
                                            endif 
                                        elseif (j.eq.2) then ! eta 
                                            if     (i .eq. 1) then 
                                                val = val + self%xstore(s, p, n, ispec) * self%dgll(p, t)
                                            elseif (i .eq. 2) then 
                                                val = val + self%ystore(s, p, n, ispec) * self%dgll(p, t)
                                            else
                                                val = val + self%zstore(s, p, n, ispec) * self%dgll(p, t)
                                            endif
                                        else                 ! zeta
                                            if     (i .eq. 1) then 
                                                val = val + self%xstore(s, t, p, ispec) * self%dgll(p, n)
                                            elseif (i .eq. 2) then 
                                                val = val + self%ystore(s, t, p, ispec) * self%dgll(p, n)
                                            else
                                                val = val + self%zstore(s, t, p, ispec) * self%dgll(p, n)
                                            endif
                                        endif 
                                    enddo              
                                    
                                    self%jac(i,j,s,t,n,ispec) = val 
                                enddo !j 
                            enddo! i 
    

                            

                            ! Compute determinant for this GLL point 
                            jl(:,:) = self%jac(:,:,s,t,n,ispec)
    
                            self%detjac(s,t,n,ispec) = jl(1,1) * ( jl(2,2)*jl(3,3) - jl(2,3)*jl(3,2))  - & 
                                                       jl(1,2) * ( jl(2,1)*jl(3,3) - jl(2,3)*jl(3,1))  + & 
                                                       jl(1,3) * ( jl(2,1)*jl(3,2) - jl(2,2)*jl(3,1)) 
                            
                            ! Compute the inverse jacobian: 
                            ! J^{-1}ij = del Xi_i/del x_j
                            self%jacinv(:,:,s,t,n,ispec) = mat_inv(jl)

                        enddo ! n 
                    enddo !t
                enddo ! s 
            enddo ! ispec


    
            if(verbose.ge.2)write(*,'(a)')'  --> done'
    
            if(save)then 
                if(verbose.ge.2)write(*,'(a)')'  --> saving jacobian to binary'
                call self%save_jacobian()
            endif 
        end subroutine compute_jacobian



        ! This subroutine recomputes the 3D Jacobian for one element
        ! based upon 125 GLL points
        ! Hejun Zhu OCT16,2009

        ! input: myrank,
        !        xstore,ystore,zstore ----- input GLL point coordinate
        !        xigll,yigll,zigll ----- GLL points position
        !        ispec,nspec       ----- element number

        ! output: xixstore, xiystore, xizstore,
        !         etaxstore,etaystore,etazstore,
        !         gammaxstore,gammaystore,gammazstore ------ parameters used to calculate Jacobian
        subroutine recalc_jacobian_gll3D(self)
            use params, only: verbose
            implicit none
            include "constants.h"
            class(SetMesh) :: self

            ! input parameter
            integer :: ispec

            ! local parameters for this subroutine
            double precision,dimension(self%NGLLX) :: hxir,hpxir
            double precision,dimension(self%NGLLY) :: hetar,hpetar
            double precision,dimension(self%NGLLZ) :: hgammar,hpgammar

            double precision :: xxi,xeta,xgamma,yxi,yeta,ygamma,zxi,zeta,zgamma
            double precision :: xi1 ,eta,gamma
            double precision :: hlagrange,hlagrange_xi,hlagrange_eta,hlagrange_gamma
            double precision :: jacobian,jacobian_inv
            double precision :: xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz
            double precision :: x,y,z

            integer :: i,j,k,i1,j1,k1
                

            if(verbose.ge.2)then
                write(*,'(/,a)')'• Setting up jacobian'
            endif 
    
            call allocate_if_unallocated(3, 3, self%ngllx, self%nglly, self%ngllz, self%nspec, self%jac)
            call allocate_if_unallocated(3, 3, self%ngllx, self%nglly, self%ngllz, self%nspec, self%jacinv)
            call allocate_if_unallocated(self%ngllx, self%nglly, self%ngllz, self%nspec, self%detjac)
    
            self%detjac = zero 
            self%jac    = zero 
    


            ! first go over all 125 GLL points
            do ispec = 1, self%nspec
                do k = 1, self%NGLLZ
                    do j = 1, self%NGLLY
                        do i = 1, self%NGLLX

                        ! Assumes ngll same in all directions
                        xi1   = self%xi(i)
                        eta   = self%xi(j)
                        gamma = self%xi(k)


                        ! calculate Lagrange polynomial and its derivative
                        call lagrange_any(xi1,    self%NGLLX, self%xi, hxir,   hpxir)
                        call lagrange_any(eta,   self%NGLLY,  self%xi, hetar,  hpetar)
                        call lagrange_any(gamma, self%NGLLZ,  self%xi, hgammar,hpgammar)

                        xxi     = ZERO
                        xeta    = ZERO
                        xgamma  = ZERO
                        yxi     = ZERO
                        yeta    = ZERO
                        ygamma  = ZERO
                        zxi     = ZERO
                        zeta    = ZERO
                        zgamma  = ZERO

                        do k1 = 1, self%NGLLZ
                            do j1 = 1, self%NGLLY
                                do i1 = 1, self%NGLLX

                                hlagrange = hxir(i1) * hetar(j1) * hgammar(k1)
                                hlagrange_xi = hpxir(i1) * hetar(j1) * hgammar(k1)
                                hlagrange_eta = hxir(i1) * hpetar(j1) * hgammar(k1)
                                hlagrange_gamma = hxir(i1) * hetar(j1) * hpgammar(k1)

                                x = self%xstore(i1,j1,k1,ispec)
                                y = self%ystore(i1,j1,k1,ispec)
                                z = self%zstore(i1,j1,k1,ispec)

                                xxi     = xxi    + x * hlagrange_xi
                                xeta    = xeta   + x * hlagrange_eta
                                xgamma  = xgamma + x * hlagrange_gamma

                                yxi     = yxi    + y * hlagrange_xi
                                yeta    = yeta   + y * hlagrange_eta
                                ygamma  = ygamma + y * hlagrange_gamma

                                zxi     = zxi    + z * hlagrange_xi
                                zeta    = zeta   + z * hlagrange_eta
                                zgamma  = zgamma + z * hlagrange_gamma
                                enddo
                            enddo
                        enddo

                        ! Determinany of the Jacobian calculation
                        self%detjac(i,j,k,ispec) = xxi*(yeta*zgamma-ygamma*zeta) - &
                                                   xeta*(yxi*zgamma-ygamma*zxi)  + &
                                                   xgamma*(yxi*zeta-yeta*zxi)

                        ! invert the relation (Fletcher p. 50 vol. 2)
                        jacobian_inv = ONE / self%detjac(i,j,k,ispec)

                        xix    = (yeta*zgamma-ygamma*zeta) * jacobian_inv
                        xiy    = (xgamma*zeta-xeta*zgamma) * jacobian_inv
                        xiz    = (xeta*ygamma-xgamma*yeta) * jacobian_inv
                        etax   = (ygamma*zxi-yxi*zgamma) * jacobian_inv
                        etay   = (xxi*zgamma-xgamma*zxi) * jacobian_inv
                        etaz   = (xgamma*yxi-xxi*ygamma) * jacobian_inv
                        gammax = (yxi*zeta-yeta*zxi) * jacobian_inv
                        gammay = (xeta*zxi-xxi*zeta) * jacobian_inv
                        gammaz = (xxi*yeta-xeta*yxi) * jacobian_inv

                        ! resave the derivatives and the Jacobian
                        ! distinguish between single and double precision for reals
                        self%jacinv(1,1,i,j,k,ispec) = real(xix,    kind=CUSTOM_REAL)
                        self%jacinv(2,1,i,j,k,ispec) = real(etax,   kind=CUSTOM_REAL)
                        self%jacinv(3,1,i,j,k,ispec) = real(gammax, kind=CUSTOM_REAL)

                        self%jacinv(1,2,i,j,k,ispec) = real(xiy,    kind=CUSTOM_REAL)
                        self%jacinv(2,2,i,j,k,ispec) = real(etay,   kind=CUSTOM_REAL)
                        self%jacinv(3,2,i,j,k,ispec) = real(gammay, kind=CUSTOM_REAL)

                        self%jacinv(1,3,i,j,k,ispec) = real(xiz,    kind=CUSTOM_REAL)
                        self%jacinv(2,3,i,j,k,ispec) = real(etaz,   kind=CUSTOM_REAL)
                        self%jacinv(3,3,i,j,k,ispec) = real(gammaz, kind=CUSTOM_REAL)

                        enddo ! i
                    enddo ! j
                enddo ! k
            enddo ! ispec

        end subroutine recalc_jacobian_gll3D









        subroutine save_get_mesh_radii_results(self)
            use params, only: datadir, IOUT
            implicit none 

            class(SetMesh) :: self

            call self%interp%save(trim(datadir)//'/store/mesh_radii_data/data_'//self%set_str)
        end subroutine save_get_mesh_radii_results


        subroutine load_get_mesh_radii_results(self)
            use params, only: datadir, IIN 
            implicit none 
            class(SetMesh) :: self

            call self%interp%load(trim(datadir)//'/store/mesh_radii_data/data_'//self%set_str)
    
        end subroutine load_get_mesh_radii_results



        subroutine compute_mode_displacement(self, m, m0de, disp)
            use ylm_plm, only: ylm_complex, ylm_deriv
            use mesh_utils, only: delta_spline
            use math, only: sinp, sqrtp
            
            implicit none
            include "constants.h"
            
            class(SetMesh) :: self
            type(Mode)     :: m0de
            integer        :: m
        
            complex(kind=SPLINE_REAL) :: disp(3, self%ngllx, self%nglly, &
                                                 self%ngllz, self%nspec)


            ! Local variables: 
            
            complex(kind=CUSTOM_REAL) :: ylm, dylm_theta, dylm_phi
            real(kind=CUSTOM_REAL)    ::  theta, phi 
            integer :: i, j, k, ispec
        
            real(kind=SPLINE_REAL)    ::  mf, tl14p, mone_l, pref, u_r, v_r, & 
                                          w_r, dm0, dd1, sinth
            complex(kind=SPLINE_REAL) :: sp_ylm, spl_dylm_theta
        
            ! Allocate the GLL displacement array
            disp = SPLINE_iZERO
        
            ! Float version of m
            mf = real(m, kind=SPLINE_REAL)
    
            tl14p  = ((SPLINE_TWO*m0de%lf + SPLINE_ONE)/(SPLINE_FOUR*SPLINE_PI))**SPLINE_HALF    ! (2l+1/4π)^1/2
            dm0    = delta_spline(m, 0)                                                          ! δ_{m0}
            dd1    = delta_spline(m, -1) - delta_spline(m, 1)                                    ! δ_{m -1} - δ_{m 1}
            mone_l = (-SPLINE_ONE)**m0de%lf                                                      ! -1 ^l 
            pref   = SPLINE_HALF * tl14p * m0de%kf                                               ! 0.5 (2l+1/4π)^1/2  (l(l+1))^1/2
        
        
            if (m0de%t.eq.'S')then 
                do ispec = 1, self%nspec
                    do i = 1, self%ngllx
                        do j = 1, self%nglly
                            do k = 1, self%ngllz
                            
                                ! Get ylm and the partial derivatives
                                theta = real(self%thetastore(i,j,k,ispec), kind=CUSTOM_REAL)
                                phi   = real(self%phistore(i,j,k,ispec),   kind=CUSTOM_REAL)
                                ylm   = ylm_complex(m0de%l, m, theta, phi)
                                call ylm_deriv(m0de%l, m, theta, phi, dylm_theta, dylm_phi)
        
                                ! u, v at at the radius of this node
                                u_r   =    m0de%u_spl(self%rad_id(i,j,k,ispec))
                                v_r   =    m0de%v_spl(self%rad_id(i,j,k,ispec)) / m0de%kf
        
                                sinth = real(sinp(theta), kind=SPLINE_REAL)
        
                                if (theta.ge.zero .and. theta.le.pole_tolerance) then 
                                    ! North pole 
                                    ! S_r     (DT98 D.8)
                                    disp(1,i,j,k,ispec) = tl14p * u_r * dm0
                                    ! S_theta (DT98 D.9)
                                    disp(2,i,j,k,ispec) =  pref * v_r * dd1
                                    ! S_phi   (DT98 D.10)
                                    disp(3,i,j,k,ispec) = pref * SPLINE_iONE * mf * v_r * dd1
        
        
                                elseif (abs(theta-PI).le.pole_tolerance) then
                                    ! South pole
                                    ! S_r     (DT98 D.11)
                                    disp(1,i,j,k,ispec) =  mone_l * tl14p * u_r * dm0
                                    ! S_theta (DT98 D.12)
                                    disp(2,i,j,k,ispec) =  mone_l * pref * v_r * dd1
                                    ! S_phi   (DT98 D.13)
                                    disp(3,i,j,k,ispec) = -mone_l * pref * SPLINE_iONE * mf * v_r * dd1
                                else 
                                    ! Convert to SPLINE_REAL precision
                                    sp_ylm         = cmplx(ylm,        kind=SPLINE_REAL)
                                    spl_dylm_theta = cmplx(dylm_theta, kind=SPLINE_REAL)
        
                                    ! S_r     (DT98 D.4)
                                    disp(1,i,j,k,ispec) = sp_ylm*u_r  
                                    ! S_theta (DT98 D.5)
                                    disp(2,i,j,k,ispec) = v_r * spl_dylm_theta
                                    ! S_phi   (DT98 D.6)
                                    disp(3,i,j,k,ispec) = SPLINE_iONE * mf * v_r * sp_ylm / sinth
                                endif 
        
                            enddo 
                        enddo 
                    enddo 
                enddo 
            elseif (m0de%t.eq.'T')then 
                do ispec = 1, self%nspec
                    do i = 1, self%ngllx
                        do j = 1, self%nglly
                            do k = 1, self%ngllz
                                ! Get ylm and the partial derivatives, and 
                                ! w at the radius of this node
                                
                                w_r   = m0de%w_spl(self%rad_id(i,j,k,ispec))/m0de%kf
        
                                theta = real(self%thetastore(i,j,k,ispec), kind=CUSTOM_REAL)
                                phi   = real(self%phistore(i,j,k,ispec), kind=CUSTOM_REAL)
                                ylm   = ylm_complex(m0de%l, m, theta, phi)
                                call ylm_deriv(m0de%l, m, theta, phi, dylm_theta, dylm_phi)
        
                                if (theta.ge.zero .and. theta.le.pole_tolerance) then 
                                    ! North pole 
                                    ! S_theta (DT98 D.9)
                                    disp(2,i,j,k,ispec) = pref * SPLINE_iONE * mf * w_r * dd1
                                    ! S_phi   (DT98 D.10)
                                    disp(3,i,j,k,ispec) = - pref * w_r * dd1
                                    
                                elseif (abs(theta-PI).le.pole_tolerance) then
                                    ! South pole
                                    ! S_theta (DT98 D.12)
                                    disp(2,i,j,k,ispec) = - mone_l * pref * SPLINE_iONE * mf * w_r * dd1
        
                                    ! S_phi   (DT98 D.13)
                                    disp(3,i,j,k,ispec) = - mone_l * pref * w_r * dd1
                                else 
        
                                    ! Convert to SPLINE_REAL precision
                                    sp_ylm         = cmplx(ylm,        kind=SPLINE_REAL)
                                    spl_dylm_theta = cmplx(dylm_theta, kind=SPLINE_REAL)
        
                                    ! S_theta (DT98 D.5)
                                    disp(2,i,j,k,ispec) = SPLINE_iONE * mf * w_r * sp_ylm / sinth
                                    ! S_phi   (DT98 D.6)
                                    disp(3,i,j,k,ispec) = -w_r * spl_dylm_theta 
                                endif 
        
                                ! No radial displacement for toroidal modes
                                disp(1,i,j,k,ispec) = SPLINE_iZERO
                            enddo 
                        enddo 
                    enddo 
                enddo 
        
            else
                write(*,*)'ERROR: Mode type must be S or T but was ', m0de%t
                stop
            endif 
        
        end subroutine compute_mode_displacement
        


        subroutine compute_mode_gradphi(self, m, m0de, gradphi)
            use ylm_plm, only: ylm_complex, ylm_deriv
            use mesh_utils, only: delta_spline
            use math, only: sinp, sqrtp
            
            implicit none
            include "constants.h"
            
            class(SetMesh) :: self
            type(Mode)     :: m0de
            integer        :: m
        
            complex(kind=SPLINE_REAL) :: gradphi(3, self%ngllx, self%nglly, &
                                                    self%ngllz, self%nspec)      

            ! Local variables: 
            complex(kind=CUSTOM_REAL) :: ylm, dylm_theta, dylm_phi
            real(kind=CUSTOM_REAL)    ::  theta, phi
            integer :: i, j, k, ispec
        
            real(kind=SPLINE_REAL)    ::  mf, tl14p, mone_l, pref, p_r, dp_r, & 
                                          w_r, dm0, dd1, sinth, unq_r
            complex(kind=SPLINE_REAL) :: sp_ylm, spl_dylm_theta
        
            ! Allocate the GLL displacement array
            gradphi = SPLINE_iZERO
        
            ! Float version of m
            mf = real(m, kind=SPLINE_REAL)
            
            if (m0de%t.eq.'S')then 
                do ispec = 1, self%nspec
                    do i = 1, self%ngllx
                        do j = 1, self%nglly
                            do k = 1, self%ngllz
                            
                                ! Get ylm and the partial derivatives
                                unq_r = real(self%rstore(i,j,k,ispec),     kind=SPLINE_REAL)
                                theta = real(self%thetastore(i,j,k,ispec), kind=CUSTOM_REAL)
                                phi   = real(self%phistore(i,j,k,ispec),   kind=CUSTOM_REAL)
                                ylm   = ylm_complex(m0de%l, m, theta, phi)
                                call ylm_deriv(m0de%l, m, theta, phi, dylm_theta, dylm_phi)
        
                                ! u, v at at the radius of this node
                                p_r   = m0de%p_spl(self%rad_id(i,j,k,ispec))
                                dp_r  = m0de%dp_spl(self%rad_id(i,j,k,ispec))
                                sinth = real(sinp(theta), kind=SPLINE_REAL)
        
                                if (theta.ge.zero .and. theta.le.pole_tolerance) then 
                                    ! North pole 
                                    gradphi(:,i,j,k,ispec) = SPLINE_iZERO
        
                                elseif (abs(theta-PI).le.pole_tolerance) then
                                    ! South pole
                                    gradphi(:,i,j,k,ispec) = SPLINE_iZERO
                                elseif (unq_r.eq.zero .or. sinth.eq.zero)then 
                                    ! Set artificially at centre
                                    gradphi(:,i,j,k,ispec) = SPLINE_iZERO
                                else 
                                    ! Convert to SPLINE_REAL precision
                                    sp_ylm         = cmplx(ylm,        kind=SPLINE_REAL)
                                    spl_dylm_theta = cmplx(dylm_theta, kind=SPLINE_REAL)
        
                                    ! radial
                                    gradphi(1,i,j,k,ispec) = sp_ylm*dp_r  
                                    ! theta 
                                    gradphi(2,i,j,k,ispec) = p_r * spl_dylm_theta / unq_r
                                    ! phi  
                                    gradphi(3,i,j,k,ispec) = SPLINE_iONE * mf * p_r * sp_ylm / (sinth*unq_r)
                                endif 
                            enddo 
                        enddo 
                    enddo 
                enddo 
            elseif (m0de%t.eq.'T')then 
                ! Since p is going to be zero so will the grad phi 
                gradphi = SPLINE_iZERO
            else
                write(*,*)'ERROR: Mode type must be S or T but was ', m0de%t
                stop
            endif 
        
        end subroutine compute_mode_gradphi




        subroutine compute_mode_grad_SR(self, m, m0de, gradSR)
            use ylm_plm, only: ylm_complex, ylm_deriv
            use mesh_utils, only: delta_spline
            use math, only: sinp, sqrtp
            
            implicit none
            include "constants.h"
            
            class(SetMesh) :: self
            type(Mode)     :: m0de
            integer        :: m
        
            complex(kind=SPLINE_REAL) :: gradSR(3, self%ngllx, self%nglly, &
                                                    self%ngllz, self%nspec)      

            ! Local variables: 
            complex(kind=CUSTOM_REAL) :: ylm, dylm_theta, dylm_phi
            real(kind=CUSTOM_REAL)    ::  theta, phi
            integer :: i, j, k, ispec
        
            real(kind=SPLINE_REAL)    ::  mf, tl14p, mone_l, pref, u_r, du_r, & 
                                          dm0, dd1, sinth, unq_r
            complex(kind=SPLINE_REAL) :: sp_ylm, spl_dylm_theta
        
            ! Allocate the GLL displacement array
            gradSR = SPLINE_iZERO
        
            ! Float version of m
            mf = real(m, kind=SPLINE_REAL)
    
            tl14p  = ((SPLINE_TWO*m0de%lf + SPLINE_ONE)/(SPLINE_FOUR*SPLINE_PI))**SPLINE_HALF    ! (2l+1/4π)^1/2
            dm0    = delta_spline(m, 0)                                                          ! δ_{m0}
            dd1    = delta_spline(m, -1) - delta_spline(m, 1)                                    ! δ_{m -1} - δ_{m 1}
            mone_l = (-SPLINE_ONE)**m0de%lf                                                      ! -1 ^l 
            pref   = SPLINE_HALF * tl14p * m0de%kf                                               ! 0.5 (2l+1/4π)^1/2  (l(l+1))^1/2
        
            if (m0de%t.eq.'S')then 
                do ispec = 1, self%nspec
                    do i = 1, self%ngllx
                        do j = 1, self%nglly
                            do k = 1, self%ngllz
                            
                                ! Get ylm and the partial derivatives
                                unq_r = real(self%rstore(i,j,k,ispec),     kind=SPLINE_REAL)
                                theta = real(self%thetastore(i,j,k,ispec), kind=CUSTOM_REAL)
                                phi   = real(self%phistore(i,j,k,ispec),   kind=CUSTOM_REAL)
                                ylm   = ylm_complex(m0de%l, m, theta, phi)
                                call ylm_deriv(m0de%l, m, theta, phi, dylm_theta, dylm_phi)
        
                                ! u, v at at the radius of this node
                                u_r   = m0de%u_spl(self%rad_id(i,j,k,ispec))
                                du_r  = m0de%du_spl(self%rad_id(i,j,k,ispec))
                                sinth = real(sinp(theta), kind=SPLINE_REAL)
        
                                if (theta.ge.zero .and. theta.le.pole_tolerance) then 
                                    ! North pole 
                                    gradSR(:,i,j,k,ispec) = SPLINE_iZERO
        
                                elseif (abs(theta-PI).le.pole_tolerance) then
                                    ! South pole
                                    gradSR(:,i,j,k,ispec) = SPLINE_iZERO
                                elseif (unq_r.eq.zero .or.sinth.eq.zero)then 
                                    gradSR(:,i,j,k,ispec) = SPLINE_iZERO
                                else 
                                    ! Convert to SPLINE_REAL precision
                                    sp_ylm         = cmplx(ylm,        kind=SPLINE_REAL)
                                    spl_dylm_theta = cmplx(dylm_theta, kind=SPLINE_REAL)
        
                                    ! radial
                                    gradSR(1,i,j,k,ispec) = sp_ylm*du_r  
                                    ! theta 
                                    gradSR(2,i,j,k,ispec) = u_r * spl_dylm_theta / unq_r
                                    ! phi  
                                    gradSR(3,i,j,k,ispec) = SPLINE_iONE * mf * u_r * sp_ylm / (unq_r*sinth)
                                endif 
                            enddo 
                        enddo 
                    enddo 
                enddo 
            elseif (m0de%t.eq.'T')then 
                ! Since sr is radial, it is going to be zero 
                gradSR = SPLINE_iZERO
            else
                write(*,*)'ERROR: Mode type must be S or T but was ', m0de%t
                stop
            endif 
        
        end subroutine compute_mode_grad_SR





        subroutine compute_mode_strain(self, m, m0de, strain)
            ! Strain ordering is as follows: 
            ! 1  =  E_rr     4  =  E_tp
            ! 2  =  E_tt     5  =  E_rp
            ! 3  =  E_pp     6  =  E_rt
            ! Note that this is differnet from what I would usually do as an order (e.g. that of a moment tensor)
            ! This is voigt notation convention
            use params, only: all_warnings
            use ylm_plm, only: ylm_complex, ylm_deriv
            use mesh_utils, only: delta_spline
            use math, only: sinp, tanp, sqrtp
        
            implicit none
            include "constants.h"
        
            ! IO variables: 
            class(SetMesh)   :: self
            integer          :: m             ! m value (order)  of mode
            type(Mode)       :: m0de          ! Mode of interest
            complex(kind=SPLINE_REAL) :: strain(6, self%ngllx, self%nglly, self%ngllz, self%nspec)
        
            ! Local variables: 
            complex(kind=CUSTOM_REAL) :: ylm, dylm_theta, dylm_phi
            real(kind=CUSTOM_REAL)    :: theta, phi
            complex(kind=SPLINE_REAL) :: sp_ylm, spl_dylm_theta
            real(kind=SPLINE_REAL)    :: sinth, tanth, unq_r, du_r, u_r, & 
                                         dv_r, v_r, w_r, dw_r, ff, &
                                         dm0, dd1, dd2
            ! GLL level variables
            real(kind=SPLINE_REAL)  :: xx_r, zz_r, mf, lf, ll1, tl14p, kr2
            integer :: rid 

            ! Loop variables 
            integer :: ispec, i, j, k, h
        

            mf  = real(m, kind=SPLINE_REAL)           ! float m

            tl14p = ((SPLINE_TWO*m0de%lf + SPLINE_ONE)/(SPLINE_FOUR*SPLINE_PI))**SPLINE_HALF    ! (2l+1/4π)^1/2
            kr2   = m0de%kf*((m0de%kf*m0de%kf - SPLINE_TWO)**SPLINE_HALF)      
        
            dm0 = delta_spline(m, 0)                         ! δ_{m0}
            dd1 = delta_spline(m, -1) - delta_spline(m, 1)   ! δ_{m -1} - δ_{m 1}
            dd2 = delta_spline(m, -2) + delta_spline(m, 2)   ! δ_{m -2} + δ_{m 2}
        
            ! Convert eigenfunctions to auxillary form
            ! DT98 D.1: 
            !    u = U     v = V/k     w = W/k     p = P
        
            ! Compute x and z auxillary variables - DT98 D.20
            ! If spheroidal mode then u,v are non zero but w, wdot are 0
            !   --> z should be 0 
            ! If toroidal mode then u,v are zero but w, wdot are not
            !   --> x should be 0
            if (m0de%t.eq.'S')then 
                ! Loop over each GLL point: 
                do ispec = 1, self%nspec
                    do i = 1, self%ngllx
                        do j = 1, self%nglly
                            do k = 1, self%ngllz
                                ! Get ylm and the partial derivatives
                                theta = self%thetastore(i,j,k,ispec)
                                phi   = self%phistore(i,j,k,ispec)
                                ylm   = ylm_complex(m0de%l, m, theta, phi)
                                call ylm_deriv(m0de%l, m, theta, phi, dylm_theta, dylm_phi)
                                
                                ! u, v, udot, vdo, x, at at the radius of this node
                                rid   = self%rad_id(i,j,k,ispec)
                                u_r   = m0de%u_spl(rid)
                                du_r  = m0de%du_spl(rid)
                                
                                v_r   = m0de%v_spl(rid)/m0de%kf
                                xx_r  = m0de%aux_x(rid)
        
                                ! r and theta at the node  
                                unq_r = real(self%rstore(i,j,k,ispec), kind=SPLINE_REAL)
                                sinth = real(sinp(theta),         kind=SPLINE_REAL)
                                tanth = real(tanp(theta),         kind=SPLINE_REAL)
        
        
                                if(unq_r.eq.0)then 
                                    if(all_warnings) write(*,*)'Warning: setting strain at centre to 0 artificially'
                                    strain(:,i,j,k,ispec) = SPLINE_iZERO
                                else
                                    if (theta.ge.zero .and. theta.le.pole_tolerance) then 
                                        ! Values at the pole 
                                        ! D.28
                                        ff = (SPLINE_TWO*u_r - m0de%kf*m0de%kf*v_r) / unq_r    
                                        
                                        ! E_rr: DT98 D.22
                                        strain(1,i,j,k,ispec) = tl14p * du_r * dm0
        
                                        ! E_tt: DT98 D.23
                                        strain(2,i,j,k,ispec) = SPLINE_HALF * tl14p * &
                                                                (ff*dm0 + & 
                                                                kr2*SPLINE_TWO*v_r * & 
                                                                dd2/(SPLINE_FOUR*unq_r))
        
                                        ! E_pp: DT98 D.24
                                        strain(3,i,j,k,ispec) = SPLINE_HALF * tl14p * &
                                                                (ff*dm0 - & 
                                                                kr2*SPLINE_TWO*v_r * & 
                                                                dd2/(SPLINE_FOUR*unq_r))
        
                                        ! E_rt: DT98 D.25
                                        strain(6,i,j,k,ispec) = tl14p * m0de%kf * xx_r * dd1/SPLINE_FOUR
                                        ! E_rp: DT98 D.26
                                        strain(5,i,j,k,ispec) = tl14p * m0de%kf * SPLINE_iONE * mf * xx_r/SPLINE_FOUR
                                        ! E_tp: DT98 D.27
                                        strain(4,i,j,k,ispec) = tl14p * kr2 * SPLINE_iONE * mf * v_r * dd2/(SPLINE_EIGHT*unq_r)
                                    else 
        
                                        ! Convert to SPLINE_REAL precision
                                        sp_ylm         = cmplx(ylm, kind=SPLINE_REAL)
                                        spl_dylm_theta = cmplx(dylm_theta, kind=SPLINE_REAL)
        
                                        ! Values away from the pole
                                        ! E_rr: DT98 D.14
                                        strain(1,i,j,k,ispec) = sp_ylm*du_r  
                                        ! E_tt: DT98 D.15
                                        strain(2,i,j,k,ispec) = (sp_ylm*u_r - v_r*(spl_dylm_theta/tanth -  &
                                                                sp_ylm*(mf/sinth)**SPLINE_TWO + m0de%kf*m0de%kf*sp_ylm))/unq_r
                                        ! E_pp: DT98 D.16
                                        strain(3,i,j,k,ispec) = (sp_ylm*u_r + v_r*(spl_dylm_theta/tanth -  &
                                                                 sp_ylm*(mf/sinth)**SPLINE_TWO))/unq_r
                                        ! E_rt: DT98 D.17
                                        strain(6,i,j,k,ispec) = SPLINE_HALF * xx_r * spl_dylm_theta
                                        ! E_rp: DT98 D.18
                                        strain(5,i,j,k,ispec) = SPLINE_HALF * SPLINE_iONE * mf * xx_r * sp_ylm/sinth
                                        ! E_tp: DT98 D.19
                                        strain(4,i,j,k,ispec) = SPLINE_iONE * mf * v_r * (spl_dylm_theta - sp_ylm/tanth) / (unq_r * sinth)
                                    endif        
                                endif !unique r
        
                                
                                do h = 1, 6
                                    if (abs(strain(h,i,j,k,ispec)).gt. 100.0)then
                                        strain(h,i,j,k,ispec) = SPLINE_ZERO
                                        if(all_warnings)write(*,*)'Setting strain values at pi to 0,'
                                        !stop 
                                    endif 
                                enddo 
        
                            enddo 
                        enddo 
                    enddo 
                enddo
        
            elseif (m0de%t.eq.'T')then 
            
                ! Loop over each GLL point: 
                do ispec = 1, self%nspec
                    do i = 1, self%ngllx
                        do j = 1, self%nglly
                            do k = 1, self%ngllz
                                ! Get ylm and the partial derivatives
                                theta = self%thetastore(i,j,k,ispec)
                                phi   = self%phistore(i,j,k,ispec)
                                ylm   = ylm_complex(m0de%l,m, theta, phi)
                                call ylm_deriv(m0de%l, m, theta, phi, dylm_theta, dylm_phi)
        
                                ! w, wdot, z, at at the radius of this node
                                rid   = self%rad_id(i,j,k,ispec)
                                w_r   = m0de%w_spl(rid)/m0de%kf
                                dw_r  = m0de%dw_spl(rid)/m0de%kf
                                zz_r  = m0de%aux_z(rid) ! No division as divided above
        
                                ! r and theta at the node  
                                unq_r = real(self%rstore(i,j,k,ispec), kind=SPLINE_REAL)
                                sinth = real(sinp(theta),              kind=SPLINE_REAL)
                                tanth = real(tanp(theta),              kind=SPLINE_REAL)
        
                                if (theta.ge.zero .and. theta.le.pole_tolerance) then 
                                    ! Values at the pole 
                                    ! E_rr: DT98 D.22
                                    strain(1,i,j,k,ispec) = SPLINE_iZERO
                                    ! E_tt: DT98 D.23
                                    strain(2,i,j,k,ispec) = SPLINE_HALF * tl14p * SPLINE_iONE * mf * w_r * dd2 * kr2/(SPLINE_FOUR * unq_r)
                                    ! E_pp: DT98 D.24
                                    strain(3,i,j,k,ispec) = -SPLINE_HALF * tl14p * SPLINE_iONE * mf * w_r * dd2 * kr2/(SPLINE_FOUR * unq_r)
                                    ! E_rt: DT98 D.25
                                    strain(6,i,j,k,ispec) = tl14p * ll1 * SPLINE_iONE * mf * zz_r * dd1 / SPLINE_FOUR
                                    ! E_rp: DT98 D.26
                                    strain(5,i,j,k,ispec) = - tl14p * ll1 * zz_r * dd1 / SPLINE_FOUR
                                    ! E_tp: DT98 D.27
                                    strain(4,i,j,k,ispec) = - tl14p * kr2 * w_r * dd2 / (unq_r * SPLINE_FOUR)
                                else
        
                                    ! Convert to SPLINE_REAL precision
                                    sp_ylm         = cmplx(ylm, kind=SPLINE_REAL)
                                    spl_dylm_theta = cmplx(dylm_theta, kind=SPLINE_REAL)
        
                                    ! E_rr: DT98 D.14
                                    strain(1,i,j,k,ispec) = SPLINE_iZERO
                                    ! E_tt: DT98 D.15
                                    strain(2,i,j,k,ispec) = SPLINE_iONE * mf * w_r * (spl_dylm_theta - sp_ylm/tanth) / (unq_r*sinth)
                                    ! E_pp: DT98 D.16
                                    strain(3,i,j,k,ispec) = - strain(2,i,j,k,ispec)
                                    ! E_rt: DT98 D.17
                                    strain(6,i,j,k,ispec) = SPLINE_HALF * SPLINE_iONE * mf * zz_r * sp_ylm / sinth
                                    ! E_rp: DT98 D.18
                                    strain(5,i,j,k,ispec) =  - SPLINE_HALF * zz_r * spl_dylm_theta
                                    ! E_tp: DT98 D.19
                                    strain(4,i,j,k,ispec) = (spl_dylm_theta/tanth + & 
                                                                 (SPLINE_HALF*m0de%kf*m0de%kf - (mf/sinth)**SPLINE_TWO)*sp_ylm & 
                                                                 )* w_r / unq_r
                                endif 
                            enddo 
                        enddo 
                    enddo 
                enddo
            else
                write(*,*)'Error: mode_type must be S or T'
                stop
            endif 
        
        
        end subroutine compute_mode_strain





        subroutine compute_mode_gradS(self, m, m0de, gradS)
            use params, only: all_warnings
            use ylm_plm, only: ylm_complex, ylm_deriv
            use mesh_utils, only: delta_spline
            use math, only: sinp, tanp, sqrtp, cosp
        
            implicit none
            include "constants.h"
        
            ! IO variables: 
            class(SetMesh)   :: self
            integer          :: m             ! m value (order)  of mode
            type(Mode)       :: m0de          ! Mode of interest
            complex(kind=SPLINE_REAL) :: gradS(3,3, self%ngllx, self%nglly, self%ngllz, self%nspec)
        
            ! Local variables: 
            complex(kind=CUSTOM_REAL) :: ylm, dylm_theta, dylm_phi
            real(kind=CUSTOM_REAL)    :: theta, phi
            complex(kind=SPLINE_REAL) :: sp_ylm, spl_dylm_theta
            real(kind=SPLINE_REAL)    :: sinth, tanth, costh, unq_r, du_r, u_r, & 
                                         dv_r, v_r, w_r, dw_r, ff, &
                                         dm0, dd1, dd2
            ! GLL level variables
            real(kind=SPLINE_REAL)  :: xx_r, zz_r, mf, lf, ll1, tl14p, kr2
            integer :: rid 

            ! Loop variables 
            integer :: ispec, i, j, k, h, q
        
            mf  = real(m, kind=SPLINE_REAL)           ! float m

            if (m0de%t.eq.'S')then 
                ! Loop over each GLL point: 
                do ispec = 1, self%nspec
                    do i = 1, self%ngllx
                        do j = 1, self%nglly
                            do k = 1, self%ngllz
                                ! Get ylm and the partial derivatives
                                theta = self%thetastore(i,j,k,ispec)
                                phi   = self%phistore(i,j,k,ispec)
                                ylm   = ylm_complex(m0de%l, m, theta, phi)
                                call ylm_deriv(m0de%l, m, theta, phi, dylm_theta, dylm_phi)
                                
                                ! u, v, udot, vdo, x, at at the radius of this node
                                rid   = self%rad_id(i,j,k,ispec)
                                u_r   = m0de%u_spl(rid)
                                du_r  = m0de%du_spl(rid)
                                
                                v_r   = m0de%v_spl(rid)/m0de%kf
                                dv_r  = m0de%dv_spl(rid)/m0de%kf
                                xx_r  = m0de%aux_x(rid)
        
                                ! r and theta at the node  
                                unq_r = real(self%rstore(i,j,k,ispec), kind=SPLINE_REAL)
                                sinth = real(sinp(theta),              kind=SPLINE_REAL)
                                tanth = real(tanp(theta),              kind=SPLINE_REAL)
                                costh = real(cosp(theta),              kind=SPLINE_REAL)
        
        
                                if(unq_r.eq.0)then 
                                    if(all_warnings) write(*,*)'Warning: setting Grad S at centre to 0 artificially'
                                    gradS(:,:,i,j,k,ispec) = SPLINE_iZERO
                                else
                                    if (theta.ge.zero .and. theta.le.pole_tolerance) then 
                                        ! Values at the pole 
                                        ! D.28
                                        !ff = (SPLINE_TWO*u_r - m0de%kf*m0de%kf*v_r) / unq_r    
                                        ! E_rr: DT98 D.22
                                        gradS(:,:,i,j,k,ispec) = ZERO
        
                                    else 
        
                                        ! Convert to SPLINE_REAL precision
                                        sp_ylm         = cmplx(ylm,        kind=SPLINE_REAL)
                                        spl_dylm_theta = cmplx(dylm_theta, kind=SPLINE_REAL)
        
                                        ! Values away from the pole
                                        ! d_r S_r - same as  E_rr: DT98 D.14
                                        gradS(1,1,i,j,k,ispec) = sp_ylm*du_r  
                                        
                                        ! d_r S_t  (4):
                                        gradS(2,1,i,j,k,ispec) = dv_r*spl_dylm_theta 
                                        
                                        ! d_r S_p  (6):
                                        gradS(3,1,i,j,k,ispec) = SPLINE_iONE*mf*dv_r*sp_ylm/sinth

                                        ! d_t S_r  (5):
                                        gradS(1,2,i,j,k,ispec) = (u_r - v_r)*spl_dylm_theta/unq_r
                                        
                                        ! d_t S_t - same as E_tt: DT98 D.15
                                        gradS(2,2,i,j,k,ispec) = (sp_ylm*u_r - v_r*(spl_dylm_theta/tanth -  &
                                                                  sp_ylm*(mf/sinth)**SPLINE_TWO + m0de%kf*m0de%kf*sp_ylm))/unq_r

                                        ! d_t S_p  (8):
                                        gradS(3,2,i,j,k,ispec) = SPLINE_iONE*mf*v_r *(spl_dylm_theta - sp_ylm*costh/sinth )/(unq_r*sinth)
                                        

                                        ! d_p S_r  (7):
                                        gradS(1,3,i,j,k,ispec) = SPLINE_iONE * mf * sp_ylm * (u_r - v_r)/(unq_r*sinth)

                                        ! d_p S_t  (9):
                                        gradS(2,3,i,j,k,ispec) =  SPLINE_iONE*mf*v_r*(spl_dylm_theta - costh*sp_ylm/sinth)/(unq_r*sinth)

                                        ! d_p S_p - same as E_pp: DT98 D.16
                                        gradS(3,3,i,j,k,ispec) = (sp_ylm*u_r + v_r*(spl_dylm_theta/tanth -  &
                                                                  sp_ylm*(mf/sinth)**SPLINE_TWO))/unq_r



                                    endif        
                                endif !unique r
                                
                                do q = 1, 3
                                    do h = 1, 3
                                        if (abs(gradS(h,q,i,j,k,ispec)).gt. 100.0)then
                                            gradS(h,q,i,j,k,ispec) = SPLINE_ZERO
                                            if(all_warnings)write(*,*)'Setting strain values at pi to 0,'
                                            !stop 
                                        endif 
                                    enddo 
                                enddo 
        
                            enddo 
                        enddo 
                    enddo 
                enddo
        
            elseif (m0de%t.eq.'T')then 
                write(*,*)'Grad U not implemented yet for Toroidal modes.'
                stop
            else
                write(*,*)'Error: mode_type must be S or T'
                stop
            endif 
        
        
        end subroutine compute_mode_gradS




        subroutine compute_strain_deviator(self, GradS, devE)
            ! Computes the deviatoric strain based on the gradient of the mode 
            ! displacement. Since d_ij = E_ij - 1/3 E_kk \delta_ij
            ! d_ij = 1/2 (del_i s_j + del_j s_i) - 1/3 del_k s_k \delta_ij

            class(SetMesh)   :: self
            complex(kind=SPLINE_REAL), dimension(3, 3, self%ngllx, & 
                                                       self%nglly, &
                                                       self%ngllz, &
                                                       self%nspec) :: GradS, devE
            complex(kind=SPLINE_REAL) :: gs(3,3), trace

            ! local variables: 
            integer :: i,j,k,ispec, p, q

            devE = SPLINE_iZERO

            do i = 1, self%ngllx
                do j = 1, self%nglly
                    do k = 1, self%ngllz
                        do ispec = 1, self%nspec

                            gs(:,:) = GradS(:,:,i,j,k,ispec)
                            trace   = (gs(1,1)+gs(2,2)+gs(3,3))/SPLINE_THREE

                            do p = 1, 3
                                do q = 1, 3
                                    devE(p,q,i,j,k,ispec) =  half*(gs(p,q)+gs(q,p))
                                    if(p.eq.q)then 
                                        devE(p,q,i,j,k,ispec) = devE(p,q,i,j,k,ispec) - trace 
                                    endif
                                enddo ! q
                            enddo ! p
                            
                        enddo ! i
                    enddo ! j
                enddo ! k
            enddo ! ispec
        end subroutine compute_strain_deviator


        subroutine compute_strain_deviator_congj(self, GradS, devE)
            ! Computes the deviatoric strain based on the gradient of the mode 
            ! displacement. Since d_ij = E_ij - 1/3 E_kk \delta_ij
            ! d_ij = 1/2 (del_i s_j + del_j s_i) - 1/3 del_k s_k \delta_ij

            class(SetMesh)   :: self
            complex(kind=SPLINE_REAL), dimension(3, 3, self%ngllx, & 
                                                       self%nglly, &
                                                       self%ngllz, &
                                                       self%nspec) :: GradS, devE
            complex(kind=SPLINE_REAL) :: gs(3,3), trace

            ! local variables: 
            integer :: i,j,k,ispec, p, q

            devE = SPLINE_iZERO

            do i = 1, self%ngllx
                do j = 1, self%nglly
                    do k = 1, self%ngllz
                        do ispec = 1, self%nspec

                            gs(:,:) = GradS(:,:,i,j,k,ispec)
                            trace   = (conjg(gs(1,1)+gs(2,2)+gs(3,3)))/SPLINE_THREE

                            do p = 1, 3
                                do q = 1, 3
                                    devE(p,q,i,j,k,ispec) =  half*(conjg(gs(p,q)+gs(q,p)))
                                    if(p.eq.q)then 
                                        devE(p,q,i,j,k,ispec) = devE(p,q,i,j,k,ispec) - trace 
                                    endif
                                enddo ! q
                            enddo ! p
                            
                        enddo ! i
                    enddo ! j
                enddo ! k
            enddo ! ispec
        end subroutine compute_strain_deviator_congj







        subroutine compute_background_g(self)
            ! Computes the background gravity acceleraiton magnitude
            !|g| at each point in g 
            use piecewise_interpolation, only: InterpPiecewise, create_PieceInterp
            use params, only: rho_spl, g_spl
            use mineos_model, only: MineosModel, mineos
            use gravitation, only: compute_background_g
            implicit none
            include "precision.h"
            
            class(SetMesh) :: self

            real(kind=CUSTOM_REAL) :: temp_r, thisrad
            integer :: temp_i, size, i, j, npoints, knot_lower, iknot, iout
            integer, dimension(self%n_unique_rad) :: idx

            real(kind=CUSTOM_REAL) :: r_lower, r_upper, store_unqr(self%n_unique_rad)
            type(InterpPiecewise) :: minInterp

            ! Store a copy
            store_unqr = self%unique_r

            ! First we want to sort the radial values into an order 
            ! but keep track of their indices
            size = self%n_unique_rad
            idx = (/ (i, i=1, size) /)


            ! Simple bubble sort (can be replaced by more efficient algorithms)
            do i = 1, size-1
                do j = 1, size-i
                    if (self%unique_r(j) > self%unique_r(j+1)) then
                        ! Swap values in arr
                        temp_r =  self%unique_r(j)
                        self%unique_r(j) = self%unique_r(j+1)
                        self%unique_r(j+1) = temp_r
                        
                        ! Swap corresponding indices
                        temp_i = idx(j)
                        idx(j) = idx(j+1)
                        idx(j+1) = temp_i
                    end if
                end do
            end do

    
            ! Now the radii are in order and we can use the idx to map back 
            ! let us now compute the gravity at 1500 points between 
            ! the centre and surface, using the mineos model: 
            npoints = 1500 
            r_lower    = zero    
            r_upper    = SCALE_R


            minInterp = create_PieceInterp(npoints)


            minInterp%radial = [((r_lower + (real(j-1)/real(npoints-1))*(r_upper-r_lower))/SCALE_R, &
                               j = 1, npoints)] 


            call minInterp%setup()
            call minInterp%create_interpolation_radial_map()


            call deallocate_if_allocated(rho_spl)
            allocate(rho_spl(npoints))


            call minInterp%interpolate_mineos_variable(real(mineos%rho_mineos, kind=SPLINE_REAL), rho_spl)
            
            ! Compute the g at these 1500 points
            allocate(g_spl(npoints))


            call compute_background_g(minInterp%radial, rho_spl, npoints, g_spl)

            ! Now let us treat these 1500 points as the function which we will spline to 
            ! get the values at each of the sorted 'unq rad' values
            ! Since then function is continuous lets have a naive interpolation method: 

            knot_lower = 0

            allocate(self%gmag_at_r(self%n_unique_rad))

            do iout = 1, self%n_unique_rad
                thisrad = self%unique_r(iout)
                do iknot = knot_lower, npoints
                    ! loop until we hit the knot above the radii we want
                    if (thisrad.lt.minInterp%radial(iknot)) exit
                enddo 

                knot_lower = iknot- 1
                
                self%gmag_at_r(idx(iout)) = g_spl(knot_lower) + & 
                                       (thisrad - minInterp%radial(knot_lower)) *& 
                                       (g_spl(knot_lower+1)-g_spl(knot_lower))
            enddo 


            ! Restore the original unique r: 
            self%unique_r = store_unqr


            deallocate(g_spl)
        end subroutine compute_background_g


















        
        

        subroutine save_mode_disp_binary(self, n, t, l, m, disp_id)
            use params, only: datadir
            implicit none 
            include "constants.h"
            
            class(SetMesh) :: self 
            integer :: m, n, l, disp_id
            character :: t
            character(len = 250) :: binname
            
            call create_mode_binary_fname(n, t, l, m, self%iset, binname)
        
            open(1, file=trim(datadir)//'/store/disp/'//trim(binname), form='unformatted')
            if(disp_id.eq.1)then 
                write(1)self%disp1
            elseif(disp_id.eq.2)then 
                write(1)self%disp2
            endif 
            close(1)
        
        end subroutine
        
        
        subroutine load_mode_disp_binary(self, n, t, l, m, disp_id)
            use params, only: datadir
            implicit none 
            include "constants.h"
            class(SetMesh) :: self 
            integer :: m, n, l, disp_id
            character :: t
            character(len = 250) ::  binname
            
            call create_mode_binary_fname(n, t, l, m, self%iset, binname)

            open(1, file=trim(datadir)//'/store/disp/'//trim(binname), form='unformatted')

            if(disp_id.eq.1)then 
                read(1)self%disp1
            elseif(disp_id.eq.2)then 
                read(1)self%disp2
            endif 

            close(1)
        
        end subroutine
        

        subroutine save_mode_strain_binary(self, n, t, l, m, mode_id)
            use params, only: datadir
            implicit none 
            include "constants.h"
            class(SetMesh) :: self 
            integer :: m, n, l, mode_id
            character :: t
            character(len = 250) ::  binname
            
            call create_mode_binary_fname(n, t, l, m, self%iset, binname)
            open(1, file=trim(datadir)//'/store/strain/'//trim(binname), form='unformatted')
            if(mode_id.eq.1)then 
                write(1)self%strain1
            elseif(mode_id.eq.2)then 
                write(1)self%strain2
            endif 
            close(1)
        end subroutine save_mode_strain_binary
        
        
        subroutine load_mode_strain_binary(self, n, t, l, m, strain_arr)
            
            use params, only: datadir
            implicit none 
            include "constants.h"
        
            class(SetMesh) :: self 
            integer     :: m, n, l, ios, mode_id
            character   :: t
            character(len = 250) :: binname
            complex(kind=SPLINE_REAL) :: strain_arr(6, self%ngllx, self%nglly, self%ngllz, self%nspec)
            
            call create_mode_binary_fname(n, t, l, m, self%iset, binname)
        
            open(1, file=trim(datadir)//'/store/strain/'//trim(binname), form='unformatted', iostat=ios)
            if(ios.ne.0)then
                write(*,*)'Could not open mode', trim(binname)
                stop
            endif 

            read(1) strain_arr
            close(1)
        
        end subroutine load_mode_strain_binary
        


        subroutine load_rho_spline(self)
            use params, only: datadir, rho_spl
            implicit none 
            class(SetMesh) :: self

            ! Local 
            integer :: stored_len
            character(len=250) :: fname 

            call create_rhospline_fname(self%iset, fname)

            open(1, file=trim(datadir)//'/store/rho/'//trim(fname), form='unformatted')

            ! Read the data from the binary file
            read(1) stored_len

            if(stored_len.ne.self%n_unique_rad)then 
            write(*,*)'Error: n_unique_rad is different from the length of rhospline stored in the binary you are trying to load: '
            write(*,*)'n_unique_rad : ', self%n_unique_rad
            write(*,*)'stored value : ', stored_len
            stop 
            endif

            call deallocate_if_allocated(rho_spl)
            allocate(rho_spl(self%n_unique_rad))

            read(1)self%unique_r
            read(1)rho_spl
            close(1)
        end subroutine load_rho_spline


        subroutine load_jacobian(self)
            use params, only: datadir, IIN
            implicit none 

            class(SetMesh)     :: self
        
            call allocate_if_unallocated(3, 3, self%ngllx, self%nglly, self%ngllz, self%nspec, self%jac)
            call allocate_if_unallocated(3, 3, self%ngllx, self%nglly, self%ngllz, self%nspec, self%jacinv)
            call allocate_if_unallocated(self%ngllx, self%nglly, self%ngllz, self%nspec, self%detjac)
        
            ! Density
            open(IIN, file=trim(datadir)//'/store/jacobian/jacdata_'//self%set_str, form='unformatted')
            read(IIN)self%jac
            read(IIN)self%detjac
            read(IIN)self%jacinv
            close(IIN)
        
        end subroutine load_jacobian

    
        subroutine save_jacobian(self)
            use params, only: datadir, IOUT
            implicit none 
            class(SetMesh)     :: self

            open(IOUT, file=trim(datadir)//'/store/jacobian/jacdata_'//self%set_str, form='unformatted')
            write(IOUT)self%jac
            write(IOUT)self%detjac
            write(IOUT)self%jacinv
            close(IOUT)
        end subroutine save_jacobian
        
        subroutine save_wglljac(self)
            use params, only: datadir, IOUT
            implicit none 
            class(SetMesh)     :: self
            open(IOUT, file=trim(datadir)//'/store/jacobian/wglljac_'//self%set_str, form='unformatted')
            write(IOUT)self%wglljac
            close(IOUT)
        end subroutine save_wglljac
        
        
        subroutine load_wglljac(self)
            use params, only: datadir, IIN
            implicit none 
            class(SetMesh)     :: self

            call deallocate_if_allocated(self%wglljac)
            allocate(self%wglljac(self%ngllx,self%nglly,self%ngllz,self%nspec))
            open(IIN, file=trim(datadir)//'/store/jacobian/wglljac_'//self%set_str, form='unformatted')
            read(IIN)self%wglljac
            close(IIN)
        end subroutine load_wglljac
        
        
        subroutine save_global_xyz(self)
            use params, only: datadir, IOUT
            implicit none         
            class(SetMesh)     :: self
            open(IOUT, file=trim(datadir)//'/store/global_xyz/coords_'//self%set_str, form='unformatted')
            write(IOUT)self%x_glob
            write(IOUT)self%y_glob
            write(IOUT)self%z_glob
            close(IOUT)
        end subroutine save_global_xyz
        
        
        subroutine load_global_xyz(self)
            use params, only: datadir, IIN
            implicit none       
            class(SetMesh)     :: self

            call allocate_if_unallocated(self%nglob, self%x_glob)
            call allocate_if_unallocated(self%nglob, self%y_glob)
            call allocate_if_unallocated(self%nglob, self%z_glob)
        
            open(IIN, file=trim(datadir)//'/store/global_xyz/coords_'//self%set_str, form='unformatted')
            read(IIN)self%x_glob
            read(IIN)self%y_glob
            read(IIN)self%z_glob
            close(IIN)
        end subroutine load_global_xyz
        
        
    
        subroutine save_elem_rtp(self)
            use params, only: datadir, IOUT
            implicit none 
            class(SetMesh) :: self

            open(IOUT, file=trim(datadir)//'/store/elemental_rtp/elem_coords_'//self%set_str, form='unformatted')
            write(IOUT)self%rstore
            write(IOUT)self%thetastore
            write(IOUT)self%phistore
            close(IOUT)
        end subroutine save_elem_rtp
        
        
        subroutine load_elem_rtp(self)
            use params, only: datadir, IIN
            implicit none 
            class(SetMesh) :: self
        
            call allocate_if_unallocated(self%ngllx, self%nglly, self%ngllz, self%nspec, self%rstore)
            call allocate_if_unallocated(self%ngllx, self%nglly, self%ngllz, self%nspec, self%thetastore)
            call allocate_if_unallocated(self%ngllx, self%nglly, self%ngllz, self%nspec, self%phistore)
        
            open(IIN, file=trim(datadir)//'/store/elemental_rtp/elem_coords_'//self%set_str, form='unformatted')
            read(IIN)self%rstore
            read(IIN)self%thetastore
            read(IIN)self%phistore
            close(IIN)
        end subroutine load_elem_rtp


        subroutine map_local_global_real_4(self, loc, glob, direction)
            ! mappings between local and global variables
            ! direction: 0   local  --> global 
            !            1   global --> local 
            implicit none 
            include "precision.h"
            class(SetMesh) :: self 

            ! I/O variables: 
            real(kind=4) :: loc(self%ngllx,self%nglly,self%ngllz,self%nspec)
            real(kind=4) :: glob(self%nglob)
            integer :: direction
    
            ! Local variables: 
            integer :: ispec, i, j, k
        
            ! Ensure we have ibool ready to use
            call self%check_ibool_is_defined()
    
            if (direction.eq.0)then 
                ! Map local to global
                do ispec = 1, self%nspec 
                    do i = 1, self%ngllx 
                        do j = 1, self%nglly
                            do k = 1, self%ngllz
                                glob(self%ibool(i,j,k,ispec)) = loc(i,j,k,ispec)
                            enddo 
                        enddo 
                    enddo 
                enddo
            elseif(direction.eq.1)then 
                ! Map global to local
                do ispec = 1,self% nspec 
                    do i = 1, self%ngllx 
                        do j = 1, self%nglly
                            do k = 1, self%ngllz
                                loc(i,j,k,ispec) = glob(self%ibool(i,j,k,ispec))
                            enddo 
                        enddo 
                    enddo 
                enddo
            else
                write(*,*)'ERROR: direction can only be 1 or 0 but has value ', direction
                stop
            endif 
    
            return 
        end subroutine map_local_global_real_4
    
        subroutine map_local_global_double_precision(self, loc, glob, direction)
            ! mappings between local and global variables
            ! direction: 0   local  --> global 
            !            1   global --> local 
            implicit none 
            include "precision.h"
            class(SetMesh) :: self 

            ! I/O variables: 
            real(kind=8) :: loc(self%ngllx,self%nglly,self%ngllz,self%nspec)
            real(kind=8) :: glob(self%nglob)
            integer :: direction
    
            ! Local variables: 
            integer :: ispec, i, j, k
        
            ! Ensure we have ibool ready to use
            call self%check_ibool_is_defined()
    
            if (direction.eq.0)then 
                ! Map local to global
                do ispec = 1, self%nspec 
                    do i = 1, self%ngllx 
                        do j = 1, self%nglly
                            do k = 1, self%ngllz
                                glob(self%ibool(i,j,k,ispec)) = loc(i,j,k,ispec)
                            enddo 
                        enddo 
                    enddo 
                enddo
            elseif(direction.eq.1)then 
                ! Map global to local
                do ispec = 1, self%nspec 
                    do i = 1, self%ngllx 
                        do j = 1, self%nglly
                            do k = 1, self%ngllz
                                loc(i,j,k,ispec) = glob(self%ibool(i,j,k,ispec))
                            enddo 
                        enddo 
                    enddo 
                enddo
            else
                write(*,*)'ERROR: direction can only be 1 or 0 but has value ', direction
                stop
            endif 
    
            return 
        end subroutine map_local_global_double_precision
    
    
    
        subroutine map_local_global_complex_4(self, loc, glob, direction)
            ! mappings between local and global variables
            ! direction: 0   local  --> global 
            !            1   global --> local 
    
            implicit none 
            include "precision.h"
            class(SetMesh) :: self
    
            ! I/O variables: 
            complex(kind=4) :: loc(self%ngllx,self%nglly,self%ngllz,self%nspec)
            complex(kind=4) :: glob(self%nglob)
            integer :: direction
    
            ! Local variables: 
            integer :: ispec, i, j, k
        
            ! Ensure we have ibool ready to use
            call self%check_ibool_is_defined()
    
            if (direction.eq.0)then 
                ! Map local to global
                do ispec = 1, self%nspec 
                    do i = 1, self%ngllx 
                        do j = 1, self%nglly
                            do k = 1, self%ngllz
                                glob(self%ibool(i,j,k,ispec)) = loc(i,j,k,ispec)
                            enddo 
                        enddo 
                    enddo 
                enddo
            elseif(direction.eq.1)then 
                ! Map global to local
                do ispec = 1, self%nspec 
                    do i = 1, self%ngllx 
                        do j = 1, self%nglly
                            do k = 1, self%ngllz
                                loc(i,j,k,ispec) = glob(self%ibool(i,j,k,ispec))
                            enddo 
                        enddo 
                    enddo 
                enddo
            else
                write(*,*)'ERROR: direction can only be 1 or 0 but has value ', direction
                stop
            endif 
            return 
        end subroutine map_local_global_complex_4

    
        subroutine map_local_global_complex_8(self, loc, glob, direction)
            ! mappings between local and global variables
            ! direction: 0   local  --> global 
            !            1   global --> local 
            implicit none 
            include "precision.h"
            class(SetMesh) :: self
    
            ! I/O variables: 
            complex(kind=8) :: loc(self%ngllx,self%nglly,self%ngllz,self%nspec)
            complex(kind=8) :: glob(self%nglob)
            integer :: direction
    
            ! Local variables: 
            integer :: ispec, i, j, k
        
            ! Ensure we have ibool ready to use
            call self%check_ibool_is_defined()
    
            if (direction.eq.0)then 
                ! Map local to global
                do ispec = 1, self%nspec 
                    do i = 1, self%ngllx 
                        do j = 1, self%nglly
                            do k = 1, self%ngllz
                                glob(self%ibool(i,j,k,ispec)) = loc(i,j,k,ispec)
                            enddo 
                        enddo 
                    enddo 
                enddo
            elseif(direction.eq.1)then 
                ! Map global to local
                do ispec = 1, self%nspec 
                    do i = 1, self%ngllx 
                        do j = 1, self%nglly
                            do k = 1, self%ngllz
                                loc(i,j,k,ispec) = glob(self%ibool(i,j,k,ispec))
                            enddo 
                        enddo 
                    enddo 
                enddo
            else
                write(*,*)'ERROR: direction can only be 1 or 0 but has value ', direction
                stop
            endif 
            return 
        end subroutine map_local_global_complex_8
    
    
        subroutine map_complex_vector_4(self, dim1, loc, glob, direction)
            ! Wrapper of map_local_global_complex that maps vector
            ! arrays from local <--> global 
            ! direction: 0   local  --> global 
            !            1   global --> local 
            implicit none 
            include "precision.h"
            class(SetMesh) :: self
    
            ! IO variables:
            integer :: dim1, direction
            complex(kind=4) :: loc(dim1, self%ngllx, self%nglly, self%ngllz, self%nspec)
            complex(kind=4) :: glob(dim1, self%nglob)
    
            ! Local: 
            complex(kind=4) :: loc_tmp(self%ngllx, self%nglly, self%ngllz, self%nspec)
            complex(kind=4) :: glob_tmp(self%nglob)
            integer         :: i
    
            if (direction.eq.0)then
                ! local -> global
                do i = 1, dim1
                    loc_tmp(:, :, :, :) = loc(i, :, :, :, :)
                    call self%map_local_global_complex_4(loc_tmp, glob_tmp, direction)
                    glob(i,:) = glob_tmp
                enddo
            elseif(direction.eq.1)then
                ! global -> local
                do i = 1, dim1
                    glob_tmp(:) = glob(i, :)
                    call self%map_local_global_complex_4(loc_tmp, glob_tmp, direction)
                    loc(i, :, :, :, :) = loc_tmp(:, :, :, :)
                enddo
            endif  
        end subroutine map_complex_vector_4
    
    
    
        subroutine map_complex_vector_8(self, dim1, loc, glob, direction)
            ! Wrapper of map_local_global_complex that maps vector
            ! arrays from local <--> global 
            ! direction: 0   local  --> global 
            !            1   global --> local 
            implicit none 
            include "precision.h"
            class(SetMesh) :: self

            ! IO variables:
            integer :: dim1, direction
            complex(kind=8) :: loc(dim1, self%ngllx, self%nglly, self%ngllz, self%nspec)
            complex(kind=8) :: glob(dim1, self%nglob)
    
            ! Local: 
            complex(kind=8) :: loc_tmp(self%ngllx, self%nglly, self%ngllz, self%nspec)
            complex(kind=8) :: glob_tmp(self%nglob)
            integer         :: i
    
            if (direction.eq.0)then
                ! local -> global
                do i = 1, dim1
                    loc_tmp(:, :, :, :) = loc(i, :, :, :, :)
                    call self%map_local_global_complex_8(loc_tmp, glob_tmp, direction)
                    glob(i,:) = glob_tmp
                enddo
            elseif(direction.eq.1)then
                ! global -> local
                do i = 1, dim1
                    glob_tmp(:) = glob(i, :)
                    call self%map_local_global_complex_8(loc_tmp, glob_tmp, direction)
                    loc(i, :, :, :, :) = loc_tmp(:, :, :, :)
                enddo
            endif  
        end subroutine map_complex_vector_8




            
        subroutine compute_wglljac(self, save)
            implicit none
            integer :: iproc 
            logical :: save 
            integer :: i, j, k, ispec
            class(SetMesh) :: self

            call deallocate_if_allocated(self%wglljac)
            allocate(self%wglljac(self%ngllx,self%nglly,self%ngllz,self%nspec))

            do ispec = 1, self%nspec
                do i = 1, self%ngllx
                    do j = 1, self%nglly
                        do k = 1, self%ngllz 
                            self%wglljac(i,j,k,ispec) = self%wgll(i) * & 
                                                        self%wgll(j) * & 
                                                        self%wgll(k) * & 
                                                        self%detjac(i,j,k,ispec)
                        enddo 
                    enddo
                enddo 
            enddo 

            if(save)then 
                call self%save_wglljac()
            endif 

        end subroutine compute_wglljac



        subroutine setup_gll(self)
            use params, only: verbose
            use gll, only: get_gll, lagrange1st, zwgljd
            implicit none 
            class(SetMesh) :: self

            if (self%ngllx.ne.self%nglly .or. self%ngllx.ne.self%ngllz .or. &
                self%nglly.ne.self%ngllz)then
                write(*,*)'Error: currently only working for case of ngllx = nglly = ngllz'
                stop 
            endif 

            call allocate_if_unallocated(self%ngllx, self%xi)
            call allocate_if_unallocated(self%ngllx, self%wgll)

            call zwgljd(self%xi,self%wgll,self%ngllx,zero,zero)


            !call get_gll(self%ngllx-1, self%xi, self%wgll)

            if (verbose.ge.5)then 
                write(*,'(/,a)')'• Setup GLL points'
                write(*,'(/,a, i1)')'  --> ngll: ', self%ngllx
                write(*,'(/,a)')'  -->  gll: '
                write(*,*) self%xi
                write(*,'(/,a)')'  -->  wgll: '
                write(*,*)self%wgll
            endif 

            ! Get derivative of lagrange polynomials
            call allocate_if_unallocated(self%ngllx, self%ngllx, self%dgll)
            call lagrange1st(self%ngllx-1, self%xi, self%dgll)
            
        end subroutine setup_gll



        subroutine setup_mesh_sem_details(self, load_from_bin, save_to_bin)
            implicit none 
            class(SetMesh) :: self 
            logical        :: load_from_bin, save_to_bin

            ! Things that need to be done for each set  
            call self%read_proc_coordinates()
            call self%load_ibool()
            call self%setup_gll()

            if(load_from_bin)then 
                call self%load_jacobian()
                call self%load_wglljac()
                call self%load_global_xyz()
                call self%load_elem_rtp()
                call self%load_get_mesh_radii_results()
            else
                call self%compute_jacobian(save_to_bin)
                call self%compute_wglljac(save_to_bin)
                call self%setup_global_coordinate_arrays(save_to_bin)
                call self%compute_rtp_from_xyz(save_to_bin)
                call self%get_unique_radii(save_to_bin)
            endif 

        end subroutine setup_mesh_sem_details



        subroutine cleanup(self)
            implicit none 
            class(SetMesh) :: self

            call deallocate_if_allocated(self%xi)
            call deallocate_if_allocated(self%wgll)
            call deallocate_if_allocated(self%dgll)
            call deallocate_if_allocated(self%wglljac)
        
            call deallocate_if_allocated(self%xstore)
            call deallocate_if_allocated(self%ystore)
            call deallocate_if_allocated(self%zstore)
            call deallocate_if_allocated(self%rstore)
            call deallocate_if_allocated(self%thetastore)
            call deallocate_if_allocated(self%phistore)

            call deallocate_if_allocated(self%x_glob)
            call deallocate_if_allocated(self%y_glob)
            call deallocate_if_allocated(self%z_glob)

            call deallocate_if_allocated(self%detjac)
            call deallocate_if_allocated(self%jacinv)
            call deallocate_if_allocated(self%jac)

            call deallocate_if_allocated(self%ibool)

            call deallocate_if_allocated(self%unique_r)
            call deallocate_if_allocated(self%rad_id)

            call self%interp%cleanup()

            call deallocate_if_allocated(self%disp1)
            call deallocate_if_allocated(self%disp2)

            call deallocate_if_allocated(self%strain1)
            call deallocate_if_allocated(self%strain2)

            call deallocate_if_allocated(self%gradS_1)
            call deallocate_if_allocated(self%gradS_2)

            call deallocate_if_allocated(self%devE_1)
            call deallocate_if_allocated(self%devE_2)

            call deallocate_if_allocated(self%gradphi_1)
            call deallocate_if_allocated(self%gradphi_2)

            call deallocate_if_allocated(self%gSR_1)
            call deallocate_if_allocated(self%gSR_2)

            call deallocate_if_allocated(self%gmag_at_r)

        end subroutine cleanup


        function create_SetMesh(iset, region) result(SM)
            implicit none 
            integer, intent(in) :: iset, region
            type(SetMesh)       :: SM
            character(len=5)    :: proc

            SM%iset = iset
            SM%region = region
            call buffer_int(proc, SM%iset)
            SM%set_str = proc
        end function 

end module specfem_mesh