module piecewise_interpolation
    use mineos_model, only: MineosModel, mineos
    use params, only: verbose, safety_checks
    use modes, only: Mode 

    implicit none 
    include "constants.h"

    type ::  InterpPiecewise
        ! Radial values to be interpolated
        integer :: n_radial

        ! Interpolation maps: 
        integer              :: min_knot_id
        integer              :: max_knot_id
        integer, allocatable :: interp_map(:)
        integer, allocatable :: sec_interp_map(:,:)

        ! Piecewise breakup of Mineos model 
        integer :: num_sections ! Number of piecewise sections covered 
        integer :: lower_disc_id
        integer :: upper_disc_id
        integer :: max_per_sec   
        integer, allocatable :: orig_id(:,:)
        integer, allocatable :: n_per_sec(:)
        real(kind=CUSTOM_REAL), allocatable :: sec_rad(:,:)
        real(kind=SPLINE_REAL), allocatable :: sec_spl(:,:)

        real(kind=CUSTOM_REAL), allocatable :: radial(:)
        real(kind=CUSTOM_REAL) :: min_r, max_r

        ! Mineos model to be pointed to: 
        type(MineosModel), POINTER :: M    

        contains 
            procedure :: setup
            procedure :: create_interpolation_radial_map
            procedure :: interpolate_mineos_variable
            procedure :: interpolate_mode_eigenfunctions
            procedure :: print
            procedure :: load 
            procedure :: save 
            procedure :: cleanup 

    end type InterpPiecewise

    ! module contains
    contains 

        subroutine setup(self)
            implicit none 
            class(InterpPiecewise) :: self
    
            integer :: idisc, section_id, ss, i
            real(kind=CUSTOM_REAL) :: rr

            ! For now we are assuming that the 'mineos' object is being
            ! pointed to 
            self%M => mineos

            ! Get minimum and maximum radial values
            self%min_r = minval(self%radial)
            self%max_r = maxval(self%radial)      
            


            if(verbose.ge.4)then 
                write(*,*)'-- interpolating mineos variable -- '
                write(*,*)'  min_r = ', self%min_r
                write(*,*)'  max_r = ', self%max_r
                write(*,*)'  ndisc = ', self%M%ndisc
            endif 
    

            ! Find out the number of 'sections' separated by discontinuities
            ! that are represented in the radial array
            do idisc = 1, self%M%ndisc
                if(self%M%rdisc_norm(idisc).le.self%min_r)then 
                    self%lower_disc_id = idisc
                endif 
                if(self%M%rdisc_norm(idisc).lt.self%max_r)then 
                    self%upper_disc_id = idisc
                endif
            enddo 

            self%num_sections = self%upper_disc_id - self%lower_disc_id + 1



            ! We will now split the radial values into their piecewise sections
            ! Three arrays: 
            !   orig_id - the position in the overal new_rad array of that value
            !   sec_rad - the radial value that will be interpolated
            !   sec_spl - the spline-interpolated value 
            ! The first index stores the index in the section spline array
            !   so that each sections spline array is continuous
            ! The index index is the section id 
            ! This will then be pieced back together to give the overall interp. 
            allocate(self%orig_id(self%n_radial, self%num_sections))
            allocate(self%sec_rad(self%n_radial, self%num_sections))
            allocate(self%sec_spl(self%n_radial, self%num_sections))
            allocate(self%n_per_sec(self%num_sections))

            ! Signals a value not filled 
            self%orig_id = -1
            self%sec_rad = -one

            ! IDs start at 1 for each section 
            self%n_per_sec(:) = 1

            ! Loop through each radial value to be interpolated to 
            ! and, for each, find the section the radius is within
            do i = 1, self%n_radial
                rr = self%radial(i)
                do idisc = self%lower_disc_id, self%upper_disc_id 
                    ! SECTIONING CONDITION: 
                    ! lies within section if it is above the lower disc radius 
                    ! up to (and including the upper disc radius)
                    ! This should be consistent with how we interpolate the cubic 
                    ! splines (e.g. how the interp map is made)
                    if(rr.gt.self%M%rdisc_norm(idisc) .and. rr.le.self%M%rdisc_norm(idisc+1))then
                        ! Point lies in this section 
                        section_id = idisc - self%lower_disc_id + 1
                        ! index in original array
                        self%orig_id(self%n_per_sec(section_id), section_id) = i 
                        ! radius value to be interpolated
                        self%sec_rad(self%n_per_sec(section_id), section_id) = self%radial(i)
                        ! update counter for this section 
                        self%n_per_sec(section_id) = self%n_per_sec(section_id) +  1
                    elseif (rr.eq.self%M%rdisc_norm(1) .and. idisc.eq.self%lower_disc_id) then 
                        ! special case when rr = centre: 
                        ! caution: will this always work? 
                        ! just give the spline the value from the first element of
                        ! the input array - since this would need to be zero 
                        section_id = 1
                        ! index in original array
                        self%orig_id(self%n_per_sec(section_id), section_id) = i 
                        ! radius value to be interpolated
                        self%sec_rad(self%n_per_sec(section_id), section_id) = self%radial(i)
                        ! update counter for this section 
                        self%n_per_sec(section_id) = self%n_per_sec(section_id) +  1
                    endif 
                    
 
                enddo ! idisc
            enddo ! i 

            ! Subtract 1 from the counter values to give the last index of the
            ! section-wise arrays. This is equivalent to the number of radial
            ! values in each section that need to be interpolated 
            self%n_per_sec = self%n_per_sec - 1
            self%max_per_sec = maxval(self%n_per_sec)
            
            if(safety_checks)then
                do i = 1, self%num_sections
                    if(self%n_per_sec(i).lt.0)then 
                    write(*,*)'Error, some counters not filled. Empty section?'
                    stop 
                    endif 
                enddo 
                ss = SUM(self%n_per_sec)
                if(ss.ne.self%n_radial)then
                    write(*,*)'Error in piecewise interpolation setup'
                    write(*,*)' Sum of counters  ', ss 
                    write(*,*)' Number of radial vals is ', self%n_radial 
                    write(*,*)'should be the same '
                    stop
                endif 
            endif


            ! Set the min and maximum knots of the model: 
            self%min_knot_id = self%M%disc(self%lower_disc_id)
            self%max_knot_id = self%M%disc(self%upper_disc_id+1)


        end subroutine setup



        subroutine create_interpolation_radial_map(self)
            implicit none 

            class(InterpPiecewise) :: self


            ! Local variables
            integer :: i_unq, i_knot, i, isec

            if(verbose.ge.2)write(*,*)'Creating interpolation map'
            allocate(self%interp_map(self%n_radial))

            if(safety_checks)then 
                if(self%min_knot_id.lt.0 .or. self%max_knot_id.lt.1)then 
                    write(*,*)'Error in create_interpolation_radial_map'
                    write(*,*)'Min knot: ', self%min_knot_id
                    write(*,*)'Max knot: ', self%max_knot_id
                    write(*,*)'Stop.'
                    stop 
                endif 
            endif

            ! Find the knot id's of the radii just below the points we will want to interpolate to 
            self%interp_map(:) = 0
            do i_unq = 1, self%n_radial
                ! For each radius we want to find the last knot that it is larger than

                do i_knot = self%min_knot_id, self%max_knot_id

                    if (self%radial(i_unq) .lt. self%M%rad_mineos(i_knot))then 
                        ! This should be the first time that the knot is above the
                        ! radius in question and so we want the i_knot - 1 to be stored
                        self%interp_map(i_unq) = i_knot - 1 

                        exit 
                    else if (self%radial(i_unq) .eq. self%M%rad_mineos(i_knot))then 
                        ! It will equal a knot at the boundaries 
                        ! If at lower boundary (iknot = min_knot_id) then 
                        ! we need to use the min_knot_id value, else can use the 
                        ! lower end (i.e. the point is at the far (upper) end)
                        ! of this piecewise continuous section 
                        ! In this case set it equal to that knot with the assumption the spline interpolation
                        ! will hold at the the boundaries of each knot
                        if (i_knot .eq. self%min_knot_id) then 
                            self%interp_map(i_unq) = i_knot 
                        else
                            self%interp_map(i_unq) = i_knot - 1 
                        endif
                    endif
                enddo 
            enddo 

            ! Check the interp_map 
            if(minval(self%interp_map).lt.self%min_knot_id .or. maxval(self%interp_map).gt.self%max_knot_id)then
                write(*,*)'Error in assigning values to map'
                write(*,*)' -- min value: ', minval(self%interp_map)
                write(*,*)' -- max value: ', maxval(self%interp_map)

                write(*,*)' ------ some useful debugging things -----'
                write(*,*)' -- min value of your radial array  : ', minval(self%radial)
                write(*,*)' -- min value of mineos radial array: ', self%M%rad_mineos(self%min_knot_id)

                write(*,*)' -- max value of your radial array  : ', maxval(self%radial)
                write(*,*)' -- max value of mineos radial array: ', self%M%rad_mineos(self%max_knot_id)
                stop
            endif 


            ! For piecewise interpolation we will want to have an interp map for each section: 
            allocate(self%sec_interp_map(self%max_per_sec, self%num_sections))
            do isec = 1, self%num_sections
                do i = 1, self%n_per_sec(isec)
                    self%sec_interp_map(i, isec) = self%interp_map(self%orig_id(i,isec))
                enddo
            enddo 
        end subroutine create_interpolation_radial_map
        



        subroutine interpolate_mineos_variable(self, in_variable, spl_out)
            ! Subroutine interpolates a mineos variable (such as the 1D density profile)
            ! To a set of new radial points  
            use allocation_module, only: allocate_if_unallocated, deallocate_if_allocated
            use spline, only: cubic_spline_interp
            implicit none

            class(InterpPiecewise) :: self

            ! IO variables
            real(kind=SPLINE_REAL)  :: in_variable(self%M%NR)
            real(kind=SPLINE_REAL)  :: spl_out(self%n_radial)
            
            ! Local variables: 
            integer :: isection, nout_sec, min_knot, max_knot, nin, i
            
            ! Initialise the output spline 
            spl_out = SPLINE_ZERO        



            ! Interpolate each piecewise section 
            do isection = 1, self%num_sections
                nout_sec = self%n_per_sec(isection)                     
                min_knot = self%M%disc(self%lower_disc_id + isection -1)+1  
                max_knot = self%M%disc(self%lower_disc_id + isection )
                nin      = max_knot - min_knot + 1


                select case(nout_sec)
                    case(1)
                        self%sec_spl(1, isection)  = in_variable(min_knot)
                    case(2)
                        self%sec_spl(1, isection)  = in_variable(min_knot)
                        self%sec_spl(2, isection)  = in_variable(max_knot)
                    case default
                    call cubic_spline_interp(nin,                                    &
                                            self%M%rad_mineos(min_knot:max_knot),    &
                                            in_variable(min_knot:max_knot),          & 
                                            nout_sec,                                &
                                            self%sec_rad(1:nout_sec, isection),      &
                                            self%sec_spl(1:nout_sec, isection),      & 
                                            self%sec_interp_map(1:nout_sec, isection), &
                                            min_knot)
                end select 


                ! Stitch it back together into a single splined array: 
                do i = 1, nout_sec
                    spl_out(self%orig_id(i,isection)) = self%sec_spl(i, isection)
                enddo 
            enddo 




        end subroutine interpolate_mineos_variable



        subroutine interpolate_mode_eigenfunctions(self, Mmode)
            ! Wrapper around interpolate_mineos_variable that uses get_mode to load a eigenfunction 
            ! before interpolating it for u, v, du, dv
            use mineos_model, only: mineos
            use allocation_module, only: allocate_if_unallocated, deallocate_if_allocated
            implicit none

            class(InterpPiecewise) :: self
            type(Mode) :: Mmode  

            if(Mmode%t.eq.'T' .or. Mmode%t.eq.'C')then 
                ! toroidal
                allocate(Mmode%w_spl(self%n_radial))
                allocate(Mmode%dw_spl(self%n_radial))
                call self%interpolate_mineos_variable(Mmode%w, Mmode%w_spl)
                call self%interpolate_mineos_variable(Mmode%dw, Mmode%dw_spl)
            else 
                ! spheroidal 
                allocate(Mmode%u_spl(self%n_radial))
                allocate(Mmode%v_spl(self%n_radial))
                allocate(Mmode%du_spl(self%n_radial))
                allocate(Mmode%dv_spl(self%n_radial))
                call self%interpolate_mineos_variable(Mmode%u, Mmode%u_spl)
                call self%interpolate_mineos_variable(Mmode%du, Mmode%du_spl)
                call self%interpolate_mineos_variable(Mmode%v, Mmode%v_spl)
                call self%interpolate_mineos_variable(Mmode%dv, Mmode%dv_spl)
            endif 

            Mmode%spl_len = self%n_radial

        end subroutine interpolate_mode_eigenfunctions


        subroutine print(self)

            class(InterpPiecewise) :: self

            write(*,*)
            write(*,*)'------- Radial array --------'
            write(*,*)'no. radial points : ', self%n_radial
            write(*,*)'min. radius       : ', self%min_r
            write(*,*)'max. radius       : ', self%max_r
            write(*,*)
            write(*,*)'----- Interpolation Map -----'
            write(*,*)'min. knot id      : ', self%min_knot_id
            write(*,*)'max. knot id      : ', self%max_knot_id
            write(*,*)
            write(*,*)'------ Piecewise setup ------'
            write(*,*)'no. sections      : ', self%num_sections
            write(*,*)'lower disc id     : ', self%lower_disc_id
            write(*,*)'upper disc id     : ', self%upper_disc_id
            write(*,*)'max. per section  : ', self%max_per_sec

        end subroutine print


        subroutine load(self, fname)
            implicit none 
            class(InterpPiecewise) :: self
            character(len=*) :: fname
            write(*,*)
        end subroutine load


        subroutine save(self, fname)
            implicit none 
            class(InterpPiecewise) :: self
            character(len=*) :: fname

            write(*,*)
        end subroutine save


        subroutine cleanup(self)
            implicit none 

            class(InterpPiecewise) :: self
            write(*,*)
        end subroutine cleanup


        ! Constructor function for the InterpPiece
        function create_PieceInterp(n) result(P_I)
            implicit none 
            integer, intent(in) :: n
            type(InterpPiecewise)   :: P_I

            P_I%n_radial = n 
            allocate(P_I%radial(n))
            
        end function create_PieceInterp
        

end module 