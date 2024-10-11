module voronoi 

    use params, only: fname_voronoi
    use variableKind, only: i32, r64

    implicit none 
    include "constants.h"

    real(r64), allocatable :: vor_x(:), vor_y(:), vor_z(:), & 
                                           vor_eta1(:), vor_eta2(:)
    real(kind=CUSTOM_REAL) :: vor_A, vor_C, vor_L, vor_N, vor_F
    integer                :: vor_len

    contains

    subroutine load_voronoi_model()
        implicit none 

        integer :: i 

        open(1, file=fname_voronoi, form='formatted')
        ! Read number of voronoi points
        read(1, *)vor_len  
        ! Read anisotropic parameters 
        read(1, *)vor_A
        read(1, *)vor_C
        read(1, *)vor_L
        read(1, *)vor_N
        read(1, *)vor_F

        ! Allocate arrays based on number of voronoi points
        allocate(vor_x(vor_len))
        allocate(vor_y(vor_len))
        allocate(vor_z(vor_len))
        allocate(vor_eta1(vor_len))
        allocate(vor_eta2(vor_len))

        ! Read in the actual values
        do i = 1, vor_len
            read(1, *)vor_x(i), vor_y(i), vor_z(i), vor_eta1(i), vor_eta2(i)
        enddo 

        close(1)

    end subroutine load_voronoi_model



    subroutine project_voroni_to_gll(tree)
        ! Projects eta1 and eta2 values stored at a cloud of Voronoi points
        ! to the GLL mesh using a K-d tree to evaluate the Voronoi cell the
        ! GLL point lies within
        use params, only: nglob, x_glob, y_glob, z_glob, glob_eta1, glob_eta2
        use dArgDynamicArray_Class, only: dArgDynamicArray
        use m_KdTree, only: KdTree, KdTreeSearch

        implicit none 
        
        ! Derived types from coretran 
        type(KdTree)           :: tree
        type(dArgDynamicArray) :: da
        type(KdTreeSearch)     :: search

        ! Local variables: 
        integer      :: i
        integer(i32) :: ind

        ! For each global coordinate we want to map its voroni cell variable value 
        do i = 1, nglob
            da = search%kNearest(tree, vor_x, vor_y, vor_z, & 
                                 xQuery = x_glob(i), &
                                 yQuery = y_glob(i), &
                                 zQuery = z_glob(i), & 
                                 k = 1)
                                 
            ! Index of the closest voronoi point 
            ind =  da%i%values(1)
            glob_eta1(i) = vor_eta1(ind)
            glob_eta2(i) = vor_eta2(ind)
        enddo 

    end subroutine project_voroni_to_gll

end module voronoi