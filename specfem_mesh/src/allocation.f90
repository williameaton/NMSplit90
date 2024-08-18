module allocation_module
    implicit none
    include "constants.h"

    interface allocate_if_unallocated
        module procedure allocate_integer_array_1d
        module procedure allocate_real_array_1d
        module procedure allocate_double_array_1d
        module procedure allocate_complex_array_1d

        module procedure allocate_integer_array_2d
        module procedure allocate_real_array_2d
        module procedure allocate_double_array_2d
        module procedure allocate_complex_array_2d

        module procedure allocate_integer_array_3d
        module procedure allocate_real_array_3d
        module procedure allocate_double_array_3d
        module procedure allocate_complex_array_3d

        module procedure allocate_integer_array_4d
        module procedure allocate_real_array_4d
        module procedure allocate_double_array_4d
        module procedure allocate_complex_array_4d

        module procedure allocate_integer_array_5d
        module procedure allocate_real_array_5d
        module procedure allocate_double_array_5d
        module procedure allocate_complex_array_5d
    end interface

    interface deallocate_if_allocated
        module procedure deallocate_integer_array_1d
        module procedure deallocate_real_array_1d
        module procedure deallocate_double_array_1d
        module procedure deallocate_complex_array_1d

        module procedure deallocate_integer_array_2d
        module procedure deallocate_real_array_2d
        module procedure deallocate_double_array_2d
        module procedure deallocate_complex_array_2d

        module procedure deallocate_integer_array_3d
        module procedure deallocate_real_array_3d
        module procedure deallocate_double_array_3d
        module procedure deallocate_complex_array_3d

        module procedure deallocate_integer_array_4d
        module procedure deallocate_real_array_4d
        module procedure deallocate_double_array_4d
        module procedure deallocate_complex_array_4d

        module procedure deallocate_integer_array_5d
        module procedure deallocate_real_array_5d
        module procedure deallocate_double_array_5d
        module procedure deallocate_complex_array_5d
    end interface



  contains
  
    ! ......................... 1D allocation .........................

    subroutine allocate_integer_array_1d(n, array)
      integer, intent(in) :: n
      integer, allocatable, intent(inout) :: array(:)
      if (.not. allocated(array)) then
        allocate(array(n))
      endif
    end subroutine allocate_integer_array_1d
  
    subroutine allocate_real_array_1d(n, array)
      integer, intent(in) :: n
      real(kind=CUSTOM_REAL), allocatable, intent(inout) :: array(:)
      if (.not. allocated(array)) then
        allocate(array(n))
      endif
    end subroutine allocate_real_array_1d
  
    subroutine allocate_double_array_1d(n, array)
      integer, intent(in) :: n
      double precision, allocatable, intent(inout) :: array(:)
      if (.not. allocated(array)) then
        allocate(array(n))
      endif
    end subroutine allocate_double_array_1d

    subroutine allocate_complex_array_1d(n, array)
        integer, intent(in) :: n
        complex(kind=CUSTOM_REAL), allocatable, intent(inout) :: array(:)
        if (.not. allocated(array)) then
          allocate(array(n))
        endif
      end subroutine allocate_complex_array_1d

    ! ......................... 2D allocation .........................

    subroutine allocate_integer_array_2d(m, n, array)
      integer, intent(in) :: m,n
      integer, allocatable, intent(inout) :: array(:,:)
      if (.not. allocated(array)) then
        allocate(array(m,n))
      endif
    end subroutine allocate_integer_array_2d
  
    subroutine allocate_real_array_2d(m, n, array)
      integer, intent(in) :: m,n
      real(kind=CUSTOM_REAL), allocatable, intent(inout) :: array(:,:)
      if (.not. allocated(array)) then
        allocate(array(m,n))
      endif
    end subroutine allocate_real_array_2d

    subroutine allocate_double_array_2d(m, n, array)
      integer, intent(in) :: m,n
      double precision, allocatable, intent(inout) :: array(:,:)
      if (.not. allocated(array)) then
        allocate(array(m,n))
      endif
    end subroutine allocate_double_array_2d
  
    subroutine allocate_complex_array_2d(m, n, array)
        integer, intent(in) :: m,n
        complex(kind=CUSTOM_REAL), allocatable, intent(inout) :: array(:,:)
        if (.not. allocated(array)) then
          allocate(array(m,n))
        endif
      end subroutine allocate_complex_array_2d

    ! ......................... 3D allocation .........................

      subroutine allocate_integer_array_3d(l, m, n, array)
        integer, intent(in) :: l,m,n
        integer, allocatable, intent(inout) :: array(:,:,:)
        if (.not. allocated(array)) then
          allocate(array(l,m,n))
        endif
      end subroutine allocate_integer_array_3d
    
      subroutine allocate_real_array_3d(l, m, n, array)
        integer, intent(in) :: l,m,n
        real(kind=CUSTOM_REAL), allocatable, intent(inout) :: array(:,:,:)
        if (.not. allocated(array)) then
          allocate(array(l,m,n))
        endif
      end subroutine allocate_real_array_3d
  
      subroutine allocate_double_array_3d(l, m, n, array)
        integer, intent(in) :: l,m,n
        double precision, allocatable, intent(inout) :: array(:,:,:)
        if (.not. allocated(array)) then
          allocate(array(l,m,n))
        endif
      end subroutine allocate_double_array_3d
    
      subroutine allocate_complex_array_3d(l, m, n, array)
          integer, intent(in) :: l,m,n
          complex(kind=CUSTOM_REAL), allocatable, intent(inout) :: array(:,:,:)
          if (.not. allocated(array)) then
            allocate(array(l,m,n))
          endif
        end subroutine allocate_complex_array_3d

    ! ......................... 4D allocation .........................

    subroutine allocate_integer_array_4d(k,l,m, n, array)
      integer, intent(in) :: k,l,m,n
      integer, allocatable, intent(inout) :: array(:,:,:,:)
      if (.not. allocated(array)) then
        allocate(array(k,l,m,n))
      endif
    end subroutine allocate_integer_array_4d
  
    subroutine allocate_real_array_4d(k,l,m, n, array)
      integer, intent(in) :: k,l,m, n
      real(kind=CUSTOM_REAL), allocatable, intent(inout) :: array(:,:,:,:)
      if (.not. allocated(array)) then
        allocate(array(k,l,m,n))
      endif
    end subroutine allocate_real_array_4d

    subroutine allocate_double_array_4d(k,l,m, n, array)
      integer, intent(in) :: k,l,m, n
      double precision, allocatable, intent(inout) :: array(:,:,:,:)
      if (.not. allocated(array)) then
        allocate(array(k,l,m,n))
      endif
    end subroutine allocate_double_array_4d
  
    subroutine allocate_complex_array_4d(k,l,m, n, array)
        integer, intent(in) :: k,l,m, n
        complex(kind=CUSTOM_REAL), allocatable, intent(inout) :: array(:,:,:,:)
        if (.not. allocated(array)) then
          allocate(array(k,l,m,n))
        endif
    end subroutine allocate_complex_array_4d 

    ! ......................... 5D allocation .........................

    subroutine allocate_integer_array_5d(j,k,l,m,n, array)
      integer, intent(in) :: j,k,l,m,n
      integer, allocatable, intent(inout) :: array(:,:,:,:,:)
      if (.not. allocated(array)) then
        allocate(array(j,k,l,m,n))
      endif
    end subroutine allocate_integer_array_5d
  
    subroutine allocate_real_array_5d(j,k,l,m,n, array)
      integer, intent(in) :: j,k,l,m, n
      real(kind=CUSTOM_REAL), allocatable, intent(inout) :: array(:,:,:,:,:)
      if (.not. allocated(array)) then
        allocate(array(j,k,l,m,n))
      endif
    end subroutine allocate_real_array_5d

    subroutine allocate_double_array_5d(j,k,l,m,n, array)
      integer, intent(in) :: j,k,l,m, n
      double precision, allocatable, intent(inout) :: array(:,:,:,:,:)
      if (.not. allocated(array)) then
        allocate(array(j,k,l,m,n))
      endif
    end subroutine allocate_double_array_5d
  
    subroutine allocate_complex_array_5d(j,k,l,m,n, array)
        integer, intent(in) :: j,k,l,m, n
        complex(kind=CUSTOM_REAL), allocatable, intent(inout) :: array(:,:,:,:,:)
        if (.not. allocated(array)) then
          allocate(array(j,k,l,m,n))
        endif
    end subroutine allocate_complex_array_5d
  
    ! ........................ 1D deallocation ........................

    subroutine deallocate_integer_array_1d(array)
      integer, allocatable, intent(inout) :: array(:)
      if (allocated(array)) then
        deallocate(array)
      endif
    end subroutine deallocate_integer_array_1d

    subroutine deallocate_real_array_1d(array)
      real(kind=CUSTOM_REAL), allocatable, intent(inout) :: array(:)
      if (allocated(array)) then
        deallocate(array)
      endif
    end subroutine deallocate_real_array_1d

    subroutine deallocate_double_array_1d(array)
      double precision, allocatable, intent(inout) :: array(:)
      if (allocated(array)) then
        deallocate(array)
      endif
    end subroutine deallocate_double_array_1d

    subroutine deallocate_complex_array_1d(array)
      complex(kind=CUSTOM_REAL), allocatable, intent(inout) :: array(:)
      if (allocated(array)) then
        deallocate(array)
      endif
    end subroutine deallocate_complex_array_1d

    ! ........................ 2D deallocation ........................

    subroutine deallocate_integer_array_2d(array)
      integer, allocatable, intent(inout) :: array(:,:)
      if (allocated(array)) then
        deallocate(array)
      endif
    end subroutine deallocate_integer_array_2d

    subroutine deallocate_real_array_2d(array)
      real(kind=CUSTOM_REAL), allocatable, intent(inout) :: array(:,:)
      if (allocated(array)) then
        deallocate(array)
      endif
    end subroutine deallocate_real_array_2d

    subroutine deallocate_double_array_2d(array)
      double precision, allocatable, intent(inout) :: array(:,:)
      if (allocated(array)) then
        deallocate(array)
      endif
    end subroutine deallocate_double_array_2d


    subroutine deallocate_complex_array_2d(array)
      complex(kind=CUSTOM_REAL), allocatable, intent(inout) :: array(:,:)
      if (allocated(array)) then
        deallocate(array)
      endif
    end subroutine deallocate_complex_array_2d

    ! ........................ 3D deallocation ........................

    subroutine deallocate_integer_array_3d(array)
      integer, allocatable, intent(inout) :: array(:,:,:)
      if (allocated(array)) then
        deallocate(array)
      endif
    end subroutine deallocate_integer_array_3d

    subroutine deallocate_real_array_3d(array)
      real(kind=CUSTOM_REAL), allocatable, intent(inout) :: array(:,:,:)
      if (allocated(array)) then
        deallocate(array)
      endif
    end subroutine deallocate_real_array_3d

    subroutine deallocate_double_array_3d(array)
      double precision, allocatable, intent(inout) :: array(:,:,:)
      if (allocated(array)) then
        deallocate(array)
      endif
    end subroutine deallocate_double_array_3d


    subroutine deallocate_complex_array_3d(array)
      complex(kind=CUSTOM_REAL), allocatable, intent(inout) :: array(:,:,:)
      if (allocated(array)) then
        deallocate(array)
      endif
    end subroutine deallocate_complex_array_3d

    ! ........................ 4D deallocation ........................

    subroutine deallocate_integer_array_4d(array)
      integer, allocatable, intent(inout) :: array(:,:,:,:)
      if (allocated(array)) then
        deallocate(array)
      endif
    end subroutine deallocate_integer_array_4d

    subroutine deallocate_real_array_4d(array)
      real(kind=CUSTOM_REAL), allocatable, intent(inout) :: array(:,:,:,:)
      if (allocated(array)) then
        deallocate(array)
      endif
    end subroutine deallocate_real_array_4d

    subroutine deallocate_double_array_4d(array)
      double precision, allocatable, intent(inout) :: array(:,:,:,:)
      if (allocated(array)) then
        deallocate(array)
      endif
    end subroutine deallocate_double_array_4d

    subroutine deallocate_complex_array_4d(array)
      complex(kind=CUSTOM_REAL), allocatable, intent(inout) :: array(:,:,:,:)
      if (allocated(array)) then
        deallocate(array)
      endif
    end subroutine deallocate_complex_array_4d

    ! ........................ 5D deallocation ........................

    subroutine deallocate_integer_array_5d(array)
      integer, allocatable, intent(inout) :: array(:,:,:,:,:)
      if (allocated(array)) then
        deallocate(array)
      endif
    end subroutine deallocate_integer_array_5d

    subroutine deallocate_real_array_5d(array)
      real(kind=CUSTOM_REAL), allocatable, intent(inout) :: array(:,:,:,:,:)
      if (allocated(array)) then
        deallocate(array)
      endif
    end subroutine deallocate_real_array_5d

    subroutine deallocate_double_array_5d(array)
      double precision, allocatable, intent(inout) :: array(:,:,:,:,:)
      if (allocated(array)) then
        deallocate(array)
      endif
    end subroutine deallocate_double_array_5d

    subroutine deallocate_complex_array_5d(array)
      complex(kind=CUSTOM_REAL), allocatable, intent(inout) :: array(:,:,:,:,:)
      if (allocated(array)) then
        deallocate(array)
      endif
    end subroutine deallocate_complex_array_5d


  end module allocation_module
  

