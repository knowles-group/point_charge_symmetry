module point_charge_symmetry
  private
  public :: SymmetryMeasureValue, SymmetryMeasureRefine, SymmetryMeasureDiscoverGroup, SymmetryMeasureOptimiseFrame
  interface
    function SymmetryMeasureValueC(groupname, atoms, coordinates, charges) bind(C, name = 'SymmetryMeasureValue')
      use iso_c_binding
      real(c_double) :: SymmetryMeasureValue
      character(c_char), dimension(*), intent(in) :: groupname
      integer(c_size_t), value, intent(in) :: atoms
      real(c_double), dimension(*), intent(in) :: coordinates
      real(c_double), dimension(*), intent(in) :: charges
    end function SymmetryMeasureValueC

    function SymmetryMeasureOptimiseFrameC(groupname, atoms, coordinates, charges) bind(C, name = 'SymmetryMeasureOptimiseFrame')
      use iso_c_binding
      integer(c_int) :: SymmetryMeasureOptimiseFrame
      character(c_char), dimension(*), intent(in) :: groupname
      integer(c_size_t), value, intent(in) :: atoms
      real(c_double), dimension(*), intent(inout) :: coordinates
      real(c_double), dimension(*), intent(in) :: charges
    end function SymmetryMeasureOptimiseFrameC

    function SymmetryMeasureDiscoverGroupC(threshold, atoms, coordinates, charges) bind(C, name = 'SymmetryMeasureDiscoverGroup')
      use iso_c_binding
      type(c_ptr) :: SymmetryMeasureDiscoverGroupC
      real(c_double), value, intent(in) :: threshold
      integer(c_size_t), value, intent(in) :: atoms
      real(c_double), dimension(*), intent(inout) :: coordinates
      real(c_double), dimension(*), intent(in) :: charges
    end function SymmetryMeasureDiscoverGroupC

    subroutine SymmetryMeasureRefineC(groupname, atoms, coordinates, charges) bind(C, name = 'SymmetryMeasureRefine')
      use iso_c_binding
      character(c_char), dimension(*), intent(in) :: groupname
      integer(c_size_t), value, intent(in) :: atoms
      real(c_double), dimension(*), intent(inout) :: coordinates
      real(c_double), dimension(*), intent(in) :: charges
    end subroutine SymmetryMeasureRefineC

  end interface
contains
  function SymmetryMeasureValue(groupname, coordinates, charges)
    use iso_c_binding, only : c_size_t, c_char
    double precision :: SymmetryMeasureValue
    character(len = *), intent(in) :: groupname
    double precision, dimension(:, :), intent(in) :: coordinates
    double precision, dimension(:), intent(in) :: charges
    SymmetryMeasureValue = SymmetryMeasureValueC(c_string_c(trim(groupname)), &
        int(ubound(coordinates, 2) - lbound(coordinates, 2) + 1, c_size_t), coordinates, charges)
  end function SymmetryMeasureValue

  function SymmetryMeasureOptimiseFrame(groupname, coordinates, charges)
    use iso_c_binding, only : c_size_t, c_char
    integer :: SymmetryMeasureOptimiseFrame
    character(len = *), intent(in) :: groupname
    double precision, dimension(:, :), intent(inout) :: coordinates
    double precision, dimension(:), intent(in) :: charges
    SymmetryMeasureOptimiseFrame = SymmetryMeasureOptimiseFrameC(c_string_c(trim(groupname)), &
        int(ubound(coordinates, 2) - lbound(coordinates, 2) + 1, c_size_t), coordinates, charges)
  end function SymmetryMeasureOptimiseFrame

  function SymmetryMeasureDiscoverGroup(threshold, coordinates, charges)
    use iso_c_binding, only : c_size_t, c_ptr, c_char, c_f_pointer, c_null_char
    CHARACTER(len = :), ALLOCATABLE :: SymmetryMeasureDiscoverGroup
    double precision, intent(in) :: threshold
    double precision, dimension(:, :), intent(inout) :: coordinates
    double precision, dimension(:), intent(in) :: charges
    character(c_char), dimension(:), pointer :: sp
    type(c_ptr) :: result
    integer :: i
    result = SymmetryMeasureDiscoverGroupC(threshold, &
        int(ubound(coordinates, 2) - lbound(coordinates, 2) + 1, c_size_t), coordinates, charges)
    call c_f_pointer(result, sp, [1])
    do i = 0, 100000000
      if (sp(i + 1).eq.c_null_char) goto 1
    end do
    stop '!'
    1 allocate(character(len = i) :: SymmetryMeasureDiscoverGroup)
    do i = 1, len(SymmetryMeasureDiscoverGroup)
      SymmetryMeasureDiscoverGroup(i:i) = sp(i)
    end do
  end function SymmetryMeasureDiscoverGroup

  subroutine SymmetryMeasureRefine(groupname, coordinates, charges)
    use iso_c_binding, only : c_size_t, c_char
    character(len = *), intent(in) :: groupname
    double precision, dimension(:, :), intent(inout) :: coordinates
    double precision, dimension(:), intent(in) :: charges
    call SymmetryMeasureRefineC(c_string_c(trim(groupname)), &
        int(ubound(coordinates, 2) - lbound(coordinates, 2) + 1, c_size_t), coordinates, charges)
  end subroutine SymmetryMeasureRefine

  !> @brief Convenience function to convert to C string from Fortran string

  FUNCTION c_string_c(fstring)
    use iso_c_binding, only : c_char, c_null_char
    CHARACTER(*), INTENT(in) :: fstring
    CHARACTER(kind = c_char), DIMENSION(:), ALLOCATABLE :: c_string_c
    INTEGER :: i
    ALLOCATE(CHARACTER(kind = c_char) :: c_string_c(len_TRIM(fstring) + 1))
    DO i = 1, len_TRIM(fstring)
      c_string_c(i) = fstring(i:i)
    END DO
    c_string_c(len_TRIM(fstring) + 1) = c_null_char
  END FUNCTION c_string_c

end module point_charge_symmetry