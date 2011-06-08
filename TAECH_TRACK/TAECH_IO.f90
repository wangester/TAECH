module TAECH_IO
  implicit none
  private
  
! math contants
  real(kind=8), parameter, public :: pi = 3.141592653589793238462643383279D0
  complex(kind=8), parameter, public :: imagj = (0.0D0,1.0D0) 

  public :: Deltam, eps, eta, deltap
  public :: t0, tend, dt, tsave
  public :: modes, smax, ds
  public :: v1, v0
  public :: sbnd, cln2
  
  public :: load_cfg
  public :: get_input, save_output

  real(kind=8) :: Deltam, eps, eta, deltap
  integer :: t0, tend, tsave     ! unit in time step
  real(kind=8) :: dt
  integer :: modes
  real(kind=8) :: smax, ds
  real(kind=8) :: v1, v0
  real(kind=8) :: sbnd, cln2

  character(len=9) :: real_pattern = "ES50.30E5"  
  character(len=24) :: complex_pattern = "ES50.30E5, 1X, ES50.30E5"
  character(len=256) :: output_pattern

interface get_input
   module procedure get_input_real_1d
   module procedure get_input_real_2d
   module procedure get_input_complex_1d
   module procedure get_input_complex_2d
end interface

interface save_output
   module procedure save_output_real_1d
   module procedure save_output_real_2d
   module procedure save_output_complex_1d
   module procedure save_output_complex_2d
end interface


contains
  subroutine load_cfg(filename)
    implicit none
    character(len=*) :: filename
    namelist /parameters/ Deltam, eps, eta, deltap, &
         modes, smax, ds, t0, tend, dt, tsave, &
         v1, v0, sbnd, cln2
    logical :: is_file_alive
    inquire(file=trim(filename), exist=is_file_alive)
    if (.not. is_file_alive) then
       write(*,*) trim(filename), " doesn't exist. Check config file in current directory!(filename?)"
       stop
    end if
    open(unit=10, file=trim(filename), action='read')
    read(10, nml=parameters)
    close(unit=10)
    return
  end subroutine load_cfg

    
  function format_pattern(dtype, col)
    implicit none
    character(len=*) :: dtype
    integer :: col
    character(len=256) :: format_pattern
    character(len=100) :: pattern
    character(len=10) :: num_col
    if (dtype=="real") then
       pattern = real_pattern
    else if (dtype=="complex") then
       pattern = complex_pattern
    else
       write(*,*)  "Specify I/O datatype!(real or complex)"
       stop
    end if
    if (col == 1) then
       format_pattern = "("//trim(pattern)//")"
    else if (col > 1) then
       write(num_col, '(I10)') col
       format_pattern = "("//trim(adjustl(num_col))//"("//trim(pattern)//", ','))"
    else
       write(*,*) "Number of I/O data should be greater than 0!(col>0)"
       stop
    end if
    return 
  end function format_pattern

  
  subroutine get_input_real_1d(filename, array1d, row_of_array1d)
    implicit none
    character(len=*) :: filename
    integer :: row_of_array1d
    real(kind=8) :: array1d(row_of_array1d)
    integer :: i
    logical :: is_file_alive
    inquire(file=trim(filename), exist=is_file_alive)
    if (.not. is_file_alive) then
       write(*,*) trim(filename), " doesn't exist! (filename?)"
       stop
    end if
    open(unit=10, file=trim(filename), action='read', status='old')
    output_pattern = format_pattern("real", 1)
    do i = 1, row_of_array1d
       read(10, trim(output_pattern)) array1d(i)
    end do
    close(unit=10)
    return
  end subroutine get_input_real_1d

 
  subroutine get_input_real_2d(filename, array2d, row_of_array2d, col_of_array2d)
    implicit none
    character(len=*) :: filename
    integer :: row_of_array2d, col_of_array2d
    real(kind=8) :: array2d(row_of_array2d, col_of_array2d)
    integer :: i, j
    logical :: is_file_alive
    inquire(file=trim(filename), exist=is_file_alive)
    if (.not. is_file_alive) then
       write(*,*) trim(filename), " doesn't exist! (filename?)"
       stop
    end if
    open(unit=10, file=trim(filename), action='read', status='old')
    output_pattern = format_pattern("real", col_of_array2d)
    do i = 1, row_of_array2d
       read(10, trim(output_pattern)) (array2d(i,j), j=1,col_of_array2d)
    end do
    close(unit=10)
    return
  end subroutine get_input_real_2d


  subroutine get_input_complex_1d(filename, zarray1d, row_of_zarray1d)
    implicit none
    character(len=*) :: filename
    integer :: row_of_zarray1d
    complex(kind=8) :: zarray1d(row_of_zarray1d)
    integer :: i
    logical :: is_file_alive
    inquire(file=trim(filename), exist=is_file_alive)
    if (.not. is_file_alive) then
       write(*,*) trim(filename), " doesn't exist! (filename?)"
       stop
    end if
    open(unit=10, file=trim(filename), action='read', status='old')
    output_pattern = format_pattern("complex", 1)
    do i = 1, row_of_zarray1d
       read(10, trim(output_pattern)) zarray1d(i)
    end do
    close(unit=10)
    return
  end subroutine get_input_complex_1d


  subroutine get_input_complex_2d(filename, zarray2d, row_of_zarray2d, col_of_zarray2d)
    implicit none
    character(len=*) :: filename
    integer :: row_of_zarray2d, col_of_zarray2d
    complex(kind=8) :: zarray2d(row_of_zarray2d, col_of_zarray2d)
    integer :: i, j
    logical :: is_file_alive
    inquire(file=trim(filename), exist=is_file_alive)
    if (.not. is_file_alive) then
       write(*,*) trim(filename), " doesn't exist! (filename?)"
       stop
    end if
    open(unit=10, file=trim(filename), action='read', status='old')
    output_pattern = format_pattern("complex", col_of_zarray2d)
    do i = 1, row_of_zarray2d
       read(10, trim(output_pattern)) (zarray2d(i,j), j=1,col_of_zarray2d)
    end do
    close(unit=10)
    return
  end subroutine get_input_complex_2d


  subroutine save_output_real_1d(filename, array1d, row_of_array1d, append)
    implicit none
    character(len=*) :: filename
    integer :: row_of_array1d
    real(kind=8) :: array1d(row_of_array1d)
    logical :: append
    integer :: i
    logical :: is_file_alive
    if (append) then
       inquire(file=trim(filename), exist=is_file_alive)
       if (.not. is_file_alive) then
          write(*,*) trim(filename), " doesn't exist! (filename?)"
          stop
       end if
       open(unit=20, file=trim(filename), action='write', position='append')
       output_pattern = format_pattern("real", 1)
       do i = 1, row_of_array1d
          write(20, trim(output_pattern)) array1d(i)
       end do
       close(unit=20)
    else
       open(unit=20, file=trim(filename), action='write')
       output_pattern = format_pattern("real", 1)
       do i = 1, row_of_array1d
          write(20, trim(output_pattern)) array1d(i)
       end do
       close(unit=20)
    end if
    return
  end subroutine save_output_real_1d

 
  subroutine save_output_real_2d(filename, array2d, row_of_array2d, col_of_array2d, append)
    implicit none
    character(len=*) :: filename
    integer :: row_of_array2d, col_of_array2d
    real(kind=8) :: array2d(row_of_array2d, col_of_array2d)
    logical :: append
    integer :: i, j
    logical :: is_file_alive
    if (append) then
       inquire(file=trim(filename), exist=is_file_alive)
       if (.not. is_file_alive) then
          write(*,*) trim(filename), " doesn't exist! (filename?)"
          stop
       end if
       open(unit=20, file=trim(filename), action='write', position='append')
       output_pattern = format_pattern("real", col_of_array2d)
       do i = 1, row_of_array2d
          write(20, trim(output_pattern)) (array2d(i,j), j=1,col_of_array2d)
       end do
       close(unit=20)
    else
       open(unit=20, file=trim(filename), action='write')
       output_pattern = format_pattern("real", col_of_array2d)
       do i = 1, row_of_array2d
          write(20, trim(output_pattern)) (array2d(i,j), j=1,col_of_array2d)
       end do
       close(unit=20)
    end if
    return
  end subroutine save_output_real_2d


  subroutine save_output_complex_1d(filename, zarray1d, row_of_zarray1d, append)
    implicit none
    character(len=*) :: filename
    integer :: row_of_zarray1d
    complex(kind=8) :: zarray1d(row_of_zarray1d)
    logical :: append
    integer :: i
    logical :: is_file_alive
    if (append) then
       inquire(file=trim(filename), exist=is_file_alive)
       if (.not. is_file_alive) then
          write(*,*) trim(filename), " doesn't exist! (filename?)"
          stop
       end if
       open(unit=20, file=trim(filename), action='write', position='append')
       output_pattern = format_pattern("complex", 1)
       do i = 1, row_of_zarray1d
          write(20, trim(output_pattern)) zarray1d(i)
       end do
       close(unit=20)
    else
       open(unit=20, file=trim(filename), action='write')
       output_pattern = format_pattern("complex", 1)
       do i = 1, row_of_zarray1d
          write(20, trim(output_pattern)) zarray1d(i)
       end do
       close(unit=20)
    end if
    return
  end subroutine save_output_complex_1d


  subroutine save_output_complex_2d(filename, zarray2d, row_of_zarray2d, col_of_zarray2d, append)
    implicit none
    character(len=*) :: filename
    integer :: row_of_zarray2d, col_of_zarray2d
    complex(kind=8) :: zarray2d(row_of_zarray2d, col_of_zarray2d)
    logical :: append
    integer :: i, j
    logical :: is_file_alive
    if (append) then
       inquire(file=trim(filename), exist=is_file_alive)
       if (.not. is_file_alive) then
          write(*,*) trim(filename), " doesn't exist! (filename?)"
          stop
       end if
       open(unit=20, file=trim(filename), action='write', position='append')
       output_pattern = format_pattern("complex", col_of_zarray2d)
       do i = 1, row_of_zarray2d
          write(20, trim(output_pattern)) (zarray2d(i,j), j=1,col_of_zarray2d)
       end do
       close(unit=20)
    else
       open(unit=20, file=trim(filename), action='write')
       output_pattern = format_pattern("complex", col_of_zarray2d)
       do i = 1, row_of_zarray2d
          write(20, trim(output_pattern)) (zarray2d(i,j), j=1,col_of_zarray2d)
       end do
       close(unit=20)
    end if
    return
  end subroutine save_output_complex_2d
end module TAECH_IO


!program io_test
!  use TAECH_IO
!  implicit none
!  complex(kind=8) :: data(10)
!  call save_output("data.dat", data, 10, .True.)
!end program io_test
