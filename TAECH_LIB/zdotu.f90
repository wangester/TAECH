module zdotu_mod
  use TAECH_IO
  implicit none
  private

  public :: zdotu

contains
  pure function zdotu(N,ZX,INCX,ZY,INCY)
    implicit none
    complex(kind=8) :: zdotu
    integer, intent(in) :: INCX,INCY,N
    complex(kind=8), intent(in) :: ZX(*),ZY(*)
    complex(kind=8) :: ZTEMP
    integer :: I,IX,IY
    ZTEMP = (0.0d0,0.0d0)
    zdotu = (0.0d0,0.0d0)
    IF (N.LE.0) RETURN
    IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
       DO I = 1,N
          ZTEMP = ZTEMP + ZX(I)*ZY(I)
       END DO
    ELSE
       IX = 1
       IY = 1
       IF (INCX.LT.0) IX = (-N+1)*INCX + 1
       IF (INCY.LT.0) IY = (-N+1)*INCY + 1
       DO I = 1,N
          ZTEMP = ZTEMP + ZX(IX)*ZY(IY)
          IX = IX + INCX
          IY = IY + INCY
       END DO
    END IF
    zdotu = ZTEMP
    RETURN
  END function zdotu
end module zdotu_mod


program test_zdotu
  use TAECH_IO
  use zdotu_mod
  implicit none
  complex(kind=8) :: a(3)
  complex(kind=8) :: b(3)
  integer :: i
  complex(kind=8) :: c(3)
  a(1) = (1.0D0, 1.0D0)
  a(2) = (3.0D0, 3.0D0)
  a(3) = (4.0D0, 4.0D0)
  b(1) = (-1.0D0, -1.0D0)
  b(2) = (-2.0D0, -2.0D0)
  b(3) = (-3.0D0, -3.0D0)
  forall(i=1:3) c(i) =  zdotu(i, a(1:i), 1, b(3:4-i:-1), -1)
  write(*,*) c
  return
end program test_zdotu
