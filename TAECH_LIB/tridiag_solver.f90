module tridiag_solver
  use TAECH_IO
  implicit none
  private
  
  public :: solve_tridiag

contains
  subroutine solve_tridiag(sub_diag, diag, sup_diag, rhs, sol, num_eqns)
    implicit none
    !      sub_diag - sub-diagonal (means it is the diagonal below the main diagonal)
    !      diag - the main diagonal
    !      sup_diag - sup-diagonal (means it is the diagonal above the main diagonal)
    !      rhs - right part
    !      sol - solution
    !      num_eqns - number of equations
    integer, intent(in) :: num_eqns
    complex(kind=8), intent(in) :: sub_diag(1:num_eqns)
    complex(kind=8), intent(in) :: diag(1:num_eqns)
    complex(kind=8), intent(in) :: sup_diag(1:num_eqns)
    complex(kind=8), intent(in) :: rhs(1:num_eqns)
    complex(kind=8), intent(out) :: sol(num_eqns)
    complex(kind=8) :: bp(1:num_eqns), vp(1:num_eqns)
    complex(kind=8) :: m
    integer :: i
    ! Make copies of the b and v variables so that they are unaltered by this sub
    bp(1) = diag(1)
    vp(1) = rhs(1)
    !The first pass (setting coefficients):
    firstpass: do i = 2, num_eqns
       m = sub_diag(i)/bp(i-1)
       bp(i) = diag(i) - m*sup_diag(i-1)
       vp(i) = rhs(i) - m*vp(i-1)
    end do firstpass
    sol(num_eqns) = vp(num_eqns)/bp(num_eqns)
    !The second pass (back-substition)
    backsub: do i = num_eqns-1, 1, -1
       sol(i) = (vp(i) - sup_diag(i)*sol(i+1))/bp(i)
    end do backsub
    return
  end subroutine solve_tridiag
end module tridiag_solver


program test_tridiag_solver
  use tridiag_solver
  implicit none
  complex(kind=8) :: sub_diag(3), diag(3), sup_diag(3)
  complex(kind=8) :: rhs(3)
  complex(kind=8) :: sol(3)
  sub_diag(2:3) = (1.0D0, 0.0D0)
  diag(1:3) = (0.01D0, 0.0001D0)
  sup_diag(1:2) = (1.0D0, 0.0D0)
  rhs(1) = (1.0D0, 0.0D0)
  rhs(2) = (2.0D0, 0.0D0)
  rhs(3) = (3.0D0, 0.0D0)
  call solve_tridiag(sub_diag, diag, sup_diag, rhs, sol, 3)
  write(*,*) sol
  return
end program test_tridiag_solver
