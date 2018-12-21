module mod_rock

abstract interface
  function func(n, y)
     integer :: n
     real(kind=kind(1.0d0)), dimension(n) :: y, func
  end function func
end interface
procedure(func), pointer :: pf => null()

contains

subroutine rock4_integration(f, tol, n, y, nt, tini, tend, info)
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  ! arguments
  interface
    function f(n, y)
       integer :: n
       real(kind=kind(1.0d0)), dimension(n) :: y, f
    end function f
  end interface
  real(kind=dp) :: tol
  integer :: n
  real(kind=dp), dimension(n) :: y
  integer :: nt
  real(kind=dp) :: tini, tend
  logical :: info
  ! local variables
  real(dp) :: dt, t
  real(dp) :: atol, rtol
  integer, dimension(12) :: iwork
  real(dp), dimension(8*n) :: work
  integer :: idid
 
  ! initialisation
  pf => f
 
  t = tini
  ! initial step size  
  dt = (tend-tini)/nt
  dt = 1.d-06

  ! initialize iwork: 
  ! iwork(1)=0  ROCK4 attempts to compute the spectral radius internaly
  ! iwork(2)=1  The Jacobian is constant 
  ! iwork(3)=0  Return and solution at tend.
  ! iwork(4)=0  Atol and rtol are scalars.
  iwork(1) = 0
  iwork(2) = 1
  iwork(3) = 0
  iwork(4) = 0

  ! required tolerance
  rtol = tol
  atol = rtol


  ! rock4
  idid = 1
  call rock4(n, t, tend, dt , y, f_rock, atol, rtol, work, iwork, idid)

  !! print statistics
  if (info) then 
    print "(a, es12.3)", "   Tolerance used: ", rtol
    write(6,*) '  Max number of stages used =', iwork(10)
    write(6,*) '  Number of steps =', iwork(6)
    write(6,*) '  Number of accepted steps =', iwork(7)
  end if

end subroutine rock4_integration

subroutine f_rock(n, t, y, dy)
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  ! parameter
  integer n
  real(kind=dp) :: t
  real(kind=dp), dimension(n) :: y, dy
  dy = pf(n, y)
  return
end subroutine f_rock

end module mod_rock
