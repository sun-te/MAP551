!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program bz_1d_main

  use mod_precision
  use mod_cartesian_grid
  use mod_utils
  use mod_bz_2eq_1d
  use mod_bz_3eq_1d
  use mod_integration

  implicit none

  integer :: neq 
  real(kind=dp) :: tini, tend
  integer :: nt
  real(kind=dp) :: xmin, xmax
  integer :: nxib
  type(cartesian_grid_type) :: grid
  real(kind=dp), allocatable, dimension(:) :: unum
  integer :: ntot = 2
  character(len=20) :: method
  real(kind=dp) :: tol
  real(kind=dp) :: norm_err
  real(kind=dp) :: loc_elapsed, elapsed

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! initialisation

  call read_data(neq, tini, tend, nt, xmin, xmax, nxib, method, tol)

  print *, "Resolution of BZ equation" 
  print '(" Model with ", i1, " species")', neq

  print *, "  neq   =", neq
  print "(a, f12.3)", "   tini = ", tini
  print "(a, f12.3)", "   tend = ", tend
  print *, "  nt   =", nt
  print "(a, f12.3)", "   xmin =", xmin
  print "(a, f12.3)", "   xmax =", xmax
  print *, "  nxib =", nxib
  print *, "  integration method : ", method
  print "(a, es12.3)", "   tolerance (for Radau5 and Rock4) =", tol
 
  call init_cartesian_grid(xmin, xmax, nxib, grid)
  
  ! init model
  if (neq==2) then
    call init_bz_2eq(grid)
  else if (neq==3) then
    call init_bz_3eq(grid)
  else 
    print *, "  The values of neq must be equal to 2 or 3, not ", neq
    print *, "  Program stopped"
    call exit(1)
  end if
 
  ! allocation
  ntot = neq*grid%nx
  allocate(unum(ntot)) 
  
  ! compute and save initial solution 
  if (neq==2) then
    call bz_2eq_init_sol(unum)
  else 
    call bz_3eq_init_sol(unum)
  end if
  call save_sol("sol_ini.dat", neq, grid, tini, unum)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! integration
  call integrate(method, tol, neq, grid%nx, unum, nt, tini, tend)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! print sol and compute error

  ! write numeric solution
  call save_sol("sol_num.dat", neq, grid, tend, unum)
 
  ! compute error
  ! use the sol_ref as the exact solution
  ! compute only for neq=2
  call compute_error(neq, grid, unum)
   
  deallocate(unum) 

end program bz_1d_main
