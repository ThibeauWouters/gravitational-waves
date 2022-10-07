module mod_cfc_parameters
  use mod_global_parameters, only: bigdouble, biginteger
  implicit none
  public 

  logical                                   :: use_cfc = .False.
  
  !-------------------------------------------------------------------!
  !          solver related parameters
  !-------------------------------------------------------------------!

  !> tolerence for 1: alp, 2: psi, 3: beta/X
  double precision, dimension(1:3)          ::  cfc_tol = (/ 1.0d-6, 1.0d-6, 1.0d-6 /)
  !> maximun iteration for each metric solver
  integer                                   ::  cfc_it_max = 1000
  !> Print out the status of the solver at every cfc_print cfc iteration
  integer                                   ::  cfc_print = 10000000
  !> N smooths will be applied for 1: cycling up, 2: cycling down
  integer, dimension(1:2)                   ::  cfc_n_cycle = (/2,2/)
  logical                                   ::  cfc_redblack = .True.
  !> use previous solution as the initial guess
  logical                                   ::  cfc_initial_guess = .True.

  !> outer boundary of beta, either robin (default) or dirichlet
  logical                                   ::  cfc_beta_robin_bc = .True.
  
  !-------------------------------------------------------------------!
  !          initialisation related parameters
  !-------------------------------------------------------------------!

  !> tolerence for initializing psi
  double precision                          ::  cfc_psi_tol_init = 1.0d-8
  !> use the input solution as the initial guess
  logical                                   ::  cfc_initial_guess_init = .False.

  !-------------------------------------------------------------------!
  !          update related parameters
  !-------------------------------------------------------------------!

  !> evolve metric or not
  logical                                   ::  cfc_evolve = .True.
  !> last metric update time
  double precision                          ::  cfc_t_last_update = 0.0d0
  integer                                   ::  cfc_it_last_update = 0
  !> solve the metric at every N steps
  integer                                   ::  cfc_dit_update = biginteger
  !> allowed smallest time step for metric solver (for dit_update only)
  double precision                          ::  cfc_smallest_dt = 0.0d0
  !> solve the metric at every dt
  double precision                          ::  cfc_dt_update = bigdouble
  !> Maximum tolerence for psi
  double precision                          ::  cfc_psi_tol_max = 1.0d-3

  !-------------------------------------------------------------------!
  !          interpolation related parameters
  !-------------------------------------------------------------------!
  !> number of points used to interpolate metric at the cell interface
  integer                                   ::  cfc_n_interpolation = 4
  !> Instead of using interpolation, metric variables can also be reconstructed by limiters, currently only wenozp5 is used.
  logical                                   ::  reconstruct_cfc = .False.

  !-------------------------------------------------------------------!
  !          debug related parameters
  !-------------------------------------------------------------------!
  !> Debug flag, test the performance of the metric solver 
  logical                                   ::  cfc_test_alp = .False.

end module mod_cfc_parameters
