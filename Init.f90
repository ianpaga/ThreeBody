module Init
use Global
use Func
use Newton
contains

  Subroutine Init_Random_Seed()
    implicit none

    integer :: ised , i , pid
    integer*8 :: t
    integer , allocatable :: sed(:)

    call random_seed( size = ised )
    allocate( sed(ised) )
    call system_clock(t)
    pid = getpid()
    t = ieor(t, int(pid, kind(t)))
    do i = 1, ised
        sed(i) = t!lcg(t)
    end do
    call random_seed( put=sed )
  End Subroutine Init_Random_Seed

  !Function lcg(s) ! Linear congruential generator
    !implicit none

    !integer :: lcg
    !integer(kind=16) :: s

    !if (s == 0) then
       !s = 104729
    !else
       !s = mod(s, 4294967296_i32)
    !end if
    !s = mod(s * 279470273, 4294967291_i32)
    !lcg = int(mod(s, int(huge(0))), kind(0))
  !End Function lcg

  subroutine Init_Parameters(m, q0, a0, x0, omega0, t0, v0, mdot)
    implicit none

    !Randomly generated
    real*8 :: m, q0, xxx
    real*8 :: a0, mdot, x0, x00, omega0, t0, v0
    integer :: iters
    logical :: debug

    !Generate m and q
    call random_number(q0)
    call random_number(xxx)
    m = xxx * 1.44 * (1 + q0) * Ms
    print*, m / Ms, q0

    !Find out the position of L1
    debug = .False.
    x00 = -0.5
    call solve(phi, phi_prime, x00, x0, q0, iters, debug)
    a0 = radius(m * q0 / (1 + q0) / Ms) * Rs / (1 / (1 + q0) + x0)
    mdot = 1e-1 * Ms / spy

    v0 = sqrt(G*m/a0)
    omega0 = sqrt(G*m/a0**3)
    t0 = 1/omega0

  end subroutine Init_Parameters


end module Init
