program Three_Body

  use Func
  use Global
  use Init
  use Newton
  implicit none
  !real*8, parameter :: RL = a*(1/(1+q) + x)
  !real*8, parameter :: RL1 = a*(q/(1+q) - x)
  !real*8, parameter :: Rc = (1+q)*RL1**4/a**3
  !Randomly generated
  real*8 :: m, q0
  real*8 :: a0, mdot, x0, x00, v0, omega0, t0
  !Other initial parameters
  real*8 :: theta, radius1, radius2
  !model parameters
  real*8 :: e, t, dt, thetac, a, q, x, thetas_p
  !Mass
  real*8 :: ma, md, m3
  !velocity & location
  real*8, dimension(2) :: vd, va, v, rd, ra, r, dra, drd, r_temp, ra_temp, rd_temp, vth
  !angular velocity
  real*8 :: omega, fd, fa
  !angular momentum
  real*8 :: jdorb, jdsp, jaorb, jasp, j, jdorb0, jdsp0, jaorb0, jasp0, jtot
  logical :: stop, debug !For Newton's method
  character(len=5) :: binary_type, file_name
  integer:: i, num_of_case
  common/ejection/ radius1
  common/accretion/ radius2

  call Init_Random_Seed()
  num_of_case = 0

  do
    num_of_case = num_of_case + 1
    write(file_name,'(i5)') num_of_case
    print*, file_name
    !output header
    open(unit=71, file="./data/R_"//trim(adjustl(file_name)))
    open(unit=72, file="./data/J_"//trim(adjustl(file_name)))
    open(unit=73, file="./data/Para_"//trim(adjustl(file_name)))
    open(unit=74, file="./data/Type_"//trim(adjustl(file_name)))

    call Init_Parameters(m, q0, a0, x0, omega0, t0, v0, mdot)

    q = q0
    a = a0
    x = x0

    !Initial condition
    !!Mass & Radius
    md = q / (1 + q)
    ma = 1 / (1 + q)
    m3 = 0
    radius1 = radius(m * md / Ms) * Rs / a0
    radius2 = radius(m * ma / Ms) * Rs / a0
    if (radius1 + radius2 >=1.0) then
      binary_type = 'MG'
      cycle
    end if
    !!Position
    rd(1) = -1 / (1 + q)
    rd(2) = 0
    ra(1) = q / (1 + q)
    ra(2) = 0
    r(1) = x
    r(2) = 0
    !!Spin
    omega = omega0
    fd = 1
    fa = 1
    !!Orbital velocity
    theta = thetas(q, x)
    vth(1) = 15731.09*sqrt(Temp)*cos(theta) / v0
    vth(2) = 15731.09*sqrt(Temp)*sin(theta) / v0
    vd(1) = 0
    vd(2) = rd(1)
    va(1) = 0
    va(2) = ra(1)
    v(1) = 0
    v(2) = fd * x

    jdorb0 = md*cross_product_2d(rd, vd)
    jdsp0  = 0.4*md*radius1**2 * fd
    jaorb0 = ma*cross_product_2d(ra, va)
    jasp0  = 0.4*ma*radius2**2 * fa
    j = 0
    jtot = jdorb0 + jdsp0 + jaorb0 + jasp0 + j

    t = 0d0
    write(71,'(a)') '#Time, R1x, R1y, R2x, R2y, R3x, R3y'
    write(71,'(999E22.12)') t, rd, ra, r
    write(72,'(a)') '#Time, Jdorb, Jdsp, Jaorb, Jasp, J, Jtot'
    write(72,'(999E22.12)') t, jdorb0, jdsp0, jaorb0, jasp0, j, jtot

    m3 = 2 * pi * t0 * mdot / m
    write(73,'(a)') '#Md (Ms), Ma (Ms), Mp (Ms), Period (min), Omega (s^-1)'
    write(73,'(999E22.12)') md * m / Ms, ma * m / Ms, m3 * m / Ms, 2*pi*t0/60, omega0
    call ejection(md, m3, rd, r, vd, v, fd, vth)

    i = 0
    thetac = 0
    dra = r - ra
    thetas_p = atan(v(2)/v(1))
    do
      jdorb = md*cross_product_2d(rd, vd)
      jdsp  = 0.4*md*radius1**2 * fd
      jaorb = ma*cross_product_2d(ra, va)
      jasp  = 0.4*ma*radius2**2 * fa
      j = m3*cross_product_2d(r, v)
      jtot = jdorb + jdsp + jaorb + jasp + j
      !a = a0/(2/sqrt(dot_product(Ra-Rd, Ra-Rd))-dot_product(va-vd, va-vd))
      !e = sqrt(1 - cross_product_2d((Ra-Rd), (va-vd))**2*a0*v0**2/G/M/a)
      !print*, a / a0, e
      if (mod(i, plot_interval) == 0) then
        write(71,'(999E22.12)') t, rd, ra, r
        write(72,'(999E22.12)') t, jdorb, jdsp, jaorb, jasp, j, jtot
      end if
      !print*, '#Time, R1x, R1y, R2x, R2y, R3x, R3y'
      !print*, vd, v, md, m3, fd
      !print*, '#Time, Jdorb, Jdsp, Jaorb, Jasp, J, Jtot'
      !print*, t, jdorb, jdsp, jaorb, jasp, j, jtot

      i = i + 1
      dt = time_interval/(max(sqrt(dot_product(v,v)),1.))
      r_temp = r
      ra_temp = ra
      rd_temp = rd
      call RK8(f, r, v, t, dt, rd_temp, ra_temp, md, ma)
      call RK8(f, rd, vd, t, dt, ra_temp, r_temp, ma, m3)
      call RK8(f, ra, va, t, dt, rd_temp, r_temp, md, m3)
      t = t + dt
      thetac = thetac + asin(cross_product_2d(r,r_temp)/sqrt(dot_product(r,r)*dot_product(r_temp,r_temp)))
      dra = r - ra
      drd = r - rd

      !Stop when direct impact occurs
      if (sqrt(dot_product(dra,dra)) < radius2) then
        binary_type = 'DI' !Direct Impact
        exit
      end if
      if ((i > 10) .and. (sqrt(dot_product(drd,drd)) < radius1)) then
        print*, sqrt(dot_product(drd,drd)), radius1
        binary_type = 'SA' !Self Accretion
        exit
      end if
      if (abs(thetac-thetas_p)>=2*pi) then
        if (dot_product(dra, dra)<1.d0) then
          binary_type = 'DF' !Disk Formation
        else
          binary_type = 'NA'
        end if
        exit
      end if
      !Stop when evolve for too long ('nstep' steps or one complete orbit
      if (i>nstep) then
        binary_type = 'NA'
        exit
      end if
    end do

    print*, binary_type
    if (binary_type == 'DI') call accretion(ma, m3, ra, r, va, v, fa)
    jdorb = md*cross_product_2d(rd, vd)
    jdsp  = 0.4*md*radius1**2 * fd
    jaorb = ma*cross_product_2d(ra, va)
    jasp  = 0.4*ma*radius2**2 * fa
    j = m3*cross_product_2d(r, v)
    jtot = jdorb + jdsp + jaorb + jasp + j
    write(71,'(999E22.12)') t, rd, ra, r
    write(72,'(999E22.12)') t, jdorb, jdsp, jaorb, jasp, j, jtot
    write(74,'(a)') binary_type
    close(unit=71)
    close(unit=72)
    close(unit=73)
    close(unit=74)
    if (num_of_case>=num_of_case_lim) exit
  end do
end program Three_Body
