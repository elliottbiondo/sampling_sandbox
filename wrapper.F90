program psuedo_mcnp
  implicit none
  logical:: analog  = .true.
  double precision :: xxx, yyy, zzz, erg, wgt
  double precision, dimension(6) :: rands
  integer:: i, j
  integer, parameter :: out_unit=20
  call mcnp_sampling_setup(analog)

  open(unit=out_unit,file="samples.out", action="write",status="replace")

  do i=1,5000
    do j=1,6
      rands(j) = RAND()
    end do
      call fparticle_birth(rands, xxx, yyy, zzz, erg, wgt)
      write(out_unit,*) xxx, yyy, zzz, erg, wgt
  end do
end program psuedo_mcnp
