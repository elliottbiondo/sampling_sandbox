program blah
  implicit none
  logical:: analog  = .false.
  double precision :: xxx, yyy, zzz, erg, wgt
  double precision, dimension(6) :: rands
  rands(1) = 0.332
  rands(2) = 0.7654
  rands(3) = 0.234
  rands(4) = 0.63453
  rands(5) = 0.234234
  rands(6) = 0.1345
  call mcnp_sampling_setup(analog)
  !call gggsampling_setup()
  call fparticle_birth(rands, xxx, yyy, zzz, erg, wgt)
  write(*,*) 'x', xxx, 'y', yyy, 'z', zzz, 'erg', erg, 'w', wgt
end program blah
