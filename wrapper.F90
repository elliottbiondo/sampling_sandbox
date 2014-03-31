program blah
  implicit none
  double precision unicorn
  double precision :: xxx, yyy, zzz, erg, wgt
  double precision, dimension(6) :: rands
  rands(1) = 0.11
  rands(2) = 0.22
  rands(3) = 0.33
  rands(4) = 0.44
  rands(5) = 0.55
  rands(6) = 0.66
  !call fsampling_setup("test.h5m", "phtn_src2", "e_bounds_file")
  call gggsampling_setup()
  call fparticle_birth(rands, xxx, yyy, zzz, erg, wgt)
  write(*,*) 'x', xxx
end program blah
