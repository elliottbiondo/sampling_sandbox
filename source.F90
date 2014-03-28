! $Id: source.F90,v 1.2 2006/10/03 02:10:16 mashnik Exp $
! Copyright LANS/LANL/DOE - see file COPYRIGHT_INFO

subroutine source
  ! dummy subroutine.  aborts job if source subroutine is missing.
  ! if nsr=0, subroutine source must be furnished by the user.
  ! at entrance, a random set of uuu,vvv,www has been defined.  the
  ! following variables must be defined within the subroutine:
  ! xxx,yyy,zzz,icl,jsu,erg,wgt,tme and possibly ipt,uuu,vvv,www.
  ! subroutine srcdx may also be needed.
  use mcnp_global
  use mcnp_debug

  implicit real(dknd) (a-h,o-z)
  logical, save :: first_run = .true.
  integer, save :: icl_tmp
  integer :: source_i
  integer :: in_out
  real(dknd), dimension(6) :: rands


  if(first_run .eqv. .true.)then
    icl_tmp = 0
    call sampling_setup("source.h5m", "phtn_src", "e_bounds_file", 1)
    source_i = 0
    do while (icl_tmp .eq. 0)
      source_i = source_i +1
      if(idum(3) .eq. ncl(source_i)) then
        icl_tmp = source_i
      endif
    enddo
    first_run = .false.
  endif

100 continue
  rands(1) = rang()
  rands(2) = rang()
  rands(3) = rang()
  rands(4) = rang()
  rands(5) = rang()
  rands(6) = rang()

  call particle_birth(rands, xxx, yyy, zzz, erg, wgt)
  call chkcel(icl_tmp,2,in_out)
  if(in_out .ne. 0)then
    goto 100
  endif

 icl = icl_temp
 tme = 0.0
 ipt = 1
 jsu = 0

 ! xxx = 1.0
 ! yyy = 1.0
 ! zzz = 1.0
 ! erg = 1.0
 ! wgt = 1.0
 ! icl = 3

  return
end subroutine source
