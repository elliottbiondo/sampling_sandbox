!+ $Id: source.F90,v 1.2 2006/10/03 02:10:16 mashnik Exp $
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

  ! 
  call expirx(0,'source','you need a source subroutine.')
  return
end subroutine source
