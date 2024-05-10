!------------------------------------------------------------------------------------------------------------
! m_aijpj - compute displacements ui from influence coefficients Aij and tractions pj
!
! Copyright 1993-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------

module m_AijPj

use m_hierarch_data
#ifdef WITH_MKLFFT
   use mkl_dfti
#endif
#if defined _OPENMP
   use omp_lib, only : omp_in_parallel
#endif

! provide for compilation with "internal parallelization" (within each contact patch) enabled or disabled
#if defined OMP_INSIDE
#define _INT_OMP $omp
#else
#define _INT_OMP
#endif

#ifdef _WIN32
   !dir$ attributes c, reference  :: mkl_set_num_threads ! on IA-32: use C calling convention
#endif

implicit none
private

public AijPj
public vecAijPj
public IgsRange
#if defined WITH_MKLFFT
public fft_makePrec
public fft_vecAijPj
public fft_cleanup
#endif

interface AijPj
   module procedure gf3_AijPj
end interface AijPj

interface vecAijPj
   module procedure gf3_VecAijPj
end interface vecAijPj

#ifdef WITH_MKLFFT

   ! internal variables for use in the FFT algorithm

   integer, parameter             :: fft_ndim = 2
   logical                        :: mkl_1st_call = .true.
   logical                        :: desc_present = .false.
   integer                        :: desc_arrsiz(fft_ndim)
   type(dfti_descriptor), pointer :: dfti_desc
   type(dfti_descriptor), pointer :: dfti_back
!$omp threadprivate(mkl_1st_call, desc_present, desc_arrsiz, dfti_desc, dfti_back)
#endif

contains

!------------------------------------------------------------------------------------------------------------

subroutine IgsRange(iigs, igs0, igs1)
!--purpose: Expand the possibly generic selection code "iigs" for elements of the potential contact area
!           into a range of concrete selection codes "igs0:igs1"
!           actual codes: ordered list Exter, Adhes, Slip, Plast
!           AllElm  ==  Exter .. Plast
!           AllExt  ==  Exter .. Exter
!           AllInt  ==  Adhes .. Plast
   implicit none
!--subroutine parameters:
   integer, intent(in)  :: iigs
   integer, intent(out) :: igs0, igs1

   if (iigs.eq.AllElm) then
      igs0 = Exter
      igs1 = Plast
   elseif (iigs.eq.AllExt) then
      igs0 = Exter
      igs1 = Exter
   elseif (iigs.eq.AllInt) then
      igs0 = Adhes
      igs1 = Plast
   elseif (iigs.ge.Exter .and. iigs.le.Plast) then
      igs0 = iigs
      igs1 = iigs
   else
      igs0 = Adhes
      igs1 = Exter
   endif
end subroutine IgsRange

!------------------------------------------------------------------------------------------------------------

function gf3_AijPj (ii, ik, p, jkarg, c)
!--purpose: Compute the displacement in one element ii in direction ik due to the tractions in direction jk,
!           with influence coefficients c. "inner product of one row of matrix A times tractions vector p"
!
!           jk == jkALL = -3 means directions 1, 2, 3, i.e. all tractions,
!           jk == jkTANG = -2 means directions 1 and 2, i.e. the tangential tractions,
   implicit none
!--subroutine parameters :
   type(t_inflcf),   intent(in) :: c
   type(t_gridfnc3), intent(in) :: p
   integer,          intent(in) :: ii, ik, jkarg
!--function result :
   real(kind=8)                 :: gf3_aijpj
!--local variables :
   integer      :: mx, my, ix, iy, jx, jy, jh2, jk0, jk1, jk
   real(kind=8) :: rowsum, partsum
!--functions called :
   intrinsic :: random_number

   gf3_aijpj = 0.0

   ! determine range for traction directions jk
   ! Note: assuming that jkXDIR==1, jkYDIR==2, jkZDIR==3
   ! Note: if .not.nt_cpl, then there is no influence of normal tractions on tangential displacements
   !       and vice versa.

   if (jkarg.eq.-3) then
      jk0 = 1
      jk1 = 3
      if (.not.c%nt_cpl .and. ik.le.2 .and. c%itypcf.eq.0) jk1 = 2
      if (.not.c%nt_cpl .and. ik.eq.3 .and. c%itypcf.eq.0) jk0 = 3
   elseif (jkarg.eq.-2) then
      jk0 = 1
      jk1 = 2
      if (.not.c%nt_cpl .and. ik.eq.3 .and. c%itypcf.eq.0) jk1 = 0
   elseif (jkarg.ge.1 .and. jkarg.le.3) then
      jk0 = jkarg
      jk1 = jkarg
      if (.not.c%nt_cpl .and. ik.le.2 .and. c%itypcf.eq.0) jk1 = min(2,jk1)
      if (.not.c%nt_cpl .and. ik.eq.3 .and. c%itypcf.eq.0) jk0 = 3
   else
      jk0 = 1
      jk1 = 0
   endif

   if (.not.associated(p%grid)) then
      call write_log('gf3_AijPj: ERROR: tractions do not contain grid info.')
      call abort_run()
   endif
   if (.not.associated(p%eldiv)) then
      call write_log('gf3_AijPj: ERROR: tractions do not contain element division.')
      call abort_run()
   endif

   mx = p%grid%nx
   my = p%grid%ny
   ix = p%grid%ix(ii)
   iy = p%grid%iy(ii)

   if (c%itypcf.eq.0 .and. my.ge.2 .and. mx.ge.35) then

      ! Standard case, one matrix of influence coefficients, Parallel

      ! Note: function gf3_AijPj may be called from outside and from within parallel regions. In the latter
      !       case, gf3_AijPj should not use parallelisation itself.

!_INT_OMP parallel if (omp_in_parallel().eq.0)                                                              &
!_INT_OMP         default(none)                                                                             &
!_INT_OMP         shared(ix, iy, ik, jk, jk0, jk1, c, p, gf3_aijpj, mx, my)                                 &
!_INT_OMP         private(jx, jy, jh2, rowsum, partsum)

      partsum = 0d0  

!_INT_OMP do schedule(static,1)
      do jy = 1, my
         do jk = jk0, jk1
            jh2 = (jy - 1) * mx
            rowsum = 0d0
            do jx = max(1,p%eldiv%row1st(jy)-1), min(mx, p%eldiv%rowlst(jy)+1)

               ! element  ii = (ix,iy) is translated to (0,0)

               rowsum = rowsum + c%cf(ix-jx, iy-jy, ik, jk) * p%val(jx+jh2,jk)
            enddo
            partsum = partsum + rowsum
         enddo
      enddo
!_INT_OMP end do nowait

!_INT_OMP atomic    
      gf3_aijpj = gf3_aijpj + partsum
!_INT_OMP end parallel

   elseif (c%itypcf.eq.0) then

      ! Standard case, one matrix of influence coefficients, Sequential

      partsum = 0d0  

      do jy = 1, my
         do jk = jk0, jk1
            jh2 = (jy - 1) * mx
            rowsum = 0d0
            do jx = max(1,p%eldiv%row1st(jy)-1), min(mx, p%eldiv%rowlst(jy)+1)

               ! element  ii = (ix,iy) is translated to (0,0)

               rowsum = rowsum + c%cf(ix-jx, iy-jy, ik, jk) * p%val(jx+jh2,jk)
            enddo
            partsum = partsum + rowsum
         enddo
      enddo
      gf3_aijpj = gf3_aijpj + partsum

   else

      ! Separate matrices of influence coefficients for each y-coordinate

      partsum = 0d0  
      do jy = 1, my
         do jk = jk0, jk1
            jh2 = (jy - 1) * mx
            rowsum = 0d0
            do jx = max(1,p%eldiv%row1st(jy)-1), min(mx, p%eldiv%rowlst(jy)+1)

               ! element  ii = (ix,iy) is translated to (0,0)

               rowsum = rowsum + c%cy(ix-jx,iy-jy, jy, ik,jk) * p%val(jx+jh2,jk)
            enddo
            partsum = partsum + rowsum
         enddo
      enddo
      gf3_aijpj = gf3_aijpj + partsum

   endif

   ! Influence coefficients need be multiplied by 1/G

   gf3_aijpj = gf3_aijpj * c%ga_inv

   ! Add contribution of normal compressible layer when needed
   !   --> material model B=1, contribution of pn on un

   if (c%use_flxz .and. ik.eq.ikZDIR .and. (jkarg.eq.ik .or. jkarg.eq.jkALL)) then
      gf3_aijpj = gf3_aijpj + c%flx_z * p%val(ii,ik)
   endif

   ! Add contribution of tangential third body layer when needed
   !   --> material model M=4, contribution of px on ux, and of py on uy
   !       Note: assuming ikX,Y <= ikZ, jkALL<=jkTANG<0

   if (c%use_3bl .and. ik.le.ikYDIR .and. (jkarg.eq.ik .or. jkarg.le.jkTANG)) then
      gf3_aijpj = gf3_aijpj + c%flx_3bl * p%val(ii,ik)
   endif

end function gf3_AijPj

!------------------------------------------------------------------------------------------------------------

subroutine gf3_VecAijPj (Igs, iigs, U, ikarg, P, jkarg, C)
!--purpose: Compute the displacements U (which is a gf3) in elements "iigs" of Igs in direction(s) "ikarg"
!           due to the tractions P (which is a gf3) in direction(s) "jkarg", with influence coefficients C.
!
!            iigs == indication of the set of elements where U is required:
!                     AllElm, AllExt, AllInt, Exter, Adhes, Slip, Plast.
!
!            ik == -3 means directions 1, 2, 3, i.e. all displacements,
!            ik == -2 means directions 1 and 2, i.e. the tang.displacements
!
!            jk == -3 means directions 1, 2, 3, i.e. due to all tractions,
!            jk == -2 means directions 1 and 2, i.e. due to tangential tractions
   implicit none
!--subroutine parameters:
   type(t_eldiv)    :: igs
   type(t_inflcf)   :: c
   type(t_gridfnc3) :: u, p
   integer, intent(in) :: iigs, ikarg, jkarg
!--local variables:
   integer                :: idebug = 0
   integer                :: ik, ik0, ik1, ii, igs0, igs1
   integer                :: jk, jk0, jk1
   logical                :: usefft, ladd

   ! determine range for element division codes igs

   if (idebug.ge.5) then
      write(bufout,*) 'Starting VecAijPj for ik,jk=',ikarg,jkarg
      call write_log(1, bufout)
   endif
   call IgsRange(iigs, igs0, igs1)

   ! determine range for displacement directions ik

   if (ikarg.eq.-3) then
      ik0 = 1
      ik1 = 3
   elseif (ikarg.eq.-2) then
      ik0 = 1
      ik1 = 2
   elseif (ikarg.ge.1 .and. ikarg.le.3) then
      ik0 = ikarg
      ik1 = ikarg
   else
      ik0 = 1
      ik1 = 0
   endif

   associate(npot => p%grid%ntot, mx => p%grid%nx, my => p%grid%ny)

   if (mx.gt.c%cf_mx .or. my.gt.c%cf_my) then
      write(*,*) 'Internal error: contact grid is larger than grid used for infl.coeff.'
      write(*,*) '                mx,my=',mx,my,', cf_mx,my=',c%cf_mx, c%cf_my
      call abort_run()
   endif

   ! Switch depending on whether the FFT routine should be used or not

#if defined WITH_MKLFFT
   usefft = c%itypcf.eq.0
#else
   usefft = .false.
#endif

   ! determine range for traction directions jk

   if (jkarg.eq.jkALL) then
      jk0 = jkXDIR
      jk1 = jkZDIR
   elseif (jkarg.eq.jkTANG) then
      jk0 = jkXDIR
      jk1 = jkYDIR
   elseif (jkarg.ge.1 .and. jkarg.le.3) then
      jk0 = jkarg
      jk1 = jkarg
   else
      jk0 = 1
      jk1 = 0
   endif
   if (idebug.ge.5) then
      write(bufout,*) 'ik=',ik0,':',ik1,', jk=',jk0,':',jk1
      call write_log(1, bufout)
   endif

   if (usefft) then

      ! calculate displacement difference u for all requested coord.dirs ik

      do ik = ik0, ik1

         ! overwrite previous u at first "jk", then add further contributions

         ladd = .false.
         do jk = jk0, jk1
            if (idebug.ge.8) then
               write(bufout,*) 'starting jk=',jk,' for ik=',ik
               call write_log(1, bufout)
               write(bufout,*) 'nt_cpl=',c%nt_cpl
               call write_log(1, bufout)
            endif
            if (.not.C%nt_cpl .and. (ik*jk.eq.3 .or. ik*jk.eq.6)) then

               ! No influence of normal tractions on tangential displacements and vice versa.

               if (idebug.ge.5) then
                  write(bufout,'(2(a,i1),a)') 'skipping contribution of p(:,',jk,') on u(:,',ik,')'
                  call write_log(1, bufout)
               endif

               ! Initialize result variable u = 0 for coord.dir ik

               if (.not.ladd) call gf3_set(iigs, 0d0, u, ik)
            else

               ! compute c^{ik,jk} * p^{jk} and copy/add to u^{ik} in selected elements

               if (idebug.ge.5) then
                  write(bufout,'(2(a,i1),a)') ' FFT for contribution of p(:,',jk,') on U(:,',ik,')'
                  call write_log(1, bufout)
               endif
               call fft_VecAijPj (igs, ladd, iigs, u, ik, p, jk, c)
               ladd = .true.

               ! write(*,*) 'ik=',ik,', jk=',jk,': result u'
               ! do ii = 1, npot
               !    write(*,*) 'ii=',ii,':',u%val(ii,ik)
               ! enddo
            endif
         enddo
      enddo

      ! Influence coefficients must be multiplied by 1/G to get true deformation

      if (idebug.ge.10) then
         write(bufout,*) 'fft: ga_inv=',c%ga_inv
         call write_log(1, bufout)
      endif
      call gf3_scal(iigs, c%ga_inv, u, ikarg)
 
      ! Add contribution of normal compressible layer when needed
      !   --> material model B=1, contribution of pn on un

      do ik = ik0, ik1
         if (c%use_flxz .and. ik.eq.ikZDIR .and. (jkarg.eq.ik .or. jkarg.eq.jkALL)) then
            do ii = 1, npot
               if (igs%el(ii).ge.igs0 .and. igs%el(ii).le.igs1) then
                  u%val(ii,ik) = u%val(ii,ik) + c%flx_z * p%val(ii,ik)
               endif
            enddo
         endif
      enddo

      ! Add contribution of tangential third body layer when needed
      !   --> material model M=4, contribution of px on ux, and of py on uy
      !       Note: assuming ikX,Y <= ikZ, jkALL<=jkTANG<0

      do ik = ik0, ik1
         if (c%use_3bl .and. ik.le.ikYDIR .and. (jkarg.eq.ik .or. jkarg.le.jkTANG)) then
            do ii = 1, npot
               if (igs%el(ii).ge.igs0 .and. igs%el(ii).le.igs1) then
                  u%val(ii,ik) = u%val(ii,ik) + c%flx_3bl * p%val(ii,ik)
               endif
            enddo
         endif
      enddo

   else

      ! original implementation using function AijPj

!_INT_OMP parallel do if (omp_in_parallel().eq.0)                                                           &
!_INT_OMP         default(none)                                                                             &
!_INT_OMP         shared(c, p, u, igs, npot, jkarg, igs0, igs1, ik0, ik1)                                   &
!_INT_OMP         private(ii, ik)

      ! loop over all elements of potential contact area
      !    check whether calculation is requested or not

      do ii = 1, npot
         if (igs%el(ii).ge.igs0 .and. igs%el(ii).le.igs1) then

            ! calculate displ.diff. u for all requested coord.dirs ik

            do ik = ik0, ik1
               u%val(ii,ik) = gf3_AijPj(ii, ik, p, jkarg, c)
            enddo
         endif
      enddo
!_INT_OMP end parallel do

   endif ! usefft

   end associate
end subroutine gf3_VecAijPj

!------------------------------------------------------------------------------------------------------------

#ifdef WITH_MKLFFT

subroutine fft_makePrec (ik, c, jk, m)
!--purpose: Compute the influence coefficients M that give the pressures needed to attain unit displacement
!           in the central element.
!
!           ik and jk must be a single coordinate direction "1".
   implicit none
!--subroutine parameters:
   type(t_inflcf)              :: c, m
   integer,         intent(in) :: ik, jk
!--local variables:
   integer                :: idebug
   logical, parameter     :: time_fft = .false.
   integer                :: fft_mx, fft_my, ix, ix0, ix1, iy, iy0, iy1
   integer                :: len_cf, ii, iof1, arrsiz(fft_ndim)
   integer                :: ierr
   integer                :: strides_spac(fft_ndim+1), strides_four(fft_ndim+1)
   integer(kind=4)        :: status
   real(kind=8)           :: scale
!dir$attributes align : 16 :: cs, cf, ms, mf
   real(kind=8),    dimension(:), allocatable :: cs, ms
   complex(kind=8), dimension(:), pointer     :: cf, mf

   ! disable use of OpenMP within MKL (fft) routines

   if (mkl_1st_call) then
      call mkl_set_num_threads(1)
      mkl_1st_call = .false.
   endif

   if (time_fft) call timer_start(itimer_ffttot)
   idebug = 0
   if (idebug.ge.1) then
      call write_log(' mkl_fft_makeprec: starting for AllElm ...')
   endif

   ! Note: ik and jk must be equal to each other (single matrix block, diagonal)

   if (ik.ne.jk) then
      write(bufout,*) 'ERROR: Prec allowed for diagonal block only: ik,jk=', ik, jk
      call write_log(1, bufout)
      call abort_run()
   endif
   m%nt_cpl = c%nt_cpl

   ! Get the selection of the input infl.coeff to be used in the Fourier transform
   ! Use full potential contact area (AllElm) 

   ix0    = 1
   ix1    = c%cf_mx
   iy0    = 1
   iy1    = c%cf_my

   ! Compute favourable data sizes for use in the FFT
   ! NOTE: Cannot be used in makePrec (?!!)

   fft_mx = ix1 - ix0 + 1
   fft_my = iy1 - iy0 + 1
   if (.false.) then
      call write_log(' mkl_fft_makeprec: WARNING: using opt_fft_size')
      fft_mx = opt_fft_size(fft_mx)
      fft_my = opt_fft_size(fft_my)
   endif

   if (idebug.ge.2) then
      write(bufout,'(4(a,i5))') ' mkl_fft_makeprec: grid size=',ix1-ix0+1,' x', &
        iy1-iy0+1,', array size=',2*fft_mx,' x',2*fft_my
      call write_log(1, bufout)
   endif

   ! Create/update descriptor for Intel MKL's FFT routines when needed

   arrsiz       = (/ 2*fft_mx, 2*fft_my /)

   if (.not.desc_present .or. any(arrsiz.ne.desc_arrsiz)) then
      if (time_fft) call timer_start(itimer_fftdsc)

      ! Destroy descriptors that were created before

      if (desc_present) then
         status = DftiFreeDescriptor(dfti_desc)
         status = DftiFreeDescriptor(dfti_back)
         if (idebug.ge.5) call write_log(' mkl_fft_makeprec: descriptors freed...')
      endif

      ! Create new forward and backward descriptors for the current grid-sizes

      desc_arrsiz  = arrsiz
      strides_spac = (/ 0, 1, 2*fft_mx /)
      strides_four = (/ 0, 1, fft_mx+1 /)

                                           ! precision, forward_domain, dimension, length
      status = DftiCreateDescriptor(dfti_desc, DFTI_DOUBLE, DFTI_REAL, fft_ndim, arrsiz)
      if (idebug.ge.5) call write_log(' mkl_fft_makeprec: forward descriptor created...')

      ! status = DftiSetValue(dfti_desc, DFTI_NUMBER_OF_TRANSFORMS,   1)
      status = DftiSetValue(dfti_desc, DFTI_PLACEMENT,      DFTI_NOT_INPLACE)
      if (idebug.ge.5) then
         write(bufout,*) 'mkl_fft_makeprec: status(1) =',status
         call write_log(1, bufout)
      endif
      ! status = DftiSetValue(dfti_desc, DFTI_PACKED_FORMAT,  DFTI_CCS_FORMAT)
      ! status = DftiSetValue(dfti_desc, DFTI_REAL_STORAGE,   DFTI_REAL_REAL)
      status = DftiSetValue(dfti_desc, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX)
      if (idebug.ge.5) then
         write(bufout,*) 'mkl_fft_makeprec: status(2) =',status
         call write_log(1, bufout)
      endif
      ! status = DftiSetValue(dfti_desc, DFTI_INPUT_STRIDES,  strides_spac)
      status = DftiSetValue(dfti_desc, DFTI_OUTPUT_STRIDES, strides_four)
      if (idebug.ge.5) then
         write(bufout,*) 'mkl_fft_makeprec: status(3) =',status
         call write_log(1, bufout)
      endif
      ! status = DftiSetValue(dfti_desc, DFTI_ORDERING, DFTI_BACKWARD_SCRAMBLED)
      if (idebug.ge.5) then
         write(bufout,*) 'mkl_fft_makeprec: status(4) =',status
         call write_log(1, bufout)
      endif

      if (idebug.ge.5) call write_log(' mkl_fft_makeprec: committing descriptor...')
      status = DftiCommitDescriptor(dfti_desc)
      if (idebug.ge.5) call write_log(' mkl_fft_makeprec: descriptor committed...')

      ! Create descriptor for the back-transform (with strides on input)

      status = DftiCreateDescriptor(dfti_back, DFTI_DOUBLE, DFTI_REAL, fft_ndim, arrsiz)
      if (idebug.ge.5) call write_log(' mkl_fft_makeprec: backward descriptor created...')
      scale  = 1d0/(4d0*fft_mx*fft_my)
      status = DftiSetValue(dfti_back, DFTI_BACKWARD_SCALE, scale)
      status = DftiSetValue(dfti_back, DFTI_PLACEMENT,      DFTI_NOT_INPLACE)
      status = DftiSetValue(dfti_back, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX)
      status = DftiSetValue(dfti_back, DFTI_INPUT_STRIDES,  strides_four)
      ! status = DftiSetValue(dfti_back, DFTI_ORDERING,       DFTI_BACKWARD_SCRAMBLED)
      ! status = DftiSetValue(dfti_back, DFTI_OUTPUT_STRIDES, strides_spac)
      status = DftiCommitDescriptor(dfti_back)
      if (idebug.ge.5) call write_log(' mkl_fft_makeprec: descriptor committed...')

      desc_present = .true.
      if (time_fft) call timer_stop(itimer_fftdsc)
   endif

   ! Check validity of Fourier transform data of the influence coefficients
   !   - make sure that the fft_cf array is present and of the correct length
   !   - invalidate all sections of fft_cf when the grid sizes have changed

   len_cf = (fft_mx+1) * 2*fft_my
   call reallocate_arr(c%fft_cf, len_cf, 3, 3)

   if (c%fft_mx.ne.fft_mx .or. c%fft_my.ne.fft_my) then
      c%fft_ok = .false.
      c%fft_mx = fft_mx
      c%fft_my = fft_my
   endif

   ! Compute the Fourier transform of the influence coefficients for coordinate
   ! directions (ik,jk) when not yet available.

   if (.not.c%fft_ok(ik,jk)) then

      ! Copy input data C to work array cs. 
      ! Note that c%cf is (-mx:mx-1,-my:my-1), and can both be smaller and larger than cs.

      allocate(cs( 2*fft_mx  * 2*fft_my))
      cs = 0d0
      do iy = -min(fft_my,c%cf_my), min(fft_my,c%cf_my)-1
         iof1 = (iy+fft_my)*2*fft_mx + fft_mx + 1
         do ix = -min(fft_mx,c%cf_mx), min(fft_mx,c%cf_mx)-1
            cs(iof1+ix) = c%cf(ix,iy,ik,jk)
         enddo
      enddo
      if (idebug.ge.2) call write_log(' mkl_fft_makeprec: done copying c...')
      if (idebug.ge.13) write(*,'(4(6f8.2,/))') cs

      ! Perform the transformation of space domain cs to Fourier domain fft_cf

      cf => c%fft_cf(:,ik,jk)
      status = DftiComputeForward(dfti_desc, cs, cf)
      if (idebug.ge.4) then
         write(bufout,*) 'mkl_fft_makeprec: status=',status,', FFT(c)=...'
         call write_log(1, bufout)
      endif
      if (idebug.ge.7) then
         do iy = 1, 2*fft_my
            write(bufout,'(4f8.2,5x,4f8.2)') (real(cf((iy-1)*(fft_mx+1)+ix)), ix=1,fft_mx+1),   &
                                             (imag(cf((iy-1)*(fft_mx+1)+ix)), ix=1,fft_mx+1)
            call write_log(1, bufout)
         enddo
      endif

      ! Destroy work-array, mark availability of column (ik,jk) in fft_cf.

      deallocate(cs)
      c%fft_ok(ik,jk) = .true.
   endif

   ! Allocate 1D work-arrays for spatial and Fourier transformed M

   if (idebug.ge.2) call write_log(' mkl_fft_makeprec: allocating ms')
   allocate(ms( 2*fft_mx  * 2*fft_my), stat=ierr)
   if (ierr.ne.0) then
      write(bufout,*) 'allocation of ms failed with ierr=',ierr
      call write_log(1, bufout)
      call abort_run()
   endif
   if (idebug.ge.2) call write_log(' mkl_fft_makeprec: allocating mf')
   allocate(mf((fft_mx+1) * 2*fft_my), stat=ierr)
   if (ierr.ne.0) then
      write(bufout,*) 'allocation of mf failed with ierr=',ierr
      call write_log(1, bufout)
      call abort_run()
   endif
   if (idebug.ge.10) call write_log(' mkl_fft_makeprec: allocate ok...')

   ! Perform the element-wise division of fft_cf in Fourier domain,
   !    scale the data w.r.t. back-transform

   do ii = 1, (fft_mx+1) * 2*fft_my
      mf(ii) = 1d0 / c%fft_cf(ii,ik,jk)
   enddo

   if (idebug.ge.4) call write_log(' mkl_fft_makeprec: computed mf = 1 ./ cf ...')
   if (idebug.ge.7) then
      do iy=1,2*fft_my
         write(bufout,'(4f8.2,5x,4f8.2)') (real(mf((iy-1)*(fft_mx+1)+ix)), ix=1,fft_mx+1),       &
                                          (imag(mf((iy-1)*(fft_mx+1)+ix)), ix=1,fft_mx+1)
         call write_log(1, bufout)
      enddo
   endif

   ! Perform the back-transform of mf to space domain ms

   status = DftiComputeBackward(dfti_back, mf, ms)
   if (idebug.ge.2) then
      write(bufout,*) 'mkl_fft_makeprec: status=',status,', IFFT(m)=...'
      call write_log(1, bufout)
   endif
   if (idebug.ge.13) write(*,'(4(6f8.2,/))') ms

   ! Copy spatial preconditioner coefficients from work array ms to M

   do iy = -min(fft_my,c%cf_my), min(fft_my,c%cf_my)-1
      iof1 = (iy+fft_my)*2*fft_mx + fft_mx + 1
      do ix = -min(fft_mx,c%cf_mx), min(fft_mx,c%cf_mx)-1
         m%cf(ix,iy,ik,jk) = c%ga**2 * ms(iof1+ix) 
      enddo
   enddo
   deallocate(ms, mf)
   if (idebug.ge.2) call write_log(' mkl_fft_makeprec: done copying m...')

   if (time_fft) call timer_stop(itimer_ffttot)
   if (idebug.ge.1) call write_log(' mkl_fft_makeprec: done, returning...')
end subroutine fft_makePrec 

!------------------------------------------------------------------------------------------------------------

subroutine fft_VecAijPj (Igs, ladd, iigs, U, ik, P, jk, C)
!--purpose: Compute the displacements U (which is a gf3) in elements "iigs" of Igs in direction(s)
!           "ik" due to the tractions P (which is a gf3 and has its own Igs) in direction(s) "jk",
!           with influence coefficients C, using the Fast Fourier Transform.
!
!            ladd == flag that indicates whether the displacements must overwrite (.false.) or
!                    be added to (.true.) the initial contents of U.
!
!            iigs == indication of the set of elements where U is required:
!                     AllElm, AllExt, AllInt, Exter, Adhes, Slip, Plast.
!
!            ik and jk must be a single coordinate direction 1, 2 or 3.
   implicit none
!--subroutine parameters:
   type(t_eldiv)               :: igs
   type(t_inflcf)              :: c
   type(t_gridfnc3)            :: u, p
   logical,         intent(in) :: ladd
   integer,         intent(in) :: iigs, ik, jk
!--local variables:
   integer, parameter     :: idebug = 0
   logical, parameter     :: time_fft = .false.
   integer                :: fft_mx, fft_my, ix0, ix1, iy0, iy1, iarea
   integer                :: ix, iy, ii, igs0, igs1
   integer                :: len_cf, iof1, iof2, arrsiz(fft_ndim)
   integer                :: strides_spac(fft_ndim+1), strides_four(fft_ndim+1)
   integer(kind=4)        :: status
   real(kind=8)           :: scale
   real(kind=8),    dimension(:), allocatable :: ps, cs, us
   complex(kind=8), dimension(:), allocatable :: pf, uf

   ! disable use of OpenMP within MKL (fft) routines

   if (mkl_1st_call) then
      call mkl_set_num_threads(1)
      mkl_1st_call = .false.
   endif

   if (time_fft) call timer_start(itimer_ffttot)
   if (idebug.ge.1) then
      if (iigs.eq.AllElm) write(*,*) 'mkl_fft_vecaijpj: starting for AllElm ...'
      if (iigs.eq.AllInt) write(*,*) 'mkl_fft_vecaijpj: starting for AllInt ...'
   endif

   associate( mx => p%grid%nx, my => p%grid%ny )

   ! Note: ik and jk must be 1, 2 or 3 (single matrix block).
   !       iigs must be AllElm or AllInt.

   if (ik.le.0 .or. ik.gt.3 .or. jk.le.0 .or. jk.gt.3) then
      write(*,*) 'ERROR: FFT allowed for single coord.direction only: ik,jk=', ik, jk
      call abort_run()
   endif
   if (iigs.ne.AllElm .and. iigs.ne.AllInt) then
      write(*,*) 'ERROR: FFT allowed only for AllElm or AllInt (',iigs,')'
      call abort_run()
   endif

   ! Get the selection of the input grid to be used in the Fourier transform
   ! Note: include one column to the left and right of the actual contact area (ixmin-1, ixmax+1)
   !       for products with dp in SteadyGS.

   if (iigs.eq.AllInt) then
      ix0    = max(1, min(igs%ixmin-1, p%eldiv%ixmin-1))
      ix1    = min(mx, max(igs%ixmax+1, p%eldiv%ixmax+1))
      iy0    = min(igs%iymin, p%eldiv%iymin)
      iy1    = max(igs%iymax, p%eldiv%iymax)
      iarea = (ix1-ix0+1)*(iy1-iy0+1)
   else
      iarea = 0
   endif

   ! Use full potential contact area when requested (AllElm) or if the encompassing rectangle is not
   ! much smaller (factor 1/2.5).
   ! Note: Norm+CG is faster with encompassing rectangle, NormGPCG with always using fullbox

   if (iigs.ne.AllInt .or. 1.1*iarea.gt.mx*my .or. .false.) then
      ix0    = 1
      ix1    = mx
      iy0    = 1
      iy1    = my
   endif

   if (ix1.lt.ix0 .or. iy1.lt.iy0) then

      ! Empty grid: set the relevant output data in array U to zero

      if (.not.ladd) call gf3_set(iigs, 0d0, u, ik)
      return
   endif

   ! Compute favourable data sizes for use in the FFT

   fft_mx = ix1 - ix0 + 1
   fft_my = iy1 - iy0 + 1
   if (.true.) then
      fft_mx = opt_fft_size(fft_mx)
      fft_my = opt_fft_size(fft_my)
   endif
   ! do iof1 = 100, 128
   !    iof2 = opt_fft_size( iof1 )
   ! enddo

   if (idebug.ge.2) write(*,'(4(a,i5))') ' mkl_fft_vecaijpj: grid size=',ix1-ix0+1,' x', iy1-iy0+1,     &
        ', array size=',2*fft_mx,' x',2*fft_my

   ! Create/update descriptor for Intel MKL's FFT routines when needed

   arrsiz       = (/ 2*fft_mx, 2*fft_my /)

   if (.not.desc_present .or. any(arrsiz.ne.desc_arrsiz)) then
      if (time_fft) call timer_start(itimer_fftdsc)

      ! Destroy descriptors that were created before

      if (desc_present) then
         status = DftiFreeDescriptor(dfti_desc)
         status = DftiFreeDescriptor(dfti_back)
         if (idebug.ge.5) write(*,*) 'mkl_fft_vecaijpj: descriptors freed...'
      endif

      ! Create new forward and backward descriptors for the current grid-sizes

      desc_arrsiz  = arrsiz
      strides_spac = (/ 0, 1, 2*fft_mx /)
      strides_four = (/ 0, 1, fft_mx+1 /)

                                   ! precision, forward_domain, dimension, length

      status = DftiCreateDescriptor(dfti_desc, DFTI_DOUBLE, DFTI_REAL, fft_ndim, arrsiz)
      if (idebug.ge.5) write(*,*) 'mkl_fft_vecaijpj: forward descriptor created...'

      ! status = DftiSetValue(dfti_desc, DFTI_NUMBER_OF_TRANSFORMS,   1)
      status = DftiSetValue(dfti_desc, DFTI_PLACEMENT,       DFTI_NOT_INPLACE)
      ! status = DftiSetValue(dfti_desc, DFTI_PACKED_FORMAT, DFTI_CCS_FORMAT)
      ! status = DftiSetValue(dfti_desc, DFTI_REAL_STORAGE,  DFTI_REAL_REAL)
      status = DftiSetValue(dfti_desc, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX)
      ! status = DftiSetValue(dfti_desc, DFTI_INPUT_STRIDES, strides_spac)
      status = DftiSetValue(dfti_desc, DFTI_OUTPUT_STRIDES,  strides_four)
      ! status = DftiSetValue(dfti_desc, DFTI_ORDERING,        DFTI_BACKWARD_SCRAMBLED)
      status = DftiCommitDescriptor(dfti_desc)
      if (idebug.ge.5) write(*,*) 'mkl_fft_vecaijpj: descriptor committed...'

      ! Create descriptor for the back-transform (with strides on input)

      status = DftiCreateDescriptor(dfti_back, DFTI_DOUBLE, DFTI_REAL, fft_ndim, arrsiz)
      if (idebug.ge.5) write(*,*) 'mkl_fft_vecaijpj: backward descriptor created...'
      scale  = 1d0/(4d0*fft_mx*fft_my)
      status = DftiSetValue(dfti_back, DFTI_BACKWARD_SCALE, scale)
      status = DftiSetValue(dfti_back, DFTI_PLACEMENT,      DFTI_NOT_INPLACE)
      status = DftiSetValue(dfti_back, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX)
      status = DftiSetValue(dfti_back, DFTI_INPUT_STRIDES,  strides_four)
      ! status = DftiSetValue(dfti_back, DFTI_ORDERING,       DFTI_BACKWARD_SCRAMBLED)
      ! status = DftiSetValue(dfti_back, DFTI_OUTPUT_STRIDES, strides_spac)
      status = DftiCommitDescriptor(dfti_back)
      if (idebug.ge.5)  write(*,*) 'mkl_fft_vecaijpj: descriptor committed...'

      desc_present = .true.
      if (time_fft) call timer_stop(itimer_fftdsc)
   endif

   ! Check validity of Fourier transform data of the influence coefficients
   !   - make sure that the fft_cf array is present and of the correct length
   !   - invalidate all sections of fft_cf when the grid sizes have changed

   len_cf = (fft_mx+1) * 2*fft_my
   call reallocate_arr(c%fft_cf, len_cf, 3, 3)

   if (c%fft_mx.ne.fft_mx .or. c%fft_my.ne.fft_my) then
      c%fft_ok = .false.
      c%fft_mx = fft_mx
      c%fft_my = fft_my
   endif

   ! Compute the Fourier transform of the influence coefficients for coordinate directions (ik,jk)
   ! when not yet available.

   if (.not.c%fft_ok(ik,jk)) then

      ! Copy input data C to work array cs. 
      ! Note that c%cf is (-mx:mx-1,-my:my-1), and can both be smaller and larger than cs.

      allocate(cs( 2*fft_mx  * 2*fft_my))
      cs = 0d0
      do iy = -min(fft_my,my), min(fft_my,my)-1
         iof1 = (iy+fft_my)*2*fft_mx + fft_mx + 1
         do ix = -min(fft_mx,mx), min(fft_mx,mx)-1
            cs(iof1+ix) = c%cf(ix,iy,ik,jk)
         enddo
      enddo
      if (idebug.ge.2) write(*,*) 'mkl_fft_vecaijpj: done copying c...'
      if (idebug.ge.3) write(*,'(4(6f8.2,/))') cs

      ! Perform the transformation of space domain cs to Fourier domain fft_cf

      status = DftiComputeForward(dfti_desc, cs, c%fft_cf(:,ik,jk))
      if (idebug.ge.4) write(*,*) 'mkl_fft_vecaijpj: status=',status,', FFT(c)=...'
      if (idebug.ge.7 .and. fft_mx.eq.8) then
         do iy=1,2*fft_my
            write(*,'(4f8.2,5x,4f8.2)') (real(c%fft_cf((iy-1)*(fft_mx+1)+ix,ik,jk)), ix=1,fft_mx+1),   &
                                        (imag(c%fft_cf((iy-1)*(fft_mx+1)+ix,ik,jk)), ix=1,fft_mx+1)
         enddo
      endif

      ! Destroy work-array, mark availability of column (ik,jk) in fft_cf.

      deallocate(cs)
      c%fft_ok(ik,jk) = .true.
   endif

   ! Allocate 1D work-arrays for spatial and Fourier transformed P and U

   allocate(ps( 2*fft_mx  * 2*fft_my))
   allocate(us( 2*fft_mx  * 2*fft_my))
   allocate(pf((fft_mx+1) * 2*fft_my))
   allocate(uf((fft_mx+1) * 2*fft_my))
   if (idebug.ge.10) write(*,*) 'mkl_fft_vecaijpj: allocate ok...'

   ! Copy input data P to work array, padded with zeros

   ps = 0d0
   do iy = iy0, iy1
      iof1 = (iy-iy0)*2*fft_mx - ix0 + 1
      iof2 = (iy-1)*mx
      do ix = ix0, ix1
         ps(iof1+ix) = p%val(iof2+ix,jk)
      enddo
   enddo
   if (idebug.ge.2) write(*,*) 'mkl_fft_vecaijpj: done copying p...'
   if (idebug.ge.3) write(*,'(4(6f8.2,/))') ps

   ! Perform the transformation of space domain ps to Fourier domain pf

   status = DftiComputeForward(dfti_desc, ps, pf)
   if (idebug.ge.4) write(*,*) 'mkl_fft_vecaijpj: status=',status,', FFT(p)=...'
   if (idebug.ge.7) then
      do iy=1,2*fft_my
         write(*,'(4f8.2,5x,4f8.2)') (real(pf((iy-1)*(fft_mx+1)+ix)), ix=1,fft_mx+1),       &
                                     (imag(pf((iy-1)*(fft_mx+1)+ix)), ix=1,fft_mx+1)
      enddo
   endif

   ! Perform the element-wise multiplication with fft_cf in Fourier domain

   do ii = 1, (fft_mx+1) * 2*fft_my
      uf(ii) = c%fft_cf(ii,ik,jk) * pf(ii)
   enddo

   if (idebug.ge.4) write(*,*) 'mkl_fft_vecaijpj: multiplied uf = cf * pf ...'
   if (idebug.ge.7) then
      do iy=1,2*fft_my
         write(*,'(4f8.2,5x,4f8.2)') (real(uf((iy-1)*(fft_mx+1)+ix)), ix=1,fft_mx+1),       &
                                     (imag(uf((iy-1)*(fft_mx+1)+ix)), ix=1,fft_mx+1)
      enddo
   endif

   ! Perform the back-transform of uf to space domain us

   status = DftiComputeBackward(dfti_back, uf, us)
   if (idebug.ge.2) write(*,*) 'mkl_fft_vecaijpj: status=',status,', IFFT(u)=...'
   if (idebug.ge.3) write(*,'(4(6f8.2,/))') us

   ! Determine range for element division codes igs

   call IgsRange(iigs, igs0, igs1)

   if (.not.ladd) then

      ! Copy the relevant output data from the work array to array U
      ! uf(f_mx+[1:f_mx], f_my+[1:f_my]), no flipping required (?!)

      do iy = iy0, iy1
         iof1 = (fft_my+iy-iy0)*2*fft_mx + fft_mx-ix0+1
         iof2 =    (iy-1)*  mx
         if (idebug.ge.8) write(*,*) 'iy=',iy,': iof2=',iof2,', iof1=',iof1
         do ix = ix0, ix1
            ii = iof2 + ix
            if (idebug.ge.8) write(*,*) '   ix=',ix,': adr2=',ii, ', adr1=', iof1+ix
            if (igs%el(ii).ge.igs0 .and. igs%el(ii).le.igs1) u%val(ii,ik) = us(iof1+ix)
         enddo
      enddo
   else

      ! Add the relevant output data from the work array to array U

      do iy = iy0, iy1
         iof1 = (fft_my+iy-iy0)*2*fft_mx + fft_mx-ix0+1
         iof2 =    (iy-1)*  mx
         do ix = ix0, ix1
            ii = iof2 + ix
            if (igs%el(ii).ge.igs0 .and. igs%el(ii).le.igs1) u%val(ii,ik) = u%val(ii,ik) + us(iof1+ix)
         enddo
      enddo
   endif
   if (idebug.ge.5) write(*,*) 'mkl_fft_vecaijpj: done copying u...'

   ! Cleanup, destroy work-arrays

   deallocate(ps, us, pf, uf)
   end associate

   if (time_fft) call timer_stop(itimer_ffttot)
   if (idebug.ge.1) write(*,*) 'mkl_fft_vecaijpj: done, returning...'
end subroutine fft_VecAijPj 

! end of version using MKL FFT
#endif

!------------------------------------------------------------------------------------------------------------

integer function opt_fft_size( fft_size )
!--purpose: Compute a favourable size for the FFT routines >= fft_size.
   implicit none
!--function parameters:
   integer,         intent(in) :: fft_size
!--local variables:
   integer, parameter :: idebug   =  0
   integer, parameter :: use_fac  =  7
   integer, parameter :: npattern = 15
   integer, parameter :: kpattern(8,npattern) = reshape( (/                    &
                -2,  1,  0,  0,  0,       4,   3,       3,                     &
                -5,  3,  0,  0,  0,      32,  27,       3,                     &
                -1, -1,  1,  0,  0,       6,   5,       5,                     &
                -4,  1,  1,  0,  0,      16,  15,       5,                     &
                 0, -3,  2,  0,  0,      27,  25,       5,                     &
                -3,  0,  0,  1,  0,       8,   7,       7,                     &
                 1, -1, -1,  1,  0,      15,  14,       7,                     &
                -1,  2, -1,  0,  0,      10,   9,       5,                     &
                -2, -1,  0,  0,  1,      12,  11,      11,                     &
                 0,  0,  0,  0,  0,       0,   0,       0,                     &
                 0,  0,  0,  0,  0,       0,   0,       0,                     &
                 0,  0,  0,  0,  0,       0,   0,       0,                     &
                 0,  0,  0,  0,  0,       0,   0,       0,                     &
                 0,  0,  0,  0,  0,       0,   0,       0,                     &
                 0,  0,  0,  0,  0,       0,   0,       0                      &
               /), (/ 8, npattern /) )
               !  -2, -2,  1,  1,  0,      36,  35,       7,                     &
               !  -1,  0, -2,  2,  0,      50,  49,       7,                     &
               !   2, -2, -1,  0,  2,      45,  44,      11,                     &
               !   0,  1, -1, -1,  1,      35,  33,      11,                     &
               !   0, -4,  0,  1,  1,      81,  77,      11,                     &
   integer :: kbase(5) = (/ 2, 3, 5, 7, 11 /)
   integer :: kfactors(5)
   integer :: k2, ipat, ifac, prod_fac
   logical :: lchanged

   ! Start with the next power of 2: opt_size = 2^k2

   k2       = ceiling( log(1d0*fft_size) / log(2d0) )
   prod_fac = 2**k2
   kfactors = (/ k2, 0, 0, 0, 0 /)

   ! Trade factors according to patterns, e.g. replace two factors 2 by
   !    one factor 3 as long as the resulting value is >= the input value

   ! Re-use the list of patterns until no changes occur anymore

   lchanged = .true.
   do while (lchanged)
      lchanged = .false.

      ! Loop over all patterns,
      !    skip patterns that are not filled in
      !    skip patterns with factors higher than maximum "use_fac"

      do ipat = 1, npattern
         if (kpattern(6,ipat).gt.0 .and. kpattern(8,ipat).le.use_fac) then

            ! Check conditions required for trading:
            !  - prod_fac (kfactors) contains enough of the required factors
            !  - prod_fac is sufficiently larger than the input size

            do while (all(kfactors+kpattern(1:5,ipat).ge.0) .and.              &
                       kpattern(7,ipat)*prod_fac.ge.kpattern(6,ipat)*fft_size)

               ! Conditions ok: apply pattern

               kfactors = kfactors + kpattern(1:5,ipat)
               prod_fac = kpattern(7,ipat)*prod_fac / kpattern(6,ipat)
               lchanged = .true.
            enddo
         endif
      enddo
   enddo

   ! Optionally print the resulting size and its factorisation

   if (idebug.ge.1) then
      write(*,'(2(a,i5),$)') 'size=',fft_size,': opt_size=',prod_fac
      lchanged = .false.
      do ifac = 1, 5
         if (kfactors(ifac).ne.0) then
            if (.not.lchanged) then
               write(*,'(2(a,i1),$)') ' = ',kbase(ifac),'^',kfactors(ifac)
            else
               write(*,'(2(a,i1),$)') ' * ',kbase(ifac),'^',kfactors(ifac)
            endif
            lchanged = .true.
         endif
      enddo
      write(*,*)
   endif

   ! Return the resulting size

   opt_fft_size = prod_fac

end function opt_fft_size

!------------------------------------------------------------------------------------------------------------
#ifdef WITH_MKLFFT
subroutine fft_cleanup()
!--purpose: Cleanup, esp. to avoid memory leak report in Valgrind 
   implicit none
!--local variables:
   integer(kind=4)        :: status

   ! Destroy descriptors that were created before

   if (desc_present) then
      status = DftiFreeDescriptor(dfti_desc)
      status = DftiFreeDescriptor(dfti_back)
   endif
end subroutine fft_cleanup
#endif

!------------------------------------------------------------------------------------------------------------

end module m_AijPj
