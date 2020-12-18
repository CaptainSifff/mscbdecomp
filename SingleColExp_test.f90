! MIT License
! 
! Copyright (c) 2020-2021 Florian Goth
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights 
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
! 
! The above copyright notice and this permission notice shall be included in 
! all copies or substantial portions of the Software.
! 
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
! OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
! THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
! FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
! DEALINGS IN THE SOFTWARE.

program singlecolexptest
    use MvG_mod
    use Exponentials_mod
    implicit none
    integer :: ndim, i, j, k, n, myl
    integer :: nredges, dn, IERR, incx, seed
    real(kind=kind(0.D0)) :: hop, r
    complex (kind=kind(0.d0)), ALLOCATABLE, DIMENSION(:,:) :: A !< the full matrix A
    complex (kind=kind(0.d0)), ALLOCATABLE, DIMENSION(:,:) :: U, Identity, M1,M2, M3, ref !< A temporary matrix
    complex (kind=kind(0.d0)), ALLOCATABLE, DIMENSION(:) :: vec, lwork, rwork, res, res2 !< the vector that we will test on
    real(kind=kind(0.D0)), allocatable, dimension(:) :: energ

    type(GraphData) :: gd
    type(EulerExp) :: ee
    real(kind=kind(0.D0)) :: dznrm2, zlange
    integer, allocatable, dimension(:) :: seedarr
    complex(kind=kind(0.D0)) :: alpha, beta

! First create some test matrix
!

    ! initialize A with some data
    hop = 0.05
! coresponds to chain with next-nearest neighbour hopping and OBC
 ndim = 6
 allocate(A(ndim, ndim))
 A=0
 do I = 1, ndim-1, 1
 A(I,I+1) = hop
 A(i+1, i) = hop
 if (i+2 <= ndim) then
    A(I,I+2) = hop
    A(i+2, i) = hop
    endif
 enddo

    allocate(U(ndim, ndim), vec(ndim), energ(ndim), M1(ndim, ndim), M2(ndim, ndim), M3(ndim,ndim), Identity(ndim, ndim))
    Allocate(ref(ndim,ndim))
    
! convert to the internal GraphData structure
!
    
    gd = mat2verts(A)
    
! perform the actual color decomposition
!
    
    call MvG_decomp(gd%verts)
    
! Output some useful information
!
    
    ! Determine the number of used colors and the number of edges
    gd%usedcolors = 0
    gd%nredges = 0
    do i = 1, gd%ndim
        gd%deltag = max(gd%deltag, gd%verts(i)%degree)
        do k = 1, gd%verts(i)%degree
            if (gd%verts(i)%nbrs(k) > i) gd%nredges = gd%nredges + 1
            if (gd%verts(i)%nbrs(k) > gd%ndim) then
                write(*,*) "invalid nbr!!!"
                STOP
            endif
            gd%usedcolors = max(gd%usedcolors, gd%verts(i)%cols(k))
        enddo
    enddo
    write (*,*) "Nr edges: ", gd%nredges
    if (gd%usedcolors == gd%deltag) then
        write(*,*) "Maximum Degree", gd%deltag, ". Found", gd%usedcolors," Families -> optimal decomposition"
    else
        write(*,*) "Maximum Degree", gd%deltag, ". Found", gd%usedcolors," Families"
    endif

! create an Exponential from the color information and the weights of the graph
!
    
    ee = createEulerExponentialfromGraphData(gd)

! Now follows some testing and the comparison to straight-forward exponentiation via lapack
!
! Let's retrieve the entries from the first color
U = 0
do i = 1, ee%singleexps(1)%nrofentries
U(ee%singleexps(1)%nodes(i)%x, ee%singleexps(1)%nodes(i)%y) = acosh(ee%singleexps(1)%nodes(i)%c)
U(ee%singleexps(1)%nodes(i)%y, ee%singleexps(1)%nodes(i)%x) = acosh(ee%singleexps(1)%nodes(i)%c)
enddo

Identity=0
do i = 1, ndim
Identity(i,i) = 1
enddo

do i=1,ndim
write (*,*) (dble(U(i,j)), j=1,ndim)
enddo

! Test 1 - check that inversion brings us close to the identity...
call ee%singleexps(1)%lmult(Identity)
call ee%singleexps(1)%lmultinv(Identity)

write (*,*) "========"
do i=1,ndim
write (*,*) (dble(Identity(i,j)), j=1,ndim)
enddo

! Test 2 - the same for rmult:
Identity=0
do i = 1, ndim
Identity(i,i) = 1
enddo
call ee%singleexps(1)%rmult(Identity)
call ee%singleexps(1)%rmultinv(Identity)

write (*,*) "========"
do i=1,ndim
write (*,*) (dble(Identity(i,j)), j=1,ndim)
enddo


! Test 3 - compare to output of lapack routines
Identity=0
do i = 1, ndim
Identity(i,i) = 1
enddo
  dn = 3*ndim
 allocate(lwork(dn), rwork(dn))
   M1 = U
   call zheev('V', 'U', ndim, M1, ndim, energ, lwork, dn, rwork, IERR)

   energ = exp(energ)
   ! apply to Identity
   alpha = 1.D0
   beta = 0.D0

   call ZGEMM('C', 'N', ndim, ndim, ndim, alpha, M1, ndim, Identity, ndim, beta, M2, ndim)
   do i = 1, ndim
   do j = 1, ndim
   M2(i, j) = M2(i, j) * energ(i)
   enddo
   enddo
   call ZGEMM('N', 'N', ndim, ndim, ndim, alpha, M1, ndim, M2, ndim, beta, Identity, ndim)

   write (*,*) "========"
do i=1,ndim
write (*,*) (dble(Identity(i,j)), j=1,ndim)
enddo

! Test 4 - test the Euler type approximations
Identity=0
do i = 1, ndim
Identity(i,i) = 1
enddo

call ee%lmult(Identity)
call ee%lmult_T(Identity)
   write (*,*) "======== lmult ===="
do i=1,ndim
write (*,*) (dble(Identity(i,j)), j=1,ndim)
enddo
! and from the right
Identity=0
do i = 1, ndim
Identity(i,i) = 1
enddo

call ee%rmult(Identity)
call ee%rmult_T(Identity)
   write (*,*) "======== rmult ===="
do i=1,ndim
write (*,*) (dble(Identity(i,j)), j=1,ndim)
enddo
M3 = Identity

! compare with lapack
Identity=0
do i = 1, ndim
Identity(i,i) = 1
enddo
M1 = 2*A

   call zheev('V', 'U', ndim, M1, ndim, energ, lwork, dn, rwork, IERR)

   energ = exp(energ)
   ! apply to Identity
   alpha = 1.D0
   beta = 0.D0

   call ZGEMM('C', 'N', ndim, ndim, ndim, alpha, M1, ndim, Identity, ndim, beta, M2, ndim)
   do i = 1, ndim
   do j = 1, ndim
   M2(i, j) = M2(i, j) * energ(i)
   enddo
   enddo
   call ZGEMM('N', 'N', ndim, ndim, ndim, alpha, M1, ndim, M2, ndim, beta, ref, ndim)

   write (*,*) "======== ref ======"
do i=1,ndim
write (*,*) (dble(ref(i,j)), j=1,ndim)
enddo
M3 = M3 - Ref
write (*,*) zlange('F', ndim, ndim, M3, ndim, lwork)


STOP 2   
   do i = 1, gd%ndim
    call gd%verts(i)%destruct()
   enddo
   call ee%dealloc()
   deallocate(U, vec, energ, M1, M2, M3, gd%verts, gd%elems)
end program singlecolexptest
