! MIT License
! 
! Copyright (c) 2018-2020 Florian Goth
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

program mscbdecomp
    use MvG_mod
    use Exponentials_mod
    use graphdata_mod
    implicit none
    integer :: ndim, i, j, k, n, myl
    integer :: nredges, dn, IERR, incx, seed
    real(kind=kind(0.D0)) :: hop, r
    complex (kind=kind(0.d0)), ALLOCATABLE, DIMENSION(:,:) :: A !< the full matrix A
    complex (kind=kind(0.d0)), ALLOCATABLE, DIMENSION(:,:) :: U, M1,M2, M3 !< A temporary matrix
    complex (kind=kind(0.d0)), ALLOCATABLE, DIMENSION(:) :: vec, lwork, rwork, res, res2 !< the vector that we will test on
    real(kind=kind(0.D0)), allocatable, dimension(:) :: energ

    type(GraphData) :: gd
    type(EulerExp) :: fe
    real(kind=kind(0.D0)) :: dznrm2, zlange
    integer, allocatable, dimension(:) :: seedarr
    complex(kind=kind(0.D0)) :: alpha, beta

! First create some test matrix
!

    ! initialize A with some data
    hop = 0.1
! coresponds to chain with next-nearest neighbour hopping and OBC
!  ndim = 50
!  allocate(A(ndim, ndim))
!  do I = 1, ndim-1, 1
!  A(I,I+1) = hop
!  A(i+1, i) = hop
! ! A(I,I+2) = hop
! ! A(i+2, i) = hop
!  enddo
! A(1,10) = hop
! A(10,1) = hop
! A(3,8) = hop
! A(8, 3) = hop
! A(9,2) = hop
! A(2, 9) = hop
! A(7,4) = hop
! A(4, 7) = hop
! A(1,5) = hop
! A(5,1) = hop
! A(10,6) = hop
! A(6,10) = hop

    myl = 50
ndim = myl*myl
allocate(A(ndim, ndim))
A = 0
do i = 1, ndim
if((mod(i + 1,myl) .ne. 0) .and. (i + 1 < ndim) ) then
A(i,i+1) = 1
A(i+1,i) = 1
endif
if((mod(i - 1,myl) .ne. 0) .and. (i-1 > 0) ) then
A(i,i-1) = 1
A(i-1,i) = 1
endif
if((i + myl < ndim) ) then
A(i,i + myl) = 1
A(i + myl,i) = 1
endif
if((i - myl > 0) ) then
A(i,i - myl) = 1
A(i - myl,i) = 1
endif
enddo


!    ndim = 100
!     !ndim=7
!     call random_seed(size = n)
!     allocate(A(ndim, ndim), seedarr(n))
! !    do seed = 1002,1010
! !    seed = 0
!     !seed = 1061
!     seed = 99
! !    write (*,*) "seed", seed
!     seedarr = seed + 37 * (/ (i-1, i=1, n) /)
!     call random_seed(put = seedarr)
!     A=0
!    nredges = 0
!    do i = 1, ndim-1
!    do j = i+1, ndim
!    call random_number(r)
!    if (r > 0.8) then ! 0.2
! !    write (*,*) i,j
!    A(i,j) = hop + r-0.8
!    A(j,i) = A(i, j)
!    nredges = nredges + 1
!    endif
!    enddo
!    enddo
!    write (*,*) "created matrix with", nredges, "edges."
    allocate(U(ndim, ndim), vec(ndim), energ(ndim), M1(ndim, ndim), M2(ndim, ndim), M3(ndim,ndim))
    
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
    
    fe = createEulerExponentialfromGraphData(gd, energ)

! Now follows some testing and the comparison to straight-forward exponentiation via lapack
!

vec = 1.D0
! ! ! !   call fe%vecmult(vec)
! ! ! ! !  write (*,*) vec
! ! ! ! !  write(*,*) "generating comparison data"
! ! ! !   res = vec
! ! ! !   vec = 1
! ! ! !   dn = 3*ndim
! ! ! !  allocate(lwork(dn), rwork(dn), res2(ndim))
! ! ! !    U = A
! ! ! !    call zheev('V', 'U', ndim, U, ndim, energ, lwork, dn, rwork, IERR)
! ! ! !    energ = exp(energ)
! ! ! !    ! apply to vec
! ! ! !    alpha = 1.D0
! ! ! !    beta = 0.D0
! ! ! !    incx = 1
! ! ! !    call ZGEMV('C', ndim, ndim, alpha, U, ndim, vec, incx, beta, res2, incx)
! ! ! !    do i = 1, ndim
! ! ! !         res2(i) = res2(i) * energ(i)
! ! ! !    enddo
! ! ! !    call ZGEMV('N', ndim, ndim, alpha, U, ndim, res2, incx, beta, vec, incx)
! ! ! ! !   write(*, *) vec
! ! ! !    res2 = res-vec
! ! ! ! !    write (*,*) res2
! ! ! !    write (*,*) "norm error: ", dznrm2(ndim, res2, incx)
! ! ! !    deallocate(lwork, rwork, res2)

!enddo ! seed loop
M1 = 1.D0
do i = 1,80
    call fe%lmult(M1)
enddo
! ! ! !   write(*,*) DBLE(M1)
! ! !    M2 = 1.D0
! ! !    call ZGEMM('C', 'N', ndim, ndim, ndim, alpha, U, ndim, M2, ndim, beta, M3, ndim)
! ! !       do i = 1, ndim
! ! !         M3(i,:) = M3(i,:) * energ(i)
! ! !    enddo
! ! !    call ZGEMM('N', 'N', ndim, ndim, ndim, alpha, U, ndim, M3, ndim, beta, M2, ndim)
! ! !    
! ! ! !   write (*,*) DBLE(M2)
! ! !    M3 = M2-M1
! ! !    write (*, *) "Difference in 1-Norm:", zlange('1', ndim, ndim, M3, ndim, lwork)
   do i = 1, gd%ndim
    call gd%verts(i)%destruct()
   enddo
   call fe%dealloc()
   deallocate(U, vec, energ, M1, M2, M3, gd%verts, gd%elems)
end program mscbdecomp
