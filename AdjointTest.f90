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
    integer :: ndim, i, j, k, n, myl, myiostat, iounit
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
ndim = 6
allocate(A(ndim, ndim), lwork(3*ndim))
OPEN(newunit=iounit, file="norms.dat", iostat=myiostat)
IF(myiostat /= 0) THEN
    WRITE(*,*) 'Failed to open norms.dat'
    STOP
ENDIF

    ! initialize A with some data
do myl = 0, 397,1
! corresponds to chain with next-nearest neighbour hopping and OBC
hop = 0.01*myl
A=0
do I = 1, ndim-1, 1
A(I,I+1) = hop
A(i+1, i) = hop
if (i+2 <= ndim) then
    A(I,I+2) = hop
    A(i+2, i) = hop
    endif
enddo

    allocate(U(ndim, ndim), vec(ndim), energ(ndim), M1(ndim, ndim), M2(ndim, ndim), M3(ndim,ndim))
    allocate(Identity(ndim, ndim))
    
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
!     write (*,*) "Nr edges: ", gd%nredges
!     if (gd%usedcolors == gd%deltag) then
!         write(*,*) "Maximum Degree", gd%deltag, ". Found", gd%usedcolors," Families -> optimal decomposition"
!     else
!         write(*,*) "Maximum Degree", gd%deltag, ". Found", gd%usedcolors," Families"
!     endif

! create an Exponential from the color information and the weights of the graph
!
    
    ee = createEulerExponentialfromGraphData(gd)

! Now follows some testing and the comparison to straight-forward exponentiation via lapack
!

Identity=0
do i = 1, ndim
Identity(i,i) = 1
enddo

write (*,*) zlange('F', ndim, ndim, M3, ndim, lwork)

! test adjoint over two
Identity=0
do i = 1, ndim
Identity(i,i) = 1
enddo
call ee%singleexps(1)%adjoint_over_two(Identity)

do i=1,ndim
write (*,*) (dble(Identity(i,j)), j=1,ndim)
enddo

Identity=0
do i = 1, ndim
Identity(i,i) = 1
enddo
write (*,*) ee%nrofcols
!call ee%singleexps(1)%lmult(Identity)
!call ee%singleexps(1)%rmultinv(Identity)
do k = 1,4
    call ee%adjoint_over_two(Identity)
enddo
do k = 1, ndim
Identity(k, k) = Identity(k, k) - 1.0
enddo
write (unit=iounit,fmt=*) hop, zlange('F', ndim, ndim, Identity, ndim, lwork)
! do i=1,ndim
! write (unit=iounit,fmt=*) (dble(Identity(i,j)), j=1,ndim)
! enddo
do i = 1, gd%ndim
    call gd%verts(i)%destruct()
enddo
call ee%dealloc()
deallocate(U, vec, energ, M1, M2, M3, gd%verts, gd%elems, Identity)
enddo
end program singlecolexptest
