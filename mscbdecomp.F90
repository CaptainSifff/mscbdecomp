! MIT License
! 
! Copyright (c) 2018 Florian Goth
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
    use vertex_mod
    implicit none
    integer :: ndim, i, j, k, l, deltag, usedcolors
    integer :: nredges, dn, IERR, incx, nbr1
    real(kind=kind(0.D0)) :: hop
    complex (kind=kind(0.d0)), ALLOCATABLE, DIMENSION(:,:) :: A !< the full matrix A
    complex (kind=kind(0.d0)), ALLOCATABLE, DIMENSION(:,:) :: U, M1,M2, M3 !< A temporary matrix
    complex (kind=kind(0.d0)), ALLOCATABLE, DIMENSION(:) :: vec, lwork, rwork, res, res2 !< the vector that we will test on
    real(kind=kind(0.D0)), allocatable, dimension(:) :: energ
    
    type(Vertex), allocatable, dimension(:) :: verts
    logical, allocatable, dimension(:) :: usedcols
    type(node), allocatable, dimension(:) :: nodes
    type(FullExp) :: fe
    real(kind=kind(0.D0)) :: dznrm2, zlange
    integer :: seed
    complex(kind=kind(0.D0)) :: alpha, beta
    ! initialize A with some data
    hop = 0.01
    nredges = 0
! coresponds to chain with next-nearest neighbour hopping and OBC
! ndim = 10
! allocate(A(ndim, ndim))
! do I = 1, ndim-1
! A(I,I+1) = hop
! A(i+1, i) = hop
! enddo
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

ndim = 400
!ndim=7
allocate(A(ndim, ndim))
do seed = 1000,1010
!seed = 4887
!seed = 1061
write (*,*) "seed", seed
call srand(seed)
A=0
nredges = 0
do i = 1, ndim-1
do j = i+1, ndim
if (rand() > 0.6) then ! 0.2
!write (*,*) i,j
A(i,j) = hop
A(j,i) = hop
nredges = nredges + 1
endif
enddo
enddo
write (*,*) "created matrix with", nredges, "edges."
    allocate(U(ndim, ndim), vec(ndim), energ(ndim), M1(ndim, ndim), M2(ndim, ndim), M3(ndim,ndim))
    verts = mat2verts(A)
    call MvG_decomp(verts)
    ! Determine the number of used colors and the number of edges
    deltag = 0
    usedcolors = 0
    nredges = 0
    do i = 1, ndim
        deltag = max(deltag, verts(i)%degree)
        do k = 1, verts(i)%degree
            if (verts(i)%nbrs(k) > i) nredges = nredges + 1
            if (verts(i)%nbrs(k) > size(verts)) then
                write(*,*) "invalid nbr!!!"
                STOP
            endif
            usedcolors = max(usedcolors, verts(i)%cols(k))
        enddo
    enddo
    write (*,*) "Nr edges: ", nredges
    if (usedcolors == deltag) then
        write(*,*) "Maximum Degree", deltag, ". Found", usedcolors," Families -> optimal decomposition"
    else
        write(*,*) "Maximum Degree", deltag, ". Found", usedcolors," Families"
    endif
! set up data in an edges based layout
k = 0
allocate( nodes(nredges), usedcols(usedcolors))
    do i = 1, ndim-1
        ! check validity of the coloring locally
        usedcols = .false.
        do l = 1, verts(i)%degree
            if(verts(i)%cols(l) == 0) then
                write (*,*) "forgotten edge found!"
                STOP
            endif
            if (usedcols(verts(i)%cols(l)) .eqv. .true. ) then
                write (*,*) "invalid coloring!!"
                STOP
            else
                usedcols(verts(i)%cols(l)) = .true.
            endif
        enddo
        do l = 1, usedcolors
            nbr1 = verts(i)%nbrs(verts(i)%nbrbycol(l))
            if (nbr1 > i) then ! nbr1 could be zero if there is no such edge
                k = k+1
                nodes(k)%x = i
                nodes(k)%y = nbr1
                nodes(k)%axy = A(i, nbr1)
                nodes(k)%col = l
            endif
        enddo
    enddo
!STOP
    call fe%init(nodes, usedcolors)
    deallocate(nodes, usedcols)
!STOP
! Now we have to return the decomposed matrices/or setup objects for multiplication with the 
! exponentiated variants.
! ! ! !     do i = 1, ndim
! ! ! !         do j = 1, verts(i)%degree
! ! ! !         write (*,*) i, "->", verts(i)%nbrs(j), " = ", verts(i)%cols(j)
! ! ! !         enddo
! ! ! !     enddo
vec = 1.D0
  call fe%vecmult(vec)
!  write (*,*) vec
!  write(*,*) "generating comparison data"
  res = vec
  vec = 1
  dn = 3*ndim
 allocate(lwork(dn), rwork(dn), res2(ndim))
! ! ! ! ! ! ! ! ! ! ! ! !    U = A
! ! ! ! ! ! ! ! ! ! ! ! !    call zheev('V', 'U', ndim, U, ndim, energ, lwork, dn, rwork, IERR)
! ! ! ! ! ! ! ! ! ! ! ! !    energ = exp(energ)
! ! ! ! ! ! ! ! ! ! ! ! !    ! apply to vec
! ! ! ! ! ! ! ! ! ! ! ! !    alpha = 1.D0
! ! ! ! ! ! ! ! ! ! ! ! !    beta = 0.D0
! ! ! ! ! ! ! ! ! ! ! ! !    incx = 1
! ! ! ! ! ! ! ! ! ! ! ! !    call ZGEMV('C', ndim, ndim, alpha, U, ndim, vec, incx, beta, res2, incx)
! ! ! ! ! ! ! ! ! ! ! ! !    do i = 1, ndim
! ! ! ! ! ! ! ! ! ! ! ! !         res2(i) = res2(i) * energ(i)
! ! ! ! ! ! ! ! ! ! ! ! !    enddo
! ! ! ! ! ! ! ! ! ! ! ! !    call ZGEMV('N', ndim, ndim, alpha, U, ndim, res2, incx, beta, vec, incx)
! ! ! ! ! ! ! ! ! ! ! ! ! !   write(*, *) vec
! ! ! ! ! ! ! ! ! ! ! ! !    res2 = res-vec
! ! ! ! ! ! ! ! ! ! ! ! ! !    write (*,*) res2
! ! ! ! ! ! ! ! ! ! ! ! !    write (*,*) "norm error: ", dznrm2(ndim, res2, incx)
   deallocate(lwork, rwork)
   do i = 1, size(verts)
    call verts(i)%destruct()
   enddo
   call fe%dealloc()
   deallocate(U, vec, energ, M1, M2, M3, verts, res2)
enddo
! do i = 1,80
!    M1 = 1.D0
!    call fe%matmult(M1)
! enddo
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
end program mscbdecomp
