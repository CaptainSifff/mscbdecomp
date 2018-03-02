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
    integer :: ndim, i, j, k, l, deltag, cnt, maxcolors, usedcolors
    integer :: availablecolor, nredges, dn, IERR, incx, fantail, fanlen, oldcol, tmpcol
    real(kind=kind(0.D0)) :: hop
    complex (kind=kind(0.d0)), ALLOCATABLE, DIMENSION(:,:) :: A !< the full matrix A
    complex (kind=kind(0.d0)), ALLOCATABLE, DIMENSION(:,:) :: U, M1,M2, M3 !< A temporary matrix
    complex (kind=kind(0.d0)), ALLOCATABLE, DIMENSION(:) :: vec, lwork, rwork, res, res2 !< the vector that we will test on
    real(kind=kind(0.D0)), allocatable, dimension(:) :: energ
    
    type(Vertex), allocatable, dimension(:) :: verts
    integer, allocatable, dimension(:) :: fan
    logical :: check
    type(node), allocatable, dimension(:) :: nodes
    type(FullExp) :: fe
    real(kind=kind(0.D0)) :: dznrm2, zlange
    integer :: seed
    complex(kind=kind(0.D0)) :: alpha, beta
    ! initialize A with some data
    hop = 0.1
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

ndim = 6
allocate(A(ndim, ndim))
! do seed = 1, 100000
! write (*,*) "seed", seed
seed = 2158
call srand(seed)
A=0
do i = 1, ndim-1
do j = i+1, ndim
if (rand() > 0.5) then
write (*,*) i,j
A(i,j) = hop
A(j,i) = hop
endif
enddo
enddo
nredges = 0
    allocate(U(ndim, ndim), vec(ndim), energ(ndim), M1(ndim, ndim), M2(ndim, ndim), M3(ndim,ndim))
! check input
    ! first check diagonal
    do i = 1, ndim
        if(A(i, i) /= 0.D0) then
            write (*, *) "the main-diagonal must be zero!"
            stop
        endif
    enddo
    ! check symmetry
!     do i = 1, ndim
!         do j = 1, ndim
!             if(A(i,j) /= conjg(A(j,i))) then
!                 write (*, *) "Non-hermitian matrix encountered!"
!                 stop
!             endif
!         enddo
!     enddo
    allocate(verts(ndim))
! calculate Vertex degree
    deltag = 0;
    do i = 1, ndim
        cnt = 0
        do j = 1, ndim
            if(A(i, j) /= 0.D0) cnt = cnt +1
        enddo
        call verts(i)%init(cnt)
        k = 1
        do j = 1, ndim
            if(A(i, j) /= 0.D0) then
                verts(i)%nbrs(k) = j
                k = k + 1
            endif
        enddo
        deltag = max(deltag, cnt)
    enddo
    write (*,*) "Delta(G) = ", deltag
    maxcolors = deltag + 1

    ! Starting Vizings algorith as outlined in https://thorehusfeldt.files.wordpress.com/2010/08/gca.pdf
    ! we obtain the edges by looking in the upper triangular part of the matrix for non-zero entries
    do i = 1, ndim-1
        do j = i+1, ndim
        if (A(i, j) /= 0.D0) then
        ! Edge found between vertex i and j
        ! Let's check wether we have free edges at every vertex
            nredges = nredges + 1
! A debugging check that the data is consistent:
 check = .false.
do k = 1, verts(i)%degree
 if (verts(i)%nbrs(k) == j) check = .true.
enddo
if(check .eqv. .false.) then
write(*,*) "inconsistent data!"
STOP
endif
            availablecolor = find_common_free_color(verts(i), verts(j), maxcolors)
            if(availablecolor == 0) then
                ! Our starting vertex is verts(i), our target vertex is verts(j)
                ! Now we need to construct a Vizing fan around verts(i)
                write (*,*) "out of colors. Trying to downshift Vizing fan"
                allocate(fan(verts(i)%degree))
                oldcol = verts(i)%find_maximal_fan(verts, j, maxcolors, fan, fanlen)
                write (*,*) "oldcol = ", oldcol, "fanlen = ", fanlen, fan
                if (oldcol .ne. 0) then
                    ! the end of the fan has a free color -> down-shifting sufficient
                    do k = 1, fanlen-1
                        tmpcol = verts(i)%get_edge_color(fan(k+1))
                        call verts(i)%set_edge_color(fan(k), tmpcol)
                        call verts(fan(k))%set_edge_color(i, tmpcol)
                    enddo
                    call verts(i)%set_edge_color(fan(fanlen), oldcol)
                    call verts(fan(fanlen))%set_edge_color(i, oldcol)
                else
                    ! We would need to inverse a path
                    write (*,*) "inversion of paths not implemented!"
                    STOP
                endif
                deallocate(fan)
            else
                write(*,*) availablecolor
                ! set that color
                call verts(i)%set_edge_color(j, availablecolor);
                call verts(j)%set_edge_color(i, availablecolor);
            endif
        endif
        enddo
    enddo
! Determine the number of used colors
    usedcolors = deltag
    do i = 1, ndim
        do j = 1, verts(i)%degree
            if (verts(i)%cols(j) == deltag + 1) usedcolors = deltag + 1
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
allocate( nodes(nredges))
    do i = 1, ndim-1
        do j = i+1, ndim
        if (A(i, j) /= 0.D0) then
            k = k + 1
            nodes(k)%x = i
            nodes(k)%y = j
            nodes(k)%axy = A(i,j)
            do l = 1, verts(i)%degree
                if(verts(i)%nbrs(l) == j) nodes(k)%col = verts(i)%cols(l)
            enddo
        endif
        enddo
    enddo
    call fe%init(nodes, usedcolors)
    deallocate(nodes)
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
   U = A
   call zheev('V', 'U', ndim, U, ndim, energ, lwork, dn, rwork, IERR)
   deallocate(lwork, rwork)
   energ = exp(energ)
   ! apply to vec
   alpha = 1.D0
   beta = 0.D0
   incx = 1
   call ZGEMV('C', ndim, ndim, alpha, U, ndim, vec, incx, beta, res2, incx)
   do i = 1, ndim
        res2(i) = res2(i) * energ(i)
   enddo
   call ZGEMV('N', ndim, ndim, alpha, U, ndim, res2, incx, beta, vec, incx)
!   write(*, *) vec
   res2 = res-vec
!    write (*,*) res2
   write (*,*) "norm error: ", dznrm2(ndim, res2, incx)
   do i = 1, size(verts)
    call verts(i)%destruct()
   enddo
   call fe%dealloc()
   deallocate(U, vec, energ, M1, M2, M3, verts, res2)
!enddo
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
