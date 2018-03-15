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
    integer :: availablecolor, nredges, dn, IERR, incx, fantail, fanlen, oldcol, tmpcol, nbr1, nbr2
    integer :: curcol, col(2), colctr
    integer :: ver, nbr
    real(kind=kind(0.D0)) :: hop
    complex (kind=kind(0.d0)), ALLOCATABLE, DIMENSION(:,:) :: A !< the full matrix A
    complex (kind=kind(0.d0)), ALLOCATABLE, DIMENSION(:,:) :: U, M1,M2, M3 !< A temporary matrix
    complex (kind=kind(0.d0)), ALLOCATABLE, DIMENSION(:) :: vec, lwork, rwork, res, res2 !< the vector that we will test on
    real(kind=kind(0.D0)), allocatable, dimension(:) :: energ
    
    type(Vertex), allocatable, dimension(:) :: verts
    integer, allocatable, dimension(:) :: fan, fannbr, fanedgecol
    logical, allocatable, dimension(:) :: usedcols
    logical :: check, stoppath
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

ndim = 2000
!ndim=7
allocate(A(ndim, ndim))
do seed = 1000,1010
!seed = 4887
!seed = 1061
write (*,*) "seed", seed
call srand(seed)
A=0
do i = 1, ndim-1
do j = i+1, ndim
if (rand() > 0.6) then ! 0.2
!write (*,*) i,j
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
            if(A(i, j) /= 0.D0) cnt = cnt + 1
        enddo
        deltag = max(deltag, cnt)
    enddo

    write (*,*) "Delta(G) = ", deltag
    maxcolors = deltag + 1
    
    do i = 1, ndim
        cnt = 0
        do j = 1, ndim
            if(A(i, j) /= 0.D0) cnt = cnt + 1
        enddo
        call verts(i)%init(cnt, maxcolors)
        k = 1
        do j = 1, ndim
            if(A(i, j) /= 0.D0) then
                verts(i)%nbrs(k) = j
                k = k + 1
            endif
        enddo
    enddo

    ! Starting Vizings algorith as outlined in https://thorehusfeldt.files.wordpress.com/2010/08/gca.pdf
    ! we obtain the edges by looking in the upper triangular part of the matrix for non-zero entries
    do i = 1, ndim-1
        do j = i+1, ndim
        if (A(i, j) /= 0.D0) then
!         write (*,*) "current state of vert", 5
!             write (*,*) verts(5)%nbrbycol
!             write (*,*) verts(5)%cols
        ! Edge found between vertex i and j
        ! Let's check wether we have free edges at every vertex
            nredges = nredges + 1
! ! ! ! A debugging check that the data is consistent:
! ! !  check = .false.
! ! ! do k = 1, verts(i)%degree
! ! !  if (verts(i)%nbrs(k) == j) check = .true.
! ! ! enddo
! ! ! if(check .eqv. .false.) then
! ! ! write(*,*) "inconsistent data!"
! ! ! STOP
! ! ! endif
            availablecolor = find_and_try_to_set_common_free_color(i, j, verts, maxcolors)
            if(availablecolor == 0) then
                ! Our starting vertex is verts(i), our target vertex is verts(j)
                ! Now we need to construct a Vizing fan around verts(i)
#ifndef NDEBUG
                write (*,*) "Out of colors. construct fan around", i
#endif
                allocate(fan(verts(i)%degree), fannbr(verts(i)%degree), fanedgecol(verts(i)%degree))
                
                fannbr = 0
                fanedgecol = 0
                do k = 1, verts(i)%degree
                    if(verts(i)%nbrs(k) == j) fannbr(1) = k
                enddo
                oldcol = verts(i)%find_maximal_fan(verts, j, maxcolors, fan, fannbr, fanedgecol, fanlen)
#ifndef NDEBUG
                write (*,*) "oldcol = ", oldcol, "fanlen = ", fanlen, fan
#endif
                if (oldcol .ne. 0) then
#ifndef NDEBUG
                write (*,*) "Free color available in fan -> downshift"
#endif
                    ! the end of the fan has a free color -> down-shifting sufficient
                    call downshift_fan_and_set_color(i, fan, fannbr, fanedgecol, fanlen, oldcol, verts)
                else
                    ! We would need to inverse a path
                    ! set up the start of the path at the fan
                    ! determine the two colors for the c-d path
                    col(1) = verts(fan(fanlen))%findfreecolor()
                    col(2) = verts(i)%findfreecolor()
                    call construct_and_invert_path(col, verts, i)
                    ! try again to obtain a color
                    availablecolor = find_and_try_to_set_common_free_color(i, j, verts, maxcolors)
                    if (availablecolor == 0) then
                        ! We have to do a proper downshifting
                        ! find w in fan()
                        k = 1
                        ver = 0
                        do while (k <= fanlen)
                            if(verts(fan(k))%iscolorfree(col(1))) then
                                write(*,*) "color found at ",k, "of ", fanlen
                                ver = k
                                k = fanlen
                            endif
                            k = k+1
                        enddo
                        write (*,*) "determined ", ver
                        ! the path inversion has modified the information contained in the auxiliary arrays
                        ! fannbr and fanedgecol -> We need to update that
                        ! The edge that has been modified had color col(1) and now has color col(2)
                        ! search the edge that has changed
                        k = 1
                        do while ((k <= fanlen) .and. (fanedgecol(k) .ne. col(1)))
                            k = k+1
                        enddo
                        if (k > fanlen) then
                            write (*,*) "changed color ", col(1), "not found"
                            STOP
                        endif
                        fanedgecol(k) = col(2)
                        ! downshift the fan until the position of w
                        call downshift_fan_and_set_color(i, fan, fannbr, fanedgecol, ver, col(1), verts)
                    endif
                endif
                deallocate(fan, fannbr, fanedgecol)
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
allocate( nodes(nredges), usedcols(maxcolors))
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
        do l = 1, maxcolors
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
