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

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> This function returns the smallest available color at a vertex
!
!> @param[in] v The vertex that we consider
!> @param[in] edgecolors the array with the coloring state
!> @param[in] maxcolors The available maximum number colors
!> @return the smallest available color
!--------------------------------------------------------------------

function findfreecolor(v, degv, w, degw edgecolors, maxcolor) result(c)
    implicit none
    integer, intent(in) :: v, degv, w, degw, maxcolor
    integer, dimension(:), intent(in) :: edgecolors
    integer, intent(out) :: c
    integer, allocatable, dimension(:) :: usedcols
    integer :: maxdeg, i
    
    maxdeg = max(degv, degw)
    allocate(usedcols(maxdeg))
    usedcols = 0
    do i = 1, degv
        if (edgecolors(v * maxcolor + i) /= 0) usedcols(edgecolors(v * maxcolor + i)) = 1
    enddo
    do i = 1, degw
        if (edgecolors(w * maxcolor + i) /= 0) usedcols(edgecolors(w * maxcolor + i)) = 1
    enddo
    ! The point where we find zeroes in this array are colors that are unused in both vertices
    
    do i = maxdeg, 1, -1
        if (usedcols[i] == 0) c = i
    enddo
    deallocate(usedcols)
end function freecolor

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> This function returns the smallest available color at a vertex
!
!> @param[in] v The vertex that we consider
!> @param[in] edgecolors the array with the coloring state
!> @param[in] maxcolors The available maximum number colors
!> @return the smallest available color
!--------------------------------------------------------------------

function findfreeedge(v, degree, edgecolors, maxcolors) result(c)
    implicit none
    integer, intent(in) :: v, degree, maxcolors
    integer, dimension(:), intent(in) :: edgecolors
    integer, intent(out) :: edge
    integer :: i
    edge = -1
    do i = degree, 1, -1
        if (edgecolors(maxcolors*i) == 0.D0) edge = i !> Not yet used edge found
    enddo

end function findfreeedge

program mscbdecomp
    implicit none
    integer :: ndim, i, j, k, deltag, cnt, vleg, wleg, maxcolors
    real(kind=kind(0.D0) :: hop
    complex (kind=kind(0.d0)), ALLOCATABLE, DIMENSION(:,:) :: A !< the full matrix A
    complex (kind=kind(0.d0)), ALLOCATABLE, DIMENSION(:,:) :: TMP !< A temporary matrix
    complex (kind=kind(0.d0)), ALLOCATABLE, DIMENSION(:) :: vec !< the vector that we will test on
    integer, allocatable, dimension(:) :: edgecolors !< the array where we store the edge colors
    integer, allocatable, dimension(:) :: degrees !< Here we store the edge degree of a particular vertex
    ! initialize A with some data
    ndim = 10
    hop = 0.5
    allocate(A(ndim, ndim), tmp(ndim, ndim), vec(ndim), degrees(ndim))
    do i = 1, ndim-1
        A(i,i+1) = hop
        A(i+1,i) = hop
    enddo
! check input
    ! first check diagonal
    do i = 1, ndim
        if(A(i, i) /= 0.D0) then
            write (*, *) "the main-diagonal must be zero!"
            stop
        endif
    enddo
    ! check symmetry
    do i = 1, ndim
        do j = 1, ndim
            if(A(i,j) /= conjg(A(j,i))) then
                write (*, *) "Non-hermitian matrix encountered!"
                stop
            endif
        enddo
    enddo
! calculate Vertex degree
    deltag = 0;
    do i = 1, ndim
        cnt = 0
        do j = 1, ndim
            if(A(i, j) /= 0.D0) cnt = cnt +1
        enddo
        degrees(i) = cnt
        deltag = max(deltag, cnt)
    enddo
    write (*,*) "Delta(G) = ", deltag
    maxcolors = deltag + 1

    ! set up storage foer the edge colors
    ! Since we know the maximum degree, we can over-estimate the use and use a linear array
    ! We encode an unused color at an edge as 0 else we use integers from [1, Delta(G) + 1]
    allocate(edgecolors((deltag+1)*ndim))
    ! Starting Vizings algorith as outlined in https://thorehusfeldt.files.wordpress.com/2010/08/gca.pdf
    ! we obtain the edges by looking in the upper triangular part of the matrix for non-zero entries
    do i = 1, ndim-1
        do j = i+1, ndim
        if (A(i, j) /= 0.D0) then
        ! Edge found between vertex i and j
        ! Let's check wether we have free edges at every vertex
            vleg = findfreeedge(i, degree(i), edgecolors, maxcolors)
            wleg = findfreeedge(j, degree(j), edgecolors, maxcolors)
            if ((vleg == -1) .or. (wleg == -1)) then
            ! One of the vertices has no free color -> we don't know what to do yet.
            endif
            
        endif
        enddo
    enddo

end program mscbdecomp
