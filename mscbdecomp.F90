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
    integer :: ndim, i, j, k, deltag, cnt, maxcolors, usedcolors
    integer :: availablecolor
    real(kind=kind(0.D0)) :: hop
    complex (kind=kind(0.d0)), ALLOCATABLE, DIMENSION(:,:) :: A !< the full matrix A
    complex (kind=kind(0.d0)), ALLOCATABLE, DIMENSION(:,:) :: TMP !< A temporary matrix
    complex (kind=kind(0.d0)), ALLOCATABLE, DIMENSION(:) :: vec !< the vector that we will test on
    type(Vertex), allocatable, dimension(:) :: verts
    logical :: check
    ! initialize A with some data
    ndim = 10
    hop = 0.5
    allocate(A(ndim, ndim), tmp(ndim, ndim), vec(ndim))
    do i = 1, ndim-1
        A(i,i+1) = hop
        A(i+1,i) = hop
    enddo
    do i = 1, ndim-2
        A(i,i+2) = hop
        A(i+2,i) = hop
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
            ! One of the vertices has no free color -> we don't know what to do yet.
            endif
            write(*,*) availablecolor
            ! set that color
            call verts(i)%set_edge_color(j, availablecolor);
            call verts(j)%set_edge_color(i, availablecolor);
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
    if (usedcolors == deltag) then
        write(*,*) "Maximum Degree", deltag, ". Found", usedcolors," Families -> optimal decomposition"
    else
        write(*,*) "Maximum Degree", deltag, ". Found", usedcolors," Families"
    endif
! Now we have to return the decomposed matrices/or setup objects for multiplication with the 
! exponentiated variants.
    do i = 1, ndim
        do j = 1, verts(i)%degree
        write (*,*) i, "->", verts(i)%nbrs(j), " = ", verts(i)%cols(j)
        enddo
    enddo
    
end program mscbdecomp
