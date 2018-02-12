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
    implicit none
    integer :: ndim, i, j, k
    complex (kind=kind(0.d0)), ALLOCATABLE, DIMENSION(:,:) :: A !< the full matrix A
    complex (kind=kind(0.d0)), ALLOCATABLE, DIMENSION(:,:) :: TMP !< A temporary matrix
    complex (kind=kind(0.d0)), ALLOCATABLE, DIMENSION(:) :: v !< the vector that we will test on
    ! initialize A with some data
    ndim = 10
    allocate(A(ndim, ndim), tmp(ndim, ndim), v(ndim))
    do i = 1, ndim-1
        A(i,i+1) = 0.5
        A(i+1,i) = 0.5
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

    write (*,*) "Hello, world."

end program mscbdecomp
