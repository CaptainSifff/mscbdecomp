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

module vertex_mod
    implicit none

    !Another helper type for the cd_X path construction
    type :: Path
        integer :: avamem
        integer :: tail
        integer, allocatable, dimension(:) :: vertices
        integer, allocatable, dimension(:) :: nbrindex
    contains
        procedure :: init => Path_init
        procedure :: dealloc => Path_dealloc
        procedure :: pushback => Path_pushback
        procedure :: at => Path_at
        procedure :: back => Path_back
        procedure :: length => Path_length
    end type Path
    
    type :: Simplenode
        integer :: x,y
        complex (kind=kind(0.d0)) :: s,c
    end type Simplenode

    type :: node
        integer :: x,y, col
        complex (kind=kind(0.d0)) :: axy
    end type node

    type :: SingleColExp
        integer :: nrofentries
        type(Simplenode), dimension(:), allocatable :: nodes
    contains
        procedure :: init => SingleColExp_init
        procedure :: dealloc => SingleColExp_dealloc
        procedure :: vecmult => SingleColExp_vecmult
        procedure :: matmult => SingleColExp_matmult
    end type

    type :: FullExp
        integer :: nrofcols
        type(SingleColExp), dimension(:), allocatable :: singleexps
    contains
        procedure :: init => FullExp_init
        procedure :: dealloc => FullExp_dealloc
        procedure :: vecmult => FullExp_vecmult
        procedure :: matmult => FullExp_matmult
    end type FullExp
    
    type :: Vertex
        integer :: degree
        integer, allocatable, dimension(:) :: nbrs
        integer, allocatable, dimension(:) :: cols
        integer, allocatable, dimension(:) :: nbrbycol
    contains
        procedure :: init => vertex_init
        procedure :: destruct => vertex_destruct
        procedure :: any_color_available => vertex_any_color_available
        procedure :: set_edge_color => vertex_set_edge_color
        procedure :: get_edge_color => vertex_get_edge_color
        procedure :: find_maximal_fan => vertex_find_maximal_fan
        procedure :: findfreecolor => vertex_findfreecolor
        procedure :: iscolorfree => vertex_iscolorfree
    end type Vertex
contains

function vertex_iscolorfree(this, col) result(r)
    class(Vertex) :: this
    integer, intent(in) :: col
    logical :: r
    integer :: k
    
    r = .true.
    do k = 1, this%degree
    ! some debug check while we are here...
        if ((this%cols(k) == col) .and. (r .eqv. .false.)) then
        write (*,*) "duplicate color found!"
        STOP
        endif
        
        if (this%cols(k) == col) r = .false.
    enddo
end function

subroutine Path_init(this)
    class(Path) :: this
    integer :: temp
    this%tail = 1
    this%avamem = 4096/SIZEOF(temp) ! allocate a page of memory
    allocate(this%vertices(this%avamem), this%nbrindex(this%avamem))
end subroutine Path_init

subroutine Path_dealloc(this)
    class(Path) :: this
    deallocate(this%vertices, this%nbrindex)
end subroutine

subroutine Path_pushback(this, vert, nbridx)
    class(Path) :: this
    integer, intent(in) :: vert, nbridx
    integer, allocatable, dimension(:) :: temp1, temp2
    integer :: i
    if (this%tail == this%avamem) then
        ! reallocate the memory
        write (*,*) "not enough space!"
        call MOVE_ALLOC(temp1, this%vertices)
        call MOVE_ALLOC(temp2, this%nbrindex)
        allocate(this%vertices(2*this%avamem), this%nbrindex(2*this%avamem))
        do i = 1, this%avamem
            this%vertices(i) = temp1(i)
            this%nbrindex(i) = temp2(i)
        enddo
        deallocate(temp1, temp2)
        this%avamem = 2*this%avamem
        STOP
    endif
    this%vertices(this%tail) = vert
    this%nbrindex(this%tail) = nbridx
    this%tail = this%tail + 1
end subroutine

subroutine Path_at(this, pos, vert, nbridx)
    class(Path) :: this
    integer, intent(in) :: pos
    integer, intent(out) :: vert, nbridx
    vert = this%vertices(pos)
    nbridx = this%nbrindex(pos)
end subroutine

subroutine Path_back(this, vert, nbridx)
    class(Path) :: this
    integer, intent(out) :: vert, nbridx
    vert = this%vertices(this%tail-1)
    nbridx = this%nbrindex(this%tail-1)
end subroutine

function Path_length(this) result(l)
    class(Path) :: this
    integer :: l
    l = this%tail-1
end function 

function vertex_findfreecolor(this, maxcols)  result(col)
    class(Vertex) :: this
    integer :: maxcols
    integer :: col
    integer :: i
    logical :: usedcols(maxcols)

    usedcols = .false.
    do i = 1, this%degree
        if(this%cols(i) /= 0) usedcols(this%cols(i)) = .true.
    enddo
    col = 0
    do i = maxcols, 1, -1
        if(usedcols(i) .eqv. .false.) col = i
    enddo
end function vertex_findfreecolor

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> This function constructs a Vizing Fan.
!
!> @param[in] this The vertex where we construct a Vizing fan around
!> @param[in] verts the total array of all vertices
!> @param[in] v0 The incident vertex that still is uncolored.
!> @param[out] f a linear array where we construct the fan
!> @param[out] fanlen The length of the entire fan
!--------------------------------------------------------------------
function vertex_find_maximal_fan(this, verts, v0, maxcols, f, fanlen) result(rescol)
    class(Vertex) :: this
    type(Vertex), dimension(:), intent(in) :: verts
    integer, intent(out) :: fanlen
    integer :: rescol
    integer, dimension(:), intent(out) :: f
    integer :: v0, maxcols
    integer :: col, i, j, ctr
    logical, allocatable, dimension(:) :: usedcols !< in usedcols we track the colors we have chosen during construction of the fan.
    
!!!!!!!    logical, allocatable, dimension(:) :: dbgcols
    
    allocate (usedcols(maxcols))
    fanlen = 1
    f = 0
    f(1) = v0
    ctr = 1
    rescol = 0
    usedcols = .false.
    do while ((ctr < this%degree) .and. (rescol == 0))
#ifdef DEBUG
    write (*,*) "try to attach to fan"
! ! ! ! ! ! ! check that the fantail has valid color information
! ! ! ! ! ! ! check that no wrong color has crept in
! ! ! ! ! ! do i = 1, verts(f(fanlen))%degree
! ! ! ! ! ! if (verts(f(fanlen))%cols(i) > maxcols) then
! ! ! ! ! ! write (*,*) "invalid color found!"
! ! ! ! ! ! STOP
! ! ! ! ! ! endif
! ! ! ! ! ! enddo
! ! ! ! ! ! allocate (dbgcols(maxcols))
! ! ! ! ! ! dbgcols = .false.
! ! ! ! ! ! do i = 1, verts(f(fanlen))%degree
! ! ! ! ! ! if (verts(f(fanlen))%cols(i) /= 0) then
! ! ! ! ! ! if (dbgcols(verts(f(fanlen))%cols(i)) .eqv. .false.) then
! ! ! ! ! ! dbgcols(verts(f(fanlen))%cols(i)) = .true.
! ! ! ! ! ! else
! ! ! ! ! ! write(*,*) "duplicate color found!"
! ! ! ! ! ! STOP
! ! ! ! ! ! endif
! ! ! ! ! ! endif
! ! ! ! ! ! enddo
! ! ! ! ! ! deallocate(dbgcols)
#endif
        ! determine free color at end of fan that has not already been used in the construction of the fan.
        col = 0
        i = 1
        do while ((i <= maxcols) .and. (col == 0))
            if (usedcols(i) .neqv. .true.) then
                j = 1
                do while ((j < verts(f(fanlen))%degree) .and. (verts(f(fanlen))%cols(j) /= i))
                    j = j + 1
                enddo
                ! special treatment for the last entry
                if((verts(f(fanlen))%cols(j) /= i)) j = j+1
                if (j > verts(f(fanlen))%degree) col = i
#ifdef DEBUG
            else
            write (*,*) "ignoring color ", i
#endif
            endif
            i = i + 1
        enddo
        if (col /= 0) then
            ! col is now a small color that is free at f(fanlen)
            ! determine incident edge that has exactly this color
            i = 1
            do while ((f(fanlen + 1) == 0) .and. (i <= this%degree))
                if (this%cols(i) == col) f(fanlen + 1) = this%nbrs(i)
            i = i + 1
            end do
            fanlen = fanlen + 1
            usedcols(col) = .true.
            ! let's see wether we can stop the fan construction since we find matching colors
            rescol = find_common_free_color(this, verts(f(fanlen)), maxcols)
        else
            ! force termination of fan construction
            write (*,*) "fan construction stops. No color available that has NOT been already used in the fan."
            ctr = this%degree
        endif
        ctr = ctr + 1
    enddo
    deallocate(usedcols)
end function vertex_find_maximal_fan

subroutine SingleColExp_vecmult(this, vec)
    class(SingleColExp) :: this
    complex(kind=kind(0.D0)), dimension(:) :: vec
    integer :: i, j
    complex(kind=kind(0.D0)) :: t1,t2
    do i = 1, this%nrofentries! for every matrix
        t1 = vec(this%nodes(i)%x)
        t2 = vec(this%nodes(i)%y)
        vec(this%nodes(i)%x) = this%nodes(i)%c * t1 + this%nodes(i)%s* t2
        vec(this%nodes(i)%y) = this%nodes(i)%c * t2 + this%nodes(i)%s* t1
    enddo
end subroutine SingleColExp_vecmult

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> Perform the multiplication of this exponential with a matrix.
!
!> @param[in] this The vertex that we consider
!> @param[inout] mat the matrix that we modify.
!--------------------------------------------------------------------
subroutine SingleColExp_matmult(this, mat)
    class(SingleColExp) :: this
    complex(kind=kind(0.D0)), dimension(:, :) :: mat
    integer :: i, j, k, ndim
    integer, parameter :: step = 2
    complex(kind=kind(0.D0)) :: t1(step), t2(step)
    
    ndim = size(mat,1)
    do j = 1, ndim, step
        do i = 1, this%nrofentries! for every matrix
            do k = 1,step
                t1(k) = mat(this%nodes(i)%x, j+k)
                t2(k) = mat(this%nodes(i)%y, j+k)
            enddo
            do k = 1, step
                mat(this%nodes(i)%x, j+k) = this%nodes(i)%c * t1(k) + this%nodes(i)%s* t2(k)
                mat(this%nodes(i)%y, j+k) = this%nodes(i)%c * t2(k) + this%nodes(i)%s* t1(k)
            enddo
        enddo
    enddo
end subroutine SingleColExp_matmult

subroutine SingleColExp_init(this, nodes, nredges)
    class(SingleColExp) :: this
    type(node), dimension(:), intent(in) :: nodes
    integer, intent(in) :: nredges
    integer :: i
    allocate(this%nodes(nredges))
    this%nrofentries = nredges
    write(*,*) "Setting up strict. sparse matrix with ", nredges, "edges"
    do i = 1, nredges
        this%nodes(i)%x = nodes(i)%x
        this%nodes(i)%y = nodes(i)%y
        this%nodes(i)%c = cosh(nodes(i)%axy)
        this%nodes(i)%s = sinh(nodes(i)%axy)
    enddo
! All nodes that we have been passed are now from a single color.
! They constitute now a strictly sparse matrix.
! Further processing of the entries could be done here.
end subroutine SingleColExp_init

subroutine SingleColExp_dealloc(this)
    class(SingleColExp) :: this
    deallocate(this%nodes)
end subroutine SingleColExp_dealloc

subroutine FullExp_dealloc(this)
    class(FullExp) :: this
    deallocate(this%singleexps)
end subroutine FullExp_dealloc

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> This function multiplies this full exponential with a vector
!
!> @param[in] this The exponential opbject
!> @param[in] vec The vector that we multiply
!--------------------------------------------------------------------
subroutine FullExp_vecmult(this, vec)
    class(FullExp) :: this
    complex(kind=kind(0.D0)), dimension(:) :: vec
    integer :: i
    do i = 1, this%nrofcols
        call this%singleexps(i)%vecmult(vec)
    enddo
end subroutine FullExp_vecmult

subroutine FullExp_matmult(this, mat)
    class(FullExp) :: this
    complex(kind=kind(0.D0)), dimension(:, :) :: mat
    integer :: i
    do i = 1, this%nrofcols
        call this%singleexps(i)%matmult(mat)
    enddo
end subroutine FullExp_matmult

subroutine FullExp_init(this, nodes, usedcolors)
    class(FullExp) :: this
    type(node), dimension(:), intent(in) :: nodes
    integer, intent(in) :: usedcolors
    integer, dimension(:), allocatable :: nredges, edgectr
    integer :: i, maxedges
    type(node), dimension(:, :), allocatable :: simplenodes
    write(*,*) "Setting up Full Checkerboard exponential."
    ! Determine the number of matrix entries in each family
    allocate (nredges(usedcolors), edgectr(usedcolors))
    nredges = 0
    this%nrofcols = usedcolors
    do i = 1, size(nodes)
        nredges(nodes(i)%col) = nredges(nodes(i)%col) + 1
    enddo
    maxedges = maxval(nredges)
    edgectr = 1
    allocate(simplenodes(usedcolors, maxedges))
    do i = 1, size(nodes)
        simplenodes(nodes(i)%col, edgectr(nodes(i)%col)) = nodes(i)
        edgectr(nodes(i)%col) = edgectr(nodes(i)%col) + 1
    enddo
!     do i = 1, usedcolors
!     write (*,*) edgectr(i), nredges(i)
!     enddo
    allocate(this%singleexps(usedcolors))
    do i = 1, usedcolors
        call this%singleexps(i)%init(simplenodes(i, :), nredges(i))
    enddo
    deallocate(simplenodes)
    deallocate(nredges, edgectr)
end subroutine FullExp_init

subroutine vertex_init(this, deg, maxcols)
    class(vertex) :: this
    integer, intent(in) :: deg
    integer :: maxcols
    this%degree = deg
    allocate(this%nbrs(deg), this%cols(deg), this%nbrbycol(maxcols))
    this%cols = 0
end subroutine vertex_init

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> get the color of an incident edge.
!
!> @param[in] this The vertex that we consider
!> @param[in] vert the index of the neighbour
!> @param[in] col the value of the color
!--------------------------------------------------------------------
function vertex_get_edge_color(this, vert) result(col)
    class(vertex) :: this
    integer, intent(in) :: vert
    integer :: col
    integer :: i
    
    do i = 1, this%degree
        if(this%nbrs(i) == vert) col = this%cols(i)
    enddo
end function vertex_get_edge_color

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> Set the color of an incident edge.
!
!> @param[in] this The vertex that we consider
!> @param[in] vert the index of the neighbour
!> @param[in] col the value of the color
!--------------------------------------------------------------------
subroutine vertex_set_edge_color(this, vert, col)
    class(vertex) :: this
    integer, intent(in) :: vert, col
    integer :: i
    
    do i = 1, this%degree
        if(this%nbrs(i) == vert) this%cols(i) = col
    enddo
end subroutine vertex_set_edge_color

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> Tidy up allocated space
!
!> @param[in] this The vertex that we consider
!--------------------------------------------------------------------
subroutine vertex_destruct(this)
    class(vertex) :: this
    deallocate(this%nbrs, this%cols, this%nbrbycol)
end subroutine vertex_destruct

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> This function checks wether a free color is available and returns
!> the index of the neighbour vertex
!
!> @param[in] this The vertex that we consider
!> @return the index of the neighbour where the edge is still colorless
!--------------------------------------------------------------------

function vertex_any_color_available(this) result(n)
    class(vertex), intent(in) :: this
    integer :: n, i
    n = 0
    do i = this%degree, 1, -1
        if (this%cols(i) == 0) n = this%nbrs(i)
    enddo
end function vertex_any_color_available

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> This function returns the smallest available color
!> that is free at both vertices
!
!> @param[in] v one vertex that we consider
!> @param[in] w the other vertex that we consider
!> @param[in] maxcolors The available maximum number colors
!> @return the smallest common available color
!--------------------------------------------------------------------
function find_common_free_color(v, w, maxcols) result(col)
    implicit none
    class(vertex), intent(in) :: v, w
    integer, intent(in) :: maxcols
    integer :: col, i
    logical, allocatable, dimension(:) :: usedcols

    col = 0
    allocate(usedcols(maxcols))
    usedcols = .false.
    do i = 1, v%degree
        if(v%cols(i) /= 0) usedcols(v%cols(i)) = .true.
    enddo
    do i = 1, w%degree
        if(w%cols(i) /= 0) usedcols(w%cols(i)) = .true.
    enddo
    ! The locations where we find zeroes in this array are colors that are unused in both vertices
    do i = maxcols, 1, -1
        if (usedcols(i) .eqv. .false.) col = i
    enddo
    deallocate(usedcols)
end function find_common_free_color

end module vertex_mod
