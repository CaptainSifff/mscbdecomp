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
        procedure :: erase_edge_color => vertex_erase_edge_color
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

function vertex_findfreecolor(this)  result(col)
    class(Vertex) :: this
    integer :: col
    integer :: i
    logical :: usedcols(size(this%nbrbycol))

!     do i = maxcols, 1, -1
!     if(this%nbrbycol(i) == 0) col = i
!     enddo
    usedcols = .false.
    do i = 1, this%degree
        if(this%cols(i) /= 0) usedcols(this%cols(i)) = .true.
    enddo
    col = 0
    do i = size(this%nbrbycol), 1, -1
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
!> @param[out] fannbr An array with the positions of the vertices at verts(fanpos)
!> @param[out] fanedgecol The colors of the edges in the fan
!> @param[out] fanlen The length of the entire fan
!> @result rescol
!--------------------------------------------------------------------
function vertex_find_maximal_fan(this, verts, v0, maxcols, f, fannbr, fanedgecol, fanlen) result(rescol)
    class(Vertex) :: this
    type(Vertex), dimension(:), intent(in) :: verts
    integer, intent(out) :: fanlen
    integer :: rescol, tmpcol
    integer, dimension(:), intent(out) :: f, fannbr, fanedgecol
    integer :: v0, maxcols
    integer :: col, i, j, ctr, nrfreecols
    logical, allocatable, dimension(:) :: usedcols !< in usedcols we track the colors we have chosen during construction of the fan.
    integer, allocatable, dimension(:) :: freecols ! determine the colors that are available at this
    
!!!!!!!    logical, allocatable, dimension(:) :: dbgcols
    
    allocate (usedcols(maxcols))
    fanlen = 1
    f = 0
    f(1) = v0
    ctr = 1
    rescol = 0
    usedcols = .false.
    nrfreecols = 0

    do i = 1, maxcols
        if (this%nbrbycol(i) == 0) nrfreecols = nrfreecols + 1
    enddo
    allocate(freecols(nrfreecols))
    j = 1
    do i = 1, maxcols
        if (this%nbrbycol(i) == 0) then
        freecols(j) = i
        j = j + 1
        endif
    enddo
    do while ((ctr < this%degree) .and. (rescol == 0))
        ! determine free color at end of fan that has not already been used in the construction of the fan.
        col = 0
        i = 1
        tmpcol = 0
        do while ((i <= maxcols) .and. (col == 0))
            if (usedcols(i) .eqv. .false.) then
#ifndef NDEBUG
           write (*,*) "checking color ", i, "fanlen = ", fanlen
           write (*,*) "looking at vertex", f(fanlen)
            write (*,*) verts(f(fanlen))%nbrbycol
            write (*,*) verts(f(fanlen))%cols
            write(*,*) "----------------------------------"
#endif
                if(verts(f(fanlen))%nbrbycol(i) == 0) then
!                write (*,*) "free:", i
                col = i
                endif
            endif
            i = i + 1
        enddo
        if (col /= 0) then
            ! col is now a small color that is free at f(fanlen)
            ! determine incident edge that has exactly this color
            fannbr(fanlen+1) = this%nbrbycol(col)
            fanedgecol(fanlen+1) = col
            f(fanlen+1) = this%nbrs(this%nbrbycol(col))
            fanlen = fanlen + 1
            usedcols(col) = .true.
            ! let's see wether we can stop the fan construction since we find matching colors
            rescol = 0
            do i = nrfreecols, 1, -1
                if (verts(f(fanlen))%nbrbycol(freecols(i)) == 0) rescol = freecols(i)
            enddo
        else
            ! force termination of fan construction
            write (*,*) "fan construction stops. No color available that has NOT been already used in the fan."
            ctr = this%degree
        endif
        ctr = ctr + 1
    enddo
    deallocate(usedcols, freecols)
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
    this%nbrbycol = 0
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
!> erases the color of an incident edge.
!
!> @param[in] this The vertex that we consider
!> @param[in] vert the index of the neighbour
!> @param[in] col the value of the color
!--------------------------------------------------------------------
subroutine vertex_erase_edge_color(this, vert)
    class(vertex) :: this
    integer, intent(in) :: vert
    integer :: i, tmp
#ifndef NDEBUG
    write (*,*) "Called erase_color!!!!!!!!!!!!!!!!!!!!!!!!!!!"
#endif
    do i = 1, this%degree
        if(this%nbrs(i) == vert) then
            tmp = this%cols(i)
#ifndef NDEBUG
            write (*,*) "erasing color ", tmp, "of pos ", i, "which is", this%nbrs(i)
#endif
            this%cols(i) = 0
            if (tmp .ne. 0) this%nbrbycol(tmp) = 0
        endif
    enddo
#ifndef NDEBUG
    ! let's check the consistency
    do i = 1, size(this%nbrbycol)
    if (this%nbrbycol(i) .ne. 0) then
        if(this%cols(this%nbrbycol(i)) .ne. i ) then
        write (*,*) "inconsistency!"
        write (*,*) this%cols
        write (*,*) this%nbrbycol
        write (*,*) this%nbrbycol(10000)
        endif
    endif
    enddo
#endif
end subroutine vertex_erase_edge_color

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
    integer :: i, up
    
    i = 1
    up = size(this%nbrs)
    do while ((this%nbrs(i) .ne. vert) .and. (i <= up ))
    i = i+1
    enddo
#ifndef NDEBUG
        write (*,*) "setting color", col, "to idx", i, "which is", vert
#endif
            this%cols(i) = col
            this%nbrbycol(col) = i
! #ifndef NDEBUG
!     ! let's check the consistency
!     do i = 1, size(this%nbrbycol)
!     if (this%nbrbycol(i) .ne. 0) then
!         if(this%cols(this%nbrbycol(i)) .ne. i ) then
!         write (*,*) "inconsistency!"
!         write (*,*) this%cols
!         write (*,*) this%nbrbycol
!         write (*,*) this%nbrbycol(10000)
!         endif
!     endif
!     enddo
! #endif
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
!> @param[in] maxcols The available maximum number colors
!> @return the smallest common available color
!--------------------------------------------------------------------
function find_common_free_color(v, w, maxcols) result(col)
    implicit none
    class(vertex), intent(in) :: v, w
    integer, intent(in) :: maxcols
    integer :: col, i

    col = 0
    i = 1
    do while (.not.((v%nbrbycol(i) == 0) .and. (w%nbrbycol(i) == 0)) .and. i <= maxcols)
        i = i + 1
    enddo
    if (i <= maxcols) col = i
!     do i = maxcols, 1, -1
!         if ((v%nbrbycol(i) == 0) .and. (w%nbrbycol(i) == 0)) col = i
!     enddo
end function find_common_free_color

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> This function returns the smallest available color
!> that is free at both vertices and tries to set that color if
!> is free.
!
!> @param[in] v one vertex that we consider
!> @param[in] w the other vertex that we consider
!> @param[in] maxcols The available maximum number colors
!> @return the smallest common available color
!--------------------------------------------------------------------
function find_and_try_to_set_common_free_color(i, j, verts, maxcols) result(col)
    implicit none
    type(vertex), dimension(:), allocatable, intent(inout) :: verts
    integer, intent(in) :: maxcols, i, j
    integer :: col, k, l

    col = 0
    k = 1
    do while (.not.((verts(i)%nbrbycol(k) == 0) .and. (verts(j)%nbrbycol(k) == 0)) .and. k <= maxcols)
        k = k + 1
    enddo
    if (k <= maxcols) then
        col = k
        l = binarysearch(verts(i)%nbrs, j)
        verts(i)%cols(l) = col
        verts(i)%nbrbycol(col) = l
        l = binarysearch(verts(j)%nbrs, i)
        verts(j)%cols(l) = col
        verts(j)%nbrbycol(col) = l
    endif
!     do i = maxcols, 1, -1
!         if ((v%nbrbycol(i) == 0) .and. (w%nbrbycol(i) == 0)) col = i
!     enddo
end function find_and_try_to_set_common_free_color

subroutine construct_and_invert_path(col, verts, start)
    type(Vertex), allocatable, dimension(:) :: verts
    integer, intent(in) :: col(2), start
    integer :: k, nbr1, ver, colctr, curver, nbr, tmpcol
    logical :: stoppath
    Type(Path) :: p

#ifndef NDEBUG
    write (*,*) "construct path of colors", col(2), col(1), "at", i, fan(fanlen)
#endif
    call p%init()
    do k = 1, verts(start)%degree
        if (verts(start)%cols(k) == col(1)) nbr1 = k
    enddo
    ! We start the path at verts(i)
    ! (i, nbr1) is the first piece of the path. nbr1 is the position within the arrays of verts(i)
    ! The edge i -> nbr1 has color col(1)
#ifndef NDEBUG
    write (*,*) "starting path to ", nbr1, "which is", verts(i)%nbrs(nbr1)
#endif
    call p%pushback(start, nbr1)
    stoppath = .false.
    colctr = 2
    do while(stoppath .eqv. .false.)
        ! current vertex
        call p%back(ver, nbr)
        curver = verts(ver)%nbrs(nbr)
        ! find nbr of current color
        nbr1 = 0
        do k = 1, verts(curver)%degree
            if (verts(curver)%cols(k) == col(colctr)) nbr1 = k
        enddo
        if (nbr1 == 0) then
            ! presumably we have reached the end of the path
            stoppath = .true.
        else
            call p%pushback(curver, nbr1)
            colctr = colctr + 1
            if (colctr > 2) colctr = 1
        endif
    enddo
#ifndef NDEBUG
    write (*, *) "Length of path", p%length()
#endif
    ! Now we switch the colors on the path
    k = 1
    colctr = 2
    ! erase colors first
    do while(k<=p%length())
        call p%at(k, ver, nbr)
        tmpcol = verts(ver)%cols(nbr)
        verts(ver)%nbrbycol(tmpcol) = 0 ! to be on the safe side and create consistent data
                        
        call verts(verts(ver)%nbrs(nbr))%erase_edge_color(ver)
        k = k + 1
    enddo
    ! Now we perform the actual setting of the new path colors
    k = 1                   
    do while(k<=p%length())
        call p%at(k, ver, nbr)
!       tmpcol = verts(ver)%cols(nbr)
        verts(ver)%cols(nbr) = col(colctr)
        verts(ver)%nbrbycol(col(colctr)) = nbr
!       verts(ver)%nbrbycol(tmpcol) = 0 ! to be on the safe side and create consistent data
!       call verts(verts(ver)%nbrs(nbr))%erase_edge_color(ver)
        call verts(verts(ver)%nbrs(nbr))%set_edge_color(ver, col(colctr))
        k = k + 1
        colctr = colctr + 1
        if (colctr > 2) colctr = 1
    enddo
    call p%dealloc()
end subroutine

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> This function performs a downshift of a fan that was previously
!> constructed around fanpos. After that it sets the new available edge 
!> to the color col.
!
!> @param[in] fanpos the position where the vertex was constructed
!> @param[in] fan the actual length of the Vizing fan
!> @param[in] fannbr the positions of the neighbours in the array associated with verts(fanpos)
!> @param[in] fanedgecol the colors of the respective edges
!> @param[in] tail the length/end of the fan
!> @param[in] col the new color
!> @param[inout] the array of vertices.
!--------------------------------------------------------------------
subroutine downshift_fan_and_set_color(fanpos, fan, fannbr, fanedgecol, tail, col, verts)
    integer, allocatable, dimension(:), intent(in) :: fan, fannbr, fanedgecol
    type(Vertex), allocatable, dimension(:) :: verts
    integer, intent(in) :: tail, fanpos, col
    integer :: k, oldpos

    ! erase colors first that are about to become invalidated    
    do k = 2, tail
        ! erase current color at fanpos
        verts(fanpos)%cols(fannbr(k)) = 0
        verts(fanpos)%nbrbycol(fanedgecol(k)) = 0
    enddo
    
    do k = 2, tail-1
        ! set to new color
        verts(fanpos)%cols(fannbr(k)) = fanedgecol(k+1)
        verts(fanpos)%nbrbycol(fanedgecol(k+1)) = fannbr(k)
        
        ! set current color from fan(k) to fanpos. erase old entres first
        oldpos = verts(fan(k))%nbrbycol(fanedgecol(k))
        verts(fan(k))%nbrbycol(fanedgecol(k)) = 0

        verts(fan(k))%cols(oldpos) = fanedgecol(k+1)
!         
         verts(fan(k))%nbrbycol(fanedgecol(k+1)) = oldpos
    enddo
    
    verts(fanpos)%cols(fannbr(1)) = fanedgecol(2)
    verts(fanpos)%nbrbycol(fanedgecol(2)) = fannbr(1)
    call verts(fan(1))%set_edge_color(fanpos, fanedgecol(2))
    
    ! erase current color at fan(tail)
    verts(fan(tail))%cols( verts(fan(tail))%nbrbycol(fanedgecol(tail)) ) = 0
    verts(fan(tail))%nbrbycol(fanedgecol(tail)) = 0
    ! set to new color
    verts(fanpos)%cols(fannbr(tail)) = col
    verts(fanpos)%nbrbycol(col) = fannbr(tail)
    call verts(fan(tail))%set_edge_color(fanpos, col)
end subroutine

recursive subroutine quicksort(a, first, last)
  implicit none
  integer, dimension(:), intent(inout) :: a
  integer :: x, t
  integer :: first, last
  integer :: i, j

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     i=i+1
     j=j-1
  end do
  if (first < i-1) call quicksort(a, first, i-1)
  if (j+1 < last)  call quicksort(a, j+1, last)
end subroutine quicksort

function binarySearch (a, value)
    integer                  :: binarySearch
    integer, intent(in), target :: a(:)
    integer, intent(in)         :: value
    integer, pointer            :: p(:)
    integer                  :: mid, offset
 
    p => a
    binarySearch = 0
    offset = 0
    do while (size(p) > 0)
        mid = size(p)/2 + 1
        if (p(mid) > value) then
            p => p(:mid-1)
        else if (p(mid) < value) then
            offset = offset + mid
            p => p(mid+1:)
        else
            binarySearch = offset + mid    ! SUCCESS!!
            return
        end if
    end do
end function binarySearch


end module vertex_mod
