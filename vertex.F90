module vertex_mod
    implicit none

        type :: simplenode
        integer :: x,y
        complex (kind=kind(0.d0)) :: axy
    end type simplenode

    type :: node
        integer :: x,y, col
        complex (kind=kind(0.d0)) :: axy
    end type node

    type :: SingleColExp
        integer :: nrofentries
        type(simplenode), dimension(:), allocatable :: nodes
    contains
        procedure :: init => SingleColExp_init
        procedure :: dealloc => SingleColExp_dealloc
        procedure :: mult => SingleColExp_mult
    end type

    type :: FullExp
        integer :: nrofcols
        type(SingleColExp), dimension(:), allocatable :: singleexps
    contains
        procedure :: init => FullExp_init
        procedure :: dealloc => FullExp_dealloc
        procedure :: mult => FullExp_mult
    end type FullExp
    
    type :: Vertex
        integer :: degree
        integer, allocatable, dimension(:) :: nbrs
        integer, allocatable, dimension(:) :: cols
    contains
        procedure :: init => vertex_init
        procedure :: destruct => vertex_destruct
        procedure :: any_color_available => vertex_any_color_available
        procedure :: set_edge_color => vertex_set_edge_color
    end type Vertex
contains

subroutine SingleColExp_mult(this, vec)
    class(SingleColExp) :: this
    complex(kind=kind(0.D0)), dimension(:) :: vec
    integer :: i, j
    complex(kind=kind(0.D0)) :: t1,t2
    do i = 1, this%nrofentries! for every matrix
        t1 = vec(this%nodes(i)%x)
        t2 = vec(this%nodes(i)%y)
        vec(this%nodes(i)%x) = cosh(this%nodes(i)%axy) * t1 + sinh(this%nodes(i)%axy)* t2
        vec(this%nodes(i)%y) = cosh(this%nodes(i)%axy) * t2 + sinh(this%nodes(i)%axy)* t1
    enddo
end subroutine SingleColExp_mult

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
        this%nodes(i)%axy = nodes(i)%axy
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

subroutine FullExp_mult(this, vec)
    class(FullExp) :: this
    complex(kind=kind(0.D0)), dimension(:) :: vec
    integer :: i
    do i = 1, this%nrofcols
        call this%singleexps(i)%mult(vec)
    enddo
end subroutine FullExp_mult

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
    do i = 1, usedcolors
    write (*,*) edgectr(i), nredges(i)
    enddo
    allocate(this%singleexps(usedcolors))
    do i = 1, usedcolors
        call this%singleexps(i)%init(simplenodes(i, :), nredges(i))
    enddo
    deallocate(simplenodes)
    deallocate(nredges, edgectr)
end subroutine FullExp_init

subroutine vertex_init(this, deg)
    class(vertex) :: this
    integer, intent(in) :: deg
    this%degree = deg
    allocate(this%nbrs(deg), this%cols(deg))
    this%cols = 0
end subroutine vertex_init

subroutine vertex_set_edge_color(this, vert, col)
    class(vertex) :: this
    integer, intent(in) :: vert, col
    integer :: i
    
    do i = 1, this%degree
        if(this%nbrs(i) == vert) this%cols(i) = col
    enddo
end subroutine vertex_set_edge_color

subroutine vertex_destruct(this)
    class(vertex) :: this
    deallocate(this%nbrs, this%cols)
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
    if (col == 0) write(*,*) "this should not happen!"
    deallocate(usedcols)
end function find_common_free_color

end module vertex_mod
