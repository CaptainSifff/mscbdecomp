module vertex_mod
    implicit none
    type :: SingleColExp
    contains
        procedure :: init => SingleColExp_init
    end type
    
    type :: FullExp
        integer :: nrofcols
    contains
        procedure :: init => FullExp_init
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

subroutine SingleColExp_init(this)
    class(SingleColExp) :: this
end subroutine SingleColExp_init

subroutine FullExp_init(this)
    class(FullExp) :: this
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
