subroutine fmin(a, b, c)
    implicit none
    real(8), intent(in) :: a, b
    real(8), intent(out) :: c

    c = min(a,b)
end subroutine fmin