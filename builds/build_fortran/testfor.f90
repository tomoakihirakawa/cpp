module testfor
    use, intrinsic :: iso_c_binding
    implicit none
contains
    subroutine testingfor(x) bind(c, name="testingfor")
        implicit none
        double precision, intent(in) :: x
 
        print *, "# calling testfor in Fortran!!"
        print *, "x = ", x 
    end subroutine
end module testfor
