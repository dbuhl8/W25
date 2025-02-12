
program test
  implicit none

  real :: x 
  complex :: z = (2.0,0.0)

  x = real(sqrt(z))

  print *, x

end program test
