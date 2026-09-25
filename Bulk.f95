program TOVsolver
  use precision, only: wp
  use Constants
  use funciones
  use NumMethods
  implicit none
real(wp), dimension(0:puntos-1) ::eninter,printer,nbs
real(wp), dimension(0:657):: densidad, deriv2, presioninterpolada

  integer :: i

    open(2, file='EoS.dat')
  ! Leer los datos del archivo

do i=0,puntos-1
   read(2, *) nbs(i) ,eninter(i), printer(i) !Tanto la energía como la presión van de la mayor a la menor en orden
enddo
 close(2)
! Coeficiente de compresibilidadd y masa efectiva
  open(3,file='bulk_mass.dat')
      call FirstDerivate(nbs,printer,puntos-1,densidad,deriv2,657-1)
      do i=0,657-3
      if (((deriv2(i)-2._wp*(printer(i))/densidad(i))*9._wp)<9000._wp) then
            write(3,*) densidad(i), (deriv2(i)-2._wp*(printer(i))/densidad(i))*9._wp
      end if
      enddo
      close(3)
      print*, 'hecho'
  end program
