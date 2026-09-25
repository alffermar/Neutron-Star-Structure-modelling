program TOVsolver
  use precision, only: wp
  use Constants
  use funciones
  use NumMethods
  implicit none
  real(wp), dimension(0:puntos-1) :: presion,radio,masa,eninter,printer,nbs
  real(wp) :: h, conv1, msolar,dummy
  integer :: i,ilast, errores, totalerror

  ! Abrir el archivo de datos
  open(2, file='EoS.dat')
  

  ! Leer los datos del archivo

do i=0,puntos-1
   read(2, *) nbs(i) ,eninter(i), printer(i) !Tanto la energía como la presión van de la mayor a la menor en orden
enddo
do i=0,puntos-1
  write(*,*) eninter(i), printer(i)
enddo
close(2)

  eninter=real(eninter, wp)
  printer=real(printer, wp)

  !do i=0,pinter-1
  !   if (((printer(i+1)-printer(i))/(eninter(i+1)-eninter(i)).ge. 1.d0).or.&
  !   ((printer(i+1)-printer(i))/(eninter(i+1)-eninter(i)).le. 0.d0)) then
  !    print*, 'Algo salio mal', (printer(i+1)-printer(i))/(eninter(i+1)-eninter(i))
  !   end if
  !enddo
  !Factores de conversión
  conv1=qe*(real(1.d-3,wp)*real(1.d0/cl,wp))**2*&
  (real(1.d-3,wp))**3*(real(1.d0/cl,wp))**2*G*real(1.d54,wp) !(km^(-2))*(fm^3/MeV) [from Mev/fm^3 to km^(-2)]

  msolar = G*real(1.d-3,wp)**3*real(1.d0/cl,wp)**2*solarmas !km*M_sun^-1
 ! Inicializar condiciones
  h =(radiomax - radiomin) / (real(puntos,wp))
  printer=printer*real(conv1,wp)
  eninter=eninter*real(conv1,wp)
  
  open(2,file='Radio-Masa.dat')

  do i=0,puntos-2
  presion(0) = printer(i)       ! Seleccionamos una presión

  masa(0) = (4.0_wp/3.0_wp) * pi * (radiomin**3) * eninter(i)   ! Masa inicial con la energía de esa presión
  radio(0) = radiomin           ! Radio inicial
  ! Llamada a la subrutina RK4
  call DES_RK4(radio,presion,masa,eninter,printer,h,ilast,errores) !DES_RK4(radio,presion,masa,denenergia,pEoS,h,ilast)
  if((radio(ilast)<24.5_wp)) then
  write(*,*) radio(ilast), masa(ilast)/msolar,presion(ilast)/conv1, ilast
  write(2,*) radio(ilast), masa(ilast)/msolar,presion(ilast)/conv1, ilast
  end if
    totalerror=totalerror+errores
  enddo
  print*, 'falladas',totalerror
  close(2)


end program TOVsolver
