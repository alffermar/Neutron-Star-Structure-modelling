module NumMethods
      use precision, only : wp
      use constants
      use funciones
      implicit none

 !Jacobian step
      real(wp), parameter  :: J_step=1.d-6
 !Parameters of Broyden's Method
      integer, parameter :: ndim=8 !This adjust the number of particle (and equations) considered in the model
      integer, parameter :: MaxIter=int(3.d6)
      real(wp), parameter :: step=1.0d-2
      real(wp), parameter :: tol=1.d-18
 !Parameters of Runge-Kutta adaptative step size
      real(wp), parameter :: tolerance=1.d-5
!Basic numerical methods that are recurrent

contains
!AUXILIAR-------------------------------------------------------------------------------------------------------------------------
      !This subroutines are auxiliar in order to simplify different codes, it can be problematic if the calculus are massive
 subroutine Eval_Fun(x,F)
      implicit none
      real(wp), intent(in), dimension(:) :: x
      real(wp), intent(out), dimension(:) :: F
      F=[f1(x),f2(x),f3(x),f4(x),f5(x),f6(x),f7(x),f8(x)]
 end subroutine

 function Norm(x)      
      !This norm is squared
      real(wp), intent(in), dimension(:) :: x
      real(wp), dimension(size(x)) :: F
      real(wp) :: Norm
      call Eval_Fun(x,F)
      Norm=sum(F**2.d0)
 end function Norm

 subroutine ProdVec(u, v, R)
      implicit none
      real(wp), intent(in) :: u(:), v(:)
      real(wp) :: R(size(u), size(v))
      integer :: i, j
    
      do i = 1, size(u)
        do j = 1, size(v)
            R(i, j) = u(i) * v(j)
        end do
      end do
 end subroutine ProdVec

 subroutine Imprime(x)
      implicit none
      real(wp), dimension(:) :: X
      print*, '----------------------------------------------------------------------------------------------------------------'
      write(*,*) 'Densidad barionica', f3(x)
      write(*,*) 'Equilibrio Beta', f1(x)
      write(*,*) 'Neutralidad de carga:', f2(x)
      write(*,*) 'Potencial quimico leptonico:', f7(x)
      write(*,*) 'Fraccion de protones Yp por densidades:', dens(p,x)/(nb/(hbc/unit)**3)
      write(*,*) 'Cociente de masa efectivas:', meff(p,x)
      print*, '----------------------------------------------------------------------------------------------------------------'
 end subroutine
! DIFFERENTIAL CALCULUS-----------------------------------------------------------------------------------------------------------
 subroutine Jacobian(x,J)
      !I cannot use functions declared in the module as argument of subroutines so I have to write one by one. !!VER CON CALMA EN EL LIBRO VERDE
           implicit none
           integer :: i
           real(wp), intent(in),dimension(:) :: x
           real(wp), intent(out), dimension(ndim,ndim) :: J
           real(wp), dimension(size(x)) :: h
           !Remember fortran n x n matrix first possition means row and second colum
           !Richardson's numerical derivate formula with r=1/2
           do i=1,ndim
                 h(i)=J_step
                 J(1,i)=-(f1(x+h)-8.d0*f1(x+h/2.d0)-f1(x-h)+8.d0*f1(x-h/2.d0))/(6.d0*J_step)
                 J(2,i)=-(f2(x+h)-8.d0*f2(x+h/2.d0)-f2(x-h)+8.d0*f2(x-h/2.d0))/(6.d0*J_step)
                 J(3,i)=-(f3(x+h)-8.d0*f3(x+h/2.d0)-f3(x-h)+8.d0*f3(x-h/2.d0))/(6.d0*J_step) 
                 J(4,i)=-(f4(x+h)-8.d0*f4(x+h/2.d0)-f4(x-h)+8.d0*f4(x-h/2.d0))/(6.d0*J_step)
                 J(5,i)=-(f5(x+h)-8.d0*f5(x+h/2.d0)-f5(x-h)+8.d0*f5(x-h/2.d0))/(6.d0*J_step)
                 J(6,i)=-(f6(x+h)-8.d0*f6(x+h/2.d0)-f6(x-h)+8.d0*f6(x-h/2.d0))/(6.d0*J_step)
                 J(7,i)=-(f7(x+h)-8.d0*f7(x+h/2.d0)-f7(x-h)+8.d0*f7(x-h/2.d0))/(6.d0*J_step)
                 J(8,i)=-(f8(x+h)-8.d0*f8(x+h/2.d0)-f8(x-h)+8.d0*f8(x-h/2.d0))/(6.d0*J_step)
                 h=0.d0
           end do
 end subroutine Jacobian

!INTERPOLATION AND FITTERS -------------------------------------------------------------------------------------------------------
 subroutine coefs(n, f, a, b, c, d, h)
    implicit none
    integer n, i
    real(wp) ::  f(0:n), a(0:n), b(0:n), c(0:n), d(0:n)
    real(wp) :: alpha(0:n-1), beta(0:n-1), r(0:n-1), h(0:n)
    ! Copia los valores de f a d
   d=f
    ! Calcular los valores intermedios de beta y r
    do i = 1, n-1
        beta(i) = 2.0_wp * (h(i) + h(i-1))
        r(i) = 3.0_wp / h(i) * (d(i+1) - d(i)) - 3.0_wp / h(i-1) * (d(i) - d(i-1))
    !    print*, h(i) * (d(i+1) - d(i))
    end do
    ! Inicialización de alpha
    do i = 2, n-1
        alpha(i) = h(i-1)
    end do
    ! Eliminación hacia adelante
    do i = 2, n-1
        beta(i) = beta(i) - (h(i-1) * alpha(i)) / beta(i-1)
        
        r(i) = r(i) - (r(i-1) * alpha(i)) / beta(i-1)
        
    end do
        ! Definir condiciones de frontera
    b(0) = 0.0_wp
    b(n) = 0.0_wp
    ! Cálculo de b(n-1)
    b(n-1) = r(n-1) / beta(n-1)
    ! Sustitución hacia atrás para obtener b(i)
    do i = n-1, 1, -1
        b(i) = (r(i) - h(i) * b(i+1)) / beta(i)
    end do

    ! Cálculo de a(i) y c(i)
    do i = 0, n-1
        a(i) = (b(i+1) - b(i)) / (3.0_wp * h(i))
        c(i) = (d(i+1) - d(i)) / h(i) - (b(i+1) + 2.0_wp * b(i)) * h(i) / 3.0_wp
       ! print*, a(i),c(i)
    end do

    !  print*, 'FIN COEFS'
    return
 end subroutine coefs


 subroutine cubicSplines(x, y, xpoint, xint, yint, intpoint)
    implicit none
    integer :: xpoint, intpoint, i, j, m
    real(wp), dimension(0:xpoint-1) :: x, y
    real(wp), dimension(0:intpoint-1) :: xint, yint
    real(wp), dimension(0:xpoint) :: aa, bb, cc, dd
    real(wp), dimension(0:xpoint-1) :: h
    real(wp) :: dx
    
    ! Cálculo de h(i) antes de llamar a coefs
    do i = 0, xpoint-1
        h(i) = x(i+1) - x(i)
        !print*, h(i)
    end do
    h(xpoint) = 0.0_wp  ! Ajuste para evitar problemas de acceso a memoria
    ! Llamada a coefs para calcular los coeficientes
    !print*, xpoint
    call coefs(xpoint, y, aa, bb, cc, dd, h) 
    dx = (x(xpoint) - x(0)) / real(intpoint, wp)
    ! Interpolación
    do i = 0, intpoint
        xint(i) = x(0) + dx * real(i, wp)
        ! Búsqueda de j con salida eficiente
        ! Ajustar la búsqueda del intervalo para incluir el último valor
        do m = 0, xpoint-1
            if (xint(i) >= x(m) .and. xint(i) < x(m+1)) then
                j = m
            endif
        end do

        ! Evaluación del spline en el punto xint(i)
        yint(i) = aa(j) * (xint(i) - x(j))**3 + bb(j) * (xint(i) - x(j))**2 + & 
                  cc(j) * (xint(i) - x(j)) + dd(j)
      !  print*, xINT(i)/1.36d-6,yint(i)/1.36d-6, i
    end do
    print*, 'Fin Cubic'
    return
 end subroutine cubicSplines

 subroutine FirstDerivate(x, y, xpoint, xint, yint, intpoint)
    implicit none
    integer :: xpoint, intpoint, i, j, m
    real(wp), dimension(0:xpoint-1) :: x, y
    real(wp), dimension(0:intpoint-1) :: xint, yint
    real(wp), dimension(0:xpoint) :: aa, bb, cc, dd
    real(wp), dimension(0:xpoint-1) :: h
    real(wp) :: dx
    ! Cálculo de h(i) antes de llamar a coefs
    do i = 0, xpoint-2
        h(i) = x(i+1) - x(i)
    end do
    h(xpoint-1) = 0.0_wp  ! Ajuste para evitar problemas de acceso a memoria
    ! Llamada a coefs para calcular los coeficientes
    call coefs(xpoint-1, y, aa, bb, cc, dd, h) 
    
    dx = (x(xpoint) - x(0)) / real(intpoint, wp)
    ! Interpolación
    do i = 0, intpoint
        xint(i) = x(0) + dx * real(i, wp)
        ! Búsqueda de j con salida eficiente
        ! Ajustar la búsqueda del intervalo para incluir el último valor
        do m = 0, xpoint-1
            if (xint(i) >= x(m) .and. xint(i) < x(m+1)) then
                j = m
            endif
        end do

        ! Evaluación del spline en el punto xint(i)
        yint(i) = 3._wp*aa(j)*(xint(i) - x(j))**2 + 2._wp*bb(j)*(xint(i)- x(j))+cc(j)
      !  print*, xINT(i)/1.36d-6,yint(i)/1.36d-6, i
    end do
    print*, 'Fin segunda primera'
    return
 end subroutine FirstDerivate

  subroutine SecondDerivate(x, y, xpoint, xint, yint, intpoint)
    implicit none
    integer :: xpoint, intpoint, i, j, m
    real(wp), dimension(0:xpoint-1) :: x, y
    real(wp), dimension(0:intpoint-1) :: xint, yint
    real(wp), dimension(0:xpoint) :: aa, bb, cc, dd
    real(wp), dimension(0:xpoint-1) :: h
    real(wp) :: dx
    ! Cálculo de h(i) antes de llamar a coefs
    do i = 0, xpoint-2
        h(i) = x(i+1) - x(i)
    end do
    h(xpoint-1) = 0.0_wp  ! Ajuste para evitar problemas de acceso a memoria
    ! Llamada a coefs para calcular los coeficientes
    call coefs(xpoint-1, y, aa, bb, cc, dd, h) 
    
    dx = (x(xpoint) - x(0)) / real(intpoint, wp)
    ! Interpolación
    do i = 0, intpoint
        xint(i) = x(0) + dx * real(i, wp)
        ! Búsqueda de j con salida eficiente
        ! Ajustar la búsqueda del intervalo para incluir el último valor
        do m = 0, xpoint-1
            if (xint(i) >= x(m) .and. xint(i) < x(m+1)) then
                j = m
            endif
        end do

        ! Evaluación del spline en el punto xint(i)
        yint(i) = 6._wp*aa(j) * (xint(i) - x(j)) + 2._wp*bb(j)
      !  print*, xINT(i)/1.36d-6,yint(i)/1.36d-6, i
    end do
    print*, 'Fin segunda derivada'
    return
 end subroutine SecondDerivate
   
 subroutine interpolate_energy(p, denenergia, pEoS, en)
      implicit none
      real(wp) :: p
      real(wp), dimension(0:puntos-1) :: denenergia, pEoS,h
      real(wp), dimension(0:puntos) :: aa,bb,cc,dd
      real(wp), intent(out) :: en
      integer :: i,m !Para hacer la interpolación utilizamos el logaritmo ya que la función es demasiado grande
        do i=0,puntos-2
            h(i)=peos(i+1)-(peos(i))
        enddo
        h(puntos-1)=0._wp
        call coefs(puntos-1,denenergia,aa,bb,cc,dd,h)

        do i = 0, puntos-1
            if (((p >= peos(i))) .and. (p < peos(i+1))) then
                m=i
                print*, m
                exit
            endif
        end do

        en = aa(m)*(p-peos(m))**3+bb(m)*(p-peos(m))**2+& 
                  cc(m)*(p-peos(m))+dd(m)

    !    if (m==1989) print*, m

    !        en = -1.0_wp  ! Valor de salida predeterminado en caso de que p esté fuera del rango de pEoS
    
            ! Búsqueda de los índices de interpolación
          !  do i = 0, pinter-2  ! Hasta puntos-2, para evitar acceder a pEoS(i+1) fuera del rango 
          !      if (p >= pEoS(i) .and. p <= pEoS(i+1)) then  !Los operadores lógicos están así por ser en orden creciente, si fuera decreciente seria al reves
          !          ! Interpolación lineal      
          !          fraction = (p - pEoS(i)) / (pEoS(i+1) - pEoS(i))
          !          !print*, 'Pendiente', fraction
          !          en = denenergia(i) + fraction * (denenergia(i+1) - denenergia(i))
          !          !print*, en/1.3237d-6
          !          return
          !      end if
             
      
        return
 end subroutine interpolate_energy

 subroutine LeastSquare_LinearFit(x, y, n, a, b)
      implicit none
      integer, intent(in) :: n
      real(wp), dimension(1:n), intent(in) :: x, y
      real(wp), intent(out) :: a, b
      real(wp) ::  sumaX,sumaY,sumaXY,sumaX2 !internal variables
      integer :: i
      sumaX=0.d0
      sumaY=0.d0
      sumaXY=0.d0
      sumaX2=0.d0
      do i=1,n
           sumaX=sumaX+x(i) 
           sumaY=sumaY+y(i) 
           sumaXY=sumaXY+x(i)*y(i)
           sumaX2=sumaX2+x(i)**2
      enddo
      a=(real(n,wp)*sumaXY-sumaX*sumaY)/(real(n,wp)*sumaX2-sumaX**2)
      b=(sumaY-a*sumaX)/(real(n,wp))
      return
 end subroutine LeastSquare_LinearFit

! MATRIX ALGEBRA------------------------------------------------------------------------------------------------------------------
      !FORMAT_NUMBER  format(3(d10.3,3x)) !This format express the matrix in normal format
      !All subroutines are made in order to the first index will be row and second will be column
      !This format write txt archives
      !     do j=1,ndim
      !           write(chanel,fortmat) (MATRIX(j,i), i=1,ndim)
      !     end do
 subroutine Matrix_FacLU(a,L,U) 
      implicit none
      real(wp) :: suma
      real(wp), intent(in), dimension(:,:) :: a
      real(wp), intent(out), dimension(:,:) :: L, U
      integer :: i, j, k
  ! Initializing matrices      
      L = 0.d0
      U = 0.d0

  ! Setting up the first row of U and the diagonal of L
      do i = 1, size(a, 1)
            L(i, i) = 1.d0
            U(1, i) = a(1, i)
      end do

  ! LU determining      
      do k=2,size(a,1)
            do i=1,size(a,1)
              suma=0.d0
               if(i.le.k-1) then
                do j=1,i-1
                suma=suma+l(k,j)*u(j,i)
                end do
                 l(k,i)=(a(k,i)-suma)/(u(i,i))
               else if(i.ge.k) then
                do j=1,k-1
                suma=suma+l(k,j)*u(j,i)
                end do
                u(k,i)=a(k,i)-suma
               end if
            end do
      end do
 end subroutine Matrix_FacLU

 subroutine Matrix_Prod(A, B, C)  !CHECK CREO QUE LA PODEMOS QUITAR
      implicit none
      real(wp), intent(in), dimension(:,:) :: A, B
      real(wp), intent(out), dimension(size(a,1), size(b,2)) :: C
      integer :: i, j, k

      ! Initialazing matrix 
      C = 0.d0

      ! Matrix product
      do i = 1, size(a,1)
          do j = 1, size(b,2)
              do k = 1, size(a,2)
                  C(i, j) = C(i, j) + A(i, k) * B(k, j)
              end do
          end do
      end do
 end subroutine Matrix_Prod

 subroutine Matrix_Inv(A,InvA) 
      implicit none
      real(wp), intent(in), dimension(:,:) :: a
      real(wp), intent(out), dimension(:,:) :: invA
      real(wp) ::  det, suma 
      integer :: i, j, n,k
      real(wp),dimension(size(a,1),size(a,1)) :: invL, L, U
      real(wp),dimension(size(a,1)) :: sinv, var, sol
      n=size(a,1)
  ! Initialazing inverse matrix
      invL=0.d0
      invA=0.d0

      call Matrix_FacLU(A,L,U)
      det=1.d0
      do i=1,N
            det=det*U(i,i)
      end do
  ! Inverse condition, if det=0 subroutine stops
      if (det.eq.0) then
            print*, det
            write(*,*) 'Inverse matrix does not extist'
            stop 
      end if

      do i=1,n
            do j=1,n !building n linear system
                if (j.eq.i) then
                 sinv(j)=1.d0
                 var(j)=sinv(j)
                else
                  sinv(j)=0.d0
                  var(j)=sinv(j)
                end if
            end do 
    
  ! Forward substitution to solve L * y = sinv for y
            do j = 1, n
                  suma = 0.d0
                  do k = 1, j - 1
                      suma = suma + L(j, k) * var(k)
                  end do
                  var(j) = (sinv(j) - suma) / L(j, j)
                  sol(j) = var(j)
                  invL(j, i) = var(j)
              end do
  
              ! Backward substitution to solve U * x = y for x
              do j = n, 1, -1
                  suma = 0.d0
                  do k = j + 1, n
                      suma = suma + U(j, k) * sol(k)
                  end do
                  sol(j) = (var(j) - suma) / U(j, j)
                  invA(j, i) = sol(j)
              end do
        end do
  
        return
 end subroutine Matrix_Inv

 function Determinant(A) !CHECK
      implicit none
      real(wp), dimension(:,:) :: A
      real(wp) :: Determinant
      real(wp), dimension(size(a,1),size(a,1)) :: L, U
      integer :: i
      call Matrix_FacLU(a,l,U)
      Determinant=1.d0
      do i=1,size(A,1)
            Determinant=Determinant*U(i,i)
      end do
 end function

!SYSTEM SOLVER--------------------------------------------------------------------------------------------------------------------
      ! ITERATIVE METHODS
 subroutine SystemSolver_GaussSeidel(a,b,x,x0)  !hay que ajustarlo y poner el condicional
      implicit none
      real(wp) suma1,suma2
      integer n,i,j
      real(wp),intent(inout) ,dimension(:) ::  b, x, x0
      real(wp), intent(in),dimension(:,:) :: a
      n=size(b)
      do i=1,n
        suma1=0.d0
        suma2=0.d0
        do j=1,i-1
        suma1=suma1+a(i,j)*x(j)
        end do
        do j=i+1,n
          suma2=suma2+a(i,j)*x0(j)
        end do
        x(i)=(b(i)-suma1-suma2)/a(i,i)
      end do
      return
 end subroutine SystemSolver_GaussSeidel
      ! DIRECT METHODS
 subroutine SystemLinearSolver_Direct(A,b,sol) !CHECK
      !this method is not usefull to solve big linear systems
      implicit none
      real(wp), intent(in), dimension(:,:) :: A
      real(wp), intent(in), dimension(:) :: b
      real(wp), intent(out), dimension(size(A,1)) :: sol
      real(wp), dimension(size(A,1) ,size(A,1)) :: L, U 
      real(wp), dimension(size(A,1)) :: var, solL
      real(wp) :: suma
      integer :: i,j

      call Matrix_FacLU(A,L,U)
   !Forward substitution to solve L * solL =b
      do i=1,size(A,1)
            suma=0.d0
            var(1)=b(1)
            if (i.eq.1) then 
            else
                  do j=1, i-1
                  suma=suma+L(i,j)*var(j)
                  end do
            end if 
            var(i)=(b(i)-suma)/(L(i,i))
            solL(i)=var(i)
      end do
      
  !Backward subtitution to solve U *sol = solL
      do i=size(A,1),1,-1
      suma=0.d0
            do j=i+1, size(A,1)
                  suma=suma+U(i,j) *solL(j)
            end do
            solL(i)=(solL(i)-suma)/U(i,i)
            sol(i)=soll(i)
      end do
 end subroutine SystemLinearSolver_Direct
      ! NON LINEAR SYSTEM
            
 subroutine SystemSolver_Broyden(x0,x) ! CHECK, x vector is the solution vector
      implicit none
      real(wp), intent(in), dimension(ndim) :: x0 !This ndim can be replaced by : and the rest by size(x0)
      real(wp), intent(out), dimension(ndim) :: x
      real(wp), dimension(ndim,ndim) :: B, Jac,Num, j0
      real(wp), dimension(ndim) :: F, dx, Fold,dF,xold, auxdx
      real(wp) ::  denominator
      integer :: i,j,k
      common j0
      i=0!We start in 0 
      if (Determinant(j0).eq.0.d0) then
            write(*,*) 'Jacobian has not inverse'
            stop
      end if
  !We start with xold and J0 (invJ0) in order to compute invvJ1 and x1
  !First iterattion 
      xold=x0
      call Eval_Fun(xold,Fold)
      call Matrix_Inv(J0,B)
      if (Determinant(B).eq.0.d0) then
            print*, B
            write(*,*) 'Jacobian has not inverse'
            stop
      end if
      x=xold-step*matmul(B,Fold)
      do j=1,5 !This is a extra condition because momentum values must be positive except for rho messon
            x(j)=abs(x(j))
      end do
            x(6)=-1.d0*abs(x(6))
            x(7)=abs(x(7))
      !Condition to low densities
       if (nb/(hbc/unit)**3.le.0.1d0) then
       x(3)=x(1)
       x(7)=1.d-29
       else 
       endif
      call Eval_Fun(x,F)
      dx=x-xold
      df=f-Fold
      Jac=J0
  !Next iterattions 
      do while (i.lt.MaxIter.and.(norm(x).gt.tol))
  !            call ProdVec((dx-matmul(B,df)),dx,Num)
                  auxdx=dx-matmul(B,df) !!!!ENTRE ESTE PUNTO
                  do k = 1, ndim
                        do j = 1, ndim
                            Num(k, j) = auxdx(k) * dx(j)
                        end do
                  end do
                  
                  denominator=dot_product(matmul(dx,B),df)
                  B=B+matmul(num,B)/(denominator) !!!!! Y ESTE PUNTO NO TOCAR NADA, SON LA EXPRESIONES BUENAS PARA LA INVERSA
            
            xold=x
            x=xold-step*matmul(B,F)
            do j=1,5 !This is a extra condition because momentum values must be positive except for rho messon
                  x(j)=abs(x(j))
            end do
                  x(7)=abs(x(7))
            !!Condition to low densities
             if (nb/(hbc/unit)**3.le.0.001d0) then
             x(3)=x(1)
             x(7)=1.d-29
             else 
             endif
            ! if (nb/(hbc/unit)**3.le.0.09) then
            !      x(2)=(nb*(pi**2)*3.)**(1./3.)
            ! end if
             fold=f ! EQUALS TO CALL EVAL_FUN(XOLD,FOLD)
             call Eval_Fun(x,F)
             dx=x-xold
             df=F-fold
             i=i+1
      end do
      
      print*, 'Norma', norm(x), 'Iteraciones', i, 'Determinante B:', Determinant(B)
    
  !    55       format(4(d20.12,3x))
      return
 end subroutine SystemSolver_Broyden

!DIFFERENTIAL EQUATION SOLVER-----------------------------------------------------------------------------------------------------

  subroutine DES_RK4(radio, presion, masa, denenergia, pEoS, h, ilast, errores)
    implicit none
    real(wp), intent(inout), dimension(0:puntos-1) :: radio, presion, masa
    real(wp), intent(in), dimension(0:puntos-1) :: denenergia, pEoS
    real(wp), intent(in) :: h
    real(wp) :: k11, k12, k21, k22, k31, k32, k41, k42, en
    integer :: i, ilast, errores
    
    errores = 0
    i = 1
    
    do while (presion(i-1) > 0.d0 .and. i < puntos-1)
        call interpolate_energy(presion(i-1), denenergia, pEoS, en)
        
        ! Coeficientes Runge-Kutta 4
        k11 = dpdr(radio(i-1), presion(i-1), masa(i-1), en)
        k12 = dmdr(radio(i-1), en)
        k21 = dpdr(radio(i-1) + h/2.d0, presion(i-1) + h * k11 / 2.d0, masa(i-1) + h * k12 / 2.d0, en)
        k22 = dmdr(radio(i-1) + h/2.d0, en)
        k31 = dpdr(radio(i-1) + h/2.d0, presion(i-1) + h * k21 / 2.d0, masa(i-1) + h * k22 / 2.d0, en)
        k32 = dmdr(radio(i-1) + h/2.d0, en)
        k41 = dpdr(radio(i-1) + h, presion(i-1) + h * k31, masa(i-1) + h * k32, en)
        k42 = dmdr(radio(i-1) + h, en)
        
        ! Actualizar variables
        radio(i) = radiomin + h*real(i,wp)
        presion(i) = presion(i-1) + h * (k11 + 2.d0*k21 + 2.d0*k31 + k41) / 6.d0
        masa(i) = masa(i-1) + h * (k12 + 2.d0*k22 + 2.d0*k32 + k42) / 6.d0
        ilast = i-1
        i = i + 1
        
    end do
    return
end subroutine DES_RK4 

end module NumMethods
