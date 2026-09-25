module funciones
      use precision, only : wp  
      use Constants
      implicit none
! x array has dynamical memory in order to simplify the addition of equations 
contains
! Usefull Functions

real(wp) function meff(part,x) !effective mass
      implicit none
      real(wp), dimension(:) :: X
      type(particle) :: part
      meff=part%mass-(abs(part%I3)*2.d0)*(gsigma*x(sigmameson%i)) ! (abs(part%I3)*2.d0)* is used to don't include leptons
      return
end function meff

real(wp) function mu(part,x) !chemical potential
      implicit none
      real(wp), dimension(:) :: x
      type(particle) :: part
      if (part%i.ne.8) then
            mu=sqrt(x(part%i)**2+meff(part,x)**2+(abs(part%I3)*2)*gpi*&
      (2.d0*x(piminus%i)**2))&
      +(abs(part%I3)*2.d0)*(gomega*x(omegameson%i)+grho*part%I3*x(rhozero%i)) ! (abs(part%I3)*2.d0)* is used to don't include leptons
      else 
            if (sqrt(x(e%i)**2+e%mass**2).ge.piminus%mass) then
                  mu=sqrt(x(e%i)**2+e%mass**2)
            else
                  mu=0
            end if
      endif
      return
end function mu

real(wp) function dens(part,x)  !density of the specie 
      implicit none
      real(wp), dimension(:) :: x
      type(particle) :: part
      dens=(x(part%i)**3)/(3.d0*pi**2)
      return
end function dens

real(wp) function I_I(part,X) !integral that appears in the equation for the sigma meson
      implicit none
      real(wp), dimension(:) :: X
      type(particle) :: part
      real(wp) :: k,m,piones
      k=X(part%i)
      m=meff(part,X)
      piones=gpi**2*(2.d0*x(piminus%i)**2)
      !!! Mathematica Integral
      I_I=0.5d0*m*(k*sqrt(k**2+m**2+piones)-(m**2+piones)*atanh((k)/(sqrt(k**2+m**2+piones))))
      !!! Glendenning Integral
      !I_I=0.5_wp*m*(k*sqrt(k**2+m**2)-m**2*asinh(k/m))
      return
end function I_I

real(wp) function I_II(part,X) !integral that appears in the equation for the energydensity and preasure
      implicit none
      real(wp), dimension(:) :: X
      type(particle) :: part
      real(wp) :: k,m,piones
      k=X(part%i)
      m=meff(part,X)
      piones=(abs(part%I3)*2.d0)*gpi**2*(2.d0*x(piminus%i)**2)
      !!! Mathematica Intergral
      I_II=1.d0/8.d0*((k*sqrt(k**2+m**2+piones))*(2.d0*k**2+m**2+piones)-(m**2+piones)**2*&
      atanh((k)/(sqrt(k**2+m**2+piones))))
      !!! Glendenning Integral
      !I_II=(1.d0/4.d0)*(k*(k**2+m**2)**1.5d0 -0.5d0*&
      !m**2*k*sqrt(k**2+m**2)-0.5d0*m**4*log((sqrt(k**2+m**2)+k)/(m))) 
      return
end function I_II

real(wp) function I_III(part,X) !integral that appears in the equation for the preasure
      implicit none
      real(wp), dimension(:) :: X
      type(particle) :: part
      real(wp) :: k,m,piones
      k=X(part%i)
      m=meff(part,X)
      piones=(abs(part%I3)*2.d0)*gpi**2*(2.d0*x(piminus%i)**2)
      !!! Mathematica Integral
      !I_III=(1.d0/8.d0)*(k*sqrt(k**2+m**2)*(2.d0*k**2-3.d0*m**2)+ &
       !   3.d0*m**4*atanh((k)/(sqrt(k**2+m**2))))
      !!! Glendenning integral
      !I_III=(1.d0/4.d0)*(k**3*sqrt(k**2+m**2)-1.5d0*&
      !m**2*k*sqrt(k**2+m**2)+1.5d0*m**4*log((sqrt(k**2+m**2)+k)/(m)))
      !!! WOLFRAM ALPHA
      I_III=1.d0/8.d0*(-1.d0*k*sqrt(piones+m**2+k**2)*(-2.d0*&
      k**2+3.d0*(m**2+piones))+3.d0*(m**2+piones)**2*asinh(k/sqrt(m**2+piones)))
      return
end function I_III

real(wp) function I_IV(part,X) !integral that appears in the equation for the pions meson
      implicit none
      real(wp), dimension(:) :: X
      type(particle) :: part
      if (((part%i).eq.9)) then
      I_IV=(I_I(p,x)/meff(p,x)-I_I(n,x)/meff(n,x))*x(part%i)
      else
            if (min(x(p%i),x(n%i))==x(p%i)) then
                  !I_IV=I_I(p,x)!/meff(p,x)*x(part%i)
                  I_IV=meff(p,x)/(2._wp*sqrt(x(p%i)**2+meff(p,x)**2)&
                  *sqrt(x(n%i)**2+meff(n,x)**2))*(dens(p,x)-dens(n,x))
            else
                  !I_IV=I_I(n,x)!/meff(n,x)*x(part%i)
                  I_IV=meff(p,x)/(2._wp*sqrt(x(p%i)**2+meff(p,x)**2)&
                  *sqrt(x(n%i)**2+meff(n,x)**2))*(dens(p,x)-dens(n,x))
                 
            endif
      end if      
      !!! Glendenning Integral
      !I_I=0.5_wp*m*(k*sqrt(k**2+m**2)-m**2*asinh(k/m))
      return
end function I_IV

subroutine population(x,popu)
      implicit none
      real(wp), dimension(:):: x
      real(wp) denstot
      real(wp), dimension(5) :: popu
      denstot=dens(p,x)+dens(n,x)
      popu(1)=dens(p,x)/denstot
      popu(2)=dens(n,x)/denstot
      popu(3)=dens(e,x)/denstot
      popu(4)=dens(muon,x)/denstot
end subroutine
!-----------------------------------------------------------------------------------------------
! Declaring of system of equation
real(wp) function f1(x)  !Beta equilibrium
      implicit none 
      real(wp), dimension(:) :: x(:)
      f1= mu(n,x)-mu(p,x)-mu(e,x)
      return
end function f1

real(wp) function f2(x)  !Charge Neutrality
      implicit none
      real(wp), dimension(:) :: x(:)
      f2= dens(p,x)-dens(e,x)-dens(muon,x)
      return
end function f2

real(wp) function f3(x)  !Baryon Density
      implicit none
      real(wp), dimension(:) :: x(:)
      f3= nb-dens(p,x)-dens(n,x)
      return
end function f3

real(wp) function f4(x)  !Sigma Meson MFA
      implicit none
      real(wp), dimension(:) :: x(:)
      f4=(sigmameson%mass**2)*X(sigmameson%i)+b*(n%mass)*gsigma*(gsigma*X(sigmameson%i))**2 &
      +c*gsigma*(gsigma*X(sigmameson%i))**3.d0-(gsigma)/(pi**2)*((I_I(p,x))+(I_I(n,x)))
end function

real(wp) function f5(x)  !Omega Meson MFA
      implicit none
      real(wp), dimension(:) :: x(:)
      f5=x(omegameson%i)-(gomega)/(omegameson%mass**2.d0)*(dens(p,x)+dens(n,x))
      return
end function

real(wp) function f6(x)  !Rho Meson MFA
      implicit none
      real(wp), dimension(:) :: x(:)
      f6=x(rhozero%i)-(grho)/(2.d0*rhozero%mass**2.d0)*(dens(p,x)-dens(n,x))
      return
end function

real(wp) function f7(x) !Leptonic equilibrium
      implicit none
      real(wp), dimension(:) :: x(:)
      f7=mu(e,x)-mu(muon,x)
      return
end function

real(wp) function f8(x)  !Pion - MFA
      implicit none
      real(wp), dimension(:) :: x(:)
f8=(piminus%mass**2)*x(piminus%i)+sqrt(2._wp)*gpi*I_IV(piminus,x)/(pi**2)
      return
end function
!Equations of state
real(wp) function energydensity(x)
      implicit none
      real(wp), dimension(:) :: X

 energydensity=(b/3.d0)*(gsigma*x(sigmameson%i))**3&
               +(c/4.d0)*(gsigma*x(sigmameson%i))**4&
               +(1.d0/2.d0)*(sigmameson%mass*x(sigmameson%i))**2&
               +(1.d0/2.d0)*(omegameson%mass*x(omegameson%i))**2&
               +(1.d0/2.d0)*(rhozero%mass*x(rhozero%i))**2&
               +abs(2.d0*(1.d0/2.d0)*(piminus%mass*x(piminus%i))**2)&
               +1.d0/(pi**2)*(I_II(p,x)+I_II(n,x)+I_II(e,x)+I_II(muon,x))
      return
end function energydensity

real(wp) function preasure(x)
      implicit none !AJUSTAR 3 POLITROPOS PARA VER SI FUNCIONA
      real(wp), dimension(:) :: X
      preasure=-1.d0*(b/3.d0)*(gsigma*x(sigmameson%i))**3&
               -1.d0*(c/4.d0)*(gsigma*x(sigmameson%i))**4&
               -1.d0*(1.d0/2.d0)*(sigmameson%mass*x(sigmameson%i))**2&
               +(1.d0/2.d0)*(omegameson%mass*x(omegameson%i))**2&
               +(1.d0/2.d0)*(rhozero%mass*x(rhozero%i))**2&
               +abs(2.d0*(1.d0/2.d0)*(piminus%mass*x(piminus%i))**2)&
               +1.d0/(3.d0*pi**2)*(I_III(p,x)+I_III(n,x)+I_III(e,x)+I_III(muon,x))
      return
end function preasure

real(wp) function dpdr(radio,presion,masa,denenergia)
      implicit none
      real(wp) :: radio,presion,masa,denenergia
      dpdr=-1.d0*((presion+denenergia)*(masa+4._wp*pi*radio**3*presion))&
      /(radio*(radio-2*masa))
      return
end function dpdr

real(wp) function dmdr(radio,denenergia)
      implicit none
      real(wp) :: radio, denenergia
      dmdr=4.d0*pi*radio**2.d0*denenergia
    !  print*, 'radio:', radio, 'dmdr:', dmdr, 'denenergia:', denenergia
      return
end function
end module funciones