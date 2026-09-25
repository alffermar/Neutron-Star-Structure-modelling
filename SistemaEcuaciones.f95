program SistemaEcuaciones
      use precision, only : wp 
      use constants 
      use funciones 
      use NumMethods
      implicit none
      integer :: i,j
      real(wp), dimension(ndim) :: x,sol
      real(wp), dimension(5) :: popu
      real(wp), dimension(pnormaldensity) :: nbs, enden, pres
      real(wp) ,dimension(ndim,ndim) :: j0
      real(wp) :: dummy
      real(wp), dimension(plowdensity) :: lownb, lowpres, lowener
      
      common J0
      
! x order: proton neutron electron sigma ||||||| omega rho
      open(43, file='EstrellaMixta.table')
      open(44, file='momentos.dat')
      open(45, file='poblaciones.dat')
      open(46, file='Energydensity.dat')
      open(47, file='Preasure.dat')
      open(48, file='EoS.dat')
      open(49, file='masas.dat')
      x=0.02d0
     
      do i=1,plowdensity
            read(43,*) dummy, lownb(i), dummy, lowpres(i), lowener(i)

            write(46,*)  lownb(i), lowener(i)
            write(47,*) lownb(i),lowpres(i)
            write(48,*) lownb(i), lowener(i), lowpres(i)
      enddo
      do i=1,pnormaldensity
            nbs(i)=0.135d0+real(i-1)*(1.5d0-0.135d0)/(real(pnormaldensity-1))
      end do      !nb     p     n    e    sma     omga   rho   mu   pion-  pion0   pion+
      write(44,*) 1.0e-7_wp,1.0e-7_wp,1.0e-7_wp,1.0e-7_wp,1.0e-7_wp,1.0e-7_wp,-1.0e-7_wp,1.0e-10_wp,-1.0e-7_wp !Imponemos momentos nulos a densidad 0
      write(45,*) 1.0e-7_wp,1.0e-7_wp,1.e0_wp,1.0e-7_wp,1.0e-10_wp,1.0e-10_wp                                  !Imponemos pure neutron matter a densidad 0
      call Jacobian(x,J0)
! Computing for each baryon density            
      do i=1,pnormaldensity
            nb=nbs(i)*(hbc/unit)**3
            call SystemSolver_Broyden(x,sol)
            call population(sol,popu)
            do j=1,ndim
                  x(j)=sol(j)
            end do
            write(*,*) 'Poblaciones', popu
            !call Imprime(x)
            write(*,*) '  Densidad barionica        Densidad de Energia            Presion'
            write(44,*) nb/(hbc/unit)**3, sol*(unit/hbc)
            write(45,*) nb/(hbc/unit)**3, popu
            !Energy density and preassure                  
            enden(i)=energydensity(sol)
            pres(i)=preasure(sol)
            write(46,*) nbs(i), energydensity(sol)*unit**4/hbc**3
            write(47,*) nbs(i), preasure(sol)*unit**4/hbc**3
            write(48,*) nbs(i),energydensity(sol)*unit**4/hbc**3, &
            preasure(sol)*unit**4/hbc**3
            write(49,*) nbs(i),meff(p,sol)
            write(*,*) nbs(i),energydensity(sol)*unit**4/hbc**3, &
                        preasure(sol)*unit**4/hbc**3
            write(*,*) '------------------------------------------------------------------------------'
      end do         
      print*, 'Hecho'

!Closing previous archives
      close(43)
      close(44) 
      close(45)
      close(46)
      close(47)
      close(48)      
end program SistemaEcuaciones