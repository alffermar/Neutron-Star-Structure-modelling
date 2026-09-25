      module Constants
     use precision, only : wp
!Here I'm going to introduce all phyisical constants and parameters that I will need in my program            
            implicit none

!Program parameters
      integer, parameter :: pnormaldensity=423
      integer, parameter :: plowdensity=234
      integer, parameter :: puntos= pnormaldensity+plowdensity !=Points of Broyden's Methods
 !TOV equations
      !Radius in kilometers
      real(wp), parameter :: radiomin=1.3_wp !
      real(wp), parameter :: radiomax=25.0_wp
!Physics            
      real(wp) , parameter :: hbc=197.33_wp ! In Mev*fm
      real(wp) , parameter :: unit=939.0_wp ! In Mev Nucleon Mass
      real(wp) , parameter :: G=6.67430_wp *real(1.d-11,wp) !In m^3/(s^2kg)
      real(wp) , parameter :: cl=299792.458_wp !In (km/s)
      real(wp) , parameter :: qe=1.602176634_wp*real(1.d-13,wp) !In (Mega C)
      real(wp) , parameter :: solarmas=1.988416_wp*real(1.d30,wp) !In (kg)

      real(wp) :: nb=0.0_wp
      ! Useful properties of all species of particles
      type particle
            real(wp) :: mass,I3,q !Mass in Gev, I3 third isospin component, q charge of the particle in units of e
            integer i    ! It labels the particle in order to use it in vectors
      end type particle
      !Data from particle data group cite REFERENCIAR !We are working adimensionaly
      type(particle) ::p = particle(939.0_wp/unit,0.5_wp,1.0_wp,1)         
      type(particle) ::n = particle(939.0_wp/unit,-0.5_wp,0_wp,2)         
      type(particle) ::e = particle(0.511_wp/unit,0_wp,-1.0_wp,3)         
      type(particle) ::sigmameson = particle(500.0_wp/unit,0._wp,0._wp,4)         
      type(particle) ::omegameson = particle(782.66_wp/unit,0._wp,0._wp,5)    
      type(particle) ::rhozero = particle(775.26_wp/unit,0.0_wp,0._wp,6)   
      type(particle) :: muon=particle(105.66_wp/unit,0._wp,-1.0_wp,7)
      type(particle) :: piminus=particle(140_wp/unit, -1_wp, -1_wp, 8)



      !Coupling constant of the model sigma omega rho TFM MANU, preguntar si hay que citar
      real(wp), parameter::gsigma=8.15458_wp !10.327982_wp
      real(wp), parameter::gomega=9.46597_wp !13.18631_wp
      real(wp), parameter::grho=8.47736_wp !11.1596_wp
      real(wp), parameter::gpi=12.4_wp !unit/92.4_wp  masa del nucleon dividido constante de decaimiento del pion
      !real(wp), parameter :: gsigma=8.7440_wp 
      !real(wp), parameter :: gomega=9.9893_wp
      !real(wp), parameter :: grho=8.7815_wp
      real(wp), parameter ::b=0.006275_wp !sigma cubic autointeraction. 
      real(wp), parameter ::c=-0.003409_wp

      !Pions Condensate Paper PRC
      type model
            real(wp) ::  b0,b1,bre,alfa !All coefficients are in pion mass units I have to divide by pion mass to obtain the adimensional
            character(4) :: seleccion
            integer :: i
      end type model

      type(model):: T=model(-0.034_wp,-0.078_wp,0.0_wp,0.0_wp,'T',1)
      type(model):: BFG=model(-0.025_wp,-0.085_wp,-0.021_wp,0.0_wp,'BFG',2)
      type(model):: SM=model(-0.027_wp,-0.12_wp,0.0_wp,0.0_wp,'SM',3)
      type(model):: ET=model(-0.02_wp,-0.0873_wp,-0.049_wp,0.0_wp,'ET',4)
      type(model):: NOG=model(-0.013_wp,-0.105_wp,0.0_wp,0.0_wp,'NOG',5)
      type(model):: KY=model(-0.0233_wp,-0.1473_wp,-0.019_wp,0.367_wp,'KY',6)
      !These constants are determined by scattering
      type(model):: FC=model(-0.009_wp,-0.114_wp,0.04_wp,0.0_wp,'FC',7)
      type(model):: FW=model(-0.009_wp,-0.081_wp,0.04_wp,0.391_wp,'FW',8)
!Mathematics
            real(wp), parameter :: pi=atan(1.0)*4.d0
      contains
            
      end module Constants