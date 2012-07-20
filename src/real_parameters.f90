MODULE real_parameters

USE precision,            ONLY : r_single

IMPLICIT NONE
!
!.....................program constants (machine precision is about 7 digits)
!
REAL(KIND=r_single), PARAMETER :: tiny    = 1.e-30_r_single
REAL(KIND=r_single), PARAMETER :: vsmall  = 1.e-20_r_single
REAL(KIND=r_single), PARAMETER :: small   = 1.e-10_r_single

REAL(KIND=r_single), PARAMETER :: real_1e2  = 100.000000000_r_single
REAL(KIND=r_single), PARAMETER :: real_1em2 = 1._r_single/real_1e2

REAL(KIND=r_single), PARAMETER :: real_1e3  = 1000.000000000_r_single
REAL(KIND=r_single), PARAMETER :: real_1em3 = 1._r_single/real_1e3

REAL(KIND=r_single), PARAMETER :: real_1e4  = 10000.000000000_r_single
REAL(KIND=r_single), PARAMETER :: real_1em4 = 1._r_single/real_1e4

REAL(KIND=r_single), PARAMETER :: real_1e5  = 100000.000000000_r_single
REAL(KIND=r_single), PARAMETER :: real_1em5 = 1._r_single/real_1e5

REAL(KIND=r_single), PARAMETER :: real_1e6  = 1000000.000000000_r_single
REAL(KIND=r_single), PARAMETER :: real_1em6 = 1._r_single/real_1e6

REAL(KIND=r_single), PARAMETER :: real_1e7  = 10000000.000000000_r_single
REAL(KIND=r_single), PARAMETER :: real_1em7 = 1._r_single/real_1e7

REAL(KIND=r_single), PARAMETER :: large   = 1.e+10_r_single
REAL(KIND=r_single), PARAMETER :: vlarge  = 1.e+20_r_single
REAL(KIND=r_single), PARAMETER :: huge    = 1.e+30_r_single

REAL(KIND=r_single), PARAMETER :: zero    = 0.000000000_r_single
REAL(KIND=r_single), PARAMETER :: nearzero = zero+tiny

REAL(KIND=r_single), PARAMETER :: minusone  = -1.000000000_r_single
REAL(KIND=r_single), PARAMETER :: one       = 1.000000000_r_single
REAL(KIND=r_single), PARAMETER :: nearone   = one+tiny
REAL(KIND=r_single), PARAMETER :: two       = 2.000000000_r_single
REAL(KIND=r_single), PARAMETER :: three     = 3.000000000_r_single
REAL(KIND=r_single), PARAMETER :: four      = 4.000000000_r_single
REAL(KIND=r_single), PARAMETER :: five      = 5.000000000_r_single
REAL(KIND=r_single), PARAMETER :: six       = 6.000000000_r_single
REAL(KIND=r_single), PARAMETER :: seven     = 7.000000000_r_single
REAL(KIND=r_single), PARAMETER :: eight     = 8.000000000_r_single
REAL(KIND=r_single), PARAMETER :: nine      = 9.000000000_r_single
REAL(KIND=r_single), PARAMETER :: ten       = 10.000000000_r_single
REAL(KIND=r_single), PARAMETER :: twelve    = 12.000000000_r_single
REAL(KIND=r_single), PARAMETER :: svntytwo  = 72.000000000_r_single

REAL(KIND=r_single), PARAMETER :: threehalf  = three/two
REAL(KIND=r_single), PARAMETER :: fivehalf   = five/two
REAL(KIND=r_single), PARAMETER :: fivethree  = five/three
REAL(KIND=r_single), PARAMETER :: fourfifth  = four/five
REAL(KIND=r_single), PARAMETER :: sevensix   = seven/six
REAL(KIND=r_single), PARAMETER :: oneandhalf= 3.000000000_r_single/2.00000000000_r_single
REAL(KIND=r_single), PARAMETER :: elevensix = 11.00000000_r_single/6.00000000000_r_single
REAL(KIND=r_single), PARAMETER :: minuselevensix = -1._r_single*elevensix
REAL(KIND=r_single), PARAMETER :: half      = 1.000000000_r_single/2.00000000000_r_single
REAL(KIND=r_single), PARAMETER :: quarter   = 1.000000000_r_single/4.00000000000_r_single
REAL(KIND=r_single), PARAMETER :: threequarter = three*quarter
REAL(KIND=r_single), PARAMETER :: third     = 1.000000000_r_single/3.00000000000_r_single
REAL(KIND=r_single), PARAMETER :: twothird  = 2.000000000_r_single/3.00000000000_r_single
REAL(KIND=r_single), PARAMETER :: fourthird = 4.000000000_r_single/3.00000000000_r_single
REAL(KIND=r_single), PARAMETER :: sixth     = 1.000000000_r_single/6.00000000000_r_single
REAL(KIND=r_single), PARAMETER :: eighth    = 1.000000000_r_single/8.00000000000_r_single
REAL(KIND=r_single), PARAMETER :: seventh   = one/seven
REAL(KIND=r_single), PARAMETER :: fifth     = one/five
REAL(KIND=r_single), PARAMETER :: tenth     = one/ten

REAL(KIND=r_single), PARAMETER :: single_machine =1.0000000e-6_r_single

REAL(KIND=r_single), PARAMETER :: pi        = 3.1415926535897932384626434_r_single
REAL(KIND=r_single), PARAMETER :: halfpi    = pi/two
REAL(KIND=r_single), PARAMETER :: sqrtpi    = 0.9189385332046727417803297_r_single
REAL(KIND=r_single), PARAMETER :: twopi     = 2._r_single*pi
REAL(KIND=r_single), PARAMETER :: pioversix = pi/6._r_single
REAL(KIND=r_single), PARAMETER :: sixoverpi = 6._r_single/pi

REAL(KIND=r_single), PARAMETER :: radtodeg = 180.00000000_r_single/pi
REAL(KIND=r_single), PARAMETER :: degtorad = one / radtodeg
REAL(KIND=r_single), PARAMETER :: sqrtthree= 1.73205080756887
REAL(KIND=r_single), PARAMETER :: sqrttwo= 1.41421356237310

REAL(KIND=r_single), PARAMETER :: svntytwosqrtthree = svntytwo * sqrtthree !--Normalization constant for
                                                                           !--tetrahedron quality
REAL(KIND=r_single), PARAMETER :: foursqrtthree = four * sqrtthree         !--Normalization constant for
                                                                           !--triangle quality
REAL(KIND=r_single), PARAMETER :: epsilonP = two*real_1e6
REAL(KIND=r_single), PARAMETER :: distfactor = three
!
! use for setting mass flux on aentwsb (update_iflux) and switching signs on fluxes in update_<<boundary>> subroutines

REAL(KIND=r_single),PARAMETER,DIMENSION(0:1) :: rew01 = (/-1._r_single,1._r_single/)     ! input=0, rew01=-1. : input=1, rew01=1.

REAL(KIND=r_single),PARAMETER :: xbig = 171.624_r_single
REAL(KIND=r_single),PARAMETER :: xminin = 2.23e-308_r_single
REAL(KIND=r_single),PARAMETER :: eps = 2.22e-16_r_single
REAL(KIND=r_single),PARAMETER :: xinf = 1.79e308_r_single
REAL(KIND=r_single),PARAMETER,DIMENSION(8) :: P = (/-1.71618513886549492533811e+0_r_single,2.47656508055759199108314e+1_r_single,&
                                                    -3.79804256470945635097577e+2_r_single,6.29331155312818442661052e+2_r_single,&
                                                    8.66966202790413211295064e+2_r_single,-3.14512729688483675254357e+4_r_single,&
                                                    -3.61444134186911729807069e+4_r_single,6.64561438202405440627855e+4_r_single/)
REAL(KIND=r_single),PARAMETER,DIMENSION(8) :: Q = (/-3.08402300119738975254353e+1_r_single,3.15350626979604161529144e+2_r_single,&
                                                    -1.01515636749021914166146e+3_r_single,-3.10777167157231109440444e+3_r_single,&
                                                    2.25381184209801510330112e+4_r_single,4.75584627752788110767815e+3_r_single,&
                                                    -1.34659959864969306392456e+5_r_single,-1.15132259675553483497211e+5_r_single/)
REAL(KIND=r_single),PARAMETER,DIMENSION(7) :: C = (/-1.910444077728e-03_r_single,8.4171387781295e-04_r_single,&
                                                    -5.952379913043012e-04_r_single,7.93650793500350248e-04_r_single,&
                                                    -2.777777777777681622553e-03_r_single,8.333333333333333331554247e-02_r_single,&
                                                    5.7083835261e-03_r_single/)
REAL(KIND=r_single),PARAMETER,DIMENSION(0:8) :: p2 = (/ 0.99999999999980993, 676.5203681218851, -1259.1392167224028, &
                                                      771.32342877765313, -176.61502916214059, 12.507343278686905, &        
                                                      -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7 /)
END MODULE real_parameters
