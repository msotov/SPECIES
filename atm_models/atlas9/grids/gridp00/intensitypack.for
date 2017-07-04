      program intensitypack
      real wave(1221),i17(17),centi(1221)
      integer inu(17,1221)
      character*80 titlea,titleb
      CHARACTER*105 TITLE
CTEFF   3500.  GRAVITY 0.00000 LTE
CTITLE SDSC GRID  [-2.0]   VTURB 2.0 KM/S    L/H 1.25
C 17 ANGLES 1.0000 0.9000 0.8000 0.7000 0.6000 0.5000 0.4000 0.3000 0.2500 0.2000
C           0.1500 0.1250 0.1000 0.0750 0.0500 0.0250 0.0100
CINTENSITY  124    92.50   3.241000E+15
C 1.111E-35 1.000E-35 8.890E-36 7.779E-36 6.669E-36 5.558E-36 4.447E-36 3.337E-36
C 2.782E-36 2.226E-36 1.671E-36 1.394E-36 1.117E-36 8.391E-37 5.619E-37 2.849E-37
C 1.192E-37
      READ(3,33)WAVE
   33 FORMAT(8f10.2)
      do 9 mod=1,500
      read(1,1,end=10)titlea,titleb
    1 format(a80)
      title=titlea(1:30)//titleb(6:80)
      READ(1,1)
      READ(1,1)
      do 3 nu=1,1221
      centi(NU)=0.
      do 3 mu=1,17
    3 inu(mu,nu)=0
      do 4 nu=1,1222
      read(1,2)i
    2 format(9x,i5)
      if(i.eq.0)go to 5
      read(1,22)i17
   22 format(8e10.3/8e10.3/e10.3)
      centi(I)=i17(1)
      i17(1)=max(i17(1),1.e-38)
      DO 4 MU=2,17
      INU(MU,I)=I17(MU)/I17(1)*100000.+.5
C      INU(MU,I)=MIN(99999,INU(MU,I))
    4 CONTINUE
    5 write(2,6)title
      write(2,13)
   13 format('   wl(nm)    Inu(ergs/cm**2/s/hz/ster) for 17 mu',
     1' in 1221 frequency intervals'/
     2'             1.000   .900  .800  .700  .600  .500  .400 ',
     3' .300  .250  .200  .150  .125  .100  .075  .050  .025  .010')
    6 FORMAT(A105)
      DO 8 NU=1,1221
    8 WRITE(2,7)wave(nu),CENTI(NU),(INU(MU,NU),MU=2,17)
    7 format(F9.2,1pe10.3,16I6)
    9 continue
   10 call exit
      end
