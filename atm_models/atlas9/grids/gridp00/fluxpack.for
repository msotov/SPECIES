      program fluxpack
      dimension hnu(1221),hcont(1221),wave(1221)
      character*80 titlea,titleb
cTEFF   3500.  GRAVITY 0.00000 LTE
cTITLE SDSC GRID  [-2.0]   VTURB 2.0 KM/S    L/H 1.25
cFLUX  124    92.50        3.241000E+15   3.4344E-36   3.1968E-32   0.00011
c
cTEFF   3500.  GRAVITY 0.00000  SDSC GRID  [-2.0]   VTURB 2.0 KM/S    L/H 1.25
      do 9 mod=1,500
      read(1,1,end=10)titlea,titleb
    1 format(a80)
      titlea=titlea(1:30)//titleb(6:55)
      do 3 nu=1,1221
      hnu(nu)=0.
    3 hcont(nu)=0.
      do 4 nu=1,1222
      read(1,2)inu,h,hc
    2 format(4x,i5,9x,20x,e13.4,e13.4)
      if(inu.eq.0)go to 5
      hnu(inu)=h
    4 hcont(inu)=hc
    5 write(2,1)titlea
      write(2,7)hnu
      write(2,7)hcont
    7 format(1p8e10.4)
    9 continue
   10 call exit
      end
