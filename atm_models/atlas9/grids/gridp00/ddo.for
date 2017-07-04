      PROGRAM DDO
      DIMENSION HLAM(1221),WAVE(1221),HNU(1221),HNUCONT(1221)                                            
      DIMENSION AMAGI(13,1221),TRANSI(13,1221),EBVI(13),AMAG(1221)
      CHARACTER*80 TITLE
      DIMENSION A(8)                                                            
      DIMENSION F(750)
      DIMENSION W35(31),FILT35(31),F35(31),WAVE35(750),S35(750)                         
      DIMENSION W38(34),FILT38(34),F38(34),WAVE38(660),S38(660)                         
      DIMENSION W41(35),FILT41(35),F41(35),WAVE41(170),S41(170)                         
      DIMENSION W42(31),FILT42(31),F42(31),WAVE42(150),S42(150)                         
      DIMENSION W45(31),FILT45(31),F45(31),WAVE45(150),S45(150)                         
      DIMENSION W48(40),FILT48(40),F48(40),WAVE48(390),S48(390)                         
      REAL MAG35,MAG38,MAG41,MAG42,MAG45,MAG48
C      DATA EBVI/0.,.1,.2,.3,.4,.5,.6,.7,.8,.9,1./
      DATA EBVI/0.,.1,.2,.3,.4,.5,.6,.8,1.,2.,3.,4.,5./
      DATA FILT35/                                                               
     1 .000, .000, .004, .020, .053, .104, .162, .225, .277, .323,
     2 .356, .385, .396, .405, .407, .405, .395, .380, .368, .330,
     3 .294, .254, .205, .152, .102, .062, .030, .010, .002, .000,
     4 .000/
      DATA W35/
     1310.0,312.5,315.0,317.5,320.0,322.5,325.0,327.5,330.0,332.5,
     2335.0,337.5,340.0,342.5,345.0,347.5,350.0,352.5,355.0,357.5,
     3360.0,362.5,365.0,367.5,370.0,372.5,375.0,377.5,380.0,382.5,
     4385.0/
      DATA FILT38/
     1 .000, .000, .019, .056, .113, .189, .273, .358, .435, .503,
     2 .558, .597, .628, .645, .652, .651, .642, .628, .601, .561,
     3 .508, .461, .395, .331, .267, .207, .147, .101, .065, .038,
     4 .019, .007, .000, .000/
      DATA W38/
     1  352., 354., 356., 358., 360., 362., 364., 366., 368., 370.,
     2  372., 374., 376., 378., 380., 382., 384., 386., 388., 390.,
     3  392., 394., 396., 398., 400., 402., 404., 406., 408., 410.,
     4  412., 414., 416., 418./
      DATA FILT41/
     1  .000, .000, .002, .006, .014, .026, .049, .092, .154, .260,
     2  .352, .437, .476, .483, .487, .489, .490, .490, .488, .485,
     3  .482, .476, .465, .407, .319, .215, .125, .072, .046, .028, 
     4  .013, .007, .003, .000, .000/
      DATA W41/
     1 408.0,408.5,409.0,409.5,410.0,410.5,411.0,411.5,412.0,412.5,
     2 413.0,413.5,414.0,414.5,415.0,415.5,416.0,416.5,417.0,417.5,
     3 418.0,418.5,419.0,419.5,420.0,420.5,421.0,421.5,422.0,422.5,
     4 423.0,423.5,424.0,424.5,425.0/
      DATA FILT42/
     1  .000, .000, .004, .014, .023, .067, .141, .265, .400, .506,
     2  .527, .522, .514, .513, .524, .531, .536, .531, .507, .455,
     3  .385, .278, .186, .102, .057, .032, .018, .008, .003, .000,
     4  .000/
      DATA W42/
     1 418.5,419.0,419.5,420.0,420.5,421.0,421.5,422.0,422.5,423.0,
     2 423.5,424.0,424.5,425.0,425.5,426.0,426.5,427.0,427.5,428.0,
     3 428.5,429.0,429.5,430.0,430.5,431.0,431.5,432.0,432.5,433.0,
     4 433.5/
      DATA FILT45/
     1  .000, .000, .006, .015, .032, .082, .185, .310, .450, .512,
     2  .524, .527, .538, .548, .551, .540, .527, .518, .518, .512,
     3  .470, .340, .195, .100, .046, .024, .014, .008, .003, .000,
     4  .000/
      DATA W45/
     1 445.0,445.5,446.0,446.5,447.0,447.5,448.0,448.5,449.0,449.5,
     2 450.0,450.5,451.0,451.5,452.0,452.5,453.0,453.5,454.0,454.5,
     3 455.0,455.5,456.0,456.5,457.0,457.5,458.0,458.5,459.0,459.5,
     4 460.0/
      DATA FILT48/
     1 .000, .000, .006, .009, .013, .020, .032, .048, .074, .116,
     2 .183, .282, .394, .488, .542, .561, .565, .563, .560, .557,
     3 .555, .552, .550, .547, .544, .537, .522, .485, .435, .346,
     4 .260, .180, .125, .087, .062, .042, .029, .010, .000, .000/
      DATA W48/
     1  468., 469., 470., 471., 472., 473., 474., 475., 476., 477.,
     2  478., 479., 480., 481., 482., 483., 484., 485., 486., 487.,
     3  488., 489., 490., 491., 492., 493., 494., 495., 496., 497.,
     4  498., 499., 500., 501., 502., 503., 504., 505., 506., 507./
C
      NFILT=31
      N35=750
      DO 11 I=1,NFILT                                                              
   11 F35(I)=FILT35(I)*ONEP21(W35(I))*AIR(W35(I),0.)*REFLCT(W35(I))          
      DO 12 I=1,N35                                                             
   12 WAVE35(I)=310.+FLOAT(I)*.1                                                 
      CALL PINTER(W35,F35,NFILT,WAVE35,S35,N35)                                      
      SNORM35=0.                                                                  
      DO 13 I=1,N35                                                             
   13 SNORM35=SNORM35+S35(I)                                                         
      SNORM35=SNORM35*.1                                                            
      S35NOMAG=-2.5*ALOG10(SNORM35)                                                 
C      
      NFILT=34
      N38=660
      DO 21 I=1,NFILT                                                              
   21 F38(I)=FILT38(I)*ONEP21(W38(I))*AIR(W38(I),0.)*REFLCT(W38(I))          
      DO 22 I=1,N38                                                             
   22 WAVE38(I)=352.+FLOAT(I)*.1                                                 
      CALL PINTER(W38,F38,NFILT,WAVE38,S38,N38)                                      
      SNORM38=0.                                                                  
      DO 23 I=1,N38                                                             
   23 SNORM38=SNORM38+S38(I)                                                         
      SNORM38=SNORM38*.1                                                            
      S38NOMAG=-2.5*ALOG10(SNORM38)                                                 
C      
      NFILT=35
      N41=170
      DO 31 I=1,NFILT                                                              
   31 F41(I)=FILT41(I)*ONEP21(W41(I))*AIR(W41(I),0.)*REFLCT(W41(I))          
      DO 32 I=1,N41                                                             
   32 WAVE41(I)=408.+FLOAT(I)*.1                                                 
      CALL PINTER(W41,F41,NFILT,WAVE41,S41,N41)                                      
      SNORM41=0.                                                                  
      DO 33 I=1,N41                                                             
   33 SNORM41=SNORM41+S41(I)                                                         
      SNORM41=SNORM41*.1                                                            
      S41NOMAG=-2.5*ALOG10(SNORM41)                                                 
C      
      NFILT=31
      N42=150
      DO 41 I=1,NFILT                                                              
   41 F42(I)=FILT42(I)*ONEP21(W42(I))*AIR(W42(I),0.)*REFLCT(W42(I))          
      DO 42 I=1,N42                                                             
   42 WAVE42(I)=418.5+FLOAT(I)*.1                                                 
      CALL PINTER(W42,F42,NFILT,WAVE42,S42,N42)                                      
      SNORM42=0.                                                                  
      DO 43 I=1,N42                                                             
   43 SNORM42=SNORM42+S42(I)                                                         
      SNORM42=SNORM42*.1                                                            
      S42NOMAG=-2.5*ALOG10(SNORM42)                                                 
C      
      NFILT=31
      N45=150
      DO 51 I=1,NFILT                                                              
   51 F45(I)=FILT45(I)*ONEP21(W45(I))*AIR(W45(I),0.)*REFLCT(W45(I))          
      DO 52 I=1,N45                                                             
   52 WAVE45(I)=445.+FLOAT(I)*.1                                                 
      CALL PINTER(W45,F45,NFILT,WAVE45,S45,N45)                                      
      SNORM45=0.                                                                  
      DO 53 I=1,N45                                                             
   53 SNORM45=SNORM45+S45(I)                                                         
      SNORM45=SNORM45*.1                                                            
      S45NOMAG=-2.5*ALOG10(SNORM45)                                                 
C      
      NFILT=40
      N48=390
      DO 61 I=1,NFILT                                                              
   61 F48(I)=FILT48(I)*ONEP21(W48(I))*AIR(W48(I),0.)*REFLCT(W48(I))          
      DO 62 I=1,N48                                                             
   62 WAVE48(I)=468.+FLOAT(I)*.1                                                 
      CALL PINTER(W48,F48,NFILT,WAVE48,S48,N48)                                      
      SNORM48=0.                                                                  
      DO 63 I=1,N48                                                             
   63 SNORM48=SNORM48+S48(I)                                                         
      SNORM48=SNORM48*.1                                                            
      S48NOMAG=-2.5*ALOG10(SNORM48)                                                 
C      
CSDSC GRID  [+0.0]  VTURB 2.0 KM/S   L/H 1.25
      READ(1,5)ABUND,VTURB,CONVEC
    5 FORMAT(12X,F4.1,8X,F4.1,11X,F5.2)
      DO 616 ISKIP=1,21
  616 READ(1,1)
C     wavelength in nm
      READ(1,1)WAVE
    1 FORMAT(8F10.2)
      RV=3.1
      EBV=.1
C     CALL REDDENING(WAVE,RV,EBV,AMAG)
      READ(2,344)
      READ(2,344)
      DO 366 NU=1,1221
  366 READ(2,344)
      READ(2,344)
      EBVI(1)=0.
      READ(2,344)(EBVI(IRED),IRED=2,13)
  344 FORMAT(10X,12F10.1)
      DO 367 NU=1,1221
      TRANSI(1,NU)=1.
      READ(2,359)(TRANSI(IRED,NU),IRED=2,13)
  359 FORMAT(13X,12E10.3)
  367 CONTINUE
      WRITE(6,6)
      WRITE(7,6)
      WRITE(8,6)
    6 FORMAT('        Teff log g  [M]  Vturb  l/H  E(B-V)',
     1'   35     38     41     42     45     48',
     2'   35-38  38-41  41-42  42-45  45-48')
      DO 1000 NMODEL=1,500                                              
      READ(1,711,END=9)TITLE
  711 FORMAT(A80)
      PRINT 713,MODEL,TITLE
  713 FORMAT(I5,1X,A80)
      READ(TITLE,'(5X,I6,10X,F8.5)')ITEFF,GLOG
C     ergs/cm**2/s/hz/ster
      READ(1,4)Hnu
      READ(1,4)HnuCONT
    4 FORMAT(8E10.4)
      NNU=1221
      DO 900 IRED=1,13
      DO 715 NU=1,1221
      FREQ=2.99792458E17/WAVE(NU)
  715 HLAM(NU)=HNU(NU)*FREQ/WAVE(NU)*TRANSI(IRED,NU)
C     PRINT 77,(WAVE(I),HLAM(I),I=1,NNU)                                        
C
      CALL LINTER(WAVE,HLAM,NNU,WAVE35,F,N35)                                    
C     PRINT 77,F                                                                
      SF35=0.                                                                     
      DO 14 I=1,N35                                                             
   14 SF35=SF35+S35(I)*F(I)                                                          
      SF35=SF35*.1                                                                  
      MAG35=-2.5*ALOG10(SF35)-S35NOMAG                                                          
C
      CALL LINTER(WAVE,HLAM,NNU,WAVE38,F,N38)                                    
C     PRINT 77,F                                                                
      SF38=0.                                                                     
      DO 24 I=1,N38                                                             
   24 SF38=SF38+S38(I)*F(I)                                                          
      SF38=SF38*.1                                                                  
      MAG38=-2.5*ALOG10(SF38)-S38NOMAG                                                          
C
      CALL LINTER(WAVE,HLAM,NNU,WAVE41,F,N41)                                    
C     PRINT 77,F                                                                
      SF41=0.                                                                     
      DO 34 I=1,N41                                                             
   34 SF41=SF41+S41(I)*F(I)                                                          
      SF41=SF41*.1                                                                  
      MAG41=-2.5*ALOG10(SF41)-S41NOMAG                                                          
C
      CALL LINTER(WAVE,HLAM,NNU,WAVE42,F,N42)                                    
C     PRINT 77,F                                                                
      SF42=0.                                                                     
      DO 44 I=1,N42                                                             
   44 SF42=SF42+S42(I)*F(I)                                                          
      SF42=SF42*.1                                                                  
      MAG42=-2.5*ALOG10(SF42)-S42NOMAG                                                          
C
      CALL LINTER(WAVE,HLAM,NNU,WAVE45,F,N45)                                    
C     PRINT 77,F                                                                
      SF45=0.                                                                     
      DO 54 I=1,N45                                                             
   54 SF45=SF45+S45(I)*F(I)                                                          
      SF45=SF45*.1                                                                  
      MAG45=-2.5*ALOG10(SF45)-S45NOMAG                                                          
C
      CALL LINTER(WAVE,HLAM,NNU,WAVE48,F,N48)                                    
C     PRINT 77,F                                                                
      SF48=0.                                                                     
      DO 64 I=1,N48                                                             
   64 SF48=SF48+S48(I)*F(I)                                                          
      SF48=SF48*.1                                                                  
      MAG48=-2.5*ALOG10(SF48)-S48NOMAG                                                          
C
      C35MIN38=MAG38-MAG35
      C38MIN41=MAG41-MAG38
      C41MIN42=MAG42-MAG41
      C42MIN45=MAG45-MAG42
      C45MIN48=MAG48-MAG45
C     NORMALIZATION TO VEGA
C     VEGA CRAWFORD AND BARNES 1970
c     1.411   1.411   0.004   0.156   1.089   1.089
c     Crawford, D.L. and Barnes, J.V. 1970, A.J.,75,978.
      XSCALE=ABUND
c     actually l/h
      xh=CONVEC
      if(ItEFF.GE.9000)xh=0.
      WRITE(6,60)NMODEL,ITEFF,GLOG,XSCALE,VTURB,XH,EBVI(IRED),
     1MAG35,MAG38,MAG41,MAG42,MAG45,MAG48,C35MIN38,C38MIN41,C41MIN42,
     2C42MIN45,C45MIN48
   60 FORMAT(I6,I6,5F6.2,6F7.3,5F7.3)
      WRITE(7,60)NMODEL,ITEFF,GLOG,XSCALE,VTURB,XH,EBVI(IRED),
     1MAG35,MAG38,MAG41,MAG42,MAG45,MAG48,C35MIN38,C38MIN41,C41MIN42,
     2C42MIN45,C45MIN48
      IF(IRED.EQ.1)
     1WRITE(8,60)NMODEL,ITEFF,GLOG,XSCALE,VTURB,XH,EBVI(IRED),
     2MAG35,MAG38,MAG41,MAG42,MAG45,MAG48,C35MIN38,C38MIN41,C41MIN42,
     3C42MIN45,C45MIN48
  900 CONTINUE                                                                  
 1000 CONTINUE                                                                  
    9 CALL EXIT                                                                 
      END                                                                       
      SUBROUTINE LINTER(XOLD,YOLD,NOLD,XNEW,YNEW,NNEW)                          
      DIMENSION XOLD(1),YOLD(1),XNEW(1),YNEW(1)                                 
C     XOLD AND XNEW INCREASING                                                  
      IOLD=2                                                                    
      DO 2 INEW=1,NNEW                                                          
    1 IF(XNEW(INEW).LT.XOLD(IOLD))GO TO 2                                       
      IF(IOLD.EQ.NOLD)GO TO 2                                                   
      IOLD=IOLD+1                                                               
      GO TO 1                                                                   
    2 YNEW(INEW)=YOLD(IOLD-1)+(YOLD(IOLD)-YOLD(IOLD-1))/                        
     1(XOLD(IOLD)-XOLD(IOLD-1))*(XNEW(INEW)-XOLD(IOLD-1))                       
      RETURN                                                                    
      END                                                                       
      SUBROUTINE PINTER(XOLD,FOLD,NOLD,XNEW,FNEW,NNEW)                          
      DIMENSION XOLD(1),FOLD(1),XNEW(1),FNEW(1)                                 
      L=2                                                                       
      LL=0                                                                      
      DO 50 K=1,NNEW                                                            
   10 IF(XNEW(K).LT.XOLD(L))GO TO 20                                            
      L=L+1                                                                     
      IF(L.GT.NOLD)GO TO 30                                                     
      GO TO 10                                                                  
   20 IF(L.EQ.LL)GO TO 50                                                       
      IF(L.EQ.2)GO TO 30                                                        
      L1=L-1                                                                    
      IF(L.GT.LL+1.OR.L.EQ.3)GO TO 21                                           
      CBAC=CFOR                                                                 
      BBAC=BFOR                                                                 
      ABAC=AFOR                                                                 
      IF(L.EQ.NOLD)GO TO 22                                                     
      GO TO 25                                                                  
   21 L2=L-2                                                                    
      D=(FOLD(L1)-FOLD(L2))/(XOLD(L1)-XOLD(L2))                                 
      CBAC=FOLD(L)/((XOLD(L)-XOLD(L1))*(XOLD(L)-XOLD(L2)))+                     
     1(FOLD(L2)/(XOLD(L)-XOLD(L2))-FOLD(L1)/(XOLD(L)-XOLD(L1)))/                
     2(XOLD(L1)-XOLD(L2))                                                       
      BBAC=D-(XOLD(L1)+XOLD(L2))*CBAC                                           
      ABAC=FOLD(L2)-XOLD(L2)*D+XOLD(L1)*XOLD(L2)*CBAC                           
      IF(L.LT.NOLD)GO TO 25                                                     
   22 C=CBAC                                                                    
      B=BBAC                                                                    
      A=ABAC                                                                    
      LL=L                                                                      
      GO TO 50                                                                  
   25 D=(FOLD(L)-FOLD(L1))/(XOLD(L)-XOLD(L1))                                   
      CFOR=FOLD(L+1)/((XOLD(L+1)-XOLD(L))*(XOLD(L+1)-XOLD(L1)))+                
     1(FOLD(L1)/(XOLD(L+1)-XOLD(L1))-FOLD(L)/(XOLD(L+1)-XOLD(L)))/              
     2(XOLD(L)-XOLD(L1))                                                        
      BFOR=D-(XOLD(L)+XOLD(L1))*CFOR                                            
      AFOR=FOLD(L1)-XOLD(L1)*D+XOLD(L)*XOLD(L1)*CFOR                            
      WT=0.                                                                     
      IF(ABS(CFOR).NE.0.)WT=ABS(CFOR)/(ABS(CFOR)+ABS(CBAC))                     
      A=AFOR+WT*(ABAC-AFOR)                                                     
      B=BFOR+WT*(BBAC-BFOR)                                                     
      C=CFOR+WT*(CBAC-CFOR)                                                     
      LL=L                                                                      
      GO TO 50                                                                  
   30 IF(L.EQ.LL)GO TO 50                                                       
      L=AMIN0(NOLD,L)                                                           
      C=0.                                                                      
      B=(FOLD(L)-FOLD(L-1))/(XOLD(L)-XOLD(L-1))                                 
      A=FOLD(L)-XOLD(L)*B                                                       
      LL=L                                                                      
   50 FNEW(K)=A+(B+C*XNEW(K))*XNEW(K)                                           
      MAP1=LL-1                                                                 
      RETURN                                                                    
      END                                                                       
      FUNCTION ONEP21(WAVE)                                                     
C     RESPONSE OF A STANDARD 1P21  RCA CURVE                                    
C     10 NM STEPS FROM 300NM TO 600NM NORMALIZED TO 400NM                       
      DIMENSION WL(31),RESPON(31)                                               
      DATA RESPON/.08,.22,.45,.68,.81,.90,.95,.97,.99,1.00,1.00,.99,.98,        
     1 .95,.91,.87,.83,.77,.71,.65,.58,.52,.46,.40,.34,.29,.24,.20,.16,         
     2 .13,.10/                                                                 
      DATA WL/300.,310.,320.,330.,340.,350.,360.,370.,380.,390.,400.,           
     1 410.,420.,430.,440.,450.,460.,470.,480.,490.,500.,510.,520.,530.,        
     2 540.,550.,560.,570.,580.,590.,600./                                      
      CALL LINTER(WL,RESPON,31,WAVE,ONEP21,1)                                   
      RETURN                                                                    
      END                                                                       
      FUNCTION AIR(WAVE,ATMOS)                                                  
C     FROM ALLEN PAGE 126                                                       
      DIMENSION COEF(10),W(10)                                                  
      DATA COEF/4.5,1.30,.84,.68,.55,.46,.31,.23,.195,.170/                     
      DATAW/300.,320.,340.,360.,380.,400.,450.,500.,550.,600./                  
      AIR=1.                                                                    
      IF(ATMOS.EQ.0.)RETURN                                                     
      CALL LINTER(W,COEF,10,WAVE,C,1)                                           
      AIR=EXP(-C*ATMOS)                                                         
      RETURN                                                                    
      END                                                                       
      FUNCTION REFLCT(WAVE)                                                     
C     FROM ALLEN                                                                
      DIMENSION ALUM(8),WL(8)                                                   
      DATA ALUM/.82,.83,.84,.85,.86,.87,.88,.89/                                
      DATA WL/300.,350.,380.,400.,450.,500.,550.,600./                          
      CALL LINTER(WL,ALUM,8,WAVE,R,1)                                           
      REFLCT=.77                                                                
      REFLCT=R**2                                                               
      RETURN                                                                    
      END                                                                       
