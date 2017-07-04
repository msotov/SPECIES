      PROGRAM UVBY
      DIMENSION HLAM(1221),WAVE(1221),HNU(1221),HNUCONT(1221)                                            
      DIMENSION AMAGI(13,1221),TRANSI(13,1221),EBVI(13),AMAG(1221)
      CHARACTER*80 TITLE
      DIMENSION A(8)                                                            
      DIMENSION SU(700),SU1(700),SU2(700),SV(700),SB(700),SY(700)               
      DIMENSION WAVEU(700),WAVEV(700),WAVEB(700),WAVEY(700)                     
      DIMENSION F(700)                                                          
      DIMENSION UWAVE(29),VWAVE(29),BWAVE(29),YWAVE(29),CAWAVE(29)                         
      DIMENSION U1(29),U2(29)                                                   
      DIMENSION WAVECA(700),SCA(700)
      REAL M1                                                                   
      DIMENSION U(29),V(29),B(29),Y(29),CA(29)                                         
      DIMENSION UFILT(29),VFILT(29),BFILT(29),YFILT(29),CAFILT(29)                         
C      DATA EBVI/0.,.1,.2,.3,.4,.5,.6,.7,.8,.9,1./
      DATA EBVI/0.,.1,.2,.3,.4,.5,.6,.8,1.,2.,3.,4.,5./
      DATA UFILT/                                                               
     1 .0000,.0060,.0592,.1180,.1780,.2380,.2912,.3332,.3608,.3812,             
     2 .3928,.3944,.3908,.3852,.3708,.3532,.3308,.3008,.2616,.2168,             
     3 .1692,.1176,.0768,.0468,.0216,.0088,.0040,.0000,.0000/                   
      DATA VFILT/                                                               
     1 .0000,.0016,.0036,.0084,.0148,.0228,.0304,.0480,.0784,.1300,             
     2 .2000,.2988,.4000,.4728,.4932,.4800,.4352,.3720,.2820,.1812,             
     3 .1112,.0672,.0400,.0268,.0200,.0140,.0072,.0036,.0000/                   
      DATA BFILT/                                                               
     1 .0000,.0040,.0096,.0160,.0232,.0360,.0500,.0800,.1240,.2000,             
     2 .3020,.4040,.4588,.4688,.4504,.3820,.2784,.1740,.1100,.0688,             
     3 .0392,.0252,.0152,.0124,.0092,.0072,.0040,.0020,.0000/                   
C     1 .0000,.0068,.0168,.0272,.4000,.0696,.1020,.1520,.2340,.3300,             
      DATA YFILT/                                                               
     1 .0000,.0068,.0168,.0272,.0400,.0696,.1020,.1520,.2340,.3300,             
     2 .4100,.4540,.4792,.5012,.5220,.5132,.4544,.3480,.2480,.1548,             
     3 .1000,.0660,.0416,.0252,.0136,.0096,.0080,.0040,.0000/                   
      DATA CAFILT/
     1 .0000,.0000,.0000,.0000,.0000,.0100,.0100,.0300,.0700,.1500,
     2 .2600,.3300,.3600,.3800,.3900,.4000,.3900,.3500,.2500,.1400,
     3 .0700,.0300,.0200,.0100,.0000,.0000,.0000,.0000,.0000/
      DO 11 I=1,29                                                              
      UWAVE(I)=312.5+FLOAT(I)*2.5                                               
      VWAVE(I)=372.5+FLOAT(I)*2.5                                               
      BWAVE(I)=432.5+FLOAT(I)*2.5                                               
      YWAVE(I)=512.5+FLOAT(I)*2.5                                               
      CAWAVE(I)=380.+FLOAT(I)
      U(I)=UFILT(I)*ONEP21(UWAVE(I))*AIR(UWAVE(I),0.)*REFLCT(UWAVE(I))          
      U1(I)=UFILT(I)*ONEP21(UWAVE(I))*AIR(UWAVE(I),1.)*REFLCT(UWAVE(I))         
      U2(I)=UFILT(I)*ONEP21(UWAVE(I))*AIR(UWAVE(I),2.)*REFLCT(UWAVE(I))         
      V(I)=VFILT(I)*ONEP21(VWAVE(I))*AIR(VWAVE(I),0.)*REFLCT(VWAVE(I))          
      B(I)=BFILT(I)*ONEP21(BWAVE(I))*AIR(BWAVE(I),0.)*REFLCT(BWAVE(I))          
      Y(I)=YFILT(I)*ONEP21(YWAVE(I))*AIR(YWAVE(I),0.)*REFLCT(YWAVE(I))          
      CA(I)=CAFILT(I)*ONEP21(CAWAVE(I))*AIR(CAWAVE(I),0.)*
     1      REFLCT(CAWAVE(I))          
   11 CONTINUE                                                                  
   77 FORMAT(10E12.4)                                                           
C     PRINT 77,U                                                                
C     PRINT 77,U1                                                               
C     PRINT 77,U2                                                               
C     PRINT 77,V                                                                
C     PRINT 77,B                                                                
C     PRINT 77,Y                                                                
C     PRINT 77,UWAVE                                                            
C     PRINT 77,VWAVE                                                            
C     PRINT 77,BWAVE                                                            
C     PRINT 77,YWAVE                                                            
      DO 12 I=1,700                                                             
      WAVEU(I)=315.+FLOAT(I)*.1                                                 
      WAVEV(I)=375.+FLOAT(I)*.1                                                 
      WAVEB(I)=435.+FLOAT(I)*.1                                                 
      WAVEY(I)=515.+FLOAT(I)*.1                                                 
      WAVECA(I)=381.+FLOAT(I)*.1
   12 CONTINUE                                                                  
      CALL PINTER(UWAVE,U,28,WAVEU,SU,675)                                      
      CALL PINTER(UWAVE,U1,28,WAVEU,SU1,675)                                    
      CALL PINTER(UWAVE,U2,28,WAVEU,SU2,675)                                    
      CALL PINTER(VWAVE,V,29,WAVEV,SV,700)                                      
      CALL PINTER(BWAVE,B,29,WAVEB,SB,700)                                      
      CALL PINTER(YWAVE,Y,29,WAVEY,SY,700)                                      
      CALL PINTER(CAWAVE,CA,29,WAVECA,SCA,700)
C     PRINT 77,SU                                                               
C     PRINT 77,SU1                                                              
C     PRINT 77,SU2                                                              
C     PRINT 77,SV                                                               
C     PRINT 77,SB                                                               
C     PRINT 77,SY                                                               
      UNORM=0.                                                                  
      DO 13 I=1,675                                                             
   13 UNORM=UNORM+SU(I)                                                         
      VNORM=0.                                                                  
      BNORM=0.                                                                  
      YNORM=0.                                                                  
      CANORM=0.
      DO 14 I=1,700                                                             
      CANORM=CANORM+SCA(I)
      VNORM=VNORM+SV(I)                                                         
      BNORM=BNORM+SB(I)                                                         
   14 YNORM=YNORM+SY(I)                                                         
      UNORM=UNORM*.1                                                            
      VNORM=VNORM*.1                                                            
      BNORM=BNORM*.1                                                            
      YNORM=YNORM*.1                                                            
      CANORM=CANORM*.1
      UNOMAG=-2.5*ALOG10(UNORM)                                                 
      VNOMAG=-2.5*ALOG10(VNORM)                                                 
      BNOMAG=-2.5*ALOG10(BNORM)                                                 
      YNOMAG=-2.5*ALOG10(YNORM)                                                 
      CANOMAG=-2.5*ALOG10(CANORM)
C     PRINT 77,UNORM,VNORM,BNORM,YNORM                                          
C     PRINT 77,UNOMAG,VNOMAG,BNOMAG,YNOMAG                                      
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
     1'   u       ca      v       b       y',
     2'       u-b     b-y      m1      c1      hk')
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
      CALL LINTER(WAVE,HLAM,NNU,WAVEU,F,675)                                    
C     PRINT 77,F                                                                
      UF=0.                                                                     
      U1F=0.                                                                    
      U2F=0.                                                                    
      DO 22 I=1,675                                                             
      UF=UF+SU(I)*F(I)                                                          
      U1F=U1F+SU1(I)*F(I)                                                       
      U2F=U2F+SU2(I)*F(I)                                                       
   22 CONTINUE                                                                  
      UF=UF*.1                                                                  
      U1F=U1F*.1                                                                
      U2F=U2F*.1                                                                
C     PRINT 77,UF                                                               
C     PRINT 77,U1F                                                              
C     PRINT 77,U2F                                                              
      CALL LINTER(WAVE,HLAM,NNU,WAVEV,F,700)                                    
C     PRINT 77,F                                                                
      VF=0.                                                                     
      DO 32 I=1,700                                                             
   32 VF=VF+SV(I)*F(I)                                                          
      VF=VF*.1                                                                  
C     PRINT 77,VF                                                               
      CALL LINTER(WAVE,HLAM,NNU,WAVEB,F,700)                                    
C     PRINT 77,F                                                                
      BF=0.                                                                     
      DO 42 I=1,700                                                             
   42 BF=BF+SB(I)*F(I)                                                          
      BF=BF*.1                                                                  
C     PRINT 77,BF                                                               
      CALL LINTER(WAVE,HLAM,NNU,WAVEY,F,700)                                    
C     PRINT 77,F                                                                
      YF=0.                                                                     
      DO 52 I=1,700                                                             
   52 YF=YF+SY(I)*F(I)                                                          
      YF=YF*.1                                                                  
      CALL LINTER(WAVE,HLAM,NNU,WAVECA,F,700)                                    
C     PRINT 77,F                                                                
      CAF=0.                                                                     
      DO 62 I=1,700                                                             
   62 CAF=CAF+SCA(I)*F(I)                                                          
      CAF=CAF*.1                                                                  
C     PRINT 77,CAF                                                               
C     UNOMAG=0.                                                                 
C     VNOMAG=0.                                                                 
C     BNOMAG=0.                                                                 
C     YNOMAG=0.                                                                 
      UMAG=-2.5*ALOG10(UF)                                                      
      U1MAG=-2.5*ALOG10(U1F)                                                    
      U2MAG=-2.5*ALOG10(U2F)                                                    
      U0MAG=2.*U1MAG-U2MAG                                                      
      VMAG=-2.5*ALOG10(VF)                                                      
      BMAG=-2.5*ALOG10(BF)                                                      
      YMAG=-2.5*ALOG10(YF)                                                      
      CAMAG=-2.5*ALOG10(CAF)
      UMAG=UMAG-UNOMAG                                                          
      U0MAG=U0MAG-UNOMAG                                                        
      VMAG=VMAG-VNOMAG                                                          
      BMAG=BMAG-BNOMAG                                                          
      YMAG=YMAG-YNOMAG                                                          
      CAMAG=CAMAG-CANOMAG
      UMINB=UMAG-BMAG                                                           
      U0MINB=U0MAG-BMAG                                                         
      BMINY=BMAG-YMAG                                                           
      VMINB=VMAG-BMAG                                                           
      M1=VMINB-BMINY                                                            
      UMINV=UMAG-VMAG                                                           
      U0MINV=U0MAG-VMAG                                                         
      C1=UMINV-VMINB                                                            
      C10=U0MINV-VMINB                                                          
      hk=(camag-bmag)-bminy
C     NORMALIZATION TO VEGA
      UMINB=UMINB+.802
      U0MINB=U0MINB+.780
      BMINY=BMINY+.500
      M1=M1-.071
      C1=C1-.058
      C10=C10-.080
C     GUESS FROM HD83373  = .194  = A1V  = 9250,4
C     RAW AP00T9250G40K2 = 0.331
      hk=hk-0.133
C     VEGA CRAWFORD AND BARNES 1970
c     1.411  0.004   0.156   1.089
c     Crawford, D.L. and Barnes, J.V. 1970, A.J.,75,978.
      XSCALE=ABUND
c     actually l/h
      xh=CONVEC
      if(ItEFF.GE.9000)xh=0.
      WRITE(6,60)NMODEL,ITEFF,GLOG,XSCALE,VTURB,XH,EBVI(IRED),UMAG,
     1CAMAG,VMAG,BMAG,YMAG,UMINB,BMINY,M1,C1,hk
   60 FORMAT(I6,I6,5F6.2,10F8.3)
      WRITE(7,60)NMODEL,ITEFF,GLOG,XSCALE,VTURB,XH,EBVI(IRED),UMAG,
     1CAMAG,VMAG,BMAG,YMAG,UMINB,BMINY,M1,C1,hk
      IF(IRED.EQ.1)
     1WRITE(8,60)NMODEL,ITEFF,GLOG,XSCALE,VTURB,XH,EBVI(IRED),UMAG,
     2CAMAG,VMAG,BMAG,YMAG,UMINB,BMINY,M1,C1,hk
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
