      PROGRAM UBV
      DIMENSION HLAM(1221),WAVE(1221),HNU(1221),HNUCONT(1221)                                            
      DIMENSION AMAGI(13,1221),TRANSI(13,1221),EBVI(13),AMAG(1221)
      CHARACTER*80 TITLE
      DIMENSION WLIT(1221)
      DIMENSION A(8)                                                            
      DIMENSION SB0(2100),SV0(2600)                                             
      DIMENSION SU1(1300),SB1(2100),SV1(2700)                                   
      DIMENSION SU2(1300)                                                       
      DIMENSION WAVEU(1300),WAVEB(2100),WAVEV(2700)                             
      DIMENSION F(2700)                                                         
      DIMENSION UWAVE(14),BWAVE(22),VWAVE(27)                                   
      DIMENSION B0(22),V0(27)                                                   
      DIMENSION U1(14),B1(22),V1(27)                                            
      DIMENSION U2(14)                                                          
      DIMENSION VAS(54),BAS(40),BPAS(40),UPAS(24),UBUSER(24)                    
      DIMENSION WVAS(54),WBAS(40),WUAS(24)                                      
      COMMON SVAS(2700),SBAS(2100),SBPAS(2100),SUPAS(1300),SUBUS(1300)          
C      DATA EBVI/0.,.1,.2,.3,.4,.5,.6,.7,.8,.9,1./
      DATA EBVI/0.,.1,.2,.3,.4,.5,.6,.8,1.,2.,3.,4.,5./
      DATA VAS/0.,.030,.084,.163,.301,.458,.630,.780,.895,.967,.997,1.,         
     1 .988,.958,.919,.877,.819,.765,.711,.657,.602,.545,.488,.434,.386,        
     2 .331,.289,.250,.214,.181,.151,.120,.093,.069,.051,.036,.027,.021,        
     3 .018,.016,.014,.012,.011,.010,.009,.008,.007,.006,.005,.004,.003,        
     4 .002,.001,.000/                                                          
      DATA BAS/0.,.006,.030,.060,.134,.302,.567,.841,.959,.983,.996,1.,         
     1 .996,.987,.974,.957,.931,.897,.849,.800,.748,.698,.648,.597,.545,        
     2 .497,.447,.397,.345,.297,.252,.207,.166,.129,.095,.069,.043,.024,        
     3 .009,.0/                                                                 
      DATA BPAS/0.,.006,.023,.045,.106,.254,.492,.752,.881,.923,.955,           
     1 .977,.990,1.,1.,.997,.984,.958,.916,.871,.820,.775,.723,.672,.617        
     2,.569,.511,.457,.402,.347,.299,.244,.199,.154,.113,.084,.051,.029,        
     3 .010,0./                                                                 
      DATA UPAS/0.,.009,.028,.071,.127,.208,.307,.415,.524,.632,.740,           
     1 .840,.929,.995,.990,.830,.613,.406,.212,.090,.033,.014,0.,0./            
      DATA UBUSER/0.,.020,.077,.135,.204,.282,.385,.493,.600,.705,.820,         
     1 .900,.959,.993,1.,.975,.850,.645,.400,.223,.125,.057,.005,0./            
      DATA WVAS/475.,480.,485.,490.,495.,500.,505.,510.,515.,520.,525.,         
     1 530.,535.,540.,545.,550.,555.,560.,565.,570.,575.,580.,585.,590.,        
     2 595.,600.,605.,610.,615.,620.,625.,630.,635.,640.,645.,650.,655.,        
     3 660.,665.,670.,675.,680.,685.,690.,695.,700.,705.,710.,715.,720.,        
     4 725.,730.,735.,740./                                                     
      DATA WBAS/360.,365,370.,375.,380.,385.,390.,395.,400.,405.,410.,          
     1 415.,420.,425.,430.,435.,440.,445.,450.,455.,460.,465.,470.,475.,        
     2 480.,485.,490.,495.,500.,505.,510.,515.,520.,525.,530.,535.,540.,        
     3 545.,550.,555./                                                          
      DATA WUAS/305.,310.,315.,320.,325.,330.,335.,340.,345.,350.,355.,         
     1 360.,365.,370.,375.,380.,385.,390.,395.,400.,405.,410.,415.,420./        
      DATA U1/0.,.025,.250,.680,1.137,1.650,2.006,2.250,2.337,1.925,            
     1 .650,.197,.070,0./                                                       
      DATA U2/0.,.025,.060,.170,.375,.675,1.000,1.250,1.390,1.125,.600,         
     1 .140,.030,.0/                                                            
      DATA B0/0.,.015,.100,.500,1.800,3.620,3.910,4.000,3.980,3.780,            
     1 3.500,3.150,2.700,2.320,1.890,1.530,1.140,.750,.500,.250,.070,0./        
      DATA B1/0.,.006,.080,.337,1.425,2.253,2.806,2.950,3.000,2.937,            
     1 2.780,2.520,2.230,1.881,1.550,1.275,.975,.695,.430,.210,.055,0./         
      DATA V0/0.,.020,.280,1.180,2.170,2.970,3.300,3.250,3.000,2.700,           
     1 2.320,2.000,1.670,1.280,.920,.650,.430,.260,.175,.125,.100,.070,         
     2 .050,.030,.020,.010,0./                                                  
      DATA V1/0.,.020,.175,.900,1.880,2.512,2.850,2.820,2.625,2.370,            
     1 2.050,1.720,1.413,1.068,.795,.567,.387,.250,.160,.110,.081,.061,         
     1 .045,.028,.017,.007,0./                                                  
      DATA UWAVE/290.,300.,310.,320.,330.,340.,350.,360.,370.,380.,             
     1 390.,400.,410.,420./                                                     
      DATA BWAVE/350.,360.,370.,380.,390.,400.,410.,420.,430.,440.,450.,        
     1 460.,470.,480.,490.,500.,510.,520.,530.,540.,550.,560./                  
      DATA VWAVE/470.,480.,490.,500.,510.,520.,530.,540.,550.,560.,570.,        
     1 580.,590.,600.,610.,620.,630.,640.,650.,660.,670.,680.,690.,700.,        
     2 710.,720.,730./                                                          
   77 FORMAT(10E12.4)                                                           
C     PRINT 77,B0                                                               
C     PRINT 77,V0                                                               
C     PRINT 77,U1                                                               
C     PRINT 77,B1                                                               
C     PRINT 77,V1                                                               
C     PRINT 77,UWAVE                                                            
C     PRINT 77,BWAVE                                                            
C     PRINT 77,VWAVE                                                            
      DO 12 I=1,1300                                                            
   12 WAVEU(I)=290.+FLOAT(I)*.1                                                 
      DO 13 I=1,2100                                                            
   13 WAVEB(I)=350.+FLOAT(I)*.1                                                 
      DO 14 I=1,2700                                                            
   14 WAVEV(I)=470.+FLOAT(I)*.1                                                 
      CALL PINTER(UWAVE,U2,14,WAVEU,SU2,1300)                                   
      CALL PINTER(UWAVE,U1,14,WAVEU,SU1,1300)                                   
      CALL PINTER(BWAVE,B1,22,WAVEB,SB1,2100)                                   
      CALL PINTER(VWAVE,V1,27,WAVEV,SV1,2700)                                   
      CALL PINTER(BWAVE,B0,22,WAVEB,SB0,2100)                                   
      CALL PINTER(VWAVE,V0,27,WAVEV,SV0,2700)                                   
      CALL PINTER(WUAS,UPAS,24,WAVEU,SUPAS,1300)                                
      CALL PINTER(WUAS,UBUSER,24,WAVEU,SUBUS,1300)                              
      CALL PINTER(WBAS,BAS,40,WAVEB,SBAS,2100)                                  
      CALL PINTER(WBAS,BPAS,40,WAVEB,SBPAS,2100)                                
      CALL PINTER(WVAS,VAS,54,WAVEV,SVAS,2700)                                  
      DO 70 I=1,150                                                             
      SUPAS(I)=0.                                                               
   70 SUBUS(I)=0.                                                               
      DO 71 I=1,124                                                             
      SBAS(I)=0.                                                                
   71 SBPAS(I)=0.                                                               
      DO 72 I=1,50                                                              
   72 SVAS(I)=0.                                                                
      DO 73 I=2601,2700                                                         
      SV0(I)=0.                                                                 
   73 SV1(I)=0.                                                                 
      DO 74 I=2051,2100                                                         
      SBPAS(I)=0.                                                               
   74 SBAS(I)=0.                                                                
      DO 75 I=1261,1300                                                         
   75 SUBUS(I)=0.                                                               
      DO 76 I=1251,1300                                                         
   76 SUPAS(I)=0.                                                               
C     PRINT 77,SUPAS                                                            
C     PRINT 77,SUBUS                                                            
C     PRINT 77,SBAS                                                             
C     PRINT 77,SBPAS                                                            
C     PRINT 77,SVAS                                                             
C     PRINT 77,SU1                                                              
C     PRINT 77,SB1                                                              
C     PRINT 77,SV1                                                              
C     PRINT 77,SB0                                                              
C     PRINT 77,SV0                                                              
      UNORM=0.                                                                  
      U1NORM=0.                                                                 
      U2NORM=0.                                                                 
      UANORM=0.                                                                 
      UBNORM=0.                                                                 
      DO 15 I=1,1300                                                            
      UANORM=UANORM+SUPAS(I)                                                    
      UBNORM=UBNORM+SUBUS(I)                                                    
      U1NORM=U1NORM+SU1(I)                                                      
   15 U2NORM=U2NORM+SU2(I)                                                      
      U1NORM=U1NORM*.1                                                          
      U2NORM=U2NORM*.1                                                          
      UANORM=UANORM*.1                                                          
      UBNORM=UBNORM*.1                                                          
      BANORM=0.                                                                 
      BAPNORM=0.                                                                
      BNORM=0.                                                                  
      DO 16 I=1,2100                                                            
      BANORM=BANORM+SBAS(I)                                                     
      BAPNORM=BAPNORM+SBPAS(I)                                                  
   16 BNORM=BNORM+SB0(I)                                                        
      BANORM=BANORM*.1                                                          
      BAPNORM=BAPNORM*.1                                                        
      BNORM=BNORM*.1                                                            
      VANORM=0.                                                                 
      VNORM=0.                                                                  
      DO 17 I=1,2700                                                            
      VANORM=VANORM+SVAS(I)                                                     
   17 VNORM=VNORM+SV0(I)                                                        
      VANORM=VANORM*.1                                                          
      VNORM=VNORM*.1                                                            
      U1NOMG=-2.5*ALOG10(U1NORM)                                                
      U2NOMG=-2.5*ALOG10(U2NORM)                                                
      BNOMAG=-2.5*ALOG10(BNORM)                                                 
      VNOMAG=-2.5*ALOG10(VNORM)                                                 
      UANOMG=-2.5*ALOG10(UANORM)                                                
      UBNOMG=-2.5*ALOG10(UBNORM)                                                
      BANOMG=-2.5*ALOG10(BANORM)                                                
      BAPNOMG=-2.5*ALOG10(BAPNORM)                                              
      VANOMG=-2.5*ALOG10(VANORM)                                                
C                                                                               
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
     1'    U       B       V      BOL     BC      U-B     B-V')
C     wavelength in nm
      DO 1000 NMODEL=1,1000                                              
C     ergs/cm**2/s/hz/ster
      READ(1,712,END=9)TITLE
  712 FORMAT(A80)
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
      CALL LINTER(WAVE,HLAM,NNU,WAVEU,F,1300)                                   
C     PRINT 77,F                                                                
      U1F=0.                                                                    
      U2F=0.                                                                    
      UAF=0.                                                                    
      UBF=0.                                                                    
      DO 22 I=1,1300                                                            
      U1F=U1F+SU1(I)*F(I)                                                       
      U2F=U2F+SU2(I)*F(I)                                                       
      UAF=UAF+SUPAS(I)*F(I)                                                     
      UBF=UBF+SUBUS(I)*F(I)                                                     
   22 CONTINUE                                                                  
      U1F=U1F*.1                                                                
      U2F=U2F*.1                                                                
      UAF=UAF*.1                                                                
      UBF=UBF*.1                                                                
C     PRINT 77,U1F                                                              
C     PRINT 77,U2F                                                              
      CALL LINTER(WAVE,HLAM,NNU,WAVEB,F,2100)                                   
C     PRINT 77,F                                                                
      B0F=0.                                                                    
      B1F=0.                                                                    
      BAF=0.                                                                    
      BAPF=0.                                                                   
      DO 32 I=1,2100                                                            
      BAF=BAF+SBAS(I)*F(I)                                                      
      BAPF=BAPF+SBPAS(I)*F(I)                                                   
      B0F=B0F+SB0(I)*F(I)                                                       
   32 B1F=B1F+SB1(I)*F(I)                                                       
      B0F=B0F*.1                                                                
      B1F=B1F*.1                                                                
      BAF=BAF*.1                                                                
      BAPF=BAPF*.1                                                              
C     PRINT 77,B0F                                                              
C     PRINT 77,B1F                                                              
      CALL LINTER(WAVE,HLAM,NNU,WAVEV,F,2700)                                   
C     PRINT 77,F                                                                
      V0F=0.                                                                    
      V1F=0.                                                                    
      VAF=0.                                                                    
      DO 42 I=1,2700                                                            
      VAF=VAF+SVAS(I)*F(I)                                                      
      V0F=V0F+SV0(I)*F(I)                                                       
   42 V1F=V1F+SV1(I)*F(I)                                                       
      VAF=VAF*.1                                                                
      V0F=V0F*.1                                                                
      V1F=V1F*.1                                                                
      U1MAG=-2.5*ALOG10(U1F)                                                    
      U2MAG=-2.5*ALOG10(U2F)                                                    
      B0MAG=-2.5*ALOG10(B0F)                                                    
      B1MAG=-2.5*ALOG10(B1F)                                                    
      V0MAG=-2.5*ALOG10(V0F)                                                    
      V1MAG=-2.5*ALOG10(V1F)                                                    
      UAMAG=-2.5*ALOG10(UAF)                                                    
      UBMAG=-2.5*ALOG10(UBF)                                                    
      BAMAG=-2.5*ALOG10(BAF)                                                    
      BAPMAG=-2.5*ALOG10(BAPF)                                                  
      VAMAG=-2.5*ALOG10(VAF)                                                    
C     PRINT 77,U1NORM,U2NORM,BNORM,VNORM                                        
C     PRINT 77,U1NOMG,U2NOMG,BNOMAG,VNOMAG                                      
C     PRINT 77,U1MAG,U2MAG,B0MAG,B1MAG,V0MAG,V1MAG                              
C     PRINT 77,UANORM,UBNORM,BANORM,BAPNORM,VANORM                              
C     PRINT 77,UANOMG,UBNOMG,BANOMG,BAPNOMG,VANOMG                              
C     PRINT 77,UAMAG,UBMAG,BAMAG,BAPMAG,VAMAG                                   
      U0MAG=2.*(U1MAG-U1NOMG)-(U2MAG-U2NOMG)                                    
      B0MAG=B0MAG-BNOMAG                                                        
      V0MAG=V0MAG-VNOMAG                                                        
      UMINB1=U1MAG-B1MAG                                                        
      BMINV1=B1MAG-V1MAG                                                        
      BMINV=1.024*BMINV1+.81                                                    
      UMINB=.921*UMINB1-1.308                                                   
      UMINBA=UAMAG-BAPMAG                                                       
      UMINBB=UBMAG-BAPMAG                                                       
      UMINBA=UMINBA-1.33                                                        
      BMINVA=BAMAG-VAMAG                                                        
      BMINVA=BMINVA+0.67                                                        
      UMINBB=UMINBB-1.093                                                       
      BMINVB=BAMAG-VAMAG                                                        
      BMINVB=BMINVB+.710                                                        
      UBMAG=UBMAG-UBNOMG                                                        
      UAMAG=UAMAG-UANOMG                                                        
      BAMAG=BAMAG-BANOMG                                                        
      BAPMAG=BAPMAG-BAPNOMG                                                     
      VAMAG=VAMAG-VANOMG                                                        
C     NORMALIZATION TO VEGA
      UMINBB=UMINBB+.011
      BMINVB=BMINVB+.012
c     VEGA  -.005 -.003
      TEFF=ITEFF
      BOL=-2.5*ALOG10(5.66956E-5/3.14159*TEFF**4)
      BC=V0MAG-BOL
      BC=BC-8.509
      BC=-BC
      BC=BC-.125
      BCA=VAMAG-BOL
      BCA=BCA-8.365
      BCA=-BCA
c     new smallest bc -.115      
c     should be 0
      BCA=BCA+.115
      XSCALE=ABUND
c     actually l/h
      xh=CONVEC
      if(ItEFF.GE.9000)xh=0.
      WRITE(6,60)NMODEL,ITEFF,GLOG,XSCALE,VTURB,XH,EBVI(IRED),
     1UBMAG,BAMAG,VAMAG,BOL,BCA,UMINBB,BMINVB        
      WRITE(7,60)NMODEL,ITEFF,GLOG,XSCALE,VTURB,XH,EBVI(IRED),
     1UBMAG,BAMAG,VAMAG,BOL,BCA,UMINBB,BMINVB        
      IF(IRED.EQ.1)
     1WRITE(8,60)NMODEL,ITEFF,GLOG,XSCALE,VTURB,XH,EBVI(IRED),
     2UBMAG,BAMAG,VAMAG,BOL,BCA,UMINBB,BMINVB        
   60 FORMAT(I6,I6,5F6.2,7F8.3)                                                    
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
