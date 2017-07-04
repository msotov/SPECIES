      PROGRAM COUSINS
      IMPLICIT REAL*4 (A-Z)
      INTEGER I,NR,NI,NC,NV,MODEL,MODEL1,MODEL2,NSKIP,NMODEL,ITEFF
      DIMENSION HLAM(1221),WAVE(1221),HNU(1221),HNUCONT(1221)                                            
      DIMENSION AMAGI(13,1221),TRANSI(13,1221),EBVI(13),AMAG(1221)
      CHARACTER*80 TITLE
      DIMENSION A(8)                                                            
      DIMENSION F(2700)                                                         
      DIMENSION RFILT(14),IFILT(16),CATHODE75(14),VFILT(54)
      DIMENSION RWAVE(14),IWAVE(16),CWAVE(14),VWAVE(54)
      DIMENSION WAVER(820),WAVEI(1200),WAVEC(800),WAVEV(2700)                             
      DIMENSION SR(820),SI(1200),SC(800),SV(2700)
C      DATA EBVI/0.,.1,.2,.3,.4,.5,.6,.7,.8,.9,1./
      DATA EBVI/0.,.1,.2,.3,.4,.5,.6,.8,1.,2.,3.,4.,5./
      DATA NR,NI,NC,NV/820,1200,800,2700/
      DATA VFILT/0.,.03,.084,.163,.301,.458,.630,.780,.895,.967,.997,1.,         
     1 .988,.958,.919,.877,.819,.765,.711,.657,.602,.545,.488,.434,.386,        
     2 .331,.289,.250,.214,.181,.151,.120,.093,.069,.051,.036,.027,.021,        
     3 .018,.016,.014,.012,.011,.010,.009,.008,.007,.006,.005,.004,.003,        
     4 .002,.001,.000/                                                          
      DATA RFILT/.00,.01,.09,.46,.71,.77,.77,.65,.45,.25,.11,.04,.01,
     1 .00/
      DATA IFILT/.00,.06,.26,.49,.68,.78,.84,.86,.85,.84,.80,.61,.19,
     1 .01,.00,.00/
      DATA CATHODE75/.55,.72,.92,.99,1.00,.99,.97,.94,.78,.47,.10,.03,
     1 .01,.00/
      DATA VWAVE/475.,480.,485.,490.,495.,500.,505.,510.,515.,520.,525.,         
     1 530.,535.,540.,545.,550.,555.,560.,565.,570.,575.,580.,585.,590.,        
     2 595.,600.,605.,610.,615.,620.,625.,630.,635.,640.,645.,650.,655.,        
     3 660.,665.,670.,675.,680.,685.,690.,695.,700.,705.,710.,715.,720.,        
     4 725.,730.,735.,740./                                                     
      DATA RWAVE/540.,550.,560.,570.,580.,590.,600.,650.,700.,750.,
     1 800.,850.,900.,950./
      DATA IWAVE/700.,710.,720.,730.,740.,750.,770.,790.,810.,850.,
     1 900.,1000.,1100.,1200.,1300.,1400./
      DATA CWAVE/500.,600.,700.,800.,810.,820.,830.,840.,850.,860.,
     1 870.,880.,890.,900./
   77 FORMAT(10E12.4)                                                           
C      PRINT 77,RFILT                                                               
C      PRINT 77,IFILT                                                               
C      PRINT 77,CATHODE75                                                               
C      PRINT 77,VFILT                                                               
C      PRINT 77,RWAVE                                                            
C      PRINT 77,IWAVE                                                            
C      PRINT 77,CWAVE                                                            
C      PRINT 77,VWAVE
      DO 14 I=1,NV                                                            
   14 WAVEV(I)=470.+FLOAT(I)*.1                                                 
      DO 12 I=1,NR                                                            
   12 WAVER(I)=540.+FLOAT(I)*.5                                                 
      DO 13 I=1,NI                                                            
   13 WAVEI(I)=700.+FLOAT(I)*.5                                                 
      DO 114 I=1,NC
  114 WAVEC(I)=500.-.5+FLOAT(I)*.5
      CALL PINTER(RWAVE,RFILT,14,WAVER,SR,NR)                                   
      CALL PINTER(IWAVE,IFILT,16,WAVEI,SI,NI)                                   
      CALL PINTER(CWAVE,CATHODE75,14,WAVEC,SC,NC)                                   
      CALL PINTER(VWAVE,VFILT,54,WAVEV,SV,NV)                                   
C      PRINT 77,SR                                                              
C      PRINT 77,SI                                                              
C      PRINT 77,SC                                                              
C      PRINT 77,SV                                                              
      RNORM=0.                                                                  
      DO 15 I=1,NR                                                            
      SR(I)=SR(I)*CATHODE75(I+81)
   15 RNORM=RNORM+SR(I)                                                      
      RNORM=RNORM*.5                                                          
      INORM=0.                                                                  
      DO 16 I=1,NI                                                            
      SI(I)=SI(I)*CATHODE75(I+401)
   16 INORM=INORM+SI(I)                                                        
      INORM=INORM*.5                                                            
      VNORM=0.                                                                  
      DO 17 I=1,NV                                                            
   17 VNORM=VNORM+SV(I)                                                        
      VNORM=VNORM*.1                                                            
      RNOMAG=-2.5*ALOG10(RNORM)                                                
      INOMAG=-2.5*ALOG10(INORM)                                                
      VNOMAG=-2.5*ALOG10(VNORM)                                                 
C                                                                               
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
     1'   V       I       R       V-R     V-I')
C                                                                               
C                                                                               
      DO 1000 NMODEL=1,1000                                              
C     ergs/cm**2/s/hz/ster
      READ(1,2,END=9)TITLE
    2 FORMAT(A80)
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
      CALL LINTER(WAVE,HLAM,NNU,WAVER,F,NR)                                   
C     PRINT 77,F                                                                
      RF=0.                                                                    
      DO 22 I=1,NR                                                            
   22 RF=RF+SR(I)*F(I)                                                       
      RF=RF*.5                                                                
C     PRINT 77,RF                                                              
      CALL LINTER(WAVE,HLAM,NNU,WAVEI,F,NI)                                   
C     PRINT 77,F                                                                
      IF=0.                                                                    
      DO 32 I=1,NI                                                            
   32 IF=IF+SI(I)*F(I)                                                       
      IF=IF*.5                                                                
C     PRINT 77,IF                                                              
      CALL LINTER(WAVE,HLAM,NNU,WAVEV,F,NV)                                   
C     PRINT 77,F                                                                
      VF=0.                                                                    
      DO 42 I=1,NV                                                            
   42 VF=VF+SV(I)*F(I)                                                       
      VF=VF*.1                                                                
      RMAG=-2.5*ALOG10(RF)                                                    
      IMAG=-2.5*ALOG10(IF)                                                    
      VMAG=-2.5*ALOG10(VF)                                                    
C      PRINT 77,RNORM,INORM,VNORM                                        
C      PRINT 77,RNOMAG,INOMAG,VNOMAG                                      
C      PRINT 77,RMAG,IMAG,VMAG                              
      RMAG=RMAG-RNOMAG                                    
      IMAG=IMAG-INOMAG                                                        
      VMAG=VMAG-VNOMAG                                                        
      VMINR=VMAG-RMAG+0.585
      VMINI=VMAG-IMAG+1.527
C     I do not know the source of these data  Bruce Carney? Ruth Peterson?
C     colors for Vega   -0.009  -0.005
      XSCALE=ABUND
c     actually l/h
      xh=CONVEC
      if(ItEFF.GE.9000)xh=0.
      WRITE(6,60)NMODEL,ITEFF,GLOG,XSCALE,VTURB,XH,EBVI(IRED),
     1VMAG,RMAG,IMAG,VMINR,VMINI
      WRITE(7,60)NMODEL,ITEFF,GLOG,XSCALE,VTURB,XH,EBVI(IRED),
     1VMAG,RMAG,IMAG,VMINR,VMINI
      IF(IRED.EQ.1)
     1WRITE(8,60)NMODEL,ITEFF,GLOG,XSCALE,VTURB,XH,EBVI(IRED),
     2VMAG,RMAG,IMAG,VMINR,VMINI
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
