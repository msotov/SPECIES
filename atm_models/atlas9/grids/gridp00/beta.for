      PROGRAM HBETA
      DIMENSION W(5000),F(5000)
      DIMENSION BIG(5000),SMALL(5000),SMALLFILT(5000),BIGFILT(5000)
      DIMENSION WAVEBETA(125),HLAM(125)
      DIMENSION FILTERBIG(51),WAVEBIG(51),FILTERSMALL(39),WAVESMALL(39)
      CHARACTER*80 TITLE,TITLEA,TITLEB
      DATA WAVEBIG/
     1 464.,465.,466.,467.,468.,469.,470.,471.,472.,473.,474.,475.,
     2 476.,477.,478.,479.,480.,481.,482.,483.,484.,485.,486.,487.,
     3 488.,489.,490.,491.,492.,493.,494.,495.,496.,497.,498.,499.,
     4 500.,501.,502.,503.,504.,505.,506.,507.,508.,509.,510.,511.,
     5 512.,513.,514./
      DATA FILTERBIG/
     1 .000,.000,.002,.005,.007,.009,.011,.014,.019,.026,.036,.049,
     2 .064,.098,.138,.195,.274,.368,.473,.551,.604,.640,.664,.683,
     3 .698,.701,.701,.683,.645,.563,.428,.325,.244,.197,.150,.109,
     4 .086,.066,.056,.051,.046,.041,.035,.030,.025,.020,.015,.010,
     5 .005,.000,.000/
      DATA WAVESMALL/
     1 476.0,476.5,477.0,477.5,478.0,478.5,479.0,479.5,480.0,480.5,
     2 481.0,481.5,482.0,482.5,483.0,483.5,484.0,484.5,485.0,485.5,
     3 486.0,486.5,487.0,487.5,488.0,488.5,489.0,489.5,490.0,490.5,
     4 491.0,491.5,492.0,492.5,493.0,493.5,494.0,494.5,495.0/
      DATA FILTERSMALL/
     1  .000, .000, .000, .001, .002, .003, .004, .006, .008, .010,
     2  .014, .019, .032, .056, .113, .173, .263, .315, .480, .578,
     3  .623, .525, .453, .285, .139, .086, .053, .032, .019, .016,
     4  .014, .011, .009, .007, .005, .002, .000, .000, .000/
      DATA WAVEBETA/
     1 466.268,467.268,468.268,469.268,470.268,471.268,472.268,473.268,
     2 474.268,475.268,476.268,477.268,478.268,479.268,480.268,480.768,
     3 481.268,481.768,482.268,482.668,483.068,483.268,483.468,483.668,
     4 483.868,484.068,484.268,484.468,484.668,484.868,485.068,485.268,
     5 485.368,485.468,485.568,485.668,485.768,485.818,485.868,485.918,
     6 485.968,486.018,486.068,486.088,486.108,486.128,486.148,486.168,
     7 486.178,486.188,486.198,486.208,486.218,486.228,486.238,486.248,
     8 486.258,    486.268,    486.278,486.288,486.298,486.308,486.318,
     9 486.328,486.338,486.348,486.358,486.368,486.388,486.408,486.428,
     A 486.448,486.468,486.518,486.568,486.618,486.668,486.718,486.768,
     1 486.868,486.968,487.068,487.168,487.268,487.468,487.668,487.868,
     2 488.068,488.268,488.468,488.668,488.868,489.068,489.268,489.468,
     3 489.868,490.268,490.758,491.258,491.758,492.258,493.258,494.258,
     4 495.258,496.258,497.258,498.258,499.258,500.258,501.258,502.258,
     5 503.258,504.258,505.258,506.258,507.258,508.258,509.258,510.258,
     6 511.258,512.258,513.258,514.258,515.258,516.258/
C     FILTERS ARE IN AIR.  FLUX IS COMPUTED IN VACUUM.  
C     SUBTRACT .136 FROM FLUX WAVE
      DO 1 NU=1,125
    1 WAVEBETA(NU)=WAVEBETA(NU)-.136
      DO 2 I=1,5000
    2 W(I)=464.+.01*I
      CALL PINTER(WAVESMALL,FILTERSMALL,39,W,SMALLFILT,5000)
      CALL PINTER(WAVEBIG,FILTERBIG,51,W,BIGFILT,5000)
      DO 3 I=1,5000
      RESPONSES=ONEP21(W(I))*AIR(W(I),1.)*REFLCT(W(I))          
      SMALL(I)=MAX(SMALLFILT(I)*RESPONSES,0.)
    3 BIG(I)=MAX(BIGFILT(I)*RESPONSES,0.)
C     [M/H],VTURB,L/H
      READ(5,5)ABUND,VTURB,CONVEC
    5 FORMAT(3F10.3)
      WRITE(6,6)
      WRITE(7,6)
    6 FORMAT('        Teff log g  [M]  Vturb  l/H',
     1'   small     big    beta')
      BIGNORM=0.                                                                  
      SMALLNORM=0.
      DO 13 I=1,5000                                                             
      BIGNORM=BIGNORM+BIG(I)                                                         
   13 SMALLNORM=SMALLNORM+SMALL(I)                                                         
      BIGNORM=BIGNORM*.01                                                            
      SMALLNORM=SMALLNORM*.01                                                            
      SMALLNOMAG=-2.5*ALOG10(SMALLNORM)                                                 
      BIGNOMAG=-2.5*ALOG10(BIGNORM)                                                 
C                                                                               
C     wavelength in nm
      DO 1000 MODEL=1,500                                              
      READ(1,11,END=9)TITLEA,TITLEB
   11 FORMAT(A80/A80)
      TITLE=TITLEA(1:30)//TITLEB(6:55)
      PRINT 713,MODEL,TITLE
  713 FORMAT(I5,1X,A80)
      READ(TITLE,'(5X,I6,10X,F8.5)')ITEFF,GLOG
C     ergs/cm**2/s/hz/ster
      DO 715 NU=1,125
      READ(1,710)HNU
  710 format(4x,5X,9x,20x,e13.4,e13.4)
C 710 format(4x,i5,9x,20x,e13.4,e13.4)
      FREQ=2.99792458E17/WAVEBETA(NU)
  715 HLAM(NU)=HNU*FREQ/WAVEBETA(NU)
      READ(1,710)
      CALL LINTER(WAVEBETA,HLAM,125,W,F,5000)                                    
      BIGF=0.                                                                     
      SMALLF=0.                                                                    
      DO 22 I=1,5000                                                             
      BIGF=BIGF+BIG(I)*F(I)                                                          
   22 SMALLF=SMALLF+SMALL(I)*F(I)                                                       
      BIGF=BIGF*.01                                                                  
      SMALLF=SMALLF*.01                                                                
      BIGMAG=-2.5*ALOG10(BIGF)                                                      
      SMALLMAG=-2.5*ALOG10(SMALLF)                                                    
      BIGMAG=BIGMAG-BIGNOMAG                                                          
      SMALLMAG=SMALLMAG-SMALLNOMAG                                                          
      BETA=SMALLMAG-BIGMAG
C     NORMALIZATION TO VEGA
C     BETA VEGA=2.903
C     
      BETA=BETA+2.613
C        Teff log g  [M]  Vturb  l/H   small     big    beta
C     1  9550  3.95 -0.50  2.00  0.00 -19.099 -19.389   0.290
      XSCALE=ABUND
c     actually l/h
      xh=CONVEC
      if(ItEFF.GE.9000)xh=0.
      WRITE(6,60)MODEL,ITEFF,GLOG,XSCALE,VTURB,XH,SMALLMAG,BIGMAG,BETA
   60 FORMAT(I6,I6,4F6.2,3F8.3)
      WRITE(7,60)MODEL,ITEFF,GLOG,XSCALE,VTURB,XH,SMALLMAG,BIGMAG,BETA
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
