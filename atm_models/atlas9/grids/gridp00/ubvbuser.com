$SET DEF KUR5:[KURUCZ.COLORSK]
$ASSIGN NEWCOLORS:REDDEN.DAT FOR002
$ASSIGN POOL10:FP00K0.PCK FOR001
$ASSIGN POOLKU1:UBVP00K0.RED FOR007
$ASSIGN POOLKU1:UBVP00K0.DAT FOR008
$RUN NEWCOLORS:UBVBUSER
$ASSIGN POOL10:FP00K1.PCK FOR001
$ASSIGN POOLKU1:UBVP00K1.RED FOR007
$ASSIGN POOLKU1:UBVP00K1.DAT FOR008
$RUN NEWCOLORS:UBVBUSER
$ASSIGN POOL9:FP00K2.PCK FOR001
$ASSIGN POOLKU1:UBVP00K2.RED FOR007
$ASSIGN POOLKU1:UBVP00K2.DAT FOR008
$RUN NEWCOLORS:UBVBUSER
$ASSIGN POOL9:FP00K4.PCK FOR001
$ASSIGN POOLKU1:UBVP00K4.RED FOR007
$ASSIGN POOLKU1:UBVP00K4.DAT FOR008
$RUN NEWCOLORS:UBVBUSER
$ASSIGN POOL9:FP00K8.PCK FOR001
$ASSIGN POOLKU1:UBVP00K8.RED FOR007
$ASSIGN POOLKU1:UBVP00K8.DAT FOR008
$RUN NEWCOLORS:UBVBUSER
