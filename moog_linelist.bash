#/bin/bash

cd ./MOOG_linelist
if [$# -ne 1]; then
  echo "You need to input the file with the EW measurements"
else
  file=$1

  python moog_linelist.py $file

fi
exit
