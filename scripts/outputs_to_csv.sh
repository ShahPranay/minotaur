OUTPUTSDIR="/home/pranay/iitd/sem8/btp/minotaur/build/outputs"
csvfile=$1

echo "Problem,Time Taken(s),QG Status" >> ${csvfile}

for filepath in ${OUTPUTSDIR}/*.txt
do
  timetaken=$(cat ${filepath} | grep 'QG: time used' | grep -Eo '[0-9]*\.?[0-9]*')
  qgstatus=$(cat ${filepath} | grep 'QG: status' | grep -Eo '=(.*)')
  echo "$(basename "${filepath}" .txt),${timetaken},${qgstatus:2}" >> ${csvfile}
done
