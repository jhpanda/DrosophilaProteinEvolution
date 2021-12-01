
for f in `ls *polar*cutoff0.*.out `;do 
    #awk '{if ($6>0 && $7<=0.05) print($0)}' \
    awk '{if ($2+$3>0 && $4+$5>0 && $2+4>0 && $3+$5>0 && $6>0 && $7<=0.05) print($0)}' ${f} > ${f}.positive.dat
    awk '{print(substr($1,0,11))}' ${f}.positive.dat > ${f}.positive.txt
done
