
for f in `ls unpolar_*-cutoff0.*.out `;do 
    awk '{if ($7<=0.05) print($0)}' \
    ${f} > ${f}.negative.dat
done
