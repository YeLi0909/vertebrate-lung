awk '$4>=0.1 && $5 < -1.65 && $6 > 5 && $ 7 > 1.65 && $8 > 1.65 && $9 > 1.65 && $10 >= 6' merge.out > merge.out.filter
# awk '$4>=0.1 && $5 < -1.65 && $6 > 5 && $10 >= 6' merge.out > merge.out.filter
awk '$2 != "NA" && $4>=0.1 && $10>=6 && $7>1.65 && $8 > 1.65 && $9 > 1.65' merge.out > important_for_lung_link.out
awk '$2=="EN"' important_for_lung_link.out | awk '$13<1e-5'| cut -f 11 | sort | uniq -c | sort -k1 -nr | les
awk '$2=="STROMAL"' important_for_lung_link.out | awk '$13<1e-5'| cut -f 11 | sort | uniq -c | sort -k1 -nr | les
awk '$2=="IMMUNE"' important_for_lung_link.out | awk '$13<1e-5'| cut -f 11 | sort | uniq -c | sort -k1 -nr | les
awk '$2=="EP"' important_for_lung_link.out | awk '$13<1e-5'| cut -f 11 | sort | uniq -c | sort -k1 -nr | les
