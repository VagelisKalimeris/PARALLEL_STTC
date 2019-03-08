make clean
make
clear

date +"%F %H:%M:%S:%N"

time ./sttc 500 9 r1r_027_3
time ./sttc 500 9 r1r_041_3

time ./sttc 500 6 r5l_033_3
time ./sttc 500 8 r5l_034_3

date +"%F %H:%M:%S:%N"

time ./sttc 500 2 FS13_SA_OGB_adult_reduced_th2
time ./sttc 500 2 PL3_SA_OGB_adult_reduced_th2
time ./sttc 500 2 PL5_SA_OGB_adult_reduced_th2
time ./sttc 500 2 PL12_SA_OGB_adult_reduced_th2

date +"%F %H:%M:%S:%N"