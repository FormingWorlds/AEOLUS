prep_spec <<EOF
sp_b20_HITRAN_a16_sw
20
16
1
2
3
4
5
6
7
8
9
10
11
12
13
23
24
25
18
1
1
2
2
6
6
7
7
13
13
23
23
2
6
2
23
2
24
6
24
7
2
7
13
13
1
13
6
13
23
13
24
23
6
23
24
0
c
1
500
500
1000
1000
1500
1500
2000
2000
2500
2500
3000
3000
4000
4000
5000
5000
6000
6000
7000
7000
8000
8000
9000
9000
10000
10000
11000
11000
14000
14000
17000
17000
20000
20000
23000
23000
26000
26000
29000
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
n
-1
EOF
echo '###############   H2O lines - HITRAN  ###############'
Ccorr_k -F pt_grid_lbl.dat -D /home/x_ryabo/spectral_files/dat_hitran/h2o_data.par -R 1 20 -c 3000.0 -i 0.1 -l 1 1.0e1 -b 1.5e-3 -s sp_b20_HITRAN_a16_sw +S /home/x_ryabo/spectral_files/stellar_spectra/TRAPPIST-1_soc.txt -lk -k -o h2o_l20_l -m h2o_l20_lm -L h2o_lbl_lwf_pt48.nc -sm h2o_l20_1-20map.nc -np 16
prep_spec <<EOF
sp_b20_HITRAN_a16_sw
a
5
y
h2o_l20_l
-1
EOF
touch temp/done_h2o_lbl
cp -rf sp_b20_HITRAN_a16_sw temp/sp_b20_HITRAN_a16_sw_1
cp -rf sp_b20_HITRAN_a16_sw_k temp/sp_b20_HITRAN_a16_sw_1_k
echo '###############   CO2 lines – HITRAN  ###############'
Ccorr_k -F pt_grid_lbl.dat -D /home/x_ryabo/spectral_files/dat_hitran/co2_data.par -R 1 20 -c 3000.0 -i 0.1 -l 2 1.0e1 -t 1.0e-2 -k -s sp_b20_HITRAN_a16_sw +S /home/x_ryabo/spectral_files/stellar_spectra/TRAPPIST-1_soc.txt -lk -o co2_o -m co2_m -L co2_lbl.nc -np 16
prep_spec <<EOF
sp_b20_HITRAN_a16_sw
a
5
y
co2_o
-1
EOF
touch temp/done_co2_lbl
rsync -av sp_b20_HITRAN_a16_sw temp/sp_b20_HITRAN_a16_sw_2
rsync -av sp_b20_HITRAN_a16_sw_k temp/sp_b20_HITRAN_a16_sw_2_k
echo '###############   O3 lines – HITRAN    ###############'
Ccorr_k -F pt_grid_lbl.dat -D /home/x_ryabo/spectral_files/dat_hitran/o3_hitran.par -R 1 20 -c 3000.0 -i 0.1 -l 3 1.0e1 -t 1.0e-2 -k -s sp_b20_HITRAN_a16_sw +S /home/x_ryabo/spectral_files/stellar_spectra/TRAPPIST-1_soc.txt -lk -o o3_o -m o3_m -L o3_lbl.nc -np 16
prep_spec <<EOF
sp_b20_HITRAN_a16_sw
a
5
y
o3_o
-1
EOF
touch temp/done_o3_lbl
rsync -av sp_b20_HITRAN_a16_sw temp/sp_b20_HITRAN_a16_sw_3
rsync -av sp_b20_HITRAN_a16_sw_k temp/sp_b20_HITRAN_a16_sw_3_k
echo '###############   N2O lines – HITRAN  ###############'
Ccorr_k -F pt_grid_lbl.dat -D /home/x_ryabo/spectral_files/dat_hitran/n2o_hitran.par -R 1 20 -c 3000.0 -i 0.1 -l 4 1.0e1 -t 1.0e-2 -k -s sp_b20_HITRAN_a16_sw +S /home/x_ryabo/spectral_files/stellar_spectra/TRAPPIST-1_soc.txt -lk -o n2o_o -m n2o_m -L n2o_lbl.nc -np 16
prep_spec <<EOF
sp_b20_HITRAN_a16_sw
a
5
y
n2o_o
-1
EOF
touch temp/done_n2o_lbl
rsync -av sp_b20_HITRAN_a16_sw temp/sp_b20_HITRAN_a16_sw_4
rsync -av sp_b20_HITRAN_a16_sw_k temp/sp_b20_HITRAN_a16_sw_4_k
echo '###############   CO lines – HITRAN    ###############'
Ccorr_k -F pt_grid_lbl.dat -D /home/x_ryabo/spectral_files/dat_hitran/co_hitran.par -R 1 20 -c 3000.0 -i 0.1 -l 5 1.0e1 -t 1.0e-2 -k -s sp_b20_HITRAN_a16_sw +S /home/x_ryabo/spectral_files/stellar_spectra/TRAPPIST-1_soc.txt -lk -o co_o -m co_m -L co_lbl.nc -np 16
prep_spec <<EOF
sp_b20_HITRAN_a16_sw
a
5
y
co_o
-1
EOF
touch temp/done_co_lbl
rsync -av sp_b20_HITRAN_a16_sw temp/sp_b20_HITRAN_a16_sw_5
rsync -av sp_b20_HITRAN_a16_sw_k temp/sp_b20_HITRAN_a16_sw_5_k
echo '###############   CH4 lines – HITRAN   ###############'
Ccorr_k -F pt_grid_lbl.dat -D /home/x_ryabo/spectral_files/dat_hitran/ch4_hitran.par -R 1 20 -c 3000.0 -i 0.1 -l 6 1.0e1 -t 1.0e-2 -k -s sp_b20_HITRAN_a16_sw +S /home/x_ryabo/spectral_files/stellar_spectra/TRAPPIST-1_soc.txt -lk -o ch4_o -m ch4_m -L ch4_lbl.nc -np 16
prep_spec <<EOF
sp_b20_HITRAN_a16_sw
a
5
y
ch4_o
-1
EOF
touch temp/done_ch4_lbl
rsync -av sp_b20_HITRAN_a16_sw temp/sp_b20_HITRAN_a16_sw_6
rsync -av sp_b20_HITRAN_a16_sw_k temp/sp_b20_HITRAN_a16_sw_6_k
echo '###############   O2 lines – HITRAN    ###############'
Ccorr_k -F pt_grid_lbl.dat -D /home/x_ryabo/spectral_files/dat_hitran/o2_hitran.par -R 1 20 -c 3000.0 -i 0.1 -l 7 1.0e1 -t 1.0e-2 -k -s sp_b20_HITRAN_a16_sw +S /home/x_ryabo/spectral_files/stellar_spectra/TRAPPIST-1_soc.txt -lk -o o2_o -m o2_m -L o2_lbl.nc -np 16
prep_spec <<EOF
sp_b20_HITRAN_a16_sw
a
5
y
o2_o
-1
EOF
touch temp/done_o2_lbl
rsync -av sp_b20_HITRAN_a16_sw temp/sp_b20_HITRAN_a16_sw_7
rsync -av sp_b20_HITRAN_a16_sw_k temp/sp_b20_HITRAN_a16_sw_7_k
echo '###############   NO lines – HITRAN   ###############'
Ccorr_k -F pt_grid_lbl.dat -D /home/x_ryabo/spectral_files/dat_hitran/no_hitran.par -R 1 20 -c 3000.0 -i 0.1 -l 8 1.0e1 -t 1.0e-2 -k -s sp_b20_HITRAN_a16_sw +S /home/x_ryabo/spectral_files/stellar_spectra/TRAPPIST-1_soc.txt -lk -o no_o -m no_m -L no_lbl.nc -np 16
prep_spec <<EOF
sp_b20_HITRAN_a16_sw
a
5
y
no_o
-1
EOF
touch temp/done_no_lbl
rsync -av sp_b20_HITRAN_a16_sw temp/sp_b20_HITRAN_a16_sw_8
rsync -av sp_b20_HITRAN_a16_sw_k temp/sp_b20_HITRAN_a16_sw_8_k
echo '###############   SO2 lines – HITRAN    ###############'
Ccorr_k -F pt_grid_lbl.dat -D /home/x_ryabo/spectral_files/dat_hitran/so2_hitran.par -R 1 20 -c 3000.0 -i 0.1 -l 9 1.0e1 -t 1.0e-2 -k -s sp_b20_HITRAN_a16_sw +S /home/x_ryabo/spectral_files/stellar_spectra/TRAPPIST-1_soc.txt -lk -o so2_o -m so2_m -L so2_lbl.nc -np 16
prep_spec <<EOF
sp_b20_HITRAN_a16_sw
a
5
y
so2_o
-1
EOF
touch temp/done_so2_lbl
rsync -av sp_b20_HITRAN_a16_sw temp/sp_b20_HITRAN_a16_sw_9
rsync -av sp_b20_HITRAN_a16_sw_k temp/sp_b20_HITRAN_a16_sw_9_k
echo '###############   NO2 lines   ###############'
Ccorr_k -F pt_grid_lbl.dat -D /home/x_ryabo/spectral_files/dat_hitran/no2_hitran.par -R 1 20 -c 3000.0 -i 0.1 -l 10 1.0e1 -t 1.0e-2 -k -s sp_b20_HITRAN_a16_sw +S /home/x_ryabo/spectral_files/stellar_spectra/TRAPPIST-1_soc.txt -lk -o no2_o -m no2_m -L no2_lbl.nc -np 16
prep_spec <<EOF
sp_b20_HITRAN_a16_sw
a
5
y
no2_o
-1
EOF
touch temp/done_no2_lbl
rsync -av sp_b20_HITRAN_a16_sw temp/sp_b20_HITRAN_a16_sw_10
rsync -av sp_b20_HITRAN_a16_sw_k temp/sp_b20_HITRAN_a16_sw_10_k
echo '###############   NH3 lines – HITRAN    ###############'
Ccorr_k -F pt_grid_lbl.dat -D /home/x_ryabo/spectral_files/dat_hitran/nh3_hitran.par -R 1 20 -c 3000.0 -i 0.1 -l 11 1.0e1 -t 1.0e-2 -k -s sp_b20_HITRAN_a16_sw +S /home/x_ryabo/spectral_files/stellar_spectra/TRAPPIST-1_soc.txt -lk -o nh3_o -m nh3_m -L nh3_lbl.nc -np 16
prep_spec <<EOF
sp_b20_HITRAN_a16_sw
a
5
y
nh3_o
-1
EOF
touch temp/done_nh3_lbl
rsync -av sp_b20_HITRAN_a16_sw temp/sp_b20_HITRAN_a16_sw_11
rsync -av sp_b20_HITRAN_a16_sw_k temp/sp_b20_HITRAN_a16_sw_11_k
echo '###############   HNO3 lines – HITRAN    ###############'
Ccorr_k -F pt_grid_lbl.dat -D /home/x_ryabo/spectral_files/dat_hitran/hno3_hitran.par -R 1 20 -c 3000.0 -i 0.1 -l 12 1.0e1 -t 1.0e-2 -k -s sp_b20_HITRAN_a16_sw +S /home/x_ryabo/spectral_files/stellar_spectra/TRAPPIST-1_soc.txt -lk -o hno3_o -m hno3_m -L hno3_lbl.nc -np 16
prep_spec <<EOF
sp_b20_HITRAN_a16_sw
a
5
y
hno3_o
-1
EOF
touch temp/done_hno3_lbl
rsync -av sp_b20_HITRAN_a16_sw temp/sp_b20_HITRAN_a16_sw_12
rsync -av sp_b20_HITRAN_a16_sw_k temp/sp_b20_HITRAN_a16_sw_12_k
echo '###############   N2 lines – HITRAN    ###############'
Ccorr_k -F pt_grid_lbl.dat -D /home/x_ryabo/spectral_files/dat_hitran/n2_hitran.par -R 1 20 -c 3000.0 -i 0.1 -l 13 1.0e1 -t 1.0e-2 -k -s sp_b20_HITRAN_a16_sw +S /home/x_ryabo/spectral_files/stellar_spectra/TRAPPIST-1_soc.txt -lk -o n2_o -m n2_m -L n2_lbl.nc -np 16
prep_spec <<EOF
sp_b20_HITRAN_a16_sw
a
5
y
n2_o
-1
EOF
touch temp/done_n2_lbl
rsync -av sp_b20_HITRAN_a16_sw temp/sp_b20_HITRAN_a16_sw_13
rsync -av sp_b20_HITRAN_a16_sw_k temp/sp_b20_HITRAN_a16_sw_13_k
echo '###############   OCS lines – HITRAN    ###############'
Ccorr_k -F pt_grid_lbl.dat -D /home/x_ryabo/spectral_files/dat_hitran/ocs_hitran.par -R 1 20 -c 3000.0 -i 0.1 -l 25 1.0e1 -t 1.0e-2 -k -s sp_b20_HITRAN_a16_sw +S /home/x_ryabo/spectral_files/stellar_spectra/TRAPPIST-1_soc.txt -lk -o ocs_o -m ocs_m -L ocs_lbl.nc -np 16
prep_spec <<EOF
sp_b20_HITRAN_a16_sw
a
5
y
ocs_o
-1
EOF
touch temp/done_ocs_lbl
rsync -av sp_b20_HITRAN_a16_sw temp/sp_b20_HITRAN_a16_sw_14
rsync -av sp_b20_HITRAN_a16_sw_k temp/sp_b20_HITRAN_a16_sw_14_k
echo '###############   N2-N2 CIA   ###############'
Ccorr_k -F pt_grid_cia.dat -CIA /home/x_ryabo/spectral_files/dat_continua/N2-N2_2018_fixed.cia -R 1 20 -i 0.1 -ct 13 13 1.0e2 -t 1.0e-2 -s sp_b20_HITRAN_a16_sw +S /home/x_ryabo/spectral_files/stellar_spectra/TRAPPIST-1_soc.txt -lk -o n2n2_o -m n2n2_m -L n2n2_lbl.nc -np 16
prep_spec <<EOF
sp_b20_HITRAN_a16_sw
a
19
y
n2n2_o
-1
EOF
sed -i -e 's/\*\*\*\*\*/0/g' sp_b20_HITRAN_a16_sw
sed -i -e 's/\*\*\*//g' sp_b20_HITRAN_a16_sw
sed -i -e 's/            NaN/0.000000000E+00/g' sp_b20_HITRAN_a16_sw
touch temp/done_n2n2_cia
rsync -av sp_b20_HITRAN_a16_sw temp/sp_b20_HITRAN_a16_sw_15
rsync -av sp_b20_HITRAN_a16_sw_k temp/sp_b20_HITRAN_a16_sw_15_k
echo '###############   O2-O2 CIA   ###############'
Ccorr_k -F pt_grid_cia.dat -CIA /home/x_ryabo/spectral_files/dat_continua/O2-O2_2018b_fixed.cia -R 1 20 -i 0.1 -ct 7 7 1.0e2 -t 1.0e-2 -s sp_b20_HITRAN_a16_sw +S /home/x_ryabo/spectral_files/stellar_spectra/TRAPPIST-1_soc.txt -lk -o o2o2_o -m o2o2_m -L o2o2_lbl.nc -np 16
prep_spec <<EOF
sp_b20_HITRAN_a16_sw
a
19
y
o2o2_o
-1
EOF
sed -i -e 's/\*\*\*\*\*/0/g' sp_b20_HITRAN_a16_sw
sed -i -e 's/\*\*\*//g' sp_b20_HITRAN_a16_sw
sed -i -e 's/            NaN/0.000000000E+00/g' sp_b20_HITRAN_a16_sw
touch temp/done_o2o2_cia
rsync -av sp_b20_HITRAN_a16_sw temp/sp_b20_HITRAN_a16_sw_16
rsync -av sp_b20_HITRAN_a16_sw_k temp/sp_b20_HITRAN_a16_sw_16_k
echo '###############   H2-H2 CIA   ###############'
Ccorr_k -F pt_grid_cia.dat -CIA /home/x_ryabo/spectral_files/dat_continua/H2-H2_CIA_Borysow_combined.cia -R 1 20 -i 0.1 -ct 23 23 1.0e2 -t 1.0e-2 -s sp_b20_HITRAN_a16_sw +S /home/x_ryabo/spectral_files/stellar_spectra/TRAPPIST-1_soc.txt -lk -o h2h2_o -m h2h2_m -L h2h2_lbl.nc -np 16
prep_spec <<EOF
sp_b20_HITRAN_a16_sw
a
19
y
h2h2_o
-1
EOF
sed -i -e 's/\*\*\*\*\*/0/g' sp_b20_HITRAN_a16_sw
sed -i -e 's/\*\*\*//g' sp_b20_HITRAN_a16_sw
sed -i -e 's/            NaN/0.000000000E+00/g' sp_b20_HITRAN_a16_sw
touch temp/done_h2h2_cia
rsync -av sp_b20_HITRAN_a16_sw temp/sp_b20_HITRAN_a16_sw_17
rsync -av sp_b20_HITRAN_a16_sw_k temp/sp_b20_HITRAN_a16_sw_17_k
echo '###############   CO2-CO2 CIA   ###############'
Ccorr_k -F pt_grid_cia.dat -CIA /home/x_ryabo/spectral_files/dat_continua/CO2-CO2_2018_fixed.cia -R 1 20 -i 0.1 -ct 2 2 1.0e2 -t 1.0e-2 -s sp_b20_HITRAN_a16_sw +S /home/x_ryabo/spectral_files/stellar_spectra/TRAPPIST-1_soc.txt -lk -o co2co2_o -m co2co2_m -L co2co2_lbl.nc -np 16
prep_spec <<EOF
sp_b20_HITRAN_a16_sw
a
19
y
co2co2_o
-1
EOF
sed -i -e 's/\*\*\*\*\*/0/g' sp_b20_HITRAN_a16_sw
sed -i -e 's/\*\*\*//g' sp_b20_HITRAN_a16_sw
sed -i -e 's/            NaN/0.000000000E+00/g' sp_b20_HITRAN_a16_sw
touch temp/done_co2co2_cia
rsync -av sp_b20_HITRAN_a16_sw temp/sp_b20_HITRAN_a16_sw_18
rsync -av sp_b20_HITRAN_a16_sw_k temp/sp_b20_HITRAN_a16_sw_18_k
echo '###############   CH4-CH4 CIA   ###############'
Ccorr_k -F pt_grid_cia.dat -CIA /home/x_ryabo/spectral_files/dat_continua/CH4-CH4_2011.cia -R 1 20 -i 0.1 -ct 6 6 1.0e2 -t 1.0e-2 -s sp_b20_HITRAN_a16_sw +S /home/x_ryabo/spectral_files/stellar_spectra/TRAPPIST-1_soc.txt -lk -o ch4ch4_o -m ch4ch4_m -L ch4ch4_lbl.nc -np 16
prep_spec <<EOF
sp_b20_HITRAN_a16_sw
a
19
y
ch4ch4_o
-1
EOF
sed -i -e 's/\*\*\*\*\*/0/g' sp_b20_HITRAN_a16_sw
sed -i -e 's/\*\*\*//g' sp_b20_HITRAN_a16_sw
sed -i -e 's/            NaN/0.000000000E+00/g' sp_b20_HITRAN_a16_sw
touch temp/done_ch4ch4_cia
rsync -av sp_b20_HITRAN_a16_sw temp/sp_b20_HITRAN_a16_sw_19
rsync -av sp_b20_HITRAN_a16_sw_k temp/sp_b20_HITRAN_a16_sw_19_k
echo '###############   H2O-H2O MT_CKD   ###############'
Ccorr_k -F pt_grid_cia.dat -R 1 20 -c 2500.0 -i 0.1 -ct 1 1 10.0 -t 1.0e-3 -e /home/x_ryabo/spectral_files/dat_continua/mt_ckd_v3.0_s296 /home/x_ryabo/spectral_files/dat_continua/mt_ckd_v3.0_s260 -k -s sp_b20_HITRAN_a16_sw +S /home/x_ryabo/spectral_files/stellar_spectra/TRAPPIST-1_soc.txt -lk -o h2o-h2o_l20_1-20c -m h2o-h2o_l20_1-20cm -L h2o-h2o_lbl_lw.nc -lw h2o_l20_1-20map.nc -np 16
prep_spec <<EOF
sp_b20_HITRAN_a16_sw
a
19
y
h2o-h2o_l20_1-20c
-1
EOF
sed -i -e 's/\*\*\*\*\*/0/g' sp_b20_HITRAN_a16_sw
sed -i -e 's/\*\*\*//g' sp_b20_HITRAN_a16_sw
sed -i -e 's/            NaN/0.000000000E+00/g' sp_b20_HITRAN_a16_sw
touch temp/done_h2o_mtckd
rsync -av sp_b20_HITRAN_a16_sw temp/sp_b20_HITRAN_a16_sw_20
rsync -av sp_b20_HITRAN_a16_sw_k temp/sp_b20_HITRAN_a16_sw_20_k
echo '###############   CO2-CH4 CIA   ###############'
Ccorr_k -F pt_grid_cia.dat -CIA /home/x_ryabo/spectral_files/dat_continua/CO2-CH4_2018.cia -R 1 20 -i 0.1 -ct 2 6 1.0e2 -t 1.0e-2 -s sp_b20_HITRAN_a16_sw +S /home/x_ryabo/spectral_files/stellar_spectra/TRAPPIST-1_soc.txt -lk -o co2ch4_o -m co2ch4_m -L co2ch4_lbl.nc -np 16
prep_spec <<EOF
sp_b20_HITRAN_a16_sw
a
19
y
co2ch4_o
-1
EOF
sed -i -e 's/\*\*\*\*\*/0/g' sp_b20_HITRAN_a16_sw
sed -i -e 's/\*\*\*//g' sp_b20_HITRAN_a16_sw
sed -i -e 's/            NaN/0.000000000E+00/g' sp_b20_HITRAN_a16_sw
sed -i -e 's/          51/           0/g' sp_b20_HITRAN_a16_sw
touch temp/done_co2ch4_cia
rsync -av sp_b20_HITRAN_a16_sw temp/sp_b20_HITRAN_a16_sw_21
rsync -av sp_b20_HITRAN_a16_sw_k temp/sp_b20_HITRAN_a16_sw_21_k
echo '###############   CO2-H2 CIA   ###############'
Ccorr_k -F pt_grid_cia.dat -CIA /home/x_ryabo/spectral_files/dat_continua/CO2-H2_2018.cia -R 1 20 -i 0.1 -ct 2 23 1.0e2 -t 1.0e-2 -s sp_b20_HITRAN_a16_sw +S /home/x_ryabo/spectral_files/stellar_spectra/TRAPPIST-1_soc.txt -lk -o co2h2_o -m co2h2_m -L co2h2_lbl.nc -np 16
prep_spec <<EOF
sp_b20_HITRAN_a16_sw
a
19
y
co2h2_o
-1
EOF
sed -i -e 's/\*\*\*\*\*/0/g' sp_b20_HITRAN_a16_sw
sed -i -e 's/\*\*\*//g' sp_b20_HITRAN_a16_sw
sed -i -e 's/            NaN/0.000000000E+00/g' sp_b20_HITRAN_a16_sw
sed -i -e 's/          51/           0/g' sp_b20_HITRAN_a16_sw
touch temp/done_co2h2_cia
rsync -av sp_b20_HITRAN_a16_sw temp/sp_b20_HITRAN_a16_sw_22
rsync -av sp_b20_HITRAN_a16_sw_k temp/sp_b20_HITRAN_a16_sw_22_k
echo '###############   CO2-He CIA   ###############'
Ccorr_k -F pt_grid_cia.dat -CIA /home/x_ryabo/spectral_files/dat_continua/CO2-He_2018.cia -R 1 20 -i 0.1 -ct 2 24 1.0e2 -t 1.0e-2 -s sp_b20_HITRAN_a16_sw +S /home/x_ryabo/spectral_files/stellar_spectra/TRAPPIST-1_soc.txt -lk -o co2he_o -m co2he_m -L co2he_lbl.nc -np 16
prep_spec <<EOF
sp_b20_HITRAN_a16_sw
a
19
y
co2he_o
-1
EOF
sed -i -e 's/\*\*\*\*\*/0/g' sp_b20_HITRAN_a16_sw
sed -i -e 's/\*\*\*//g' sp_b20_HITRAN_a16_sw
sed -i -e 's/            NaN/0.000000000E+00/g' sp_b20_HITRAN_a16_sw
sed -i -e 's/          51/           0/g' sp_b20_HITRAN_a16_sw
touch temp/done_co2he_cia
rsync -av sp_b20_HITRAN_a16_sw temp/sp_b20_HITRAN_a16_sw_23
rsync -av sp_b20_HITRAN_a16_sw_k temp/sp_b20_HITRAN_a16_sw_23_k
echo '###############   CH4-He CIA   ###############'
Ccorr_k -F pt_grid_cia.dat -CIA /home/x_ryabo/spectral_files/dat_continua/CH4-He_2018.cia -R 1 20 -i 0.1 -ct 6 24 1.0e2 -t 1.0e-2 -s sp_b20_HITRAN_a16_sw +S /home/x_ryabo/spectral_files/stellar_spectra/TRAPPIST-1_soc.txt -lk -o ch4he_o -m ch4he_m -L ch4he_lbl.nc -np 16
prep_spec <<EOF
sp_b20_HITRAN_a16_sw
a
19
y
ch4he_o
-1
EOF
sed -i -e 's/\*\*\*\*\*/0/g' sp_b20_HITRAN_a16_sw
sed -i -e 's/\*\*\*//g' sp_b20_HITRAN_a16_sw
sed -i -e 's/            NaN/0.000000000E+00/g' sp_b20_HITRAN_a16_sw
sed -i -e 's/          51/           0/g' sp_b20_HITRAN_a16_sw
touch temp/done_ch4he_cia
rsync -av sp_b20_HITRAN_a16_sw temp/sp_b20_HITRAN_a16_sw_24
rsync -av sp_b20_HITRAN_a16_sw_k temp/sp_b20_HITRAN_a16_sw_24_k
echo '###############   O2-CO2 CIA   ###############'
Ccorr_k -F pt_grid_cia.dat -CIA /home/x_ryabo/spectral_files/dat_continua/O2-CO2_2011.cia -R 1 20 -i 0.1 -ct 7 2 1.0e2 -t 1.0e-2 -s sp_b20_HITRAN_a16_sw +S /home/x_ryabo/spectral_files/stellar_spectra/TRAPPIST-1_soc.txt -lk -o o2co2_o -m o2co2_m -L o2co2_lbl.nc -np 16
prep_spec <<EOF
sp_b20_HITRAN_a16_sw
a
19
y
o2co2_o
-1
EOF
sed -i -e 's/\*\*\*\*\*/0/g' sp_b20_HITRAN_a16_sw
sed -i -e 's/\*\*\*//g' sp_b20_HITRAN_a16_sw
sed -i -e 's/            NaN/0.000000000E+00/g' sp_b20_HITRAN_a16_sw
touch temp/done_o2co2_cia
rsync -av sp_b20_HITRAN_a16_sw temp/sp_b20_HITRAN_a16_sw_25
rsync -av sp_b20_HITRAN_a16_sw_k temp/sp_b20_HITRAN_a16_sw_25_k
echo '###############   O2-N2 CIA   ###############'
Ccorr_k -F pt_grid_cia.dat -CIA /home/x_ryabo/spectral_files/dat_continua/O2-N2_2018_fixed.cia -R 1 20 -i 0.1 -ct 7 13 1.0e2 -t 1.0e-2 -s sp_b20_HITRAN_a16_sw +S /home/x_ryabo/spectral_files/stellar_spectra/TRAPPIST-1_soc.txt -lk -o o2n2_o -m o2n2_m -L o2n2_lbl.nc -np 16
prep_spec <<EOF
sp_b20_HITRAN_a16_sw
a
19
y
o2n2_o
-1
EOF
sed -i -e 's/\*\*\*\*\*/0/g' sp_b20_HITRAN_a16_sw
sed -i -e 's/\*\*\*//g' sp_b20_HITRAN_a16_sw
sed -i -e 's/            NaN/0.000000000E+00/g' sp_b20_HITRAN_a16_sw
touch temp/done_o2n2_cia
rsync -av sp_b20_HITRAN_a16_sw temp/sp_b20_HITRAN_a16_sw_26
rsync -av sp_b20_HITRAN_a16_sw_k temp/sp_b20_HITRAN_a16_sw_26_k
echo '###############   N2-H2O CIA   ###############'
Ccorr_k -F pt_grid_cia.dat -CIA /home/x_ryabo/spectral_files/dat_continua/N2-H2O_2018.cia -R 1 20 -i 0.1 -ct 13 1 1.0e2 -t 1.0e-2 -s sp_b20_HITRAN_a16_sw +S /home/x_ryabo/spectral_files/stellar_spectra/TRAPPIST-1_soc.txt -lk -o n2h2o_o -m n2h2o_m -L n2h2o_lbl.nc -np 16
prep_spec <<EOF
sp_b20_HITRAN_a16_sw
a
19
y
n2h2o_o
-1
EOF
sed -i -e 's/\*\*\*\*\*/0/g' sp_b20_HITRAN_a16_sw
sed -i -e 's/\*\*\*//g' sp_b20_HITRAN_a16_sw
sed -i -e 's/            NaN/0.000000000E+00/g' sp_b20_HITRAN_a16_sw
touch temp/done_n2h2o_cia
rsync -av sp_b20_HITRAN_a16_sw temp/sp_b20_HITRAN_a16_sw_27
rsync -av sp_b20_HITRAN_a16_sw_k temp/sp_b20_HITRAN_a16_sw_27_k
echo '###############   N2-CH4 CIA   ###############'
Ccorr_k -F pt_grid_cia.dat -CIA /home/x_ryabo/spectral_files/dat_continua/N2-CH4_2011.cia -R 1 20 -i 0.1 -ct 13 6 1.0e2 -t 1.0e-2 -s sp_b20_HITRAN_a16_sw +S /home/x_ryabo/spectral_files/stellar_spectra/TRAPPIST-1_soc.txt -lk -o n2ch4_o -m n2ch4_m -L n2ch4_lbl.nc -np 16
prep_spec <<EOF
sp_b20_HITRAN_a16_sw
a
19
y
n2ch4_o
-1
EOF
sed -i -e 's/\*\*\*\*\*/0/g' sp_b20_HITRAN_a16_sw
sed -i -e 's/\*\*\*//g' sp_b20_HITRAN_a16_sw
sed -i -e 's/            NaN/0.000000000E+00/g' sp_b20_HITRAN_a16_sw
touch temp/done_n2ch4_cia
rsync -av sp_b20_HITRAN_a16_sw temp/sp_b20_HITRAN_a16_sw_28
rsync -av sp_b20_HITRAN_a16_sw_k temp/sp_b20_HITRAN_a16_sw_28_k
echo '###############   N2-H2 CIA   ###############'
Ccorr_k -F pt_grid_cia.dat -CIA /home/x_ryabo/spectral_files/dat_continua/N2-H2_2011.cia -R 1 20 -i 0.1 -ct 13 23 1.0e2 -t 1.0e-2 -s sp_b20_HITRAN_a16_sw +S /home/x_ryabo/spectral_files/stellar_spectra/TRAPPIST-1_soc.txt -lk -o n2h2_o -m n2h2_m -L n2h2_lbl.nc -np 16
prep_spec <<EOF
sp_b20_HITRAN_a16_sw
a
19
y
n2h2_o
-1
EOF
sed -i -e 's/\*\*\*\*\*/0/g' sp_b20_HITRAN_a16_sw
sed -i -e 's/\*\*\*//g' sp_b20_HITRAN_a16_sw
sed -i -e 's/            NaN/0.000000000E+00/g' sp_b20_HITRAN_a16_sw
touch temp/done_n2h2_cia
rsync -av sp_b20_HITRAN_a16_sw temp/sp_b20_HITRAN_a16_sw_29
rsync -av sp_b20_HITRAN_a16_sw_k temp/sp_b20_HITRAN_a16_sw_29_k
echo '###############   N2-He CIA   ###############'
Ccorr_k -F pt_grid_cia.dat -CIA /home/x_ryabo/spectral_files/dat_continua/N2-He_2018.cia -R 1 20 -i 0.1 -ct 13 24 1.0e2 -t 1.0e-2 -s sp_b20_HITRAN_a16_sw +S /home/x_ryabo/spectral_files/stellar_spectra/TRAPPIST-1_soc.txt -lk -o n2he_o -m n2he_m -L n2he_lbl.nc -np 16
prep_spec <<EOF
sp_b20_HITRAN_a16_sw
a
19
y
n2he_o
-1
EOF
sed -i -e 's/\*\*\*\*\*/0/g' sp_b20_HITRAN_a16_sw
sed -i -e 's/\*\*\*//g' sp_b20_HITRAN_a16_sw
sed -i -e 's/            NaN/0.000000000E+00/g' sp_b20_HITRAN_a16_sw
touch temp/done_n2he_cia
rsync -av sp_b20_HITRAN_a16_sw temp/sp_b20_HITRAN_a16_sw_30
rsync -av sp_b20_HITRAN_a16_sw_k temp/sp_b20_HITRAN_a16_sw_30_k
echo '###############   H2-CH4 CIA   ###############'
Ccorr_k -F pt_grid_cia.dat -CIA /home/x_ryabo/spectral_files/dat_continua/H2-CH4_eq_2011.cia -R 1 20 -i 0.1 -ct 23 6 1.0e2 -t 1.0e-2 -s sp_b20_HITRAN_a16_sw +S /home/x_ryabo/spectral_files/stellar_spectra/TRAPPIST-1_soc.txt -lk -o h2ch4_o -m h2ch4_m -L h2ch4_lbl.nc -np 16
prep_spec <<EOF
sp_b20_HITRAN_a16_sw
a
19
y
h2ch4_o
-1
EOF
sed -i -e 's/\*\*\*\*\*/0/g' sp_b20_HITRAN_a16_sw
sed -i -e 's/\*\*\*//g' sp_b20_HITRAN_a16_sw
sed -i -e 's/            NaN/0.000000000E+00/g' sp_b20_HITRAN_a16_sw
touch temp/done_h2ch4_cia
rsync -av sp_b20_HITRAN_a16_sw temp/sp_b20_HITRAN_a16_sw_31
rsync -av sp_b20_HITRAN_a16_sw_k temp/sp_b20_HITRAN_a16_sw_31_k
echo '###############   H2-He CIA   ###############'
Ccorr_k -F pt_grid_cia.dat -CIA /home/x_ryabo/spectral_files/dat_continua/H2-He_2011.cia -R 1 20 -i 0.1 -ct 23 24 1.0e2 -t 1.0e-2 -s sp_b20_HITRAN_a16_sw +S /home/x_ryabo/spectral_files/stellar_spectra/TRAPPIST-1_soc.txt -lk -o h2he_o -m h2he_m -L h2he_lbl.nc -np 16
prep_spec <<EOF
sp_b20_HITRAN_a16_sw
a
19
y
h2he_o
-1
EOF
sed -i -e 's/\*\*\*\*\*/0/g' sp_b20_HITRAN_a16_sw
sed -i -e 's/\*\*\*//g' sp_b20_HITRAN_a16_sw
sed -i -e 's/            NaN/0.000000000E+00/g' sp_b20_HITRAN_a16_sw
touch temp/done_h2he_cia
rsync -av sp_b20_HITRAN_a16_sw temp/sp_b20_HITRAN_a16_sw_32
rsync -av sp_b20_HITRAN_a16_sw_k temp/sp_b20_HITRAN_a16_sw_32_k
