#!/bin/bash
prep_spec <<EOF
test_spec
2
1
1
1
1
1
0
c
10000
20000
20000
30000
1
1
1
1
n
-1
EOF
Ccorr_k -F pt_grid_lbl.dat -D /home/x_ryabo/spectral_files/dat_hitran/h2o_data.par -R 1 2 -c 3000.0 -i 1.0 -l 1 1.0e1 -b 1.5e-3 -s test_spec +p -lk -k -o h2o_l2_l -m h2o_l2_lm -L h2o_lbl_lwf_pt48.nc -sm h2o_l2_1-2map.nc -np 16
prep_spec <<EOF
test_spec
a
5
h2o_l2_l
-1
EOF
Ccorr_k -F pt_grid_cia.dat -R 1 2 -c 2500.0 -i 1.0 -ct 1 1 10.0 -t 1.0e-3 -e mt_ckd_v3.0_s296 mt_ckd_v3.0_s260 -k -s test_spec +p -lk -o h2o-h2o_l2_1-2c -m h2o-h2o_l2_1-2cm -L h2o-h2o_lbl_lw.nc -lw h2o_l2_1-2map.nc -np 16
prep_spec <<EOF
test_spec
a
19
h2o-h2o_l2_1-2c
-1
EOF
#sed -i -e 's/\*\*\*\*\*/0/g' test_spec
#sed -i -e 's/\*\*\*//g' test_spec
#sed -i -e 's/            NaN/0.000000000E+00/g' test_spec
prep_spec <<EOF
test_spec
a
6
n
T
100. 4000.
2
2
n
lean_12
y
-1
EOF
