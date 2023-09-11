-----------------------------------
sp_gen_script.py
-----------------------------------

Line 58: download_trans is the first step. trans_xsec is the most time-consuming step by far. These two steps are not currently done in the bash exec file.

Line 97: new band spacing to span the entire exomol spectral range.

Line 538: ExoMol H2O processing starts. 

Line 568: Call to the external function to use ExoCross. 

-----------------------------------
exocross_loop.py
-----------------------------------

The function xcross_compute_xsec takes T_grid, P_grid, bands, and i_int from sp_gen_script.py.

The input/output paths are hardcoded for now.

It builds input files for ExoCross and produces a .xsec file for each PT point.

The call to xcross_compute_xsec for i_int = 0.1 took about 4 days with my student's reduced PT grid. With the full grids, I completed two PT points in 4 days. It might take about 42 days to produce all the .xsec files. These will need to be stored somewhere but they are not heavy. 14Mb each.