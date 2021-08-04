This repository contains the collection of MATLAB scripts needed to run the ice strain inversion described in Osman et al. (2021), which generates a range of Nuussuaq Ice Cap, west Greenland ice cap mass balance histories ranked in order of their “complexity” (see Supplementary Information of Osman et al., 2021). <br>

The main “driving function” to run is “invert_NU_ice_cap.m”.  Running this script requires the following dependencies: <br>
1. refine_time_scale_1D_int.m <br>
2. est_timescale_u.m <br>
3. chi2_of_core_int.m <br>
4. index_struct.m <br>
5. TableS1.xlsx <br>
6. A folder named “output” in the MATLAB home directory <br>

Running this script generates a .mat output file (~100 mB) named “Nuus_AgeDepth_post169_<date>.mat”, where <date> is the date the output file was generated.  To visualize the results contained in “Nuus_AgeDepth_post169_<date>.mat”, run the script “visualize_NU_inversion.m”.  Note this script requires the folder “cbrewer” (from https://www.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab) in the MATLAB home directory. <br>

Last successfully tested in MATLAB 2020b on July 24, 2021.  Takes about 30 minutes to run on my 2018 MacBook Pro (2.7 GHz Quad-Core Intel Core i7). <br>

*Reference*: <br>
Osman, M.B., Smith, B.E., Trusel, L.D. *et al*. Abrupt Common Era hydroclimate shifts drive west Greenland ice cap change, *accepted*, *Nature Geoscience*, 2021.
