% This script inverts an ice strain model using age-depth constraints to
% generate a series of accumulation histories for the Nuussuaq Peninsula Ice Cap, 
% west Greenland.  See Osman et al. (2021, Nature Geoscience) for details. 
% The output is saved to the folder "output" as a .mat file.  Please see
% the README.txt file for details on how to interpret the output, or run
% "visualize_NU_inversion.m" to visualize the output. 
% 
% Written by M.B. Osman (MIT/WHOI and University of Arizona; mattosman@arizona.edu) 
% and B.E. Smith (Applied Physics Lab and University of Washington; besmith@uw.edu)
% Spring 2018
% 
% To run, requires the following dependencies:
%   1. refine_time_scale_1D_int.m
%   2. est_timescale_u.m
%   3. chi2_of_core_int.m
%   4. index_struct.m
%   5. TableS1.xlsx
%   6. A folder in the home directory called "output"
clear

%% Changeable variables:

min_year = 169; % AD (or -300)
H0 = 140;  % Reference ice cap thickness (m ice)
dz = [-25:0.5:20]; % difference from baseline thickness guess
b_vals = [0.2:0.005:0.45];  % Range of accumulation rates (m ice eq) to examine

%%  Upload tie points

% extract tie point data
data = xlsread('TableS1.xlsx');
year        = data(:,1); 
water_depth = data(:,2); 
sigma_age   = data(:,3);

% edit out tie points before "min_year"
index_years = year > min_year;
    output_file = ['Nuus_AgeDepth_post_',num2str(min_year),'_',date];
    year          = year(index_years); % post 264 AD
    water_depth   = water_depth(index_years);
    sigma_age     = sigma_age(index_years);
    
%% Pre condition variables

% convert water depth to ice depth
ice_depth = water_depth / .917;
% delta age
age = 2015.5 - year; % difference from top of core age
% Assign relevant variables to structure "core"
core.depth = ice_depth;
core.age   = age;
core.sigma_age = sigma_age;

core.sigma_depth = core.depth*NaN; % preallocate sigma depth field - this is
    core.sigma_depth(1) = diff(core.depth(1:2)) ./ diff(core.age(1:2)) * core.sigma_age(1);
    core.sigma_depth(end) = diff(core.depth(end-1:end)) ./ ...
        diff(core.age(end-1:end))*core.sigma_age(end);
    for k = 2:length(core.depth)-1
        core.sigma_depth(k) = mean(diff(core.depth(k-1:k+1)) ./ ...
            diff(core.age(k-1:k+1)))*core.sigma_age(k);
    end
core = index_struct(core, core.sigma_age~=0); % remove surface values

%--------------------------------- optimize model fit over a range of ice thickness and accumulation levels 

clear all_fit_hist
all_sum_M2 = nan(length(dz), length(b_vals));
h = waitbar(0,['Testing all model combinations. Please wait...']);
for i = 1:length(dz)
    waitbar(i/ length(dz));
    clearvars fit_hist
    
    for j = 1:length(b_vals)
        this_H = H0 - dz(i);
        this_depth = 0:0.1:this_H-0.1;
        
        TS_constant(i,j).depth = this_depth;
        % estimate a reference-age using characteristic (steady state) accumulatiomn / ice sheet thickness
        TS_constant(i,j).age   = est_timescale_u('age_for_ref_model', ... 
            this_H - this_depth, this_H, b_vals(j)); % this_H - this_depth = vector w/ height above bed (z, H, bdot_ref)
        TS_constant(i,j).accum = b_vals(j);
        [fit_hist(j).ConstraintMisfit, fit_hist(j).R2, fit_hist(j).m] = ...
            refine_time_scale_1D_int(TS_constant(i, j), b_vals(j), core);
        fit_hist(j).Bdot0 = b_vals(j);
        if ~isempty(fit_hist(j).m)
            all_sum_M2(i, j) = fit_hist(j).ConstraintMisfit;
        end
    end
    
    [MinMisfit(i), ii] = min([fit_hist.ConstraintMisfit]);
    Accum_best(i)      = b_vals(ii);
    all_fit_hist{i}    = fit_hist;
    disp([num2str(round(i/length(dz)*100,0)),'% complete'])
    
[rr,cc] = find(all_sum_M2 == min(all_sum_M2(:)));
best_model.H        = H0-dz(rr); 
best_model.bdot_ref = b_vals(cc);
best_model.m        = all_fit_hist{rr}(cc).m;
best_model.R2       = all_sum_M2(rr, cc);
H_vals              = H0 - dz;
end
close(h);

clearvars -except ...
    best_model...
    all_fit_hist...
    H_vals...
    b_vals...
    all_sum_M2...
    core ...
    output_file

%% Save output
cd output
    save([output_file,'.mat'],'best_model','all_fit_hist','H_vals','b_vals','all_sum_M2','core');
cd ../
