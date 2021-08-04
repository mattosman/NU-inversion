%% Produces a few diagnostic plots of the NU age-scale inversion
% to run, first requires one to run "invert_NU_ice_cap.m"
clear

cd output
	[file,~] = uigetfile(['*','.mat'],'Select output file.');
	if isequal(file,0)
        disp('User selected cancel; quitting program.');
        return;
    else
        load(file)
    end
cd ../

cd cbrewer
	colors = cbrewer('div','RdBu',301);
cd ../

%%  Collect all R2 vals 
all_sum_R2 = nan(size(all_sum_M2));
for i = 1:size(all_sum_M2,1)
    for j = 1:size(all_sum_M2,2)
    all_sum_R2(i,j) = all_fit_hist{i}(j).R2;
    end
end

% count up the length of each accumulation period:
num_intervals = length(unique(best_model.m.scale)) + 1;
curr_scale = best_model.m.scale(1); 
k = 1; a = 1;
for i = 1:length(best_model.m.scale)
    if best_model.m.scale(i) == curr_scale
        k = k+1;
    else
        num_years(a,1) = k;
        k = 1; a = a+1;
        curr_scale = best_model.m.scale(i);
    end
end
num_years(a,1) = k;
k = 1; a = a+1;
curr_scale = best_model.m.scale(i);

mean(num_years)
max(num_years)
min(num_years)
prctile(num_years,25)
prctile(num_years,75)
% figure; histogram(num_years,9);

%% Fig

% define complexity (~uncertainty~) range (e.g., 10% more complex, etc.)
unc = 0.10;

fig1 = figure; clf;  % figure: age_depth_and_errors
    set(fig1,'PaperPositionMode','auto');         
    set(fig1,'PaperOrientation','landscape');
    set(fig1,'Position',[50 50 400 550]); 
        ax = gca; ax.Visible = 'off'; % turn off annoying auto-axis 
        
ax1 = axes('position', [0.15 0.55 0.70 0.30],'color','none'); hold on; box on;
    contourf(b_vals, H_vals, log(all_sum_M2), size(colors,1),'Linecolor','None'); shading flat; 
    cb = colorbar; set(cb,'Linewidth',1.5)
    colormap(ax1, colors)
        xlabel('Accumulation_{Ref.} (m_{ice} yr^{-1})')
        ylabel('Ice Cap Thickness (m)')
        cb.Label.String = 'log(Complexity)';
    set(ax1,'Linewidth',1.5,'Fontsize',12)
    % plot maximum point
    [~, index] = nanmin(all_sum_M2(:));
    [I_row, I_col] = ind2sub(size(all_sum_M2),index); % extract x and y indices
    scatter(b_vals(I_col),H_vals(I_row),200,'p','MarkerEdgeColor','None', ...
    	'MarkerFaceColor',[0 0 0], 'LineWidth',1)
    % 10% more complicated threshold
    contourf(b_vals,H_vals, log(all_sum_M2), log((1 + unc).*nanmin(all_sum_M2(:))).*[1 1] ,...
        'linecolor', 'k', 'Linewidth',1.5, 'LineStyle','--', 'Fill', 'off')
        
ax2 = axes('position', [0.15 0.10 0.70 0.30],'color','none'); hold on; box on;
    [r,c] = find(all_sum_M2 < (1 + unc)*min(all_sum_M2(:)));
    b_range = [min(b_vals(c)) max(b_vals(c))]; % range(b_vals(cc)); 
    [StackedFitHist.age, StackedFitHist.depth, StackedFitHist.SMB,StackedFitHist.sigma ] = ...
        deal(NaN(length(all_fit_hist{1}(1).m.depth), length(r)));
    cd cbrewer
        colors2 = cbrewer('seq','Reds',length(unique(c)));
        color_index = c - (min(c) - 1);
    cd ../
    for i = 1:length(r)
        M = all_fit_hist{r(i)}(c(i)).m;
        plot(2015.5 - M.age, M.scale*b_vals(c(i)),'color', colors2(color_index(i),:));
    end 
    caxis(ax2,b_range);  cb = colorbar('eastoutside'); set(cb,'Linewidth',1.5) 
    colormap(ax2,colors2)
    % set(cb,'position', get(cb,'position').*[1 1 .5 .5] + [.3 0 0 0]);
    xlabel('Year (CE)')
    ylabel('Accumulation_{Ref.} (m_{ice} yr^{-1})')
    cb.Label.String = 'Accumulation_{Ref.} (m_{ice} yr^{-1})';
    set(ax2,'Linewidth',1.5,'Fontsize',12)
    xlim([min(2015.5 - M.age), 2016])
    % plot the least complicated accumulation model:
    plot(2015.5 - best_model.m.age, best_model.m.scale*best_model.bdot_ref,... % best accumulation
        'Color',[0.2 0.2 0.2],'LineStyle','-','LineWidth',2);
    line([400 600],[0.27 0.27],'Color',[0.2 0.2 0.2],'LineStyle','-','LineWidth',1.5)
    text(650,0.27,'Optimum Model','Fontsize',12,'Color',[0.2 0.2 0.2])
    
    % save confidence range
    conf_range.year = 2015.5 - best_model.m.age;
    conf_range.lower_conf = b_vals(c(1)) * all_fit_hist{r(1)}(c(1)).m.scale;
    conf_range.middle_conf = best_model.m.scale*best_model.bdot_ref;
    conf_range.upper_conf = b_vals(c(end)) * all_fit_hist{r(end)}(c(end)).m.scale;
        
%% Fig 

fig2 = figure; clf;  % figure: age_depth_and_errors
    set(fig2,'PaperPositionMode','auto');         
    set(fig2,'PaperOrientation','landscape');
    set(fig2,'Position',[50 50 400 550]); 
        ax = gca; ax.Visible = 'off'; % turn off annoying auto-axis 

ax1 = axes('position', [0.20 0.55 0.70 0.30],'color','none'); hold on; box on;
    plot(2015.5 - best_model.m.age, best_model.m.depth,...
        'Color',[0.8 0.3 0.1],'LineStyle','-','LineWidth',2);
    set(gca','YDir','Reverse','XAxisLocation','Top')
    % errorbar(2015.5-core.age, core.depth, core.sigma_depth,'k.','Linewidth',3)
    scatter(2015.5-core.age, core.depth, 25,'MarkerEdgeColor',[0.3 0.3 0.3], ...
    	'MarkerFaceColor',[0.2 0.3 0.8], 'LineWidth',0.2)
    xlabel('Year (CE)')
    ylabel('Depth (m_{ice})')
    set(ax1,'Linewidth',1.5,'Fontsize',12)
    xlim([min(2015.5 - M.age), 2016]); ylim([0 130])
        
ax2 = axes('position', [0.20 0.20 0.70 0.30],'color','none'); hold on; box on;
    line([min(2015.5 - M.age), 2016],[0 0],'Color',[0.7 0.2 0.2],'LineStyle','-','LineWidth',1.0)
    errorbar(2015.5-core.age, core.age-interp1(best_model.m.depth, best_model.m.age, core.depth,'spline'), core.sigma_age,'k.');
    scatter(2015.5-core.age, core.age-interp1(best_model.m.depth, best_model.m.age, core.depth,'spline'), 40, ...
        'MarkerEdgeColor',[0.3 0.3 0.3], 'MarkerFaceColor',[0.2 0.3 0.8], 'LineWidth',0.2)
    xlabel('Year (CE)')
    ylabel({'Age Offset from',' Optimum Model (yr)'})
    set(ax2,'Linewidth',1.5,'Fontsize',12)
    xlim([min(2015.5 - M.age), 2016])
    ylim([-25 25])
        
%% Fig

fig3 = figure; clf;  % figure: age_depth_and_errors
    set(fig3,'PaperPositionMode','auto');         
    set(fig3,'PaperOrientation','landscape');
    set(fig3,'Position',[50 50 400 550]); 
        ax = gca; ax.Visible = 'off'; % turn off annoying auto-axis 
        
    thickness_range = 10;

    % identify upper and lower thicknesses using reference accumulation rate
    upper_thickness.H = best_model.H + thickness_range;
    lower_thickness.H = best_model.H - thickness_range;
    % reference age-scale for constant accumulation / ice thickness
        upper_thickness.ref_depth =  0:0.1:upper_thickness.H-0.1;
        upper_thickness.ref_age   = 2015.5 - est_timescale_u('age_for_ref_model', ...
            upper_thickness.H - upper_thickness.ref_depth,  upper_thickness.H, best_model.bdot_ref); 
        lower_thickness.ref_depth =  0:0.1:lower_thickness.H-0.1;
        lower_thickness.ref_age   = 2015.5 - est_timescale_u('age_for_ref_model', ...
            lower_thickness.H - lower_thickness.ref_depth,  lower_thickness.H, best_model.bdot_ref); 
    % inverted accumulation scale
    upper_thickness.age = all_fit_hist{1,find(H_vals == upper_thickness.H)}(I_col).m.age;
    upper_thickness.accum = all_fit_hist{1,find(H_vals == upper_thickness.H)}(I_col).m.scale .* best_model.bdot_ref;
    lower_thickness.age = all_fit_hist{1,find(H_vals == lower_thickness.H)}(I_col).m.age;
    lower_thickness.accum = all_fit_hist{1,find(H_vals == lower_thickness.H)}(I_col).m.scale .* best_model.bdot_ref;
    % calculate the reference best model
        best_model.ref_depth =  0:0.1:best_model.H-0.1;
        best_model.ref_age   = 2015.5 - est_timescale_u('age_for_ref_model', ...
            best_model.H - best_model.ref_depth,  best_model.H, best_model.bdot_ref); 
    
ax1 = axes('position', [0.15 0.50 0.70 0.35],'color','none'); hold on; box on;
    % plot the best reference model 
    p1 = plot(best_model.ref_age, best_model.ref_depth,...
        'Color',colors2(floor(size(colors2,1)/2),:),'LineStyle','-','LineWidth',2);
    % now plot 2 upper-lower ice thickness estimates
    p2 = plot(upper_thickness.ref_age, upper_thickness.ref_depth,...
        'Color',colors2(end-2,:),'LineStyle','-','LineWidth',2);
    p3 = plot(lower_thickness.ref_age, lower_thickness.ref_depth,...
        'Color',colors2(2,:),'LineStyle','-','LineWidth',2);    
    set(gca','YDir','Reverse','XAxisLocation','Top')
    % verrorbar(2015.5-core.age, core.depth, core.sigma_depth,'k.','Linewidth',3)
    s = scatter(2015.5-core.age, core.depth, 25,'MarkerEdgeColor',[0.3 0.3 0.3], ...
    	'MarkerFaceColor',[0.2 0.3 0.8], 'LineWidth',0.2);
    xlabel('Year (CE)')
    ylabel('Depth (m_{ice})')
    set(ax1,'Linewidth',1.5,'Fontsize',12)
    xlim([min(2015.5 - M.age), 2016]); ylim([0 150])
    l = legend([p2 p1 p3 s],['Ice cap thickness: ',num2str(upper_thickness.H)], ...
        ['Ice cap thickness: ',num2str(best_model.H)], ...
        ['Ice cap thickness: ',num2str(lower_thickness.H)], ...
        ['Layer Picks'],'Location','Best'); set(l,'Box','off')
    
ax2 = axes('position', [0.15 0.15 0.70 0.30],'color','none'); hold on; box on;
    % plot the least complicated accumulation model:
    plot(2015.5 - best_model.m.age, best_model.m.scale*best_model.bdot_ref,... % best accumulation
        'Color',colors2(floor(size(colors2,1)/2),:),'LineStyle','-','LineWidth',2);
    plot(2015.5 - upper_thickness.age, upper_thickness.accum,... % best accumulation
        'Color',colors2(end-2,:),'LineStyle','-','LineWidth',2);
    plot(2015.5 - lower_thickness.age, lower_thickness.accum,... % best accumulation
        'Color',colors2(2,:),'LineStyle','-','LineWidth',2);
    xlabel('Year (CE)')
    ylabel('Accumulation_{Ref.} (m_{ice} yr^{-1})')
    set(ax2,'Linewidth',1.5,'Fontsize',12)
    xlim([min(2015.5 - M.age), 2016])
