close all
clear all




%This changes data that is a) too shallow 
% or b) identified as outlier to nan. 
% The shallow limit is user configurable. 
% To identify outliers, a the interquarile range is computed for the data, with a moving window of 3 days by 3 meters that progresses on a timestep of of the user's choosing (you should probably stick to something between 3 hours and 1.5 days).
% Outliers are identified using with the upper and lower limits that are
% some multiple of the interquartile range. This multiple (outlier multiple) is user
% configurable. I have it set to a very conservative value of 10 right now.
% 
% This means that when the set of 3 m by 3 days range of data is looked at, an outlier is only flagged
% if it is 10 times the interquartile range of that data in that window.

% This code takes about 20 - 30 minutes to do a 2 month glider mission. If
% you want faster completion times, you can increase the time step by which
% the window moves, e.g. 1.5 day time steps does same job done in about 10 minutes. Don't go over 3 days or you will miss data altogether. 
% With a timesteps of 3 hours, the code runs about 3 times as fast using parallel
% processing, so the parallel processing is the default implementation. With a timestep of 1.5 days, the parallel processing is still about twice as fast as the
% forloop. The code can be still be run using a normal
% for loop, but I'd suggest using that only when you want to troubleshoot something, like exploring if you want to use 
% a different outlier multiple, as parloops are 
% difficult to troublshoot.The control to use the parloop or the forloop
% is in the user configurable settings.

% Parallel processing requires the Parallel computing toolbox. Some plots
% are output that show the results of the quality control. Also fair warning: for long
% glider deployments (like our 4 month summer one) this code wil likely run
% out of memory if you are doing a 3 hour timestep and you are working on
% other processes on your computer (like typing a dissertation). You can
% always just save out the tempAll variable for however far it got, adjust
% variables and missions and start again. My suggestion: run the code on a
% 1.5 day timestep setting while you're at lunch. run the code on a 3 hour timestep setting (if you think you need it) when you leave for the day. 

% If you want to adjust the moving window, the variables you are looking
% for in the code are: dz and dt_window. The window setting right now,
% assuming a profile occurs every 1.5 hours, provides a ~150 measurements
% total, 50 measurements per depth. I haven't played around with
% larger window sizes (say like 6 days), but it's probably worth doing.


%Once the quality control is complete the data is then gridded to 1 m bins for 0 to 300 by however many profiles are in
%the dataset. I think this gridding could be done in a way that is a little
%more flexible to future changes in variables, but it does the job for now.
%Can always improve later. I also think a rolling window for depth might be worth it.
% Other improvements that could be made that I would prioritize before those: (i)
%correcting for different sensor time-responses of the thermistor,
%conductivity and pressure sensors. (ii) thermal lag correction. (iii). 
%This is a good source for best practices: https://dx.doi.org/10.26198/5c997b5fdc9bd
%Author Isaac Reister 4/4/2025.

%variables loaded:
   % AllGliderVariables{1}=lats;
    % AllGliderVariables{2}=lons;
    % AllGliderVariables{3}=mydatenums;
    % AllGliderVariables{4}=chl;
    % AllGliderVariables{5}=pres;
    % AllGliderVariables{6}=cond;
    % AllGliderVariables{7}=temp;
    % AllGliderVariables{8}=par;
    % AllGliderVariables{9}=scatter;
    % AllGliderVariables{10}=do_c;
    % AllGliderVariables{11}=cdom;
    % AllGliderVariables{12}=Dfliu;
    % AllGliderVariables{13}=Dfliv;


%--------------notes-------
%v2 Rewrote everything to work with a moving window.IMporved output plots. made the comments for UAF gliders group. Added parallelization. 
% Moved the QC procedure away from a cell based matrix (which was too slow and
% clunky) and made everything run on 1 dimensional datasets, so that
% logical indices could be used more effectively. Gridding was also changed
% from a cell matrix to a a double matix, as we could now actually combine
% measurements within 1 meter bins, without worrying about outliers. File
% was renamed from DataProcessingGlider_PrelimGridDataCell_and_prelimQCv2_Parallel.m
% to DataProcessingGlider_QCglider_and_verticalgridding
% Finally this code was uploaded to github. All future changes will be documented there. 

%Zero Step: Partition our dataset down to just the mooring missions of
%interest for this chapter:


% NOTE NOW we are leaving the transecting missions behind!!! This code
% should still be able to run on them, but it's not fully tested.

% NOTE 9/4/2025, IR: suggestions for future changes:
%a) implement a flag system instead of just removing the data points.
%b) implement a manual check of the points prior to their removal.
%c) QC takes into account the depth, and only compares data within the same
%depth level. So when the glider is deep in the water column, 50 meter chl
%measurments are compared to other 50 meter chl measurments. The IQR for
%this depth is very small, but the sensor values still bounces around when
%it reaches the limits of its sensativity. I think this is why 
%we get more outliers at depth than at the surface for chl, but it might be
%worth investigating more.

%Right now the parallel loop and the normal loop provide the same result.
%Be aware that parallel loops run different from normal loop (main rule
%being each parallel loop has to be independent of the next). So there are
%slight differences in the how the data is stored in the parallel loop,
%versus the normal loop. Important to keep in mind when editing the code.

%% user input


%Important! Any depth measurement that registers
% less than 'tooshallowlimit' is replaced with a nan along with all associated variables.
%this is intended to ensure that the glider is actually underwater for our measurments. 
tooshallowlimt=0.25; %m  
om=10; % outlier multiplier.
useparallel=1;
steptime_in_seconds=10800;%the time window (spanning 3 days) moves forward in time 3 hours (10800 seconds) on every iteration.
%steptime_in_seconds= 129600; %different option,moves forward in time 1.5 days on every iteration.

enter_a_bad_value=0; %if equal to 1, this enters a bad value of 20 in the data at 10 meters
% somewhere within the timeseries. Useful to see if QC is actually working. Intended for Chlorophyll-a dataset. 

if useparallel == 1
    looptype = 'normal';
else
    looptype = 'parallel';
end

%QCfigurepath='C:\Users\funkb\Documents\MATLAB\Research\Figures\Chapter 3\QC\';
filenameappend=[char(string(looptype)),'loopv10']; %filenames already include mission number and short variables name. This appends whatever you want to that.
QCfigurepath=['C:\Users\funkb\Documents\MATLAB\Research\Figures\Chapter 3\QC\',char(string(looptype)),'loop\'];
QCdatasavepath='C:\Users\funkb\Documents\MATLAB\Research\data\Chapter3\'; %output directory for vertically gridded casts.
datafilename='QC_GC_v10'; %What you want to call the datafile.
yes_plot_QCfigure=1;



Stationkeepingsetall=logical([1,1,0,1,0,0,1,1]); %0s for missions you don't want to process. 
%Stationkeepingsetall=logical([0,1,0,0,0,0,0,0]); %for demo purposes. This
%is the shortest station keeping glider deployment.



%Adjust the variable names if needed. Note that processing more than ten
%variables may require modifiying how the processed data is saved, because
%the datafile can get too large.
    AllGliderVariablesnames{1}='chl';
    AllGliderVariablesnames{2}='pres';
    AllGliderVariablesnames{3}='cond';
    AllGliderVariablesnames{4}='temp';
    AllGliderVariablesnames{5}='par';
    AllGliderVariablesnames{6}='scatter';
    AllGliderVariablesnames{7}='do_c';
    AllGliderVariablesnames{8}='cdom';
    AllGliderVariablesnames{9}='Dfliu';
    AllGliderVariablesnames{10}='Dfliv';

%%************end of user input******************%

%identify first mission to process.
idx = find(Stationkeepingsetall);  % indices where it's 1
n = numel(Stationkeepingsetall);

for k = 1:numel(idx)
    Stationkeepingset = false(1,n);
    Stationkeepingset(idx(k)) = 1;
    

load('C:\Users\funkb\Documents\MATLAB\Research\data\Chapter3\DataProcessingGlider_AllGliderProfileSplit_v4.mat',"AllGliderVariables");

missionnumbers=1:length(Stationkeepingset);

for eachvariable=1:length(AllGliderVariables)
    

    SK_GliderVariables{eachvariable}=AllGliderVariables{eachvariable}(Stationkeepingset);
end
    clear AllGliderVariables
for eachvariable=1:length(SK_GliderVariables) %if a variable doesn't exist for a particular mission, replace with a nan string.
    for eachmission=1:length(SK_GliderVariables{5})
        if isempty(SK_GliderVariables{eachvariable}{eachmission})
            SK_GliderVariables{eachvariable}{eachmission}=nan(1,length(SK_GliderVariables{5}{eachmission}));
            
        end
    end

end
% For reference:
   % AllGliderVariables{1}=lats;
    % AllGliderVariables{2}=lons;
    % AllGliderVariables{3}=mydatenums;
    % AllGliderVariables{4}=chl;
    % AllGliderVariables{5}=pres;
    % AllGliderVariables{6}=cond;
    % AllGliderVariables{7}=temp;
    % AllGliderVariables{8}=par;
    % AllGliderVariables{9}=scatter;
    % AllGliderVariables{10}=do_c;
    % AllGliderVariables{11}=cdom;
    % AllGliderVariables{12}=Dfliu;
    % AllGliderVariables{13}=Dfliv;
for eachmission=1:sum(Stationkeepingset) %These are the variables we want to run through the very intensive QC.
    SK_GliderVariables2{1}{eachmission}=SK_GliderVariables{4}{eachmission};%chl;
    SK_GliderVariables2{2}{eachmission}=SK_GliderVariables{5}{eachmission};%pres;
    SK_GliderVariables2{3}{eachmission}=SK_GliderVariables{6}{eachmission};%cond;
    SK_GliderVariables2{4}{eachmission}=SK_GliderVariables{7}{eachmission};%temp;
    SK_GliderVariables2{5}{eachmission}=SK_GliderVariables{8}{eachmission};%par
    SK_GliderVariables2{6}{eachmission}=SK_GliderVariables{9}{eachmission};%scatter;
    SK_GliderVariables2{7}{eachmission}=SK_GliderVariables{10}{eachmission};%do_c;
    SK_GliderVariables2{8}{eachmission}=SK_GliderVariables{11}{eachmission};%cdom;
    SK_GliderVariables2{9}{eachmission}=SK_GliderVariables{12}{eachmission};%Dfliu;
    SK_GliderVariables2{10}{eachmission}=SK_GliderVariables{13}{eachmission};%Dfliv;
end


%remove all variables in SK_GliderVariables2 taht are less than tooshallowlimt

for eachmission=1:sum(Stationkeepingset) %These are the variables we want to run through the very intensive QC.
    missionalldepths=-1.*gsw_z_from_p(10*SK_GliderVariables2{2}{eachmission},60); % assume NGA lats.AllGliderVariables{5};
    tooshallow=missionalldepths<tooshallowlimt;
    SK_GliderVariables2{1}{eachmission}(tooshallow)=nan;
    SK_GliderVariables2{2}{eachmission}(tooshallow)=nan;
    SK_GliderVariables2{3}{eachmission}(tooshallow)=nan;
    SK_GliderVariables2{4}{eachmission}(tooshallow)=nan;
    SK_GliderVariables2{5}{eachmission}(tooshallow)=nan;
    SK_GliderVariables2{6}{eachmission}(tooshallow)=nan;
    SK_GliderVariables2{7}{eachmission}(tooshallow)=nan;
    SK_GliderVariables2{8}{eachmission}(tooshallow)=nan;
    SK_GliderVariables2{9}{eachmission}(tooshallow)=nan;
    SK_GliderVariables2{10}{eachmission}(tooshallow)=nan;
end

nVariables = size(SK_GliderVariables2,2);
nMissions  = sum(Stationkeepingset);
%preallocate
outlierslogic = cell(nVariables,1);   % outer cell, one per variable
for v = 1:nVariables
    outlierslogic{v} = cell(nMissions,1);  % inner cell, one per mission
end
            
temp_eachmission = cell(nVariables,1);   % outer cell, one per variable
for v = 1:nVariables
    temp_eachmission{v} = cell(nMissions,1);  % inner cell, one per mission
end


for eachmission=1:sum(Stationkeepingset) %for each mission...
   
    for eachvariable=1:size(SK_GliderVariables2,2) %for each variable within each mission...
    
        num_points = length(SK_GliderVariables{3}{eachmission}); %extract number of points..
        max_points = 3e6; % parallel processing will have memory issues if a glider deployment is too long (i.e. if the number of points is around 7e6). This defines the max points to be processed, and can be adjusted to accomodate computers with lower processing power.
        temp=SK_GliderVariables2{eachvariable}{eachmission}; %current variable dataset.
        time=SK_GliderVariables{3}{eachmission}; %extract time stamps.
        depth=-1.*gsw_z_from_p(10*SK_GliderVariables2{2}{eachmission},60); %obtain depth in meters.
        if sum(isnan(temp)) > round((9/10) * length(temp)) %if the variable doesn't exist: enter nans as place holder.
            outlierslogic{eachvariable}{eachmission} = temp;
            temp_eachmission{eachvariable}{eachmission} = temp;
            continue
        end
        if num_points > max_points
            % Split into chunks
            num_chunks = ceil(num_points / max_points);
            chunk_edges = round(linspace(1, num_points+1, num_chunks+1));
        

            chunking_control = 1;
            for c = 1:num_chunks
                idx_start = chunk_edges(c);
                idx_end   = chunk_edges(c+1)-1;
                idx_range = idx_start:idx_end;
        
                temp_chunk  = temp(idx_range);
                time_chunk  = time(idx_range);
                depth_chunk = depth(idx_range);
        
                % Run QC on just this chunk
                [outliers_chunk_output,tempAll_chunk_output] = run_glider_QC(temp_chunk,time_chunk,depth_chunk, eachmission, eachvariable,...
                    Stationkeepingset, missionnumbers, AllGliderVariablesnames, ...
                    steptime_in_seconds, om, ...
                    enter_a_bad_value, yes_plot_QCfigure, QCfigurepath, filenameappend, ...
                    useparallel,c);
        
                % Store into combined result
                outlierslogic{eachvariable}{eachmission} = [outlierslogic{eachmission}; outliers_chunk_output(:)];
                temp_eachmission{eachvariable}{eachmission} = [temp_eachmission{eachvariable}{eachmission}; tempAll_chunk_output(:)];
            end
        else
            % Small enough: just run once
            temp=SK_GliderVariables2{eachvariable}{eachmission}; %current variable dataset.
            time=SK_GliderVariables{3}{eachmission}; %extract time stamps.
            depth=-1.*gsw_z_from_p(10*SK_GliderVariables2{2}{eachmission},60); %obtain depth in meters.
            [outliers_output,temp_output] = run_glider_QC(temp,time,depth, eachmission, eachvariable,...
                Stationkeepingset, missionnumbers, AllGliderVariablesnames, ...
                steptime_in_seconds, om, ...
                enter_a_bad_value, yes_plot_QCfigure, QCfigurepath, filenameappend, ...
                useparallel);
            outlierslogic{eachvariable}{eachmission}=outliers_output(:);
            temp_eachmission{eachvariable}{eachmission}=temp_output(:);
        end

    %     vars = whos;  % Get info about all variables
    %     [~, idx_mem] = max([vars.bytes]);  % Find the index of the largest variable
    %     largestVar = vars(idx_mem);
    %     fprintf('Variable using the most memory: %s (%.2f MB)\n', ...
    % largestVar.name, largestVar.bytes/1e6);
    selectmissionnumbers=missionnumbers(Stationkeepingset);
    disp([char(string((eachvariable/size(SK_GliderVariables2,2))*100)),'% done with mission ',char(string(selectmissionnumbers))])
    end
end
    selectmissionnumbers=missionnumbers(Stationkeepingset);
    save([QCdatasavepath,'Tempmission',char(string(selectmissionnumbers)),'.mat'],"temp_eachmission","outlierslogic") %temporary save for conserving memory.
    
    disp(['Finished QC of mission ',char(string(selectmissionnumbers))])
end



%%
%recombining


load('C:\Users\funkb\Documents\MATLAB\Research\data\Chapter3\DataProcessingGlider_AllGliderProfileSplit_v4.mat',"downcasts","upcasts","AllGliderVariables");
SK_downcasts=downcasts(Stationkeepingsetall);
SK_upcasts=upcasts(Stationkeepingsetall);
for eachvariable=1:3
     GPSandtime_GliderVariables{eachvariable}=AllGliderVariables{eachvariable}(Stationkeepingsetall);
end
clear AllGliderVariables
nMissionsQC  = sum(Stationkeepingsetall);
tempAll = cell(nVariables,1);   % outer cell, one per variable
for v = 1:nVariables
    tempAll{v} = cell(nMissionsQC,1);  % inner cell, one per mission
end

processedmissions=missionnumbers(Stationkeepingsetall);
for eachmission=1:sum(Stationkeepingsetall)
    load([QCdatasavepath,'Tempmission',char(string(processedmissions(eachmission))),'.mat'])
     for eachvariable = 1:nVariables
        
        tempAll{eachvariable}{eachmission}=temp_eachmission{eachvariable}{1};
     
     end
end
%if needed, save the variables specified in the load command below. These
%are all the variables needed to continue.
%save([QCdatasavepath,'temporarysave_QC.mat'],'tempAll','SK_GliderVariables','SK_upcasts','SK_downcasts','Stationkeepingset')
for eachvariable=1:length(SK_GliderVariables) %if a variable doesn't exist for a particular mission, replace with a nan string.
    for eachmission=1:sum(Stationkeepingsetall)
        
        if eachvariable<=3
        SK_GliderVariablesQC{eachvariable}{eachmission}=GPSandtime_GliderVariables{eachvariable}{eachmission};
        else
        SK_GliderVariablesQC{eachvariable}{eachmission}=tempAll{eachvariable-3}{eachmission};
        end
   
    end

end


%First Step: 
%Grid all the data to a 1 m for 0 to 300 by however many profiles are in
%the dataset.





    for eachmission=1:sum(Stationkeepingsetall)

       

        mission_variable_downcast_index=SK_downcasts{eachmission};
        mission_variable_upcast_index=SK_upcasts{eachmission};
        tempgrid_d=nan(300,length(mission_variable_upcast_index),length(SK_GliderVariablesQC));
        tempgrid_u=nan(300,length(mission_variable_downcast_index),length(SK_GliderVariablesQC));

        %splitting out variables
        mission_depths=-1.*gsw_z_from_p(10*SK_GliderVariablesQC{5}{eachmission},60);
        mission_lats=SK_GliderVariablesQC{1}{eachmission};
        mission_lons=SK_GliderVariablesQC{2}{eachmission};
        mission_mydates=SK_GliderVariablesQC{3}{eachmission};
        mission_chl=SK_GliderVariablesQC{4}{eachmission};
        mission_cond=SK_GliderVariablesQC{6}{eachmission};
        mission_temp=SK_GliderVariablesQC{7}{eachmission};
        mission_par=SK_GliderVariablesQC{8}{eachmission};
        mission_scatter=SK_GliderVariablesQC{9}{eachmission};
        mission_do_c=SK_GliderVariablesQC{10}{eachmission};
        mission_cdom=SK_GliderVariablesQC{11}{eachmission};
        mission_u=SK_GliderVariablesQC{12}{eachmission};
        mission_v=SK_GliderVariablesQC{13}{eachmission};

        for eachprofile=1:length(mission_variable_downcast_index)
            bin_edges = [0:300]; 
            profiledowncastdepths=mission_depths(mission_variable_downcast_index(eachprofile,1):mission_variable_downcast_index(eachprofile,2));
            profiledowncastlats=mission_lats(mission_variable_downcast_index(eachprofile,1):mission_variable_downcast_index(eachprofile,2));
            profiledowncastlons=mission_lons(mission_variable_downcast_index(eachprofile,1):mission_variable_downcast_index(eachprofile,2));
            profiledowncastmydates=mission_mydates(mission_variable_downcast_index(eachprofile,1):mission_variable_downcast_index(eachprofile,2));
            profiledowncastchl=mission_chl(mission_variable_downcast_index(eachprofile,1):mission_variable_downcast_index(eachprofile,2));
            profiledowncastcond=mission_cond(mission_variable_downcast_index(eachprofile,1):mission_variable_downcast_index(eachprofile,2));
            profiledowncasttemp=mission_temp(mission_variable_downcast_index(eachprofile,1):mission_variable_downcast_index(eachprofile,2));
            profiledowncastpar=mission_par(mission_variable_downcast_index(eachprofile,1):mission_variable_downcast_index(eachprofile,2));
            profiledowncastscatter=mission_scatter(mission_variable_downcast_index(eachprofile,1):mission_variable_downcast_index(eachprofile,2));
            profiledowncastdo_c=mission_do_c(mission_variable_downcast_index(eachprofile,1):mission_variable_downcast_index(eachprofile,2));
            profiledowncastcdom=mission_cdom(mission_variable_downcast_index(eachprofile,1):mission_variable_downcast_index(eachprofile,2));
            profiledowncastu=mission_u(mission_variable_downcast_index(eachprofile,1):mission_variable_downcast_index(eachprofile,2));
            profiledowncastv=mission_v(mission_variable_downcast_index(eachprofile,1):mission_variable_downcast_index(eachprofile,2));


            profileupcastdepths=mission_depths(mission_variable_upcast_index(eachprofile,1):mission_variable_upcast_index(eachprofile,2));
            profileupcastlats=mission_lats(mission_variable_upcast_index(eachprofile,1):mission_variable_upcast_index(eachprofile,2));
            profileupcastlons=mission_lons(mission_variable_upcast_index(eachprofile,1):mission_variable_upcast_index(eachprofile,2));
            profileupcastmydates=mission_mydates(mission_variable_upcast_index(eachprofile,1):mission_variable_upcast_index(eachprofile,2));
            profileupcastchl=mission_chl(mission_variable_upcast_index(eachprofile,1):mission_variable_upcast_index(eachprofile,2));
            profileupcastcond=mission_cond(mission_variable_upcast_index(eachprofile,1):mission_variable_upcast_index(eachprofile,2));
            profileupcasttemp=mission_temp(mission_variable_upcast_index(eachprofile,1):mission_variable_upcast_index(eachprofile,2));
            profileupcastpar=mission_par(mission_variable_upcast_index(eachprofile,1):mission_variable_upcast_index(eachprofile,2));
            profileupcastscatter=mission_scatter(mission_variable_upcast_index(eachprofile,1):mission_variable_upcast_index(eachprofile,2));
            profileupcastdo_c=mission_do_c(mission_variable_upcast_index(eachprofile,1):mission_variable_upcast_index(eachprofile,2));
            profileupcastcdom=mission_cdom(mission_variable_upcast_index(eachprofile,1):mission_variable_upcast_index(eachprofile,2));
            profileupcastu=mission_u(mission_variable_upcast_index(eachprofile,1):mission_variable_upcast_index(eachprofile,2));
            profileupcastv=mission_v(mission_variable_upcast_index(eachprofile,1):mission_variable_upcast_index(eachprofile,2));

           
            


            [d_binned_values, d_bin_edges, d_Bin_I] = histcounts(profiledowncastdepths, bin_edges);
            [u_binned_values, u_bin_edges, u_Bin_I] = histcounts(profileupcastdepths, bin_edges);

            %Assign every variable to a large temporay grid.
            for eachdepthI=1:max(bin_edges) % eachdepthI=1 refers to depth bin 0 to 1, eachdepth=300 refers to depth bin 299 to 300

                d_thisdepthset=d_Bin_I==eachdepthI;
                u_thisdepthset=u_Bin_I==eachdepthI;
                if sum(d_thisdepthset)>0 %%upcasts. If there are things to bin, do so.
                    tempgrid_d(eachdepthI,eachprofile,5)=mean(profiledowncastdepths(d_thisdepthset),'omitmissing');
                    tempgrid_d(eachdepthI,eachprofile,1)=mean(profiledowncastlats(d_thisdepthset),'omitmissing');
                    tempgrid_d(eachdepthI,eachprofile,2)=mean(profiledowncastlons(d_thisdepthset),'omitmissing');
                    tempgrid_d(eachdepthI,eachprofile,3)=mean(profiledowncastmydates(d_thisdepthset),'omitmissing');
                    tempgrid_d(eachdepthI,eachprofile,4)=mean(profiledowncastchl(d_thisdepthset),'omitmissing');
                    tempgrid_d(eachdepthI,eachprofile,6)=mean(profiledowncastcond(d_thisdepthset),'omitmissing');
                    tempgrid_d(eachdepthI,eachprofile,7)=mean(profiledowncasttemp(d_thisdepthset),'omitmissing');
                    tempgrid_d(eachdepthI,eachprofile,8)=mean(profiledowncastpar(d_thisdepthset),'omitmissing');
                    tempgrid_d(eachdepthI,eachprofile,9)=mean(profiledowncastscatter(d_thisdepthset),'omitmissing');
                    tempgrid_d(eachdepthI,eachprofile,10)=mean(profiledowncastdo_c(d_thisdepthset),'omitmissing');
                    tempgrid_d(eachdepthI,eachprofile,11)=mean(profiledowncastcdom(d_thisdepthset),'omitmissing');
                    tempgrid_d(eachdepthI,eachprofile,12)=mean(profiledowncastu(d_thisdepthset),'omitmissing');
                    tempgrid_d(eachdepthI,eachprofile,13)=mean(profiledowncastv(d_thisdepthset),'omitmissing');
                    
                    
                end
                if sum(u_thisdepthset)>0 %upcasts
                    tempgrid_u(eachdepthI,eachprofile,5)=mean(profileupcastdepths(u_thisdepthset),'omitmissing');
                    tempgrid_u(eachdepthI,eachprofile,1)=mean(profileupcastlats(u_thisdepthset),'omitmissing');
                    tempgrid_u(eachdepthI,eachprofile,2)=mean(profileupcastlons(u_thisdepthset),'omitmissing');
                    tempgrid_u(eachdepthI,eachprofile,3)=mean(profileupcastmydates(u_thisdepthset),'omitmissing');
                    tempgrid_u(eachdepthI,eachprofile,4)=mean(profileupcastchl(u_thisdepthset),'omitmissing');
                    tempgrid_u(eachdepthI,eachprofile,6)=mean(profileupcastcond(u_thisdepthset),'omitmissing');
                    tempgrid_u(eachdepthI,eachprofile,7)=mean(profileupcasttemp(u_thisdepthset),'omitmissing');
                    tempgrid_u(eachdepthI,eachprofile,8)=mean(profileupcastpar(u_thisdepthset),'omitmissing');
                    tempgrid_u(eachdepthI,eachprofile,9)=mean(profileupcastscatter(u_thisdepthset),'omitmissing');
                    tempgrid_u(eachdepthI,eachprofile,10)=mean(profileupcastdo_c(u_thisdepthset),'omitmissing');
                    tempgrid_u(eachdepthI,eachprofile,11)=mean(profileupcastcdom(u_thisdepthset),'omitmissing');
                    tempgrid_u(eachdepthI,eachprofile,12)=mean(profileupcastu(u_thisdepthset),'omitmissing');
                    tempgrid_u(eachdepthI,eachprofile,13)=mean(profileupcastv(u_thisdepthset),'omitmissing');
                end

            end
           

        end
       
        tempgrid=nan(300,length(mission_variable_upcast_index)*2,length(SK_GliderVariablesQC));
        %Place up and down casts in order in the grid.
        tempgrid(:,1:2:end-1,:)=tempgrid_d;
        tempgrid(:,2:2:end,:)=tempgrid_u;
        QC_GC{eachmission}=tempgrid; %quality controlled glider profiles. rows are depths, columns are profiles, stacks are variables.
    end
%
   AllGliderVariables_forplot{1}='lats';
    AllGliderVariables_forplot{2}='lons';
    AllGliderVariables_forplot{3}='mydatenums';
    AllGliderVariables_forplot{4}='chl';
    AllGliderVariables_forplot{5}='pres';
    AllGliderVariables_forplot{6}='cond';
    AllGliderVariables_forplot{7}='temp';
    AllGliderVariables_forplot{8}='par';
    AllGliderVariables_forplot{9}='scatter';
    AllGliderVariables_forplot{10}='do_c';
    AllGliderVariables_forplot{11}='cdom';
    AllGliderVariables_forplot{12}='Dfliu';
    AllGliderVariables_forplot{13}='Dfliv';
%Quick plot of dataset.
 for eachmission=1:sum(Stationkeepingsetall)
     for eachvar=1:length(SK_GliderVariablesQC)
        figure(1)
        pcolor(QC_GC{eachmission}(:,:,eachvar))
        shading flat
        colorbar
        title([AllGliderVariables_forplot{eachvar},':Mission ',char(string(processedmissions(eachmission)))])
        set(gca, 'YDir', 'reverse')
     pause(1)
     end
    
 end


%save out.
 save([QCdatasavepath,datafilename,'.mat'],'QC_GC')



%variable organization saved:
   % AllGliderVariables{1}=lats;
    % AllGliderVariables{2}=lons;
    % AllGliderVariables{3}=mydatenums;
    % AllGliderVariables{4}=chl;
    % AllGliderVariables{5}=pres;
    % AllGliderVariables{6}=cond;
    % AllGliderVariables{7}=temp;
    % AllGliderVariables{8}=par;
    % AllGliderVariables{9}=scatter;
    % AllGliderVariables{10}=do_c;
    % AllGliderVariables{11}=cdom;
    % AllGliderVariables{12}=Dfliu;
    % AllGliderVariables{13}=Dfliv;