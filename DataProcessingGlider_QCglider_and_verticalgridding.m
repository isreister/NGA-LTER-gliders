close all
clear all


load('C:\Users\funkb\Documents\MATLAB\Research\data\Chapter3\DataProcessingGlider_AllGliderProfileSplit_v4.mat',"downcasts","upcasts","AllGliderVariables");

%This changes data that is a) too shallow 
% or b) identified as outlier to nan. 
% The shallow limit is user configurable. 
% To identify outliers, a the interquarile range is computed for the data, with a moving window of 3 days by 3 meters that progresses on a timestep of of the user's choosing (you should probably stick to something between 3 hours and 1.5 days).
% Outliers are identified using with the upper and lower limits that are
% some multiple of the interquartile range. This multiple (outlier multiple) is user
% configurable. I have it set to a very conservative value of 4 right now.
% 
% This means that when the set of 3 m by 3 days range of data is looked at, an outlier is only flagged
% if it is 4 times the interquartile range of that dataset.

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

%% user input


%Important! Any depth measurement that registers
% less than 'tooshallowlimit' is replaced with a nan along with all associated variables.
%this is intended to ensure that the glider is actually underwater for our measurments. 
tooshallowlimt=0.25; %m  
om=10; % outlier multiplier.
useparallel=0;
steptime_in_seconds=10800;%the time window (spanning 3 days) moves forward in time 3 hours (10800 seconds) on every iteration.
%steptime_in_seconds= 129600; %different option,moves forward in time 1.5 days on every iteration.

enter_a_bad_value=0; %if equal to 1, this enters a bad value of 20 in the data at 10 meters
% somewhere within the timeseries. Useful to see if QC is actually working. Intended for Chlorophyll-a dataset. 


%QCfigurepath='C:\Users\funkb\Documents\MATLAB\Research\Figures\Chapter 3\QC\';
filenameappend='normalloopv10'; %filenames already include mission number and short variables name. This appends whatever you want to that.
QCfigurepath='C:\Users\funkb\Documents\MATLAB\Research\Figures\Chapter 3\QC\normalloop\';
QCdatasavepath='C:\Users\funkb\Documents\MATLAB\Research\data\Chapter3\'; %output directory for vertically gridded casts.
datafilename='QC_GC_v10'; %What you want to call the datafile.
yes_plot_QCfigure=1;



%Stationkeepingset=logical([1,1,0,1,0,0,1,1]); %0s for missions you don't want to process. 
Stationkeepingset=logical([1,1,0,0,0,0,0,0]); %for demo purposes
%%

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




SK_downcasts=downcasts(Stationkeepingset);
SK_upcasts=upcasts(Stationkeepingset);
missionnumbers=1:length(Stationkeepingset);

for eachvariable=1:length(AllGliderVariables)
    SK_GliderVariables{eachvariable}=AllGliderVariables{eachvariable}(Stationkeepingset);
end

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


if useparallel==1





    for eachmission=1:sum(Stationkeepingset)
    
     for eachvariable=1:size(SK_GliderVariables2,2)
     
      
        time=SK_GliderVariables{3}{eachmission};
        time1=SK_GliderVariables{3}{eachmission}(~isnan(SK_GliderVariables{3}{eachmission}));
        time2=datenum(convrtdt(time1(1)):seconds(1):convrtdt(time1(end)));
        depth=-1.*gsw_z_from_p(10*SK_GliderVariables2{2}{eachmission},60);
        
        % Initialize
        temp=SK_GliderVariables2{eachvariable}{eachmission};
    
        if enter_a_bad_value==1
            %enter a bad value in place of a real 10 meter depth
            dindex=1:length(depth);
            somerealpoints=~isnan(depth) & (depth>10 & depth<11);
            realdepths=dindex(somerealpoints);
            temp(realdepths(1000))=20;% insert a bad point to test.

        end
    
    
        if sum(isnan(temp))==length(temp)
            tempAll{eachvariable}{eachmission}=nan(1,length(time));
            continue
        end
        dz = 3; 
    
        N = length(temp);
        outliersfinal=ones(N,1); %assume all are bad to start.
        
        myabsindex=1:N;
        mystep=steptime_in_seconds;%my step is 3 hours.
        dt_window=3; % 3 day window is specified here.
        i_list = 1:mystep:length(time2);%i_list is a vector of forward shifts in time, determined by a user selected timestep. (Usually 3 hours or up to 2 days)
    num_points = length(time);
    outliersfinal = ones(num_points, 1); % Initialize once outside

    if ~license('test','Distrib_Computing_Toolbox')
        error('Parallel Computing Toolbox is required to run this script with parfor.');
    end
    
    parfor i_idx = 1:numel(i_list)
        i = i_list(i_idx);
        t0 = time2(i);
        in_window = (time >= t0 - dt_window/2) & (time <= t0 + dt_window/2);
    
        if sum(in_window) < 100
            continue
        end
    
        tempinwindow = temp(in_window);
        windowed_abs_index = myabsindex(in_window);
        local_outliers_per_depth = zeros(num_points, 1);
    
        z_edges = 0:dz:ceil(max(depth(in_window)));%looking at all depths here.
        [~, z_bin] = histc(depth(in_window), z_edges);
    
        local_outliersfinal = ones(num_points, 1);
        
        tempinwindow(z_bin == 0)=nan; %excludes shallow depths (which were already desingated as nans earlier) and any other measurement assoicated with a nan depth.
                
         for zb = 1:max(z_bin) %z_bin is not a depth, it is an index of the bins in z_edges. So 27 refers to 27th bin in zedges which is 78 meters.           
    
            valsforIQR = tempinwindow(z_bin == zb);
            
            localindex_thiddepth = windowed_abs_index(z_bin == zb);
    
            if sum(~isnan(valsforIQR)) < 100
                continue
            end
    
            if eachvariable == 1
                min_threshold = 0.1;
                mask_clipped = (valsforIQR <= min_threshold);
                valsforIQR = valsforIQR(~mask_clipped);
            elseif eachvariable == 5
        
               min_threshold = 1;
               mask_clipped = (valsforIQR <= min_threshold);
               valsforIQR = valsforIQR(~mask_clipped);
            else
                valsforIQR = valsforIQR;
            end
    
            Q1 = quantile(valsforIQR, 0.25);
            Q3 = quantile(valsforIQR, 0.75);
            IQR = Q3 - Q1;
            lower = Q1 - om * IQR;
            upper = Q3 + om * IQR;
            is_out = (valsforIQR < lower) | (valsforIQR > upper);
    
            local_outliers_per_depth(localindex_thiddepth(is_out)) = 1;
    
            shiftingindex = ones(num_points, 1);
            shiftingindex(in_window) = 0;
            shiftinglogic = logical(shiftingindex);
            local_outliers_per_depth(shiftinglogic) = 1;
    
            local_outliersfinal = logical(local_outliersfinal) & logical(local_outliers_per_depth);
        end
    
        % Combine with outer loop
        outliersfinal = outliersfinal & double(local_outliersfinal);
    end
    
    
        outlierslogic{eachvariable}{eachmission}=logical(outliersfinal);
        z_edges = 0:dz:ceil(max(depth));%looking at all depths here.
        [~, z_bin] = histc(depth, z_edges);
        nobin=z_bin==0; %this identifies all the depths that are nan, which includes shallow depths.
        good= ~nobin & ~ outlierslogic{eachvariable}{eachmission};
        if yes_plot_QCfigure == 1
            figure(1);
            plot(convrtdt(SK_GliderVariables{3}{eachmission}),temp,'.');
            hold on
            plot(convrtdt(SK_GliderVariables{3}{eachmission}(good)),temp(good),'o','Color',[0.9294    0.6941    0.1255])
            plot(convrtdt(SK_GliderVariables{3}{eachmission}(outlierslogic{eachvariable}{eachmission})),temp(outlierslogic{eachvariable}{eachmission}),'o','Color','r', 'MarkerFaceColor','r')
            selectmissionnumbers=missionnumbers(Stationkeepingset);
            title(['Mission ',char(string(selectmissionnumbers(eachmission))),' ' AllGliderVariablesnames{eachvariable}])
            
            if isempty(temp(outlierslogic{eachvariable}{eachmission}))
                legend('Original Data','Good (not shallow, no outliers detected)')
            else
               legend('Original Data','Good (not shallow, not outliers)','Outliers')
            end
            saveas(figure(1),[QCfigurepath,char(string(eachmission)),AllGliderVariablesnames{eachvariable}, filenameappend,'.png'])
            saveas(figure(1),[QCfigurepath,char(string(eachmission)),AllGliderVariablesnames{eachvariable}, filenameappend])
            clf %enter a break point here if you want a quick view of the plots
        end 
        temp(~good)=nan;
        tempAll{eachvariable}{eachmission}=temp;
        
    
    
     end
    
    end
else

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%   Normal Forloop. Not parallel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



     for eachmission=1:sum(Stationkeepingset) %for each mission...
    
     for eachvariable=1:size(SK_GliderVariables2,2) %for each variable within each mission...
     
      
        time=SK_GliderVariables{3}{eachmission}; %extract time stamps.
        time1=SK_GliderVariables{3}{eachmission}(~isnan(SK_GliderVariables{3}{eachmission})); %non-nan times.
        time2=datenum(convrtdt(time1(1)):seconds(1):convrtdt(time1(end))); %create a regularly spaced timeseries, in seconds.
        depth=-1.*gsw_z_from_p(10*SK_GliderVariables2{2}{eachmission},60); %obtain depth in meters.
        
        % Initialize
        temp=SK_GliderVariables2{eachvariable}{eachmission}; %current variable.
    
        if enter_a_bad_value==1
            %enter a bad value in place of a real 10 meter depth
            dindex=1:length(depth);
            somerealpoints=~isnan(depth) & (depth>10 & depth<11);
            realdepths=dindex(somerealpoints);
            temp(realdepths(1000))=20;% insert a bad point to test. 20 is a "bad" value for early spring chlorophyll-a. Value can be changed

        end
    
    
        if sum(isnan(temp))==length(temp) %if temp is nan for the entire timeseries, we skip processing it.
            tempAll{eachvariable}{eachmission}=nan(1,length(time));
            continue
        end
        dz = 3; 
    
        N = length(temp);
        
        
        myabsindex=1:N; %the absolute index. This index only changes between missions.
        mystep=steptime_in_seconds;%my step is 3 hours.
        dt_window=3; % 3 day window is specified here.
        i_list = time2(1):(mystep/(24*60*60)):time2(end);%i_list is a vector of forward shifts in time, determined by a user selected timestep. The unit is datenum. (Usually 3 hours or up to 2 days)
    num_points = length(time);
    outliersfinal = logical(ones(num_points, 1)); % %This is our outlier record. We assume all values are outliers to start.


    
    for i_idx = 1:numel(i_list) %i_list is a time vector. i_idx is the number of chunks.
        t0 = i_list(i_idx);
        
        in_window = (time >= t0 - dt_window/2) & (time <= t0 + dt_window/2); %Produces logic vector for "measurements within this window of 3 days (1.5 days on either side of t0)".
    
        if sum(in_window) < 100 %if there are less than 100 measurements within the window skip as this is not a sufficiently large sample size.
            continue
        end
    
        tempinwindow = temp(in_window); %actual values.
        windowed_abs_index = myabsindex(in_window); %absolute index within the window.
        
    
        z_edges = 0:dz:ceil(max(depth(in_window))); %create a regularly spaced depth vector.
        [~, z_bin] = histc(depth(in_window), z_edges);%bin according to those depths
    
        local_outliersfinal = logical(ones(num_points, 1)); %this will assume all points are outliers for the start of the "in_window" set.
        
        tempinwindow(z_bin == 0)=nan; %excludes shallow depths (which were already desingated as nans earlier) and any other measurement assoicated with a nan depth.
        outliercount = 0;
        allnans=isnan(depth);
        outlierindex =[];
         for zb = 1:max(z_bin) %z_bin is not a depth, it is an index of the bins in z_edges. So 4 refers to 4th bin in zedges, the bin between 9 and 12 meters.           
            local_outliers_per_depth = ones(num_points, 1); %this will assume all values are outliers to start
            valsforIQR = tempinwindow(z_bin == zb); %aquire values for the 3 meter sets.
            
            localindex_nans=isnan(valsforIQR);
            localindex_thiddepth = windowed_abs_index(z_bin == zb); %aquire index for 3 meter sets.
          

            if sum(~isnan(valsforIQR)) < 100

                local_outliers_per_depth(localindex_thiddepth) = 0; %if there is not enough points to assess outlier status, we assume the data is good.
                local_outliersfinal(localindex_thiddepth) = local_outliersfinal(localindex_thiddepth) & local_outliers_per_depth(localindex_thiddepth);
                continue
            end
                  
            Q1 = quantile(valsforIQR, 0.25);
            Q3 = quantile(valsforIQR, 0.75);
            IQR = Q3 - Q1;
            lower = Q1 - om * IQR;
            upper = Q3 + om * IQR;
            is_out = (valsforIQR < lower) | (valsforIQR > upper); %outlier marked
            
            local_outliers_per_depth(localindex_thiddepth(~is_out)) = 0; %outliers remain as one in the timeseries.

            % remove later. Troubleshooting. If the outlier detector is
            % working, 311429 should always be identified as an outlier,
            % because for each interation through that depth bin, the 20
            % value should be flagged.


            local_outliersfinal(localindex_thiddepth) = local_outliersfinal(localindex_thiddepth) & local_outliers_per_depth(localindex_thiddepth);
             %As we iterate through all the depths, local_outliersfinal holds the outliers thoughout all depth iterations. local outliers carries the outliers for each depth.
             

         end
    
        % Combine with outer loop
        outliersfinal = outliersfinal & local_outliersfinal; %once all depths are analyzed, outliers are recorded.
        outliersfinal(allnans) = 0;
             %% This is a helpful plot for seeing how the iteration goes while
            %% chlorophyll-a is being QCed.
            figure(1);
            clf
            
            plot(time,temp)
            xlim([t0 - dt_window/2 - dt_window/4,t0 + dt_window/2 + dt_window/4]) %this plot is our 3 day window plus a little extra on either side. This is shifted by 3 hours (which causeing the red markers to disappear as points are verified to not be outliers)
            hold on
            plot(time(outliersfinal),temp(outliersfinal),'o','Color','r')
            plot([t0 - dt_window/2,t0 - dt_window/2],[max(temp),min(temp)],'--k') %this marks the trailing edge of our 3 day window.
            plot([t0 + dt_window/2,t0 + dt_window/2],[max(temp),min(temp)],'--k') %this marks the leading edge of our 3 day window.
            selectmissionnumbers=missionnumbers(Stationkeepingset);
            title(['Mission ',char(string(selectmissionnumbers(eachmission))),' ' AllGliderVariablesnames{eachvariable}])
            ylabel(AllGliderVariablesnames{eachvariable})
            xlabel('Datenum')
    end
    
    
        outlierslogic{eachvariable}{eachmission}=outliersfinal;
        z_edges = 0:dz:ceil(max(depth));%looking at all depths here.
        [~, z_bin] = histc(depth, z_edges);
        nobin=z_bin==0; %this identifies all the depths that are nan, which includes shallow depths.
        good= ~nobin & ~ outlierslogic{eachvariable}{eachmission};
        if yes_plot_QCfigure == 1
            figure(2);
            clf
            plot(convrtdt(SK_GliderVariables{3}{eachmission}),temp,'.');
            hold on
            plot(convrtdt(SK_GliderVariables{3}{eachmission}(good)),temp(good),'o','Color',[0.9294    0.6941    0.1255])
            plot(convrtdt(SK_GliderVariables{3}{eachmission}(outlierslogic{eachvariable}{eachmission})),temp(outlierslogic{eachvariable}{eachmission}),'o','Color','r', 'MarkerFaceColor','r')
            selectmissionnumbers=missionnumbers(Stationkeepingset);
            title(['Mission ',char(string(selectmissionnumbers(eachmission))),' ' AllGliderVariablesnames{eachvariable}])
            
            if isempty(temp(outlierslogic{eachvariable}{eachmission}))
                legend('Original Data','Good (not shallow, no outliers detected)')
            else
               legend('Original Data','Good (not shallow, not outliers)','Outliers')
            end
            saveas(figure(2),[QCfigurepath,char(string(selectmissionnumbers(eachmission))),AllGliderVariablesnames{eachvariable}, filenameappend,'.png'])
            saveas(figure(2),[QCfigurepath,char(string(selectmissionnumbers(eachmission))),AllGliderVariablesnames{eachvariable}, filenameappend])
            clf %enter a break point here if you want a quick view of the plots
        end 
        temp(~good)=nan;
        tempAll{eachvariable}{eachmission}=temp;
        
    
    
     end
    
    end

end

%%
%recombining
%if needed, save the variables specified in the load command below. These
%are all the variables needed to continue.
%save([QCdatasavepath,'temporarysave_QC.mat'],'tempAll','SK_GliderVariables','SK_upcasts','SK_downcasts','Stationkeepingset')
for eachvariable=1:length(SK_GliderVariables) %if a variable doesn't exist for a particular mission, replace with a nan string.
    for eachmission=1:length(SK_GliderVariables{5})
        
        if eachvariable<=3
        SK_GliderVariablesQC{eachvariable}{eachmission}=SK_GliderVariables{eachvariable}{eachmission};
        else
        SK_GliderVariablesQC{eachvariable}{eachmission}=tempAll{eachvariable-3}{eachmission};
        end
   
    end

end


%First Step: 
%Grid all the data to a 1 m for 0 to 300 by however many profiles are in
%the dataset.





    for eachmission=1:sum(Stationkeepingset)

       

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

            %all other variables:
            


            [d_binned_values, d_bin_edges, d_Bin_I] = histcounts(profiledowncastdepths, bin_edges);
            [u_binned_values, u_bin_edges, u_Bin_I] = histcounts(profileupcastdepths, bin_edges);
            for eachdepthI=1:max(bin_edges) % eachdepthI=1 refers to depth bin 0 to 1, eachdepth=300 refers to depth bin 299 to 300

                d_thisdepthset=d_Bin_I==eachdepthI;
                u_thisdepthset=u_Bin_I==eachdepthI;
                if sum(d_thisdepthset)>0 %if there are things to bin, do so.
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
                if sum(u_thisdepthset)>0
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
      
        tempgrid(:,1:2:end-1,:)=tempgrid_d;
        tempgrid(:,2:2:end,:)=tempgrid_u;
        QC_GC{eachmission}=tempgrid; %quality controlled glider profiles. rows are depths, columns are profiles, stacks are variables.
    end
%%
 for eachmission=1:sum(Stationkeepingset)
     for eachvar=1:length(SK_GliderVariablesQC)
        figure(1)
        pcolor(QC_GC{eachmission}(:,:,eachvar))
        shading flat
        colorbar
     pause(1)
     end
    
 end
%empty out any sheets that are just nans.


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