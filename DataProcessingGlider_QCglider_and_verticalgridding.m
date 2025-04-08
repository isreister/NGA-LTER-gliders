close all
clear all


load('C:\Users\funkb\Documents\MATLAB\Research\data\Chapter3\DataProcessingGlider_AllGliderProfileSplit_v4.mat',"downcasts","upcasts","AllGliderVariables");

%This changes data that is too shallow or identified as outlier to nan. 
% The shallow limit is user configurable. 
% To identify outliers, a the interquarile range is computed for the data, with a moving window of 3 days by 3 meters that progresses on a timestep of 3 hours.
% Outliers are identified using with the upper and lower limits that are
% some multiple of the interquartile range. This multiple (outlier multiple) is user
% configurable. I have it set to a very conservative value of 4 right now.
% This code takes about 20 - 30 minutes to do a 2 month glider mission. If
% you want faster completion times, you can increase the time step by which
% the window moves, e.g. 2 day time steps does same job done in about 5 minutes. Don't go over 3 days or you will miss data altogether. 
% With a timesteps of 3 hours, the code runs about 3 times as fast using parallel
% processing, so that is the default implementation. With a timestep of 2 days, the parallel processing is still about twice as fast as the
% forloop. The code can be still be run using a normal
% for loop, but I'd suggest using that only when you want to troubleshoot something, like exploring if you want to use 
% a different outlier multiple, as parloops are 
% difficult to troublshoot. The control to use the parloop or the forloop
% is in the user configurable settings.
% Parallel processing requires the Parallel computing toolbox. Some plots
% are output that show the results of the quality control. Also fair warning: for long
% glider deployments (like our 4 month summer one) this code wil likely run
% out of memory if you are doing a 3 hour timestep and you are working on
% other processes on your computer (like typing a dissertation). YOu can
% always just save out the tempAll variable for however far it got, adjust
% variables and missions and start again. Or better yet run this when you
% leave work.


%Once the quality control is complete the data is then gridded to 1 m bins for 0 to 300 by however many profiles are in
%the dataset. I think this gridding could be done in a way that is a little
%more flexible to future changes in variables, but it does the job for now.
%Can always improve later. Other improvements that could be made that I would prioritize: (i)
%correcting for different sensor time-responses of the thermistor,
%conductivity and pressure sensors. (ii) thermal lag correction. (iii). 
%PAR data should only be assessed on the upcast
%because there is a tilt to the PAR sensor that makes the downcast par
%inaccurate. These are just some of the best practices from here: https://dx.doi.org/10.26198/5c997b5fdc9bd
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
om=4; % outlier multiplier.
useparallel=1;
steptime_in_seconds=10800;%the time window (spanning 3 days) moves forward in time 3 hours (10800 secons) on every iteration.
%steptime_in_seconds= 172800; %different option,moves forwad in time 2 days.

enter_a_bad_value=0; %if equal to 1, this enters a bad value of 20 in the data at 10 meters
% somewhere within the timeseries. Useful to see if QC is actually working. Intended for Chlorophyll-a dataset. 
%note:  For some reaon PWS has a hard time picking out the erroneous datapoint. Not sure why.

QCfigurepath='C:\Users\funkb\Documents\MATLAB\Research\Figures\Chapter 3\QC\';
%QCfigurepath='C:\Users\funkb\Documents\MATLAB\Research\Figures\Chapter 3\QC\bigtimestep\';
QCdatasavepath='C:\Users\funkb\Documents\MATLAB\Research\data\Chapter3\'; %output directory for vertically gridded casts.
yes_plot_QCfigure=1;

Stationkeepingset=logical([1,0,0,1,0,0,1,1]); %0s for missions you don't want to process.
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
        dt_window=3;
        i_list = 1:mystep:length(time2)-mystep;
    num_points = length(time);
    outliersfinal = ones(num_points, 1); % Initialize once outside
    
    parfor i_idx = 1:numel(i_list)
        i = i_list(i_idx);
        t0 = time2(i);
        in_window = (time >= t0 - dt_window/2) & (time <= t0 + dt_window/2);
    
        if sum(in_window) < 100
            continue
        end
    
        tempinwindow = temp(in_window);
        windowed_abs_index = myabsindex(in_window);
        localoutliers = zeros(num_points, 1);
    
        z_edges = 0:dz:ceil(max(depth(in_window)));%looking at all depths here.
        [~, z_bin] = histc(depth(in_window), z_edges);
    
        local_outliersfinal = ones(num_points, 1);
        
        tempinwindow(z_bin == 0)=nan; %excludes shallow depths (which were already desingated as nans earlier) and any other measurement assoicated with a nan depth.
                
         for zb = 1:max(z_bin) %z_bin is not a depth, it is an index of the bins in z_edges. So 27 refers to 27th bin in zedges which is 78 meters.           
    
            vals = tempinwindow(z_bin == zb);
            
            localindex_thiddepth = windowed_abs_index(z_bin == zb);
    
            if sum(~isnan(vals)) < 100
                continue
            end
    
            if eachvariable == 1
                min_threshold = 0.1;
                mask_clipped = (vals <= min_threshold);
                valsforIQR = vals(~mask_clipped);
            elseif eachvariable == 5
               min_threshold = 1;
               mask_clipped = (vals <= min_threshold);
               valsforIQR = vals(~mask_clipped);
            else
                valsforIQR = vals;
            end
    
            Q1 = quantile(valsforIQR, 0.25);
            Q3 = quantile(valsforIQR, 0.75);
            IQR = Q3 - Q1;
            lower = Q1 - om * IQR;
            upper = Q3 + om * IQR;
            is_out = (vals < lower) | (vals > upper);
    
            localoutliers(localindex_thiddepth(is_out)) = 1;
    
            shiftingindex = ones(num_points, 1);
            shiftingindex(in_window) = 0;
            shiftinglogic = logical(shiftingindex);
            localoutliers(shiftinglogic) = 1;
    
            local_outliersfinal = logical(local_outliersfinal) & logical(localoutliers);
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
            title(['Mission ',char(string(eachmission)),' ' AllGliderVariablesnames{eachvariable}])
            
            legend('Original Data','Good (not shallow, not outliers)','Outliers')
            saveas(figure(1),[QCfigurepath,char(string(eachmission)),AllGliderVariablesnames{eachvariable},'parallel.png'])
            saveas(figure(1),[QCfigurepath,char(string(eachmission)),AllGliderVariablesnames{eachvariable},'parallel'])
            clf 
        end
        temp(~good)=nan;
        tempAll{eachvariable}{eachmission}=temp;
        
    
    
     end
    
    end
else

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%   Normal Forloop. Not parallel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



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
            somerealpoints=~isnan(temp) & (depth>10 & depth<11);
            realdepths=dindex(somerealpoints);
            temp(realdepths(100))=20;% insert a bad point to test.
    
            if eachmission==1
    
            end
        end
    
    
        if sum(isnan(temp))==length(temp)
            tempAll{eachvariable}{eachmission}=nan(1,length(time));
            continue
        end
        dz = 3; 
    
        N = length(temp);
        outliersfinal=ones(N,1); %assume all are bad to start.
        
        myabsindex=1:N;
        %mystep=259200; %my step should be less than the window if we want to overlap
        mystep=steptime_in_seconds;
        dt_window=3;
    
        
        for i = 1:mystep:length(time2)-round(mystep/2) %i operates within the 3 day chunk
                % Find window
                t0 = time2(i);
                display(string(convrtdt(t0)))
                %convrtdt(time(2))
                in_window = (time >= t0 - dt_window/2) & (time <= t0 + dt_window/2);
                %in_window = (time >= time(i)) & (time <= time(i+mystep-1));
                %timechunk=time(in_window);
                if sum(in_window)<100
                    continue
                end
                tempinwindow=temp(in_window);
                localindex = myabsindex(in_window);
                
                outliers = zeros(N,1);
                localoutliers = outliers(in_window);
                
                z_edges = 0:dz:ceil(max(depth(in_window))); %looking at all depths here.
    
                [~, z_bin] = histc(depth(in_window), z_edges);
    
                
                
                tempinwindow(z_bin == 0)=nan; %excludes shallow depths (which were already desingated as nans earlier) and any other measurement assoicated with a nan depth.
            for zb = 1:max(z_bin) %z_bin is not a depth, it is an index of the bins in z_edges. So 27 refers to 27th bin in zedges which is 78 meters.
                        
    
                        
                        vals = tempinwindow(z_bin==zb);
                        
                        localindex_thiddepth=localindex(z_bin==zb); %1,2,  56,57, and so on
                        if sum(~isnan(vals))<100
                            continue
                        end
                        % Compute IQR and outliers
                        % chlorohyll-a noise floor
                        if eachvariable==1
                            min_threshold = 0.1;
                            
                            % Values hitting the lower limit
                            mask_clipped = (vals <= min_threshold);
                            
                            % Use only non-clipped data to compute IQR
                            valsforIQR = vals(~mask_clipped);
                        elseif eachvariable == 5
                            min_threshold = 1;
                            mask_clipped = (vals <= min_threshold);
                            valsforIQR = vals(~mask_clipped);
            
                        else
                            valsforIQR=vals;
                        end
                        Q1 = quantile(valsforIQR, 0.25);
                        Q3 = quantile(valsforIQR, 0.75);
                        IQR = Q3 - Q1;
                        lower = Q1 - om * IQR;
                        upper = Q3 + om * IQR;
                        is_out = (vals < lower) | (vals > upper);
                        
                        %This is a helpful plot for seeing how QC works on an
                        %individual depth.
                        % figure(2)
                        %  clf
                        % blah=1:length(vals);
                        % plot(blah,vals)
                        % hold on 
                        % plot(blah(is_out),vals(is_out),'o','Color','r')
                       
                        outliers(localindex_thiddepth(is_out))=1; %as each depth is iterated through, this records the outliers for this timechunk
                        shiftingindex=ones(N,1);
                        shiftingindex(in_window)=0;
                        shiftinglogic=logical(shiftingindex); %all ones except for our current window
                        outliers(shiftinglogic)=1; %we just want to make the points that are not part of the time window be 1.
                        outliersfinal=logical(outliersfinal) & logical(outliers); %if outlierfinal(n) is 1 and outlier(n) is 0, the result is 0. 1s in outlier final persist only if outlierfinal(n) is 1 and outlier(n) is 1. 
                        
                        if sum(outliersfinal)==0
                            keyboard
                        end
                        outliersfinal=double(outliersfinal);
        
            end
            %% This is a helpful plot for seeing how the iteration goes while
            %% chlorophyll-a is being QCed.
            % figure(1);
            % clf
            % tempindex=1:length(temp);
            % plot(tempindex,temp)
            % xlim([1,16*10^4])
            % hold on
            % plot(tempindex(logical(outliersfinal)),temp(logical(outliersfinal)),'o','Color','r')
        
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
            title(['Mission ',char(string(eachmission)),' ' AllGliderVariablesnames{eachvariable}])
            
            legend('Original Data','Good (not shallow, not outliers)','Outliers')
            saveas(figure(1),[QCfigurepath,char(string(eachmission)),AllGliderVariablesnames{eachvariable},'normalloop.png'])
            saveas(figure(1),[QCfigurepath,char(string(eachmission)),AllGliderVariablesnames{eachvariable},'normalloop'])
            clf 
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
%load('C:\Users\funkb\Documents\MATLAB\Research\data\Chapter3\temporarysave_QC.mat','tempAll','SK_GliderVariables','SK_upcasts','SK_downcasts','Stationkeepingset')
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
%
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


 save([QCdatasavepath,'QC_GC_v8.mat'],'QC_GC')



