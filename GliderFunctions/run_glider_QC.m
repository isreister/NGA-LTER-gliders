function [outliersfinal, temp] = run_glider_QC(temp,time,depth, eachmission, eachvariable,...
    Stationkeepingset, missionnumbers, AllGliderVariablesnames, ...
    steptime_in_seconds, om, ...
    enter_a_bad_value, yes_plot_QCfigure, QCfigurepath, filenameappend, ...
    useparallel,c)

% This function performs outlier detection and quality control (QC) on 
% glider mission data. The logic follows the original script structure.
%
% Inputs:
%   SK_GliderVariables       - Cell array with time information
%   SK_GliderVariables2      - Cell array with data variables
%   Stationkeepingset        - Indices of missions to process
%   missionnumbers           - Identifiers for missions
%   AllGliderVariablesnames  - Names of variables for labeling
%   steptime_in_seconds      - Step size for QC time windows
%   om                       - Multiplier for IQR outlier detection
%   enter_a_bad_value        - Flag for inserting artificial bad value
%   yes_plot_QCfigure        - Flag for summary plotting
%   QCfigurepath             - Path to save plots
%   filenameappend           - Suffix for saved filenames
%   useparallel              - 1 = parfor, 0 = regular for loop
%
% Outputs:
%   outlierslogic            - Logical array of outlier flags
%   tempAll                  - QC’d data for each mission and variable



if nargin > 15 
    whatpart = ['part',char(string(c))];
else
    whatpart ='';
end


if useparallel==1

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%   Parallel Forloop
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


     
      
        
        time1=time(~isnan(time)); %non-nan times.
        time2=datenum(convrtdt(time1(1)):seconds(1):convrtdt(time1(end))); %create a regularly spaced timeseries, in seconds.
        
        

    
        if enter_a_bad_value==1
            %enter a bad value in place of a real 10 meter depth
            dindex=1:length(depth);
            somerealpoints=~isnan(depth) & (depth>10 & depth<11);
            realdepths=dindex(somerealpoints);
            temp(realdepths(1000))=20;% insert a bad point to test. 20 is a "bad" value for early spring chlorophyll-a. Value can be changed

        end
    
    
        if sum(isnan(temp))==length(temp) %if temp is nan for the entire timeseries, we skip processing it.
            tempAll{eachvariable}{eachmission}=nan(1,length(time));
            return
        end
        dz = 3; 
    
        N = length(temp);
        
        
        myabsindex=1:N; %the absolute index. This index only changes between missions.
        mystep=steptime_in_seconds;%my step is 3 hours.
        dt_window=3; % 3 day window is specified here.
        i_list = time2(1):(mystep/(24*60*60)):time2(end);%i_list is a vector of forward shifts in time, determined by a user selected timestep. The unit is datenum. (Usually 3 hours or up to 2 days)
    num_points = length(time);
    
    outliersMatrix = true(num_points, numel(i_list)); %Different from normal loop. Allows storing within parfor.
    if ~license('test','Distrib_Computing_Toolbox')
        error('Parallel Computing Toolbox is required to run this script with parfor.');
    end
    allnans=isnan(depth);
    %     vars = whos;  % Get info about all variables
    %     [~, idx] = max([vars.bytes]);  % Find the index of the largest variable
    %     largestVar = vars(idx);
    %     fprintf('Variable using the most memory in run_glider_QC: %s (%.2f MB)\n', ...
    % largestVar.name, largestVar.bytes/1e6);
    parfor i_idx = 1:numel(i_list) %i_list is a time vector. i_idx is the number of chunks.
        t0 = i_list(i_idx);
        
        in_window = (time >= t0 - dt_window/2) & (time <= t0 + dt_window/2); %Produces logic vector for "measurements within this window of 3 days (1.5 days on either side of t0)".
    
        if sum(in_window) < 100 %if there are less than 100 measurements within the window skip as this is not a sufficiently large sample size.
            continue
        end
    
        tempinwindow = temp(in_window); %actual values.
        windowed_abs_index = myabsindex(in_window); %absolute index within the window.
        
    
        z_edges = 0:dz:ceil(max(depth(in_window))); %create a regularly spaced depth vector.
        [~, z_bin] = histc(depth(in_window), z_edges);%bin according to those depths
    
        local_outliersfinal = true(num_points, 1); %this will assume all points are outliers for the start of the "in_window" set.
        
        tempinwindow(z_bin == 0)=nan; %excludes shallow depths (which were already desingated as nans earlier) and any other measurement assoicated with a nan depth.

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
            %local_outliersfinal(allnans)=0; %let nans pass through.

         end
    
        % Combine with outer loop
        outliersMatrix(:,i_idx) = local_outliersfinal; %.%once all depths are analyzed, outliers are recorded.store each result of inner loop. %Structure is different from normal loop. Allows storing within parfor.
        
        

    end
    
        outliersfinal = all(outliersMatrix, 2);   % equivalent to cumulative AND
        outliersfinal(allnans) = 0;
        
        z_edges = 0:dz:ceil(max(depth));%looking at all depths here.
        [~, z_bin] = histc(depth, z_edges);
        nobin=z_bin==0; %this identifies all the depths that are nan, which includes shallow depths.
        good= ~nobin & ~ outliersfinal;
        if yes_plot_QCfigure == 1
            fig = figure(2);
            fig.Visible = 'off';
            clf
            plot(convrtdt(time),temp,'.');
            hold on
            plot(convrtdt(time(good)),temp(good),'o','Color',[0.9294    0.6941    0.1255])
            plot(convrtdt(time(outliersfinal)),temp(outliersfinal),'o','Color','r', 'MarkerFaceColor','r')
            selectmissionnumbers=missionnumbers(Stationkeepingset);
            title(['Mission ',char(string(selectmissionnumbers(eachmission))),' ',whatpart,' ', AllGliderVariablesnames{eachvariable}])
            
            if isempty(temp(outliersfinal))
                legend('Original Data','Good (not shallow, no outliers detected)')
            else
               legend('Original Data','Good (not shallow, not outliers)','Outliers')
            end
            saveas(fig,[QCfigurepath,char(string(selectmissionnumbers(eachmission))),AllGliderVariablesnames{eachvariable}, filenameappend,whatpart,'.png'])
            saveas(fig,[QCfigurepath,char(string(selectmissionnumbers(eachmission))),AllGliderVariablesnames{eachvariable}, filenameappend,whatpart])
            clf %enter a break point here if you want a quick view of the plots
        end 
        temp(~good)=nan; %applies QC
        
        
    
    


else

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%   Normal Forloop. Not parallel
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
      

    
        if enter_a_bad_value==1
            %enter a bad value in place of a real 10 meter depth
            dindex=1:length(depth);
            somerealpoints=~isnan(depth) & (depth>10 & depth<11);
            realdepths=dindex(somerealpoints);
            temp(realdepths(1000))=20;% insert a bad point to test. 20 is a "bad" value for early spring chlorophyll-a. Value can be changed

        end
    
    
        if sum(isnan(temp))==length(temp) %if temp is nan for the entire timeseries, we skip processing it.
            tempAll{eachvariable}{eachmission}=nan(1,length(time));
            return
        end
        dz = 3; 
    
        N = length(temp);
        
        
        myabsindex=1:N; %the absolute index. This index only changes between missions.
        mystep=steptime_in_seconds;%my step is 3 hours.
        dt_window=3; % 3 day window is specified here.
        i_list = time2(1):(mystep/(24*60*60)):time2(end);%i_list is a vector of forward shifts in time, determined by a user selected timestep. The unit is datenum. (Usually 3 hours or up to 2 days)
    num_points = length(time);
    outliersfinal = true(num_points, 1); % %This is our outlier record. We assume all values are outliers to start.


    allnans=isnan(depth);
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
    
        local_outliersfinal = true(num_points, 1); %this will assume all points are outliers for the start of the "in_window" set.
        
        tempinwindow(z_bin == 0)=nan; %excludes shallow depths (which were already desingated as nans earlier) and any other measurement assoicated with a nan depth.
        
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
        outliersfinal = outliersfinal & local_outliersfinal; %once all depths are analyzed, outliers are recorded.%Structure is different from parallel loop. Allows for rolling QC plot.
        outliersfinal(allnans) = 0;
             %% This is a helpful plot for seeing how the iteration goes while
            %% chlorophyll-a is being QCed.
            figure(1);
            clf
            points_for_plotting_logic = (time >= t0 - dt_window/2 - dt_window/4) & (time <= t0 + dt_window/2 + dt_window/4);
            plot(time(points_for_plotting_logic),temp(points_for_plotting_logic))
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
    
    
        
        z_edges = 0:dz:ceil(max(depth));%looking at all depths here.
        [~, z_bin] = histc(depth, z_edges);
        nobin=z_bin==0; %this identifies all the depths that are nan, which includes shallow depths.
        good= ~nobin & ~ outliersfinal;
        if yes_plot_QCfigure == 1 %plots and save a summary QC plot.
            fig = figure(2);
            fig.Visible = 'off';
            clf
            plot(convrtdt(time),temp,'.');
            hold on
            plot(convrtdt(time(good)),temp(good),'o','Color',[0.9294    0.6941    0.1255])
            plot(convrtdt(time(outliersfinal)),temp(outliersfinal),'o','Color','r', 'MarkerFaceColor','r')
            selectmissionnumbers=missionnumbers(Stationkeepingset);
            title(['Mission ',char(string(selectmissionnumbers(eachmission))),' ',whatpart,' ', AllGliderVariablesnames{eachvariable}])
            
            if isempty(temp(outliersfinal))
                legend('Original Data','Good (not shallow, no outliers detected)')
            else
               legend('Original Data','Good (not shallow, not outliers)','Outliers')
            end
            saveas(fig,[QCfigurepath,char(string(selectmissionnumbers(eachmission))),AllGliderVariablesnames{eachvariable}, filenameappend,whatpart,'.png'])
            saveas(fig,[QCfigurepath,char(string(selectmissionnumbers(eachmission))),AllGliderVariablesnames{eachvariable}, filenameappend,whatpart])
            clf %enter a break point here if you want a quick view of the plots
        end 
        temp(~good)=nan; %applies QC
        
        
    
    

end

end
