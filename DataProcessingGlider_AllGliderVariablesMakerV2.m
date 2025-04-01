%GliderData Processing
%This prepares high resolution data for processing by standardizing the glider names.
% Removes extremely bad data points (like data points over land or far far away).
% Produces some VERY basic positional plots. And saves out the slightly
% cleaned up data as flight and science variables. Here are the notes on
% the save out files:


% Here is the cell organization:
%Level 1 (top level)
    % AllGliderVariables{1}=lats;
    % AllGliderVariables{2}=lons;
    % AllGliderVariables{3}=mydatenums;
    % AllGliderVariables{4}=chl;
    % AllGliderVariables{5}=dep;
    % AllGliderVariables{6}=sal;
    % AllGliderVariables{7}=temp;
    % AllGliderVariables{8}=par;
    % AllGliderVariables{9}=scatter;
    % AllGliderVariables{10}=do_c;
    % AllGliderVariables{11}=cdom;
    % %flight variables.
    % AllGliderFliVariables{1}=Dflidates;
    % AllGliderFliVariables{2}=Dflidepth;
    % AllGliderFliVariables{3}=Dfliinflections;
% Level 2
    %Each variable is a cell array wherein the a number of cells that correspond to
    %however many missions were processed in 'DataProcessingGlider_AllGliderVariablesMakerV2.m
% Level 3
    %Each of those cells then contains a column that is the actual data.
    %So for example if I wanted to see all the latitudes for the first
    %mission I would execute: AllGliderVariables{1}{1} 
    %this will produce a very long column vector.

% Author: Isaac Reister 3/31/2025
%V2 cleaned up comments for UAF glider group
%V1 intial creation
clear all
close all

%% User Inputs


%Here we load up .dat matlab files that were downloaded from given to me
%from hank:
%https://drive.google.com/drive/folders/1jiGXYKbKBLgEyQCFpeUC9TjxYIhIFIuV
%which is the GOA folder he shared with my isaac.reister@gmail.com account.
%Question: That GOA data is fully resoleved and have been post processed to
%make calculations. 
% Here is what Hank said in the email: 
% I'll keep adding the data from our ongoing missions in the Gulf to this directory after we recover the gliders and have access to the full resolution data sets.
% All data sets have been post-processed using the TEOA-2010 (https://www.teos-10.org/) to compute absolute salinity, practical salinity, in situ density and potential
% density anomaly with a reference pressure of 0 dbar (Sigma-Theta). [I'm going to assume you are already familiar with all of the variable names, but just in case, all ' ...
% 'of the sci.NAME variables are related to science and the fli.NAME variables are derived from the flight computer. Time is in UTC for all data files and I use the matlab ' ...
% 'datenum function throughout the processing. This data is provided as a raw timeseries, so you have all of the complete yo']s with the highest resolution data we 
% collected during the glider missions.

%Note to self: QC has been not been done on this data.

%Notes on gliders:
%191 has PAR, 507 does not have PAR. 
datarepo{1}='g191_2021_herringexpt_data'; % Total Mission Timeline: 01/20 - 05/27/2021, Deployment: PWS, Recovery: GAK1, Route Description: PWS and a shoreline transect back to GAK1. In sound unitl May 5.
datarepo{2}='g507_2022_Seward_Line_data'; % Total Mission Timeline: 02/14 - 04/11/2022, Deployment: GAK2, Recovery: GAK1, Route Description: GAK2 to GAK7. Station Keeping at GAK 7 from 02/22-02/25. Transects out to GAK 15 reaching GAK 15 on 03/2/2022. Then begins a wide loop east of the Seward Line reaching GAK1 04/19/2022, station keeping until recovery. 
datarepo{3}='g507_2022_Yakutat_to_Seward_data'; %Total Mission Timeline: 06/04 - 07/15/2022, Deployment: Yaktat Bay, Recovery: GAK14, Route Description: Zig zagging along the Alaskan shelf break.
datarepo{4}='g191_2023_Seward_Line_data'; % Total Mission Timeline: 02/15 - 05/10/2023, Deployment: Res 2.5, Recovery: GEO, Route Description: Transected from Res 2.5 to GAK1. Station kept at GAK1 from 02/17-02/22. Then transected to GEO mooring reaching GEO on 03/02, station keeping until recovery.
datarepo{5}='g507_2023_Seward_Line_data'; % Total Mission Timeline: 03/21 - 05/11/2023, Deployment: Res 2.5, Recovery: GAK14, Route Description: Transected from Res 2.5 to GAK1.Station kept at GAK1 from 03/22-03/24. Then transected to GEO mooring reaching GEO on 04/01. Station kept at GEO until 04/04. Transected to GAK 15 then returned to GAK 14 for recovery.
datarepo{6}='g507_2024_Seward_Line_data'; %  Total Mission Timeline: 03/12 - 04/28/2024, Deployment: Res 2.5, Recovery: GAK14, Route Description: Transected from Res 2.5 to GAK1.Station kept at GAK1 from 03/14-03/23.Then transected to GAK 14 where it was recovered.
datarepo{7}='g191_2024_Seward_Line_data_v3'; % Total Mission Timeline: 02/23 - 06/22/2024, Deployment: Res 2.5, Recovery: GAK1, Route Description: Transected from Res 2.5 to GAK1. Station kept at GAK1 from 03/10-03/13. Then transected to GEO mooring reaching GEO on 03/20. Station kept at GEO until 06/18. Transected to GAK 1 for recovery on along a pathway east of the Seward Line.
datarepo{8}='g1152_2024_Seward_Line_data'; %  Glider 1152's 2024 transect + GEO Mooring station keeping from June - September

myroot='C:\Users\funkb\Documents\MATLAB\Research\data\Chapter3\'; %this is where your data lives.
myoutput='C:\Users\funkb\Documents\MATLAB\Research\data\Chapter3\'; %this is where you want your slightly more processed data to live.
%load('C:\Users\funkb\Documents\MATLAB\Research\data\CopperRiverPaper\NGAlatlons.mat');
missioncolors=jet(length(datarepo));

    ngaminLAT=57.5;
    ngamaxLAT=61;
    ngaminLON=-150+360;
    ngamaxLON=-143+360;

    geominLAT=58.75;
    geomaxLAT=59.25;
    geominLON=211;
    geomaxLON=211.5;
   


    %% Load data, Standardize variable names, and create cell arrays for variables.
    for ix=1:length(datarepo)
        load([myroot,datarepo{ix}]);
      

    
        if isfield(sci, 'sci_time') %corrects the sci field name for time to be consistent throughout all deployments.
            sci.time=sci.sci_time;
            sci=rmfield(sci, 'sci_time');
        end
    
        if isfield(sci, 'pres') %corrects the sci field name for time to be consistent throughout all deployments.
            sci.press=sci.pres;
            sci=rmfield(sci, 'pres');
        end
    
        if isfield(sci, 'sci_depth') %corrects the sci field name for depth to be consistent throughout all deployments.
            sci.depth=sci.sci_depth;
            sci=rmfield(sci, 'sci_depth');
        end

        if isfield(sci, 'sal') %corrects the sci field name for time to be consistent throughout all deployments.
            sci.psal=sci.sal;
            sci=rmfield(sci, 'sal');
        end

        if isfield(sci,'do')
            sci.do_c=sci.do;
            sci=rmfield(sci,'do');
        end
        
        if isfield(sci,'do_conc')
            sci.do_c=sci.do_conc;
            sci=rmfield(sci,'do_conc');
        end

        %replace all depth with recalculated depth from pressure for
        %mission 7 since that one had this problem: 
        %Hank: So, this is a pretty good example of what transpired on the 2024 Spring Shackleton deployment. This was unknown to us during deployment, but we had a folder on the flight computer that was a kind of placeholder called "./Trash" and in it was all of the 2023 flight data we collected. The flight control CF card was backed up and then it's contents were cut and copied to the computer's "trash" folder, however, the computer's trash folder was not emptied while the CF card was still in the CF card reader. The computer created that ."Trash" folder as a backup just in case that delete was unintentional. Fast Forward to 2024, I started seeing that the flight computer's CF card was filling up fast, much faster than anticipated and this was causing problems because the flight computer had no where to write data to. So, long story short, we made space and that involved corrupting all of the raw, DBD files onboard the glider up to that day (April 24, 2024 00:21:43 UTC). However, the science computer has no such issues, so it was recording data at full resolution the entire deployment.
        %There are several data variables that are used by both the Flight and Science computers onboard the glider, as you so rightly point out in the attached snippet of code, Depth is one of those variables. There's two places where depth is measured on the glider, one is from a pressure port on the aft tail cone and the second is from the CTD's pressure channel. The raw pressure (and hence depth) data measured from the tail cone sensor was mostly lost due to the DBD data purge mid-flight. We do have the SBD data files though, so that saved the pressure timeseries from the tail-cone sensor, albeit at a much reduced volume because that is the data we sent via iridium. In my re-processing of the data, I merged the flight SBD data prior to April 24 with the flight DBD data post April 24 and that is why when you plot sci.depth, it looks like there is a mismatch on the sample rate in the earlier part of the deployment. I would try and stick to using the CTD pressure channel for computing the vehicle's depth. If you want to be fancy, use the function "gsw_z_from_p" from https://www.teos-10.org/. This variable is actually already computed and included in the matlab distrubition, it's called sci.pdepth, short for "potential depth" (edited) 
        if ix==7
            if isfield(sci, 'depth')
               sci.depth=[];
               sci.depth=-1.*gsw_z_from_p(10*sci.press,60);
            end
        end

%science variables
        Dlat{ix}=sci.lat;
        Dlon{ix}=sci.lon;
        Dscidates{ix}=sci.time;
        Ddatemax{ix}=max(Dscidates{ix});
        Ddatemin{ix}=min(Dscidates{ix});
        Dchl{ix}=sci.chl;
        Dpres{ix}=sci.press;
        Dsal{ix}=sci.psal;
        Dtemp{ix}=sci.temp;
        Ddepth{ix}=sci.depth;
        Dbb{ix}=sci.bb;
%flight variables    
        Dflidates{ix}=fli.gtime;
        Dflidepth{ix}=fli.depth;
        Dfliinflections{ix}=fli.inflecting;
        
   
       
%note that the sci depth has a higher resolution than the fli depth.


        % Vars that only some gliders have: par,  %do_c or do (dissolved oxygen) and cdom, 
        if isfield(sci,'par') & length(sci.par)>1
            Dpar{ix}=sci.par;
        else
            Dpar{ix}=[];
        end

        if isfield(sci,'do_c') & length(sci.do_c)>1% this seperates do_c from do_p (which is percentage)
            %if the field is empty length is zero, if the field is a nan, length is 1.
            Ddo_c{ix}=sci.do_c;
            
        else
            Ddo_c{ix}=[];
        end

        if isfield(sci,'cdom') & length(sci.cdom)>1
            Dcdom{ix}=sci.cdom;
        else
            Dcdom{ix}=[];
        end
    
        if sum(~isnan(Dpres{ix}))==0
            keyboard 
            %there is no pressure data here
        end
    end
    

    m_proj('lambert','long',[ngaminLON ngamaxLON],'lat', [ngaminLAT ngamaxLAT]); 
    %uncomment if you want to see some land. takes more time to plot.
    % m_gshhs_h('patch',[ 1.0000    1.0000    0.8196]);%mellow yellow
    % m_grid('xlabeldir','end','fontsize',10,'YaxisLocation','left','FontSize',14);
    
for ix=1:length(datarepo)

    outofrangelat=(Dlat{ix}>61.5) | (Dlat{ix}<40);
    outofrangelon=(Dlon{ix}>-130) | (Dlon{ix}<-150);

    outofrange=outofrangelat+outofrangelon;
    outofrange(outofrange==2)=1; %makeing it so true+true=2;
    outofrange=logical(outofrange);


    Dlat{ix}(outofrange)=nan;
    Dlon{ix}(outofrange)=nan;

    datemaxes{ix}=datetime(Ddatemax{ix},'ConvertFrom','datenum');
    datemins{ix}=datetime(Ddatemin{ix},'ConvertFrom','datenum');

    % Format the datetime to character strings
    dateRange{ix} = strcat(datestr(datemins{ix}, 'mm/dd'), ' - ', datestr(datemaxes{ix}, ' mm/dd/yyyy'));
    mysize=10;
    trancolor=missioncolors(ix,:);
    figure(1)
    hold on
    % m_scatter(goldgaklonlat(:,1)+360,goldgaklonlat(:,2),20,'Marker','o','MarkerFacecolor','y')
    % m_scatter(Dlon{ix}+360,Dlat{ix},mysize,'Marker','o','MarkerFacecolor',trancolor,'MarkerEdgecolor',trancolor,'MarkerFaceAlpha',.3,'MarkerEdgeAlph',.3);
    hold off
    figure(2)
    hold on
    mhandle{ix}=plot(datetime(Dscidates{ix},'ConvertFrom','datenum'),Dlat{ix},'Color',trancolor,'LineWidth',2);
    hold off
    display(dateRange{ix})
 
end
ylim([58.5,60])
legend([mhandle{1:end}], dateRange{1:end})



%% Narrow down our box some more:

for ix=1:length(datarepo)

    %remove non-NGA data


    outofrangelat=(Dlat{ix}>90) | (Dlat{ix}<40);
    outofrangelon=(Dlon{ix}>-130) | (Dlon{ix}<-150);



    outofrange=outofrangelat+outofrangelon;
    outofrange(outofrange==2)=1;
    outofrange=logical(outofrange);

    latemp=Dlat{ix};
    lotemp=Dlon{ix};
    scidatetemp=Dscidates{ix};
    chltemp=Dchl{ix};
    deptemp=Ddepth{ix};
    saltemp=Dsal{ix};
    temptemp=Dtemp{ix};
    partemp=Dpar{ix};
    bbtemp=Dbb{ix};
    dotemp=Ddo_c{ix};
    cdomtemp=Dcdom{ix};

    latemp(outofrangelat)=nan;
    lotemp(outofrangelat)=nan;
    scidatetemp(outofrangelat)=nan;
    chltemp(outofrange)=nan;
    deptemp(outofrange)=nan;
    saltemp(outofrange)=nan;
    temptemp(outofrange)=nan;
    bbtemp(outofrange)=nan;

    if ~isempty(partemp)
        partemp(outofrange)=nan;
    end
    if ~isempty(dotemp)
        dotemp(outofrange)=nan;
    end
    if ~isempty(cdomtemp)
        cdomtemp(outofrange)=nan;
    end

    datemaxes{ix}=datetime(Ddatemax{ix},'ConvertFrom','datenum');
    datemins{ix}=datetime(Ddatemin{ix},'ConvertFrom','datenum');

    % Format the datetime to character strings
    dateRange{ix} = strcat(datestr(datemins{ix}, 'mm/dd'), ' - ', datestr(datemaxes{ix}, ' mm/dd/yyyy'));

    lats{ix}= latemp;
    lons{ix}= lotemp;
    mydatenums{ix}= scidatetemp;
    chl{ix}= chltemp;
    dep{ix}= deptemp;
    sal{ix}= saltemp;
    temp{ix}= temptemp;
    par{ix}= partemp;
    scatter{ix} = bbtemp;
    do_c{ix} = dotemp;
    cdom{ix} = cdomtemp;
    
 
end
%sci variables.
AllGliderVariables{1}=lats;
AllGliderVariables{2}=lons;
AllGliderVariables{3}=mydatenums;
AllGliderVariables{4}=chl;
AllGliderVariables{5}=dep;
AllGliderVariables{6}=sal;
AllGliderVariables{7}=temp;
AllGliderVariables{8}=par;
AllGliderVariables{9}=scatter;
AllGliderVariables{10}=do_c;
AllGliderVariables{11}=cdom;
%flight variables.
AllGliderFliVariables{1}=Dflidates;
AllGliderFliVariables{2}=Dflidepth;
AllGliderFliVariables{3}=Dfliinflections;
save([myoutput,'AllgliderDataV2.mat'],'AllGliderVariables')
save([myoutput,'AllgliderFli_DataV2.mat'],'AllGliderFliVariables')



