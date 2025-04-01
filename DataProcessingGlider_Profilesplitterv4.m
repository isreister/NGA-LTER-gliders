% THis code takes in a cell array of glider data missions 'DataProcessingGlider_AllGliderVariablesMakerV2.m
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


%Thats the input. Now what does this code do? It separetes the glider data
%into upcasts and downcasts when applicable (Some missions have both,
%others have none, it all depends on what the mission was). We want to
%standardize all of this. Assumptions:
%1) Surface turning points occur in the upper 10 m of the water column
%2) Turning points are seperated by more than 100 indices. Hopefully this
%prevents false surface turning points due to data drops.
%) Once a pair of surface turning points are identified, there will be one deep turning point between
%them.

%Those are the main asumptions. Comments in the code pretain to special
%cases.

%Author: Isaac Reister 4/1/2025

%Record:   
%V4 Cleaning up the comments for release to UAF glider group.

%V3 2-25-2025 For some reason, mission 8 (the summer stationkeeping,
%wouldn't process with the AllgliderDataV2. This was fixed when I resaved
%that data file. Bug fixed.

%V3 V2 was not effective at obatinaing peaks

%V2 Since we determined that the transecting gliders can't really be
%corrected (there was no hope of getting workable regression out of those)
%This code works for everything


%Notes: 
% 1 has a strange up down sequence. Its like every the sci depth was not
% recroded in a pattern.

% 6 sometimes spends extra time on surface and generates a mini cast.

% 8 sometimes spends extra time on surface and generates a mini cast.

% Step one in correcting the glider data so I can use it. 

% In an earlier attempt at identifying turning points I used the ready made inflection
% points. That didn't work, because infelction points are a for a lower
% resolution dataset.

clear all
close all

load('C:\Users\funkb\Documents\MATLAB\Research\data\Chapter3\AllgliderDataV2.mat','AllGliderVariables')
load('C:\Users\funkb\Documents\MATLAB\Research\data\Chapter3\AllgliderFli_DataV2.mat','AllGliderFliVariables')

%% user inputs
missionset=[1:8];

%% Setting up our time limits
%The limits below partition the glider missions into transecting and
%stationkeeping. Some of the transecting date limmits include breif (3 days
%or less) periods of stationkeeping, which may be useful for verification
%of data at another point. The main intention here is to divide up our
%glider dataset into workable chunks. The station keeping periods include
%no transecting (so for example track out to GEO and track back from GEO is
%is not part of these date limits), though some stationkeeping gliders
%(particuarly mission 7 got blown around a lot and will need some further
%editing using a radius from the GEO mooring. 

Timelimits{1}=[datenum(datetime(2021,01,20,00,00,00)),datenum(datetime(2021,05,05,00,00,00))]; %PWS stationkeeping
Timelimits{2}=[datenum(datetime(2022,02,14,00,00,00)),datenum(datetime(2022,03,02,00,00,00))]; %transecting
Timelimits{3}=[datenum(datetime(2022,06,04,00,00,00)),datenum(datetime(2022,07,15,00,00,00))]; %transecting (Yakatat)
Timelimits{4}=[datenum(datetime(2023,03,02,00,00,00)),datenum(datetime(2023,05,04,00,00,00))]; %stationkeeping
Timelimits{5}=[datenum(datetime(2023,03,24,00,00,00)),datenum(datetime(2023,04,14,00,00,00))]; %transecting
Timelimits{6}=[datenum(datetime(2024,03,23,00,00,00)),datenum(datetime(2024,04,05,00,00,00))]; %transecting
Timelimits{7}=[datenum(datetime(2024,03,20,00,00,00)),datenum(datetime(2024,06,18,00,00,00))]; %stationkeeping
Timelimits{8}=[datenum(datetime(2024,06,20,20,00,00)),datenum(datetime(2024,09,21,20,00,00))]; %stationkeeping


mytitles{1}=[char(datetime(2021,01,20,00,00,00)),' ',char(datetime(2021,05,05,00,00,00)),':PWS Station Keeping']; %PWS stationkeeping
mytitles{2}=[char(datetime(2022,02,14,00,00,00)),' ',char(datetime(2022,03,02,00,00,00)),':Seward Transect']; %transecting
mytitles{3}=[char(datetime(2022,06,04,00,00,00)),' ',char(datetime(2022,07,15,00,00,00)),':Yakatat Transect']; %transecting (Yakatat)
mytitles{4}=[char(datetime(2023,03,02,00,00,00)),' ',char(datetime(2023,05,04,00,00,00)),'GEO Station Keeping']; %stationkeeping
mytitles{5}=[char(datetime(2023,03,24,00,00,00)),' ',char(datetime(2023,04,14,00,00,00)),'Seward Transect']; %transecting
mytitles{6}=[char(datetime(2024,03,23,00,00,00)),' ',char(datetime(2024,04,05,00,00,00)),'Seward Transect']; %transecting
mytitles{7}=[char(datetime(2024,03,20,00,00,00)),' ',char(datetime(2024,06,18,00,00,00)),'GEO Station Keeping']; %stationkeeping
mytitles{8}=[char(datetime(2024,06,20,20,00,00)),' ',char(datetime(2024,09,21,20,00,00)),'GEO Station Keeping']; %stationkeeping

Stationkeepingset=logical([1,0,0,1,0,0,1,1]);

%%

%%

T=AllGliderVariables{3}; %times for each mission.

% Apply mission limits to all variables
for eachmission=1:length(T)
outsideoftimeframe{eachmission}=T{eachmission}<Timelimits{eachmission}(1) | T{eachmission}>Timelimits{eachmission}(2);
end


for eachvar=1:length(AllGliderVariables)

    for eachmission=1:length(T)
        if ~isempty(AllGliderVariables{eachvar}{eachmission})
            AllGliderVariables{eachvar}{eachmission}(outsideoftimeframe{eachmission})=[];
        end
    end
  
end

T=AllGliderVariables{3}; %redefine T now that limits have been applied.


%% Turing point identifier.


scidepths=AllGliderVariables{5};
scitimes=AllGliderVariables{3};

for eachmission=1:length(T)


    
    missionalldepths=scidepths{eachmission}; %for this mission
    alldepthsindex=1:length(missionalldepths);

    missionalltimes=scitimes{eachmission};

    nonan_alldepthsindex=alldepthsindex(~isnan(missionalldepths)); %removenans
    nonan_alldepths=missionalldepths(~isnan(missionalldepths));

    %remove depth values greater than 10 for plotting.
    nonan_surfacedepths=nonan_alldepths;
    nonan_surfacedepths(nonan_alldepths>10)=[];
    nonan_sufacedepthsindex=nonan_alldepthsindex;
    nonan_sufacedepthsindex(nonan_alldepths>10)=[]; %so this is shorter than, nonan_alldepthsindex, but the indexes match between the two.

    figure;
    plot(nonan_alldepthsindex,nonan_alldepths,'Marker','.','LineStyle','none')
    hold on
    plot(nonan_sufacedepthsindex,nonan_surfacedepths)

    %split up surface sections
    diffs=diff(nonan_sufacedepthsindex); %consecutive numbers are grouped.
 
    % Identify the start of a new group (where the difference is greater than 1)
    group_start = [1, find(diffs > 100) + 1]; %a gap of 100 in the indicies means data drops are allowed within the surface, but the gap between different casts is much greater than 100.
    % this find is in terms of nonan_sufacedepthsindex
    group_end = [find(diffs > 100), length(nonan_sufacedepthsindex)];
   %  testing=1:length(nonan_sufacedepthsindex);
   % figure;plot(testing,nonan_sufacedepthsindex)
   %  hold on
   %  plot(testing(group_start),nonan_sufacedepthsindex(group_start),'o')
   %  plot(testing(group_end),nonan_sufacedepthsindex(group_start),'+')
    %testing
    
    %test the position
    % figure(1)
    % hold on
    % plot(nonan_sufacedepthsindex(group_start(1):group_end(1)),nonan_surfacedepths(group_start(1):group_end(1)),'o')
    % plot(nonan_sufacedepthsindex(group_start(2):group_end(2)),nonan_surfacedepths(group_start(2):group_end(2)),'o')
    
    % Group the numbers
    groups = arrayfun(@(s, e) alldepthsindex(s:e), group_start, group_end, 'UniformOutput', false);
    
    % figure;
    % 
    % plot(nonan_alldepthsindex,nonan_alldepths,'.')
    % hold on
    % for ix=1:100
    %     plot(nonan_sufacedepthsindex(groups{ix}),nonan_surfacedepths(groups{ix}),'o')
    %     hold on
    % end


% This does a good job of delineating surfacing periods. Note that surfacings can be rather jagged. Sometimes with a rise to say 2m then descend to 10 m then rise again to 1 m.
% now find the minimum of each group.

    for eachgroup=1:length(groups)
        [surfacevalue(eachgroup),thisgroupsurfaceI]=min(nonan_surfacedepths(groups{eachgroup}));
        tempindex=nonan_sufacedepthsindex(groups{eachgroup});
        surfaceI(eachgroup)=tempindex(thisgroupsurfaceI);
    end


%now find the deepest depth between each surfacing.
    for eachsurfacepair=1:length(surfaceI)-1;
    
        [deepvalue(eachsurfacepair),thispairdeepI] = max(missionalldepths(surfaceI(eachsurfacepair):surfaceI(eachsurfacepair+1)));
        deepI(eachsurfacepair)=surfaceI(eachsurfacepair)+thispairdeepI-1;
        % tempindexdeep=alldepthsindex(surfaceI(eachsurfacepair):surfaceI(eachsurfacepair+1));
        % tempdepths=missionalldepths(surfaceI(eachsurfacepair):surfaceI(eachsurfacepair+1));
        % deepI(eachsurfacepair)= tempindex(thispairdeepI);
    end

%If there is a long pause between two surface points, we assume the glider
%was at the surface for an extended period of time. This can result in an
%the deepvalue equalling the deepest surface point. Clearly this is in error since the deep point in these cases is above 10.
% In these cases
% There is always going to be one less deep point than surface point with
% this method.

%The index of the deep point is always going to lead the surfac point with
%this method. Find deeppoints less than 10


    deepI(deepvalue<10)=nan;
    fakedeep=isnan(deepI);
     %remove the fake deep values
    deepvalue(isnan(deepI))=[];
    deepI(isnan(deepI))=[];

    if sum(fakedeep)>0 %if any fake deeps exist.
        for eachdeeppair=1:length(deepI)-1
            possiblybadsurface=find((surfaceI>=deepI(eachdeeppair)) & (surfaceI<=deepI(eachdeeppair+1))); %find all surfaceI inbetween deep pairs
            if length(possiblybadsurface)>1
                %select shallowest surface and call it a turing point. Even if
                %the Glider hangs out at the surface for several hours, that
                %woks, we can fine tune the profiles later to remove excess
                %surface time.
                [~,bestsurface]=min(surfacevalue(possiblybadsurface));
                
                good=logical(zeros(1,length(possiblybadsurface)));
                good(bestsurface)=logical(1);
                bad=~good;
    
                badI=possiblybadsurface(bad);
                surfaceI(badI)=nan;
                surfacevalue(badI)=nan;
    
            end
        end
    end

surfacevalue(isnan(surfaceI))=[];
surfaceI(isnan(surfaceI))=[];

downcasts{eachmission}=[surfaceI(1:end-1).',deepI.'];
upcasts{eachmission}=[deepI.',surfaceI(2:end).'];

%


%sanity check plot. depths less than 10 are the orange line. Surface turning points are
%idenfied as yellow dots, and purples points are deep turning
%points. Note that this code does not set a required 'full depth profile'
%because often missions have staggered depths (200 meters for 3 dives, 40
%meters for another). We don't want to delete a well resolved surface
%layer. Also in a shallow areas like PWS that method wont work at all.
figure(1)
hold on
plot(surfaceI,surfacevalue,'Marker','o','LineStyle','none','MarkerFaceColor',[ 0.9294    0.6941    0.1255]) %this looks really good.
plot(deepI,deepvalue,'o') %this looks really good.
set(gca,'YDir','reverse')
ylabel('Depth (m)')
xlabel('Casts')
set(gcf,"Position",[25.6667 193 1.6373e+03 420])
saveas(figure(1),['C:\Users\funkb\Documents\MATLAB\Research\Figures\Chapter 3\prelim\Mission ', char(string(eachmission)),'.png'])
hold off
close all
clear deepI
clear deepvalue
clear surfaceI
clear surfacevalue
end



save('C:\Users\funkb\Documents\MATLAB\Research\data\Chapter3\DataProcessingGlider_AllGliderProfileSplit_v3.mat',"downcasts","upcasts","AllGliderVariables");

%the start and stop of ascending and descending profiles has now been
%obtained.




