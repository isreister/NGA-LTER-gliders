%This code implements more quality control measures, namely it removes
% any vertical bins that have salininty measurments less than 1. it also
% calculates practical salinity and depth using gibbs seawater toolkit.

%Once all of these bins are replaced with nans, the code also
%diagnoses if a cast contains data for at least 75% of it's full cast
%depth. This is an effective way to remove false casts

%Metrics like how many QC'd profiles are available and such are then printed out.

%Data is then gridded to 3 hours.

% Preliminary plots are created.

%Author Isaac Reister 4/8/2025


%Notes-----------------

%v4 add more comments. Adjusted for new variables for QC_GC. Added a better
%display of the raw metrics. Uploaded to github. This produces a usable
%QC'd fully gridded glider data product. Future changes will be documented
%in Github.

%v3 added some more comments and added some metrics for pre-interpolated
%data.

%2-25-2025 running v3 which has the corrected depths for mission7.
close all
clear all

%% User config
load('C:\Users\funkb\Documents\MATLAB\Research\data\Chapter3\QC_GC_1dayhalf_v9.mat')
%load('C:\Users\funkb\Documents\MATLAB\Research\data\Chapter3\QC_GC_v7.mat','QC_GC') %load wherever the initial rolling window QC was saved.

%variable organization as loaded:
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

placesnames{1}='PWS';
placesnames{2}='GEO';
placesnames{3}='GEO';
placesnames{4}='GEO';



filepacenames{1}='PWS2021_spring';
filepacenames{2}='GEO2023_spring';
filepacenames{3}='GEO2024_spring';
filepacenames{4}='GEO2024_summer';



QCdatasavepath='C:\Users\funkb\Documents\MATLAB\Research\data\Chapter3\';
datafilename = 'QC_GCgridded_v4testing';

Prelimplotfilesavepath= 'C:\Users\funkb\Documents\MATLAB\Research\Figures\Chapter 3\Gridded_no_interp yet\';


%% housekeeping


    
    % housekeeping, to make code work with downstream code.
    % calculate practical salinity, put in place of conductivityu, move conductivity to the end of the data set.
    % calculate depth, put in place of pressure. move pressure to end of dataset.
      for eachmission=1: length(QC_GC)
         currentSP = gsw_SP_from_C(QC_GC{eachmission}(:,:,6),QC_GC{eachmission}(:,:,7),QC_GC{eachmission}(:,:,5));
         currentD = - gsw_z_from_p(QC_GC{eachmission}(:,:,5),60);
    
         QC_GC{eachmission}(:,:,14)=QC_GC{eachmission}(:,:,5); %pressure
         QC_GC{eachmission}(:,:,15)=QC_GC{eachmission}(:,:,6); %conductivity
         QC_GC{eachmission}(:,:,5)=currentD;
         QC_GC{eachmission}(:,:,6)=currentSP;
      end
    
    %use salinity to clear our any data that has bad surface values. Looks like
    %our new QC procedure is working well, there are no bad salinities for the
    %datasets I worked with.
        for eachmission=1: length(QC_GC)
            
                %use salinity as a sanity check
                mysalinity=QC_GC{eachmission}(:,:,6);
        
                lessthan1= mysalinity<1;
                display(['Number of vertical bins with salinities that are <1 in ' ,placesnames{eachmission}, ' Mission ',char(string(eachmission)),' : ',char(string(sum(lessthan1,'all')))]) 
            for eachvar=1:size(QC_GC{1},3)
                temp=QC_GC{eachmission}(:,:,eachvar);
                temp(lessthan1)=nan;
                QC_GC{eachmission}(:,:,eachvar)=temp;
        
            end
        end
    
    %% some casts have very little data. Make a requirement that a cast must add up to at least 75% of the full depth to qualify as a cast.
    % some data is more consistent than others, but if salinity doesn't make up
    % a full cast, we assume the rest do not either. The issue crops up a lot (40-50%) in Mission 3 and Mission 4 . 
    %The main cause for this is due to sample scheme I think. Mission 4 and
    %Mission 3 collected science data during the upcast only. The way the
    %splitter works means a little bit of the data is registered after deep
    %turning point, and a little bit is registered before a surface turning
    %point. These then get temporarily identified as a 'a profile', which is
    %corrected here.
    %
    % Another issue corrected here: Mission four has some casts that only have
    % nans for depths. Not entirely clear why those exist, but clearly they
    % need to be removed from the data set.
        for eachmission=1: length(QC_GC)
            
                
                mydepth=max(QC_GC{eachmission}(:,:,5),[],'omitnan');
                tempindex=ceil(mydepth);
                
                depthcheck=isnan(tempindex);
                if sum(depthcheck)>0
                    mydepth(depthcheck)=[];
                    tempindex(depthcheck)=[];
                    QC_GC{eachmission}(:,depthcheck,:)=[];
                end
                
                nansalinity=isnan(QC_GC{eachmission}(:,:,6));
    
                for eachprofile=1:size(nansalinity,2)
                   whats_spotty(eachprofile)=sum(nansalinity(1:tempindex(eachprofile),eachprofile))>(.75*tempindex(eachprofile));  
                end
                percentbad=(sum(whats_spotty)/size(QC_GC{eachmission},2))*100;
                display([char(string(percentbad)),' % of ' ,placesnames{eachmission}, ' Mission ',char(string(eachmission)),' were identified as vertical profiles with large gaps that were too large. These were removed.'])
                QC_GC{eachmission}(:,whats_spotty,:)=[]; %note whole cast is removed.
     
            clear whats_spotty
        end
    
    
    
    
    
    
    
    %% here we get some "raw" metrics for the quality controlled data.
    for eachmission=1:length(QC_GC)
        Totalcasts=length(QC_GC{eachmission});
        
        %Centerpoint
        lat_temp=QC_GC{eachmission}(1,:,1);
        lon_temp=QC_GC{eachmission}(1,:,2);
        
        nolat=isnan(lat_temp);
        nolon=isnan(lon_temp);
        nogps=nolat & nolon; %note these are not bad casts, just casts wherein the glider didn't surface.
        lat_temp(nogps)=[];
        lon_temp(nogps)=[];
        
        eachcenterpointlat{eachmission}=mean(lat_temp);
        eachcenterpointlon{eachmission}=mean(lon_temp);
    
        distances = distance(eachcenterpointlat{eachmission}, eachcenterpointlon{eachmission}, lat_temp, lon_temp, wgs84Ellipsoid('km'));
        mean(distances);
        std(distances);
        display(['total number of profiles for ' ,placesnames{eachmission}, ' mission ', char(string(eachmission)) ,' was ', char(string(Totalcasts))])
    
    
        display(['lat lon centerpoint for' ,placesnames{eachmission}, ' mission ', char(string(eachmission)) ,' was ', char(string(round(eachcenterpointlat{eachmission},3))),'  , ',char(string(round(eachcenterpointlon{eachmission},3)))])
    
        distances = distance(eachcenterpointlat{eachmission}, eachcenterpointlon{eachmission}, lat_temp, lon_temp, wgs84Ellipsoid('km'));
        mean(distances);
        std(distances);
        display(['mean and std from that centerpoint for' ,placesnames{eachmission}, ' mission ', char(string(eachmission)) ,' was ', char(string(round(mean(distances),1))),'  , ',char(string(round(std(distances),1)))])
    
    
    
    end
    
    %all geo together:
    
        
        %Centerpoint
        lat_temp=[QC_GC{2}(1,:,1),QC_GC{3}(1,:,1),QC_GC{4}(1,:,1)];
        lon_temp=[QC_GC{2}(1,:,2),QC_GC{3}(1,:,2),QC_GC{4}(1,:,2)];
        Totalcasts=length(lat_temp);
        nolat=isnan(lat_temp);
        nolon=isnan(lon_temp);
        nogps=nolat & nolon; %note these are not bad casts, just casts wherein the glider didn't surface.
        lat_temp(nogps)=[];
        lon_temp(nogps)=[];
        
        centerpointlat=mean(lat_temp);
        centerpointlon=mean(lon_temp);
    
        display(['lat lon centerpoint for all' ,placesnames{2}, 'deployments was ', char(string(round(centerpointlat,3))),'   , ',char(string(round(centerpointlon,3)))])
    
        distances = distance(centerpointlat, centerpointlon, lat_temp, lon_temp, wgs84Ellipsoid('km'));
        mean(distances);
        std(distances);
        display(['mean and std from that centerpoint for all' ,placesnames{2}, 'deployments was ', char(string(round(mean(distances),1))),'  , ',char(string(round(std(distances),1)))])
    
    
        
        timestamps=[QC_GC{1}(1,:,3)];
        timestamps2=[QC_GC{2}(1,:,3),QC_GC{3}(1,:,3),QC_GC{4}(1,:,3)];
    
        difftimes=diff((convrtdt(timestamps)));
        difftimes2=diff((convrtdt(timestamps2)));
    
        difftimes(difftimes>hours(24))=[];
        difftimes2(difftimes2>hours(24))=[];
    
        PWStemporoalresolution=mean(difftimes,'omitnan');
        PWSstdtimersol=std(difftimes,'omitnan');
    
       display(['mean and std temporal resolution for ' ,placesnames{1}, ' deployment was ', char(string(PWStemporoalresolution)),' ',char(string(PWSstdtimersol))])
    
    
    
        GEOtemporoalresolution=mean(difftimes2,'omitnan');
        GEOstdtimersol=std(difftimes2,'omitnan');
    
       display(['mean and std temporal resolution for ' ,placesnames{1}, ' deployment was ', char(string(GEOtemporoalresolution)),' ',char(string(GEOstdtimersol))])
    
        depthmax=max([QC_GC{1}(:,:,5)]);
        depthmax2=max([QC_GC{2}(:,:,5),QC_GC{3}(:,:,5),QC_GC{4}(:,:,5)]);
        
        remove=difftimes>hours(24);
        remove2=difftimes>hours(24);
    
        PWSmeandepth=mean(depthmax(~remove),"omitmissing");
        PWSstddepth=std(depthmax(~remove),"omitmissing");
    
        display(['mean and std depth for ' ,placesnames{1}, ' deployment was ', char(string(round(PWSmeandepth,1))),' ',char(string(round(PWSstddepth,1)))])
    
    
        GEOmeandepth=mean(depthmax2(~remove),"omitmissing");
        GEOstdepth=std(depthmax2(~remove),"omitmissing");
        display(['mean and std depth for ' ,placesnames{2}, ' deployment was ', char(string(round(GEOmeandepth,1))),' ',char(string(round(GEOstdepth,1)))])
    
        % why does PWS have such temporally short dives??? Answer: The data looks
        % OK so it must have been a much steeper dive angle. Looking at the
        % data, the longer dives down to 120 meters took about 30 minutes to
        % complete. A profile is half of a dive, so that 15 minutes. Add in the
        % occasional short dive and you get about 15 minutes for a dive.
    
        % figure;
        % 
        % PWStrtime=QC_GC{1}(:,:,3);
        % PWStrdepth=QC_GC{1}(:,:,5);
        % PWStrtime=PWStrtime(:);
        % PWStrdepth=PWStrdepth(:);
        % scatter(convrtdt(PWStrtime),PWStrdepth,'o')
    
    
        
        
        
        
        
        
        %% grid through time by 3 hours.
     for eachmission=1: length(QC_GC)
        mytimes=QC_GC{eachmission}(:,:,3);
        mindate=char(datetime(min(QC_GC{eachmission}(:,:,3),[],'all'),'ConvertFrom','datenum','Format','uuuu-MM-dd'));
        maxdate=char(datetime(max(QC_GC{eachmission}(:,:,3),[],'all'),'ConvertFrom','datenum','Format','uuuu-MM-dd'));
        pre_mytimebin=datetime(mindate):hours(3):datetime(maxdate)+hours(24);
        pre_mytimebin(pre_mytimebin<datetime(min(QC_GC{eachmission}(:,:,3),[],'all'),'ConvertFrom','datenum'))=[];
        pre_mytimebin(pre_mytimebin>datetime(max(QC_GC{eachmission}(:,:,3),[],'all'),'ConvertFrom','datenum'))=[];

        bin_edges{eachmission}=pre_mytimebin;
        tempgrid=nan(300,length(bin_edges{eachmission})-1);



        for eachvar=1:size(QC_GC{eachmission},3)

                for eachdepth=1:300

                        timesforthisdepth=QC_GC{eachmission}(eachdepth,:,3);
                        [~, ~, Bin_I] = histcounts(timesforthisdepth, datenum(bin_edges{eachmission}));

                        for each3hourchunk=1:length(bin_edges{eachmission})-1
                            thistimeset=Bin_I==each3hourchunk;
                            thisvar=QC_GC{eachmission}(eachdepth,:,eachvar);

                            tempgrid(eachdepth,each3hourchunk,eachvar)=mean(thisvar(thistimeset),'omitnan');
                        end



                end
        end

        QC_GCgridded{eachmission}=tempgrid;

     end   
     save([QCdatasavepath,datafilename,'.mat'],"QC_GCgridded",'bin_edges','eachcenterpointlat','eachcenterpointlon')  

%load([QCdatasavepath,datafilename,'.mat'],"QC_GCgridded",'bin_edges','eachcenterpointlat','eachcenterpointlon')

%% Making some preliminary plots.
%Current variable order
varsnames{1}='Latitude';
varsnames{2}='Longitude';
varsnames{3}='Time';
varsnames{4}=['Chlorophyll-a (',char(181),'g L^{-1})'];
varsnames{5}='Depth (m)';
varsnames{6}='Salinity';
varsnames{7}=['Temperature (',char(176),'C)'];
varsnames{8}='PAR';
varsnames{9}='Scatter (sr m^{-1})';
varsnames{10}=['Dissolved Oxygen (',char(181),'mol L^{-1})'];
varsnames{11}='CDOM (ppb)';
varsnames{12}='U (m/s)';
varsnames{13}='V (m/s)';
varsnames{14}='Pressure (dbar)';
varsnames{15}='Conductivity (mS/cm)';

filenamevarsnames{1}='Latitude';
filenamevarsnames{2}='Longitude';
filenamevarsnames{3}='Time';
filenamevarsnames{4}='Chlorophyll';
filenamevarsnames{5}='Depth';
filenamevarsnames{6}='Salinity';
filenamevarsnames{7}='Temperature';
filenamevarsnames{8}='PAR';
filenamevarsnames{9}='Scatter';
filenamevarsnames{10}='Dissolved Oxygen';
filenamevarsnames{11}='CDOM';
filenamevarsnames{12}='currentsU';
filenamevarsnames{13}='currentsV';
filenamevarsnames{14}='Pressure';
filenamevarsnames{15}='Conductivity';



for eachmission=1:length(QC_GCgridded)
    for eachvar=[4,6,7,8,9,10,11,12,13,14]
        [gridtimes, griddepths]=meshgrid((bin_edges{eachmission}(1:end-1)),0:299);

        
        if sum(isnan(QC_GCgridded{eachmission}(:,:,eachvar)),"all")~=size(QC_GCgridded{eachmission},1)*size(QC_GCgridded{eachmission},2)
            
            pcolor(gridtimes,griddepths,QC_GCgridded{eachmission}(:,:,eachvar));shading flat

            cb=colorbar;

            set(gca,'YDir','reverse');
            ylabel('Depth(m)')
            ylim([1,250]);
            ylabel(cb,varsnames{eachvar})
            mindate=char(datetime(min(QC_GCgridded{eachmission}(:,:,3),[],'all'),'ConvertFrom','datenum','Format','uuuu-MM-dd'));
            maxdate=char(datetime(max(QC_GCgridded{eachmission}(:,:,3),[],'all'),'ConvertFrom','datenum','Format','uuuu-MM-dd'));
            title([mindate,':',maxdate,':',placesnames{eachmission}])


            saveas(gcf,[Prelimplotfilesavepath,filepacenames{eachmission},filenamevarsnames{eachvar},'.png'])

        end
    end

end


   