%GliderData Processing
%This prepares high resolution data for processing by standardizing the glider names.
% Removes extremely bad data points (like data points over land or far far away).
% Produces some VERY basic positional plots. And saves out the slightly
% cleaned up data as flight and science variables. Here are the notes on
% the save out files:


% Here is the cell organization:
%Level 1 (top level)
%sci variables. Some of these variables are grabbed now for troubleshooting
%but dropped later on.
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
    % 
    % 
    % %flight variables.
    % AllGliderFliVariables{1}=Dflidates;
    % AllGliderFliVariables{2}=Dfliu;
    % AllGliderFliVariables{3}=Dfliv;

% Level 2
    %Each variable is a cell array wherein the a number of cells that correspond to
    %however many missions were processed in 'DataProcessingGlider_AllGliderVariablesMakerV2.m
% Level 3
    %Each of those cells then contains a column that is the actual data.
    %So for example if I wanted to see all the latitudes for the first
    %mission I would execute: AllGliderVariables{1}{1} 
    %this will produce a very long column vector.

% Author: Isaac Reister 3/31/2025
%V2for_github -- made the code extract glider currents in addition to the other variables.
    % also moved the sanity check plot to the end and improved the
    % limiting of the data a little bit. Remaining changes are documented
    % in github.
%V2 -- cleaned up comments for UAF glider group
%V1 -- intial creation
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
%191 has PAR, 507 does not have PAR. 1152 does not have PAR
%191 doesn't have CDOM, 507 does have CDOM. 1152 does have CDOM
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
   
    %Enter here where your desired spatial limits for the glider.
    mymaxlat=90;
    mymaxlon=-130;
    myminlat=50;
    myminlon=-150;

    % use this to make a very simple plot of the latitudes.
   sanitycheckplot=0; %1 yes make this plot. 0 no.

    %IMPORTANT***********************************
    %data names are inconsistent between gliders, and like all things in life, there is no easy code fix.
    %but we can use code to make our lives easier. Before that happens
    %though TWO rules must apply to your dataset:

    %1) The variables you are interested in MUST have units, and those
    %units MUST contain the text the specified in Standarnames_units{:,2}

    %2) For the variables you are interested in: The fieldnames listed
    %under sci.units MUST match the fieldnames listed under sci. Also make
    %sure that given your list of variables you want to extract, make sure
    % you don't have any false positives in under sci.units. E.g. mission 8
    % has par listed under sci.units, BUT that glider (1152) does not have a
    % PAR sensor, or at least AOOS isn't reporting it if it does.
    
apply_handediting=1; %this controls manual changes to the data itself.
Ineedtomakemanualchanges=0; %this controls manual changes to fieldnames
    %Manual fixes to the current datset to make these things true.
if Ineedtomakemanualchanges==1
    load([myroot,datarepo{8}]); %loads up our sci structure.
  
             sci.press=sci.pres; %mismatch name.
             sci=rmfield(sci,'pres');
             sci.units.cdom = 'CDOM (ppb)'; %this glider has CDOM
             sci.units=rmfield(sci.units,'par');
    save([myroot,datarepo{8},'IRfix']); %saves corrected sci structure.       
    

    load([myroot,datarepo{7}]); %loads up our sci structure.
             sci.press=sci.pres; %mismatch name.
             sci=rmfield(sci,'pres');
    save([myroot,datarepo{7},'IRfix']); %saves corrected sci structure.

    load([myroot,datarepo{6}]); %loads up our sci structure.
             sci.press=sci.pres; %mismatch name.
             sci=rmfield(sci,'pres');
    save([myroot,datarepo{6},'IRfix']); %saves corrected sci structure.

    load([myroot,datarepo{5}]); %loads up our sci structure.
             sci.press=sci.pres; %mismatch name.
             sci=rmfield(sci,'pres');
    save([myroot,datarepo{5},'IRfix']); %saves corrected sci structure.

    load([myroot,datarepo{4}]); %loads up our sci structure.
             sci.press=sci.pres; %mismatch name.
             sci=rmfield(sci,'pres');
    save([myroot,datarepo{4},'IRfix']); %saves corrected sci structure.

    load([myroot,datarepo{3}]); %loads up our sci structure.
             sci.press=sci.pres; %mismatch name.
             sci=rmfield(sci,'pres');
    save([myroot,datarepo{3},'IRfix']); %saves corrected sci structure.

    load([myroot,datarepo{2}]); %loads up our sci structure.
             sci.press=sci.pres; %mismatch name.
             sci=rmfield(sci,'pres');
    save([myroot,datarepo{2},'IRfix']); %saves corrected sci structure.

    load([myroot,datarepo{1}]); %loads up our sci structure.
             sci.psal=sci.sal; %mismatch name.
             sci=rmfield(sci,'sal');
             sci.units.psal = 'Salinity (unitless)'; %mismatch name
             sci.units=rmfield(sci.units,'sal');

    save([myroot,datarepo{1},'IRfix']); %saves corrected sci structure.
end

    %now we resave all of the datarepo cells so that they will upload the
    %fixed data.
    datarepo2{1}= [datarepo{1},'IRfix'];
    datarepo2{2}= [datarepo{2},'IRfix'];
    datarepo2{3}= [datarepo{3},'IRfix'];
    datarepo2{4}= [datarepo{4},'IRfix'];
    datarepo2{5}= [datarepo{5},'IRfix'];
    datarepo2{6}= [datarepo{6},'IRfix'];
    datarepo2{7}= [datarepo{7},'IRfix'];
    datarepo2{8}= [datarepo{8},'IRfix'];
    

    %% Load data, Standardize variable names, using the sci units which are assumed to be consistent, if they aren't this code will break so you will know you need to fix the units.
    
    
    
    % units we are looking for: 
    Standardnames_units=cell(12,2);
    Standardnames_units{1,1}='time';Standardnames_units{1,2}='matlab datenum format';
    Standardnames_units{2,1}='press';Standardnames_units{2,2}='Pressure (bar)';
    Standardnames_units{3,1}='cond';Standardnames_units{3,2}='Conductivity (milliSiemens per centimeter';
    Standardnames_units{4,1}='temp';Standardnames_units{4,2}='Temperature (degrees Centigrade)';
    Standardnames_units{5,1}='chl';Standardnames_units{5,2}='Chlorophyll-a Concentration (micrograms per liter)';
    Standardnames_units{6,1}='cdom';Standardnames_units{6,2}='CDOM (ppb)';
    Standardnames_units{7,1}='bb';Standardnames_units{7,2}='Backscatter (per steridian per meter)';
    Standardnames_units{8,1}='psal';Standardnames_units{8,2}='Salinity (unitless)';
    Standardnames_units{9,1}='do_c';Standardnames_units{9,2}='Dissolved Oxygen Concentration (micromols per liter)';
    Standardnames_units{10,1}='lat';Standardnames_units{10,2}='Latitude (decimal degrees)';
    Standardnames_units{11,1}='lon';Standardnames_units{11,2}='Longitude (decimal degrees)';
    Standardnames_units{12,1}='par';Standardnames_units{12,2}='Photosynthetically Active Radiation (microEinstein per square meter per second)';


    %we need to search through each sci.units, find out which units match,
    %it if matches port that value over to our cell array for the unit, if
    %to doesn't match cell array is empty.
 

    for ix=1:length(datarepo2)
        load([myroot,datarepo2{ix}]); %loads up our sci structure.

        for eachvarIwant = 1: size(Standardnames_units,1);
        % String to search for
        searchStr = Standardnames_units{eachvarIwant,2};
        
        % Initialize a flag to track if the search string is found
        containsSubstring = false;
        
            % Loop through each field of sci.units
            fieldNames = fieldnames(sci.units);
            for i = 1:length(fieldNames)
                % Get the value of the current field
                fieldValue = sci.units.(fieldNames{i});
                
                % Check if the field value contains the search string
                if contains(fieldValue, searchStr)
                    containsSubstring = true;
                    unreliable_fieldname = fieldNames{i};
                    
                    break;  % Exit the loop once a match is found
                end
            end
        
            unreliable_fieldnames=fieldnames(sci);


            % Sort the data into seperate cells. For some we are more
            % concerned if the variable is not there than others
            if strcmp(Standardnames_units{eachvarIwant,1},'time')
                if true(containsSubstring) & (length(sci.(unreliable_fieldname))>1) %the length greater than one checks aganist a single nan entry for the field.

                    Dscidates{ix}=sci.(unreliable_fieldname);
                else
                    Dscidates{ix}=[];
                    display(['ERROR:No Time for mission ',char(string(ix))])
                    return
                end
            end
            if strcmp(Standardnames_units{eachvarIwant,1},'press')
                if true(containsSubstring) & (length(sci.(unreliable_fieldname))>1)
                    Dpres{ix}=sci.(unreliable_fieldname);
                    else
                    Dpres{ix}=[];
                    display(['ERROR:No Pressure for mission ',char(string(ix))])
                    return
                end
            end
            if strcmp(Standardnames_units{eachvarIwant,1}, 'cond')
                if true(containsSubstring) & (length(sci.(unreliable_fieldname))>1) 
                    Dcond{ix}=sci.(unreliable_fieldname);
                else
                    Dcond{ix}=[];
                    display(['ERROR:No Conductivity for mission ',char(string(ix))])
                    return
                end
            end

            if strcmp(Standardnames_units{eachvarIwant,1}, 'temp')
                if true(containsSubstring) & (length(sci.(unreliable_fieldname))>1)
                    Dtemp{ix}=sci.(unreliable_fieldname);
                else
                    Dtemp{ix}=[];
                    display(['ERROR:No Temperature for mission ',char(string(ix))])
                    return
                end
            end

            if strcmp(Standardnames_units{eachvarIwant,1}, 'chl')
                if true(containsSubstring) & (length(sci.(unreliable_fieldname))>1)
                    Dchl{ix}=sci.(unreliable_fieldname);
                else
                    Dchl{ix}=[];
                    display(['Warning:No Chorophyll-a for mission ',char(string(ix))])
                    
                end
            end
            if strcmp(Standardnames_units{eachvarIwant,1}, 'cdom')
                if true(containsSubstring) & (length(sci.(unreliable_fieldname))>1) 
                    Dcdom{ix}=sci.(unreliable_fieldname);
                else
                    Dcdom{ix}=[];
                    display(['Warning:No CDOM for mission ',char(string(ix))])
                    
                end
            end
            if strcmp(Standardnames_units{eachvarIwant,1}, 'bb')
                if true(containsSubstring) & (length(sci.(unreliable_fieldname))>1) 
                    Dbb{ix}=sci.(unreliable_fieldname);
                else
                    Dbb{ix}=[];
                    display(['Warning:No backscatter for mission ',char(string(ix))])
                    
                end
            end
            if strcmp(Standardnames_units{eachvarIwant,1}, 'psal')
                if true(containsSubstring) & (length(sci.(unreliable_fieldname))>1)
                    Dsal{ix}=sci.(unreliable_fieldname);
                else
                    Dsal{ix}=[];
                    display(['Warning:No salinity for mission ',char(string(ix))])
                    
                    
                end
            end
            if strcmp(Standardnames_units{eachvarIwant,1}, 'do_c')
                if true(containsSubstring) & (length(sci.(unreliable_fieldname))>1)
                    Ddo_c{ix}=sci.(unreliable_fieldname);
                else
                    Ddo_c{ix}=[];
                    display(['Warning:No dissolved oxygen concentration for mission ',char(string(ix))])
                    
                    
                end
            end

            if strcmp(Standardnames_units{eachvarIwant,1}, 'lat')
                if true(containsSubstring) & (length(sci.(unreliable_fieldname))>1)
                    Dlat{ix}=sci.(unreliable_fieldname);
                else
                    Dlat{ix}=[];
                    display(['ERROR:No latitude for mission ',char(string(ix))])
                    return
                    
                end
            end

            if strcmp(Standardnames_units{eachvarIwant,1}, 'lon')
                if true(containsSubstring) & (length(sci.(unreliable_fieldname))>1)
                    Dlon{ix}=sci.(unreliable_fieldname);
                else
                    Dlon{ix}=[];
                    display(['ERROR:No longitude for mission ',char(string(ix))])
                    return
                    
                end
            end

            if strcmp(Standardnames_units{eachvarIwant,1}, 'par')
                if true(containsSubstring) & (length(sci.(unreliable_fieldname))>1)
                    Dpar{ix}=sci.(unreliable_fieldname);
                else
                    Dpar{ix}=[];
                    display(['Warning:No Par for mission ',char(string(ix))])
                    
                    
                end
            end

        Dflidates{ix}=fli.gtime;
        Dflidepth{ix}=fli.depth;
        Dfliinflections{ix}=fli.inflecting;
        Dfliu{ix}=fli.u;
        Dfliv{ix}=fli.v;

        Ddatemax{ix}=max(Dscidates{ix});
        Ddatemin{ix}=min(Dscidates{ix});
        end
    end


%% Narrow down our box some more:

for ix=1:length(datarepo2)

    %remove non-NGA data


    datemaxes{ix}=datetime(Ddatemax{ix},'ConvertFrom','datenum');
    datemins{ix}=datetime(Ddatemin{ix},'ConvertFrom','datenum');

    % Format the datetime to character strings
    dateRange{ix} = strcat(datestr(datemins{ix}, 'mm/dd'), ' - ', datestr(datemaxes{ix}, ' mm/dd/yyyy'));

    outofrangelat=(Dlat{ix}>mymaxlat) | (Dlat{ix}<myminlat);
    outofrangelon=(Dlon{ix}>mymaxlon) | (Dlon{ix}<myminlon);

    display(['Processing file ' char(string(ix)),' ',dateRange{ix}])

    outofrange=outofrangelat+outofrangelon;
    outofrange(outofrange==2)=1;
    outofrange=logical(outofrange);

    latemp=Dlat{ix};
    lotemp=Dlon{ix};
    scidatetemp=Dscidates{ix};
    chltemp=Dchl{ix};
    saltemp=Dsal{ix};
    temptemp=Dtemp{ix};
    partemp=Dpar{ix};
    bbtemp=Dbb{ix};
    dotemp=Ddo_c{ix};
    cdomtemp=Dcdom{ix};
    condtemp=Dcond{ix};
    prestemp=Dpres{ix};

    latemp(outofrangelat)=nan;
    lotemp(outofrangelat)=nan;
    scidatetemp(outofrangelat)=nan;
    chltemp(outofrange)=nan;
    saltemp(outofrange)=nan;
    temptemp(outofrange)=nan;
    bbtemp(outofrange)=nan;
    condtemp(outofrange)=nan;
    prestemp(outofrange)=nan;

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
    sal{ix}= saltemp;
    temp{ix}= temptemp;
    par{ix}= partemp;
    scatter{ix} = bbtemp;
    do_c{ix} = dotemp;
    cdom{ix} = cdomtemp;
    cond{ix} = condtemp;
    pres{ix} = prestemp;
 
end

% sanity check plot of just latitudes.
if sanitycheckplot==1
    for ix=1:length(datarepo2)
    
        DdatemaxNEW{ix}=max(mydatenums{ix});
        DdateminNEW{ix}=min(mydatenums{ix});
    
        datemaxes{ix}=datetime(DdatemaxNEW{ix},'ConvertFrom','datenum');
        datemins{ix}=datetime(DdateminNEW{ix},'ConvertFrom','datenum');
    
        % Format the datetime to character strings
        dateRange{ix} = strcat(datestr(datemins{ix}, 'mm/dd'), ' - ', datestr(datemaxes{ix}, ' mm/dd/yyyy'));
        mysize=10;
        trancolor=missioncolors(ix,:);
        figure(1)
        hold on
        mhandle{ix}=plot(datetime(Dscidates{ix},'ConvertFrom','datenum'),lats{ix},'Color',trancolor,'LineWidth',2);
        hold off
        display(dateRange{ix})
     
    end
%ylim([58.5,60])
legend([mhandle{1:end}], dateRange{1:end})
end

%range table

CNDC=[0,60]; %s m^-1
TEMP=[-2.5, 40]; %C
PRES=[-5,1100]; %dbar
BBP=[0.0001,1]; %m^-1
PSAL=[2,41]; 
DOX1=[0,650]; %umol l^-1
CPHL=[0,100]; %mg m%-3
CDOM=[0,400]; 
IRRAD=[0,5000]; %(microEinstein per square meter per second)';
UCUR=[-10,10]; %m/s
VCUR=[-10,10]; %m/s

%%


% Manual corrections
%% Shackleton (glider 191) PAR before 2024 has a very small negative offset. To correct I simply add the minimum value for found in the 2021 (-0.9721) and 2023 (-3.3380). 
if apply_handediting==1
par{1}=abs(min(par{1}))+par{1};
par{4}=abs(min(par{4}))+par{4};
end
%basic range based QC to remove impossible values
for ix=1:length(datarepo2)


 cond{ix}(cond{ix}<CNDC(1) | cond{ix}>CNDC(2))=nan;
 temp{ix}(temp{ix}<CNDC(1) | temp{ix}>CNDC(2))=nan;
 pres{ix}(pres{ix}<PRES(1) | pres{ix}>PRES(2))=nan;
 scatter{ix}(scatter{ix}<BBP(1) | scatter{ix}>BBP(2))=nan;
 sal{ix}(sal{ix}<PSAL(1) | sal{ix}>PSAL(2))=nan;
 do_c{ix}(do_c{ix}<DOX1(1) | do_c{ix}>DOX1(2));
 chl{ix}(chl{ix}<CPHL(1) | chl{ix}>CPHL(2))=nan; 
 cdom{ix}(cdom{ix}<CDOM(1) | cdom{ix}>CDOM(2))=nan;
 par{ix}(par{ix}<IRRAD(1) | par{ix}>IRRAD(2))=nan;
 Dfliu{ix}(Dfliu{ix}<UCUR(1) | Dfliu{ix}>UCUR(2))=nan;
 Dfliv{ix}(Dfliv{ix}<VCUR(1) | Dfliv{ix}>VCUR(2))=nan;
 
 
end



%sci variables.(any more than this and MATLAB can't save it as a single
%file.
AllGliderVariables{1}=lats;
AllGliderVariables{2}=lons;
AllGliderVariables{3}=mydatenums;
AllGliderVariables{4}=chl;
AllGliderVariables{5}=pres;
AllGliderVariables{6}=cond;
AllGliderVariables{7}=temp;
AllGliderVariables{8}=par;
AllGliderVariables{9}=scatter;
AllGliderVariables{10}=do_c;
AllGliderVariables{11}=cdom;


%flight variables.
AllGliderFliVariables{1}=Dflidates;
AllGliderFliVariables{2}=Dfliu;
AllGliderFliVariables{3}=Dfliv;
save([myoutput,'AllgliderDataV3.mat'],'AllGliderVariables')
save([myoutput,'AllgliderFli_DataV3.mat'],'AllGliderFliVariables')



