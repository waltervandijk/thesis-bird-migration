%%% Theoretical bird migration model
% Walter van Dijk
% Bachelor Thesis
% 2018

% University of Amsterdam

%% Notes
% making the shapefile consist of less parts will improve code performance




%% Initialization %%
    %% Cleaning environment
    clear 
    close all
    clc

    %% Performance check
    tic
    profile on

    %% Loading data and assigning parameters
    year=2013;

        %% Loading data
        % Load basemap
        load('ShapeOfEurope.mat');
        load('Europemap.mat');

        % Load weather data
        Meteo=['MeteoMatrix',num2str(year), '_925'];
        load (Meteo);

        %% Control variables
        nDays=10;
        nTracks=1000;  
        nSteps=24;
        dt=0.5;         % timestep in hours
        as=10;          % airspeed [m s-1]
        winf=1;         % wind influence multiplier [0-1]
        TakeOffDeg=25;  % maximum difference in degrees between groundspeed and airspeed at takeoff, NaN = off
        WThreshold=0.5*as; % Threshold under which windspeed has no influence on take off
        CompAmount=0.5;   % flight compensation multiplier [0-1]
        MaxComp=90;     % maximum compensation angle [0-90 deg] (true max comp = CompAmount * MaxComp)
        randHD=3;       % randomness in endogenous direction [0-x deg] (+-3)
        DayS=232;       % startday
        ExhDays=1;      % amount of days the birds can keep flying after a night flying
        DispAmount=100;  % amount of individuals to visualize per group
            % note: nTrack/nDays should be more than DispAmount
        WriteVideo=0;   % write video? 0=no 1=yes
        
        nExtraSteps=(24-nSteps*dt)/dt;      % the amount of extra steps to fill 24h
        sizeGroups=floor(nTracks/nDays);        % quantity of birds per group every day
        
        % Destination locations
        DestLong(1:nTracks,1)=-6+4*rand(nTracks,1)';
        DestLat(1:nTracks,1)=30+1*rand(nTracks,1)';
        
        %% Initiation of matrices
        Time(1:nTracks,1:nSteps+nExtraSteps,1:nDays)=NaN;
        Lat(1:nTracks,1:nSteps+nExtraSteps,1:nDays)=NaN;
        Long(1:nTracks,1:nSteps+nExtraSteps,1:nDays)=NaN;
        wu(1:nTracks,1:nSteps+nExtraSteps,1:nDays)=NaN;
        wv(1:nTracks,1:nSteps+nExtraSteps,1:nDays)=NaN;
        gu(1:nTracks,1:nSteps+nExtraSteps,1:nDays)=NaN;
        gv(1:nTracks,1:nSteps+nExtraSteps,1:nDays)=NaN;
        au(1:nTracks,1:nSteps+nExtraSteps,1:nDays)=NaN;
        av(1:nTracks,1:nSteps+nExtraSteps,1:nDays)=NaN;
        GS(1:nTracks,1:nSteps+nExtraSteps,1:nDays)=NaN;
        TDist(1:nTracks,1:nSteps+nExtraSteps,1:nDays)=0;
        
        HD(1:nTracks,1:nDays)=NaN;          % Endogenous direction
          
        Distu_unit(1:nTracks,1:nSteps+nExtraSteps,1:nDays)=NaN;
        Distv_unit(1:nTracks,1:nSteps+nExtraSteps,1:nDays)=NaN;
        LongKM(1:nTracks,1:nDays)=NaN;
        LatKM(1:nTracks,1:nDays)=NaN;
        ExtraTravel=zeros(nTracks,1);
        TiredBirds=zeros(nTracks,1);
        ArrivedBirds=zeros(nTracks,1);
        AboveLand=zeros(nTracks,1);
        LastStep(1:nTracks,1:nDays)=nSteps;     % at least nSteps are made

        % Suncycle interpolation
        SSLong= -16:10:34 ; SSLong(2,:)=SSLong(1,:);
        SSLat=[65 65 65 65 65 65 ; 30 30 30 30 30 30];
        SSTime(2,6)=NaN;



        %% Visualisation
        % image width and height
        imgWidthPix = 1900;
        imgHeightPix = 1050;
        % set the image sizes, resolution and units; also see print option -r100
        imgResolution = 100; % 100 dpi
        imgWidth = imgWidthPix / imgResolution;
        imgHeight = imgHeightPix / imgResolution;
        % this line is necessary for rendering without openGL drivers/physical screen
        set(0, 'DefaultFigureRenderer', 'zbuffer');
            scrsz = get(0,'ScreenSize');
        
        % create vector of the tracks that are plotted
        plotTracks(1:DispAmount*nDays)=NaN;
        plotTracks(1:DispAmount)=1:DispAmount;
        for i=1:nDays-1
            plotTracks(i*DispAmount+1:i*DispAmount+DispAmount)= (i*sizeGroups+1):(i*sizeGroups+DispAmount);
        end
        [~,plotDelete]=size(plotTracks);      
        
        % Weather data bounding box
        minLong = min(Mlong(1,:,1))/100;
        maxLong = max(Mlong(1,:,1))/100;
        minLat = min(Mlat(:,1,1))/100;
        maxLat = max(Mlat(:,1,1))/100;

        % Video
        if WriteVideo ==1
            v = VideoWriter(['Day', num2str(DayS), '_winf', num2str(winf),'_TOdeg',num2str(TakeOffDeg),'_Comp',num2str(CompAmount) '_as', num2str(as), '.mp4'], 'MPEG-4');
            v.FrameRate=10;
            open(v)
        end
%% Dynamic calculations %%
%% Days
for ii=1:nDays
    
    % Setting bird starting locations 
        for i=sizeGroups*ii-sizeGroups+1:sizeGroups*ii
            Long(i,1,ii)=5+20*rand;
            Lat(i,1,ii)=56+5*rand;
            while ~inpolygon (Long(i,1,ii), Lat(i,1,ii), CoastLon(6).X, CoastLat(6).Y)
                Long(i,1,ii)=5+25*rand;
                Lat(i,1,ii)=56+5*rand;
            end
        end
    
    
    % Title time
        Tdum=DayS-2+ii+datenum([num2str(year),'-01-01 00:00:00']);
        TT=datestr(Tdum, 'yyyy mm dd'); TT(5)=[]; TT(7)=[];
 
    % Setting the location at first timestep
        if ii ~= 1  % days
            Long(1:sizeGroups*ii-sizeGroups,1,ii)=Long(1:sizeGroups*ii-sizeGroups,nSteps+nExtraSteps,ii-1);
            Lat(1:sizeGroups*ii-sizeGroups,1,ii)=Lat(1:sizeGroups*ii-sizeGroups,nSteps+nExtraSteps,ii-1);  
        end
    
    % Removing unrelevant birds
        % Find birds that are exhausted above sea and remove them
        if ii ~= 1  % days
            TiredBirds=TiredBirds+ExtraTravel(1:nTracks)>ExhDays;
            Lat(TiredBirds,1:nSteps,1:nDays)=NaN;
            Long(TiredBirds,1:nSteps,1:nDays)=NaN;
        end
        % Find birds that have arrived and remove them
        if ii ~= 1  % days
            ArrivedBirds=ArrivedBirds+Lat(1:nTracks,nSteps+nExtraSteps,ii-1)<35;
            Lat(ArrivedBirds,1:nSteps,1:nDays)=NaN;
            Long(ArrivedBirds,1:nSteps,1:nDays)=NaN;
        end
    
    % Calculate HD - Great circle azimuth    in function gcazimuth
    HD(1:nTracks,ii)=gcazimuth([DestLat,DestLong],[Lat(1:nTracks,1,ii),Long(1:nTracks,1,ii)]);
    HD(1:nTracks,ii)=randHD*randn(nTracks,1)+HD(1:nTracks,ii);
    StepHD(1:nTracks)=HD(1:nTracks,ii);
    
    meanHD=nanmean(HD(:,ii));
    
    %% Steps
    for i=1:nSteps      %24 steps = 0.5 day ; 60 steps max
        %% Setting the time        
        if i==1   % starttime     
            for ic=1:6
                for ir=1:2  % suncycle strting points
                   DUM(1:2)=suncycle(SSLat(ir,ic), SSLong(ir,ic), [str2double(TT(1:4)), str2double(TT(5:6)), str2double(TT(7:8))], 2880)/24; 
                   SSTime(ir,ic)=DUM(2);
                end
            end 
            DumTime=interpn(SSLat,SSLong,SSTime,Lat(1:nTracks,1,ii),Long(1:nTracks,1,ii),'linear', -999);
            Time(1:nTracks,1,ii)=DayS-1+ii+DumTime; 
%             Time(1:nTracks,i,ii)=DayS-1+ii+1800/2400;     % start at 18.00 instead of sunset
        else        % add timestep to time
           Time(1:nTracks,i,ii)=Time(1:nTracks,i-1,ii)+dt/24;      % [day]
        end
        
        %% Calculating new bird locations
        % Calculating wind speed in m s^-1
        wu(1:nTracks,i,ii)=interpn(Mlat/100,Mlong/100,MTime, Mu, Lat(1:nTracks,i,ii),Long(1:nTracks,i,ii),Time(1:nTracks,i,ii),'linear', 0);
        wv(1:nTracks,i,ii)=interpn(Mlat/100,Mlong/100,MTime, Mv, Lat(1:nTracks,i,ii),Long(1:nTracks,i,ii),Time(1:nTracks,i,ii),'linear', 0);
        % Compensation in flight behaviour
        if i~=1
            % find direction and strength of wind
            Wstr=winf*(wu(1:nTracks,i,ii).*wu(1:nTracks,i,ii)+wv(1:nTracks,i,ii).*wv(1:nTracks,i,ii)).^0.5;
            Wdir=atan2(wu(1:nTracks,i,ii),wv(1:nTracks,i,ii));
            NegDum=find(Wdir<0);
            HighDum=find(Wdir>2*pi);
            Wdir(NegDum)=Wdir(NegDum)+2*pi;
            Wdir(HighDum)=Wdir(HighDum)-2*pi;
            Wdir=Wdir/(2*pi)*360;
            % caclulate the needed compensation to keep direction at HD
            AngleA=Wdir-HD(:,ii);
            DriftFromHD=sind(AngleA).*Wstr;
            CalcAng=find(abs(DriftFromHD/as) < 1-cosd(MaxComp)); 
            NotCalcAng=find(abs(DriftFromHD/as) >= 1-cosd(MaxComp));
            AngleB=zeros(nTracks,1); % initialization
            AngleB(CalcAng)=acosd(abs(DriftFromHD(CalcAng)/as));
            CompAngleNeeded=90-AngleB;
            CompAngleNeeded(NotCalcAng)=MaxComp;
            
            CompDumAngle=find(CompAngleNeeded>MaxComp);
            Compensation=CompAngleNeeded;
            Compensation(CompDumAngle)=MaxComp;
            % set negative values
            NegAngle=find(HD(:,ii)<Wdir & HD(:,ii)-Wdir<180 | HD(:,ii)+180<Wdir);
            Compensation(NegAngle)=Compensation(NegAngle)*-1;
            % update HD
            StepHD=HD(:,ii)+Compensation*CompAmount;
        end
        % calculate airspeed and groundspeed in m s^-1
        au(1:nTracks,i,ii)=as*sin(StepHD/180*pi);
        av(1:nTracks,i,ii)=as*cos(StepHD/180*pi); 
        gu(1:nTracks,i,ii)=au(1:nTracks,i,ii)+winf*wu(1:nTracks,i,ii); % groundspeed u
        gv(1:nTracks,i,ii)=av(1:nTracks,i,ii)+winf*wv(1:nTracks,i,ii); % groundspeed v
        GS(1:nTracks,i,ii)=(gu(1:nTracks,i,ii).*gu(1:nTracks,i,ii)+gv(1:nTracks,i,ii).*gv(1:nTracks,i,ii)).^0.5; %ground speed [m s^-1]

        
        % Add distance traveled to total distance
        if i==1
           TDist(1:nTracks,i,ii)=GS(1:nTracks,i,ii)*3600*dt/1000;   %[km] 
        else
           TDist(1:nTracks,i,ii)=TDist(1:nTracks,i-1,ii)+GS(1:nTracks,i,ii)*3600*dt/1000; % total route distance [km]
        end

        % km per degree
        Distu_unit(1:nTracks,i,ii)=distWBvector([Lat(1:nTracks,i,ii) Long(1:nTracks,i,ii)-0.5] ,[Lat(1:nTracks,i,ii) Long(1:nTracks,i,ii)+0.5]);
        Distv_unit(1:nTracks,i,ii)=distWBvector([Lat(1:nTracks,i,ii)-0.5 Long(1:nTracks,i,ii)] ,[Lat(1:nTracks,i,ii)+0.5 Long(1:nTracks,i,ii)]);
       
        % Write new bird locations
       % if i<nSteps
            Lat(1:nTracks,i+1,ii)=Lat(1:nTracks,i,ii)+gv(1:nTracks,i,ii).*3600*dt./1000./Distv_unit(1:nTracks,i,ii);
            Long(1:nTracks,i+1,ii)=Long(1:nTracks,i,ii)+gu(1:nTracks,i,ii).*3600*dt./1000./Distu_unit(1:nTracks,i,ii);
       % end
        
        %% Reset location if bird needs to rest or add time to ExtraTravel
        AboveLand=zeros(nTracks,1); % reset checking parameter
        for j=1:6 % find birds that are above land
            AboveLandDum = inpolygon(Long(1:nTracks,i,ii), Lat(1:nTracks,i,ii), CoastLon(j).X, CoastLat(j).Y);
%             AboveLandDum = inpolygon(Long(1:nTracks,i,ii), Lat(1:nTracks,i,ii), ...
%                 ShapeOfEurope(j).X,ShapeOfEurope(j).Y);    % Other shapefile
            AboveLand = AboveLandDum + AboveLand;
        end

        RestInd = find(ExtraTravel~=0);
        for j=RestInd'
            if AboveLand(j) ~= 0    % Reset location and subtract dt from rest time
                Lat(j,i+1,ii)=Lat(j,i,ii);
                Long(j,i+1,ii)=Long(j,i,ii);
                if i~=1
                    TDist(j,i,ii)=TDist(j,i-1,ii);
                else
                    TDist(j,i,ii)=0;
                end
                ExtraTravel(j) = ExtraTravel(j)-dt/24;
                if ExtraTravel(j)<0
                    ExtraTravel(j)=0;
                end
            elseif AboveLand(j) == 0    % Add dt to rest time
                ExtraTravel(j) = ExtraTravel(j)+dt/24;
            end
        end
        
        %% Reset location if weather prevents takeoff
                % Birds only depart when wind is in the right direction
        if i==1
            BadWeatherDum(1:nTracks,1:4)=NaN;
            for j=1:4           % calculate wind speeds in past 2 hours
                Takeoffwu=winf*interpn(Mlat/100,Mlong/100,MTime, Mu, Lat(1:nTracks,i,ii),Long(1:nTracks,i,ii),Time(1:nTracks,i,ii)-j*dt/24,'linear', 0);
                Takeoffwv=winf*interpn(Mlat/100,Mlong/100,MTime, Mv, Lat(1:nTracks,i,ii),Long(1:nTracks,i,ii),Time(1:nTracks,i,ii)-j*dt/24,'linear', 0);
            
                TakeoffWstr=(Takeoffwu.*Takeoffwu+Takeoffwv.*Takeoffwv).^0.5;
            
                TakeoffWdir=atan2(Takeoffwu,Takeoffwv);
                NegDum=find(TakeoffWdir<0);
                HighDum=find(TakeoffWdir>2*pi);
                TakeoffWdir(NegDum)=TakeoffWdir(NegDum)+2*pi;
                TakeoffWdir(HighDum)=TakeoffWdir(HighDum)-2*pi;
                TakeoffWdir=TakeoffWdir/(2*pi)*360;
            
                BadWeatherDum(1:nTracks,j) = (abs(HD(:,ii)-TakeoffWdir) > TakeOffDeg & TakeoffWstr > WThreshold);
            end
            BadWeather=find(sum(BadWeatherDum,2)>2);
        end
        
        for j=BadWeather'
            if AboveLand(j) ~= 0    % Reset location and subtract dt from rest time
                Lat(j,i+1,ii)=Lat(j,i,ii);
                Long(j,i+1,ii)=Long(j,i,ii);
                if i~=1
                    TDist(j,i,ii)=TDist(j,i-1,ii);
                else
                    TDist(j,i,ii)=0;
                end
                ExtraTravel(j) = ExtraTravel(j)-dt/24;
                if ExtraTravel(j)<0
                    ExtraTravel(j)=0;
                end
            elseif AboveLand(j) == 0    % Add dt to rest time
                ExtraTravel(j) = ExtraTravel(j)+dt/24;
            end
        end
        
        
    end %for i=1:nSteps
    
    
%% Land and water
    % To keep birds flying when above water
    for i=nSteps:nSteps+nExtraSteps
        %% Find indexes of birds that are above sea
        AboveLand=zeros(nTracks,1); % reset checking parameter
        for j=1:6
            AboveLandDum = inpolygon(Long(1:nTracks,i,ii), Lat(1:nTracks,i,ii), CoastLon(j).X, CoastLat(j).Y);
%                 AboveLandDum = inpolygon(Long(1:nTracks,i,ii), Lat(1:nTracks,i,ii), ...
%                     ShapeOfEurope(j).X,ShapeOfEurope(j).Y);     % Other shapefile
            AboveLand = AboveLandDum + AboveLand;
            AboveLand(isnan(Long(1:nTracks,i,ii)))=1;
        end
        AboveSeaInd = find(AboveLand == 0);
        AboveLandInd = find(AboveLand ~= 0);

        %% Birds that are above land stay there
        if i<nSteps+nExtraSteps
            Lat(AboveLandInd,i+1,ii)=Lat(AboveLandInd,i,ii);
            Long(AboveLandInd,i+1,ii)=Long(AboveLandInd,i,ii);
        end
            TDist(AboveLandInd,i,ii)=TDist(AboveLandInd,i-1,ii);
            
        %% Calculate locations of sea birds at next timestep            
        if min(AboveLand) == 0      % if there are any birds above the water
            Time(AboveSeaInd,i,ii)=Time(AboveSeaInd,i-1,ii)+dt/24;      % [day]
            
            % Calculating wind speed in m s^-1            
            wu(AboveSeaInd,i,ii)=interpn(Mlat/100,Mlong/100,MTime, Mu, Lat(AboveSeaInd,i,ii),Long(AboveSeaInd,i,ii),Time(AboveSeaInd,i,ii),'linear', 0);
            wv(AboveSeaInd,i,ii)=interpn(Mlat/100,Mlong/100,MTime, Mv, Lat(AboveSeaInd,i,ii),Long(AboveSeaInd,i,ii),Time(AboveSeaInd,i,ii),'linear', 0);

            % Compensation in flight behaviour
            if i~=1
                % find direction and strength of wind
                Wstr=winf*(wu(1:nTracks,i,ii).*wu(1:nTracks,i,ii)+wv(1:nTracks,i,ii).*wv(1:nTracks,i,ii)).^0.5;
                Wdir=atan2(wu(1:nTracks,i,ii),wv(1:nTracks,i,ii));
                NegDum=find(Wdir<0);
                HighDum=find(Wdir>2*pi);
                Wdir(NegDum)=Wdir(NegDum)+2*pi;
                Wdir(HighDum)=Wdir(HighDum)-2*pi;
                Wdir=Wdir/(2*pi)*360;

                % caclulate the needed compensation to keep direction at HD
                AngleA=Wdir-HD(:,ii);
                DriftFromHD=sind(AngleA).*Wstr;
                CalcAng=find(abs(DriftFromHD/as) < 1-cosd(MaxComp)); 
                NotCalcAng=find(abs(DriftFromHD/as) >= 1-cosd(MaxComp));
                AngleB=zeros(nTracks,1); % initialization
                AngleB(CalcAng)=acosd(abs(DriftFromHD(CalcAng)/as));
                CompAngleNeeded=90-AngleB;
                CompAngleNeeded(NotCalcAng)=MaxComp;
                CompDumAngle=find(CompAngleNeeded>MaxComp);
                Compensation=CompAngleNeeded;
                Compensation(CompDumAngle)=MaxComp;

                % set negative values
                NegAngle=find(HD(:,ii)<Wdir & HD(:,ii)-Wdir<180 | HD(:,ii)+180<Wdir);
                Compensation(NegAngle)=Compensation(NegAngle)*-1;

                % update HD
                StepHD=HD(:,ii)+Compensation*CompAmount;
            end
            
            % calculate airspeed and groundspeed in m s^-1
            au(AboveSeaInd,i,ii)=as*sin(StepHD(AboveSeaInd)/180*pi);
            av(AboveSeaInd,i,ii)=as*cos(StepHD(AboveSeaInd)/180*pi);
            gu(AboveSeaInd,i,ii)=au(AboveSeaInd,i,ii)+winf*wu(AboveSeaInd,i,ii); % groundspeed u
            gv(AboveSeaInd,i,ii)=av(AboveSeaInd,i,ii)+winf*wv(AboveSeaInd,i,ii); % groundspeed v
            GS(AboveSeaInd,i,ii)=(gu(AboveSeaInd,i,ii).*gu(AboveSeaInd,i,ii)+gv(AboveSeaInd,i,ii).*gv(AboveSeaInd,i,ii)).^0.5; %ground speed [m s^-1]

            % Add distance traveled to total distance      
            TDist(AboveSeaInd,i,ii)=TDist(AboveSeaInd,i-1,ii)+GS(AboveSeaInd,i,ii)*3600*dt/1000;
            
            % km per degree
            Distu_unit(AboveSeaInd,i,ii)=distWBvector([Lat(AboveSeaInd,i,ii) Long(AboveSeaInd,i,ii)-0.5] ,[Lat(AboveSeaInd,i,ii) Long(AboveSeaInd,i,ii)+0.5]); %km per degree
            Distv_unit(AboveSeaInd,i,ii)=distWBvector([Lat(AboveSeaInd,i,ii)-0.5 Long(AboveSeaInd,i,ii)] ,[Lat(AboveSeaInd,i,ii)+0.5 Long(AboveSeaInd,i,ii)]); %km per degree
            
            % Write new bird locations
            if i<nSteps+nExtraSteps
                Lat(AboveSeaInd,i+1,ii)=Lat(AboveSeaInd,i,ii)+gv(AboveSeaInd,i,ii).*3600*dt./1000./Distv_unit(AboveSeaInd,i,ii);
                Long(AboveSeaInd,i+1,ii)=Long(AboveSeaInd,i,ii)+gu(AboveSeaInd,i,ii).*3600*dt./1000./Distu_unit(AboveSeaInd,i,ii);
                ExtraTravel(AboveSeaInd)=ExtraTravel(AboveSeaInd)+dt/24;    %store additional time travelled [day]
                LastStep(AboveSeaInd,ii)=LastStep(AboveSeaInd,ii)+1;    
            end
        end % if there are birds above water        
            
        %% Prevent movement if birds decided not to fly this day
        
            for j=BadWeather
                if AboveLand(j) ~= 0    % Reset location and subtract dt from rest time
                    if i<nSteps+nExtraSteps
                        Lat(j,i+1,ii)=Lat(j,i,ii);
                        Long(j,i+1,ii)=Long(j,i,ii);
                    end
                    if i~=1
                        TDist(j,i,ii)=TDist(j,i-1,ii);
                    else
                        TDist(j,i,ii)=0;
                    end
                    ExtraTravel(j) = ExtraTravel(j)-dt/24;
                    if ExtraTravel(j)<0
                        ExtraTravel(j)=0;
                    end
                elseif AboveLand(j) == 0    % Add dt to rest time
                    ExtraTravel(j) = ExtraTravel(j)+dt/24;
                end
            end
        
        
    end % for i=nSteps+1:nSteps+nExtraSteps

%% Removing unrelevant birds
    % Find birds that are exhausted above sea and remove them
    TiredBirds=TiredBirds+ExtraTravel(1:nTracks)>ExhDays;
    Lat(TiredBirds,1:nSteps,1:nDays)=NaN;
    Long(TiredBirds,1:nSteps,1:nDays)=NaN;

    % find birds that have arrived and remove them
    ArrivedBirds=ArrivedBirds+Lat(1:nTracks,i,ii)<35;
    Lat(ArrivedBirds,1:nSteps,1:nDays)=NaN;
    Long(ArrivedBirds,1:nSteps,1:nDays)=NaN;
    
    
%% Visualization
    %% Make figure
    if ii==1     
        scrsz = get(0,'ScreenSize');
        figure('Position',[30 scrsz(4)/20 scrsz(3)/2 scrsz(4)/1.2]); 
        xlim([-15 35]); ylim([30 65]) 
        xlabel('Longitude')
        ylabel('Latitude')
        hold on
    end

    %% Background
    if ii==1    % draw map of Europe
        for i=1:6
            h1=plot(CoastLon(i).X,CoastLat(i).Y,'k');
        end
%         for j=1:62    % other shapefile
%             h1=plot(ShapeOfEurope(j).X,ShapeOfEurope(j).Y,'k')
%         end
             % draw bounding box weather data
        BoundBox = [minLat,minLong; minLat,maxLong; ...
            maxLat,maxLong; maxLat,minLong; minLat,minLong];
        h2=plot(BoundBox(:,2),BoundBox(:,1),'r-');
    end
    
    if ii~=1
        children = get(gca, 'children');
        delete(children(1:plotDelete*3));
    end
    %% Make birds fly
    for j=plotTracks
        h3=plot(Long(j,1:LastStep(j,ii),ii),Lat(j,1:LastStep(j,ii),ii),'--','color',[0.5 0.5 0.5]);      % tracks
        
        h4=plot(Long(j,1,ii),Lat(j,1,ii),'.','color',[0.3 0.3 0.3]);         % beginpoints every day
        h5=plot(Long(j,LastStep(j,ii),ii),Lat(j,LastStep(j,ii),ii),'g.');    % endpoints every day

    end
    
    %% Formatting
     if ii==1    % add legend
         legend([h1 h2 h3 h4 h5],{'Coastline','BoundingBox','Tracks','Startpoints','Endpoints'})
     end
    
    title([TT, ', AS=',num2str(as), ' Winf', num2str(winf),' TOdeg',num2str(TakeOffDeg),' Comp',num2str(CompAmount)...
        ,' EndedInSea=',num2str(sum(TiredBirds))])
    drawnow
       
    %% Video
    if WriteVideo==1
        for k=1:8
            writeVideo(v,getframe(gcf))
        end
    end
end %for ii= nDays

if WriteVideo==1
    close(v)
end
%% Analysis
%     %% Distance travelled
%     SummedDist(1,1:nTracks)=0; %initialization
%     for j=1:nDays
%         SummedDist=nansum([SummedDist;TDist(:,48,j)']);
%     end
% 
%     
%     %% Radars
%     RadarArea=100;  % n*n km    n will be halved to obtain radar reach
%     RadarLat=[55*ones(1,16),53*ones(1,16),51*ones(1,16),49*ones(1,16),...
%         47*ones(1,16),45*ones(1,16),43*ones(1,16),41*ones(1,16),39*ones(1,16),37*ones(1,16),35*ones(1,16)];
%     RadarLong=[-5:2:25,-5:2:25,-5:2:25,-5:2:25,-5:2:25,-5:2:25,-5:2:25,...
%         -5:2:25,-5:2:25,-5:2:25,-5:2:25];
%     [~,AmountRadars]=size(RadarLat);
%     RadarOcc(1:AmountRadars)=0;
%     
%     % km degree-1
%     RadarDistu_unit=distWBvector([RadarLat(:) RadarLong(:)-0.5] ,[RadarLat(:) RadarLong(:)+0.5]);
%     RadarDistv_unit=distWBvector([RadarLat(:)-0.5 RadarLong(:)] ,[RadarLat(:)+0.5 RadarLong(:)]);
%     
%     % radar reach in lat long degrees
%     DegLatRadar=0.5*RadarArea./RadarDistv_unit;
%     DegLongRadar=0.5*RadarArea./RadarDistu_unit;
%     
%     % RadarOcc contains the amount of birds that are within reach
%     for j=1:AmountRadars
%         RadarOcc(j) = sum(Lat(:,i,ii)>RadarLat(j)-DegLatRadar(i) & Lat(:,i,ii)<RadarLat(j)+ DegLatRadar(j) ...
%             & Long(:,i,ii)>RadarLong(j)-DegLongRadar(j)  & Long(:,i,ii)<RadarLong(j)+DegLongRadar(j));
%     end
%     
%     % Show locations of radars
%         plot(RadarLong,RadarLat,'b+')
%     
%% Saving
    % Figure
    Filename=['M' TT '_' num2str(as) '.png'];
    print ('-dpng', '-r100', Filename)
    
    % Analysis
%     save(['RadarData',num2str(year),'_day',num2str(DayS+ii-1),'_winf',num2str(winf),'_as',num2str(as),'_todeg',num2str(TakeOffDeg)], 'RadarOcc');
%     save(['AllResults',num2str(year)]);

%% Performance check
    toc
    profile viewer

