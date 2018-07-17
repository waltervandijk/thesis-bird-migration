function [Dir] = gcazimuth(DestLoc,OriLoc)

% Calculate the direction when leaving from OriLoc to DestLoc in degrees
% Great circle azimuth



% Destination
% Loc1 = [52.3826 -4.637; 52.3826 -4.637; 52.3826 -4.637; 52.3826 -4.637];

% Origen
% Loc2 = [50.9235 5.7812; 50.9235 5.7812; 50.9235 5.7812; 50.9235 5.7812];

% Loc=[Lat Long]

DestLong=DestLoc(:,2);
DestLat=DestLoc(:,1);
OriLong=OriLoc(:,2);
OriLat=OriLoc(:,1);

    % Calculate Dirction in deg - Great circle azimuth
    A=sind(DestLong-OriLong);
    B=cosd(OriLat).*tand(DestLat)-sind(OriLat).*cosd(DestLat-OriLat);
    Dir=rad2deg(atan2(A,B));
    
    NegInd=find(Dir<0);
    Dir(NegInd)=Dir(NegInd)+360;