function [km] = distWBvector(Loc1, Loc2)
% Loc1 = [52.3826 -4.637];
% Loc2 = [50.9235 5.7812];
% Loc=[Lat Long]

l1=deg2rad(Loc1);
l2=deg2rad(Loc2);
lat=[l1(:,1) l2(:,1)];
lon=[l1(:,2) l2(:,2)];
% Example     
% distWBhaversine([long1 lat1], [long2 lat2]) returns distance in [km]
% distWBhaversine([53.1472 -1.8494], [52.2044 0.1406]) returns 170.2563
%
%   Inputs
%       LOC must be a 2-valued numeric array specifying the
%       location in decimal degrees.  
%
%   Notes
%       The Haversine formula is used to calculate the great-circle
%       distance between two points, which is the shortest distance over
%       the earth's surface.
%
%       This program was created using equations found on the website
%       http://www.movable-type.co.uk/scripts/latlong.html

% delta_lat = locs{2}(1) - locs{1}(1);        % difference in latitude
% delta_lon = locs{2}(2) - locs{1}(2);        % difference in longitude
% a = sin(delta_lat/2)^2 + cos(locs{1}(1)) * cos(locs{2}(1)) * ...
%     sin(delta_lon/2)^2;
%% Begin calculation

R = 6371;                                   % Earth's radius in km
delta_lat = lat(:,2) - lat(:,1);        % difference in latitude
delta_lon = lon(:,2) - lon(:,1);        % difference in longitude
a = sin(delta_lat/2).^2 + cos(lat(:,1)).*cos(lat(:,2)).*sin(delta_lon/2).^2;
%c=2*asin(min(1,sqrt(a)));
b = 2 .* atan2(sqrt(a), sqrt(1-a));
km = R * b ;                                % distance in km
