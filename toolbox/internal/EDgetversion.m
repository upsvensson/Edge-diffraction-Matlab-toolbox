function [versionnumber,lastchangedate,lastchangetime]= EDgetversion
% EDgetversion - Returns the current version of the EDtoolbox.
%
% Output parameters:
%   versionnumber       A number, like 0.101
%   lastchangedate      A string with the date of the last change,
%                       on the form '15Jan2018'
%   lastchangetime      A string with the time of the last change, 
%                       on the form '15h35m11'
%
% Peter Svensson (peter.svensson@ntnu.no) 15 Oct. 2026
%
% [versionnumber,lastchangedate,lastchangetime]= EDgetversion;

% 15 Jan. 2018 First version
% 16 Jan. 2018 Changed name to EDgetversion to avoid confusion with the
% variable name
% After 16 Jan 2018: Updated date and/or time (and/or version)
% 15 Oct. 2026 The header was updated (by Claude, prompted by PS).

versionnumber = 0.61;
lastchangedate = '14Oct2026';
lastchangetime = '22h27';
