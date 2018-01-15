function [versionnumber,lastchangedate,lastchangetime]= EDversion
% This function returns the current version of the EDtoolbox
%
% Output parameters:
%   versionnumber       A number, like 0.1
%   lastchangedate      A string with the date of the last change,
%                       on the form '15Jan2018'
%   lastchangetime      A string with the time of the last change, 
%                       on the form '15h35m11'
%
% Peter Svensson 15 Jan. 2018 (peter.svensson@ntnu.no)
%
% [versionnumber,lastchangedate,lastchangetime]= EDversion;

% 15 Jan. 2018 First version

versionnumber = 0.1;
lastchangedate = '15Jan2018';
lastchangetime = '15h35m11';
