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
% Peter Svensson (peter.svensson@ntnu.no) 21 Nov. 2026
%
% [versionnumber,lastchangedate,lastchangetime]= EDgetversion;

% 15 Jan. 2018 First version
% 16 Jan. 2018 Changed name to EDgetversion to avoid confusion with the
% variable name
% After 16 Jan 2018: Updated date and/or time (and/or version)
% 15 Oct. 2026 The header was updated (by Claude, prompted by PS).
% 15 Oct. 2026 Version 0.611.
% 15 Oct. 2026 Version 0.612.
% 16 Oct. 2026 Version 0.613.
% 18 Oct. 2026 Version 0.614.
% 20 Oct. 2026 Version 0.615.
%  9 Nov. 2026 version 0.62
% 20 Nov. 2026 Version 0.63: first step towards non-convex models.
% 21 Nov. 2026 Version 0.631: partly visible edges and edge pairs
% (non-convex models).
% 21 Nov. 2026 Version 0.632: a plane wave along a face (the source term
% 'exactpole', the new default).
% 21 Nov. 2026 Version 0.633: a V-groove (non-convex) fix; per-node angles
% and factors for edges partly across an indent.
% 21 Nov. 2026 Version 0.634: no NaN in first-order diffraction for a
% receiver within about 1e-8 rad of a zone boundary.
% 21 Nov. 2026 Version 0.635: controlparameters.specorder = 'auto' (the new
% default for non-convex models) and .specordermax.

versionnumber = 0.635;
lastchangedate = '21Nov2026';
lastchangetime = '22h30';
