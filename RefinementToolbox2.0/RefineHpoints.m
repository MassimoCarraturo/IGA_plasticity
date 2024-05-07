function [hpoints] = RefineHpoints(hspace, hpoints, nref)
% Update the hspace after refining the selected functions. 
%
% Input:
%   basis     - data structure containing information on parametric basis
%   hspace    - data structure containing information on hierarchical basis
%   refFunc   - list of function indices of functions to be refined
%
% Output:
%   hspace    - Updated with the refined function structure
%
%
% Created by Pieter van Zuijlen

if ~exist('nref','var')
    nref = numel(hpoints);
end

K        = RefinementOperator(hspace,nref);

hpoints{nref+1} = K * hpoints{nref}; 

end