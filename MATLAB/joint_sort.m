% JOINT_SORT
%
%   Sorts variable a, a union of two vectors.
%
%   Parameters:
%      a      union of two vectors
%      ia     union of indices of two vectors
% 
%   Returns: 
%      asort  sorted values from a
%      iasort sorted indices from ia   

function [asort,iasort] = joint_sort(a,ia)
  [asort,inds] = sort(a);
  iasort = ia(inds);
end
