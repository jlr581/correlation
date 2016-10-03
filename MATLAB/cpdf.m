% CPDF
%
%   Cumulative probability density function of variable x
%
%   Parameters:
%      x     vector (N x 1)
% 
%   Returns:
%      cx    CPDF vector (N x 1) of x

function cx = cpdf(x)
  if x>=0
    cx=1-0.5/((1+0.196854*x+0.115194*x.^2+0.000344*x.^3+0.019527*x.^4).^4);
    return
  else
    cx=0.5/((1-0.196854*x+0.115194*x.^2-0.000344*x.^3+0.019527*x.^4).^4);
    return
  end
end
