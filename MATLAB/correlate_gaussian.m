% CORRELATE_GAUSSIAN
%
%   Following method by Rehfeld et al. (Nonlin. Processes Geophys. 2011)
%
%   Parameters:
%      x     timeseries (N x 1)
%      y     timeseries (M x 1)
%      tx    time intervals associated with x
%      ty    time intervals associated with y
%
%   Returns:
%      corrgauss

function corrgauss=correlate_gaussian(x,y,tx,ty)
  n = length(x);
  m = length(y);

  % calculate mean for x,y series
  xmean = nanmean(x);
  txmean = (tx(end)-tx(1)) / (n-1);

  ymean = nanmean(y);
  tymean = (ty(end)-ty(1)) / (m-1);

  delta_t = max([txmean,tymean]);
  h = delta_t / 4;

  rx = x(1:end, ones(m,1));
  ry = y(1:end, ones(n,1));
  rtx = tx(1:end, ones(m,1));
  rty = ty(1:end, ones(n,1));

  d = rty - rtx';
  b = exp(-1*(d.^2)/(2*(h^2))) / sqrt(2*pi*h);
  num = (rx'-xmean).*(ry-ymean).*b;
  sdx = b.*(rx'-xmean).^2;
  sdy = b.*(ry-ymean).^2;

  corrgauss = sum(num(:)) / sqrt(sum(sdx(:))*sum(sdy(:)));
end
