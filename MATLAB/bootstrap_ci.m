% BOOSTRAP_CI
%
%   Parameters:
%      x     timeseries (N x 1)
%      y     timeseries (M x 1)
%      tx    time intervals associated with x
%      ty    time intervals associated with y
%
%   Returns:
%      bci1, bci2   upper and lower confidence intervals

function [bci1, bci2] = bootstrap_ci(x, y, tx, ty)
  max_iter = 2000;
  n = length(x);
  m = length(y);

  temp = zeros(m+n,1);
  indx = zeros(m+n,2);
  seq = zeros(m+n+1,1);
  tstar = zeros(m+n,1);
  inds = zeros(m+n,1);
  nx = zeros(m+n,1);
  ny = zeros(m+n,1);
  ntx = zeros(m+n,1);
  nty = zeros(m+n,1);
  bse = zeros(max_iter,1);
  med_min = zeros(25,1);
  med_max = zeros(25,1);

  orig = correlate_gaussian(x, y, tx, ty);

  % estimate persistence times
  txmean = (tx(n)-tx(1)) / (n-1);
  tymean = (ty(m)-ty(1)) / (m-1);
  tmean = max([txmean,tymean]);

  for i = 1:n/16
    tmp = tx(1:n) + i*txmean;
    taux = correlate_gaussian(x, x, tx, tmp(1:n));
    if abs(taux) < 0.368
      break
    end
  end
  taux = i*txmean;

  for j = 1:m/16
    tmp = ty(1:m) + j*tymean;
    tauy = correlate_gaussian(y, y, ty, tmp(1:m));
    if abs(tauy) < 0.368
      break
    end
  end
  tauy = j*tymean;

  p = 1.0 - tmean / (4*max([taux,tauy]));
  logp = log(p);

  % build integrated time-list
  temp = vertcat(tx(1:n),ty(1:m));

  % index elements for series x are positive integers, and negative for y
  tmp_indx = horzcat([1:n,-1:-1:-m]);

  [temp, tmp_indx] = joint_sort(temp, tmp_indx);
  j = 1;
  for i = 1:n+m
    tstar(j) = temp(i);
    if tmp_indx(i) > 0
      indx(j,1) = tmp_indx(i);
    else
      indx(j,2) = tmp_indx(i);
    end
    if i < n+m
      if temp(i) ~= temp(i+1)
	j = j + 1;
      end
    end
  end
  nm = j;

  % do the bootstraps
  for outer = 1:25
    n_low = 0;
    for iter = 1:max_iter % inner iteration
      j = 1;
      init = [0 0];
      rnd = rand;
      % random starting point
      i = floor((nm-1)*rnd) + 1;
      if indx(i,1) ~= 0
        seq(j) = indx(i,1);
	j = j + 1;
	init(1) = 1;
      end
      if indx(i,2) ~= 0
	seq(j) = indx(i,2);
	j = j + 1;
	init(2) = 1;
      end

      % loop for rest of points to generate seq of indices
      while j <= n+m
        rnd = rand;
	if i > 1
	  dt = tstar(i) - tstar(i-1);
	else
	  dt = tstar(2) - tstar(1);
	end
	if rnd > exp(dt/tmean*logp)
	  if ~all(init)
	    j = 1;
	    init = [0 0];
	  end
	  rnd = rand;
	  i = floor((nm-1)*rnd) + 1;
	  if indx(i,1) ~= 0
	    seq(j) = indx(i,1);
	    j = j + 1;
	  end
	  if indx(i,2) ~= 0
	    seq(j) = indx(i,2);
	    j = j + 1;
	  end
	else
	  i = i + 1;
	  if i > nm
	    i = 1;
	  end
	  if indx(i,1) ~= 0
	    seq(j) = indx(i,1);
	    j = j + 1;
	    init(1) = 1;
	  end
	  if indx(i,2) ~= 0
	    seq(j) = indx(i,2);
	    j = j + 1;
	    init(2) = 1;
	  end
	end
      end

      % build up two series and corresponding time bases
      % find start times for each series
      init = [0 0];
      j = 1;
      k = 1;
      for i = 1:n+m
        if seq(i) > 0
          nx(j) = x(seq(i));
          if init(1) % already have time base
            if seq(i) > 1
              dt = tx(seq(i)) - tx(seq(i)-1);
            else
              dt = tx(2) - tx(1);
            end
            ntx(j) = ntx(j-1) + dt;
          else
            init(1) = 1;
            if (init(2))
              if ty(-seq(i-1)) > tx(seq(i))
                if seq(i) > 1
                  dt = tx(seq(i)) - tx(seq(i)-1);
                else
                  dt = tx(2) - tx(1);
                end
                ntx(j) = nty(k-1) + dt;
	      else
                ntx(j) = nty(k-1) + tx(seq(i)) - ty(-seq(i-1));
              end
            else
              ntx(j) = tx(seq(i));
            end
          end
          j = j + 1;
        else
          ny(k) = y(-seq(i));
          if (init(2)) % already have time base
            if -seq(i) > 1
              dt = ty(-seq(i)) - ty(-seq(i)-1);
            else
              dt = ty(2) - ty(1);
            end
            nty(k) = nty(k-1) + dt;
          else
            init(2) = 1;
            if (init(1))
              if tx(seq(i-1)) > ty(-seq(i))
                if -seq(i) > 1
                  dt = ty(-seq(i)) - ty(-seq(i)-1);
                else
                  dt = ty(2) - ty(1);
                end
                nty(k) = ntx(j-1) + dt;
              else
                nty(k) = ntx(j-1) + ty(-seq(i)) - tx(seq(i-1));
              end
            else
              nty(k) = ty(-seq(i));
            end
          end
          k = k + 1;
        end
      end
      j = j - 1;
      k = k - 1;
      bse(iter) = correlate_gaussian(nx(1:j), ny(1:k), ntx(1:j), nty(1:k));

      if bse(iter) < orig
        n_low = n_low + 1;
      end
    end

    % call shell(bse(1:max_iter),max_iter)
    % calculate bias parameter
    % call PNI(real(n_low)/max_iter,1.0-real(n_low)/max_iter,real(n_low)/max_iter-0.5,z0,ierr)
    z0 = icpdf(real(n_low)/max_iter);

    % calculate acceleration parameter
    asum = sum(bse);
    astar = (sum(asum-bse)/(max_iter-1)) / max_iter;
    anum = sum((astar-(asum-bse)/(max_iter-1)).^3);
    aden = sum((astar-(asum-bse)/(max_iter-1)).^2);
    a = anum / (6.0*exp(1.5*log(aden)));

    % calculate interval points
    b0 = z0 + (z0-1.960) / (1-a*(z0-1.960));
    b1 = z0 + (z0+1.960) / (1-a*(z0+1.960));

    ci1 = nanmin([nanmax([floor(max_iter*cpdf(b0)),1]),max_iter]);
    ci2 = nanmax([nanmin([ceil(max_iter*cpdf(b1)),max_iter]),1]);

    bse = sort(bse);
    med_min(outer) = bse(ci1);
    med_max(outer) = bse(ci2);
    fprintf('%f, %f\n', bse(ci1), bse(ci2));
  end
  bci1 = median(med_min);
  bci2 = median(med_max);
end

