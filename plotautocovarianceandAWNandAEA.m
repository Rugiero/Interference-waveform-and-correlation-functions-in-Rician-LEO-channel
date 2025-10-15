function [AEAsLOS AEAsRICE AWNtotalsLOS AWNtotalsRICE] ...
	 = plotautocovarianceandAWNandAEA()
  ## Simulate and plot the total interference average envelope amplitude
  ## (AEA) and the corresponding interference waveform. Requires the 
  ## Octave  signal  processing and statistics packages:
  ## https://packages.octave.org
  ## Output:
  ## AEAsLOS: A vector of the AEA in the LoS channel
  ## AEAsRICE: A vector of the AEA in the Rician channel
  ## AWNtotalsLOS: A vector of the interference waveform in the LoS channel
  ## AWNtotalsRICE: A vector of the waveform in the Rician channel 
  
  ## tkappa (tildekappa) implicitly encompasses the terrestrial
  ## interferer density, LEO BS altitude, beamwidth, and elevation angle. 
  ## kappa = log(2) * tildekappa is the average number of transmitters
  ## inside the -3 dB power footprint---further information in the thesis
  tkappa = 1;
  fc = 1 * 10 ^ 6; # Carrier frequency
  cohtime =  4 * pi * 10 ^ 5 / (5 * fc); # Coherence time (s) (tau_c)
  orbitalspeed = 7.4 * 10 ^ 3; # The LEO BS orbital speed
  h = 200 * 10 ^ 3; # The LEO BS altitude
  ## The constant D_{h, epsilon, varphi} for epsilon = 90 deg and
  ## varphi = 1.6 degrees
  dd = 1 / (h ^ 2 * deg2rad(1.6) ^ 2);
  ## Limit the number of fading blocks, set to frac1 = 1 for more accurate
  ## autocovariance estimate
  frac1 = 1 / 1; 
  ## The number of fading blocks is defined such that the LEO BS does not
  ## move outside the rectangular area in planarHPPrefs(.)
  scaling = sqrt(10 ^ 5 / dd); # Scales the PPP region [-0.5 0.5] ^ 2
  fadingblocksN = ceil(frac1 * (scaling / orbitalspeed) / cohtime);
  
  symbolN = 1; # The symbol sample length 
  ## Number of samples per fading block (must have symbolN as a factor). 
  ## Each transmitter is in the same block (with independent fading gains)
  signallength = 1 * symbolN; 
  K = 1000; # The Rice-K fading parameter
  ## Generate a Rice-K distribution of normalized power
  pd = RicianDistribution(sqrt(K / (1 + K)), sqrt(1 / (2 * (1 + K))));
  fadingmean = mean(pd); # The mean of the amplitude fading gain
  fadingvar = var(pd); # The variance of the amplitude fading gain
  
  tic
  ## Simulate the signals
  [AEAsLOS AEAsRICE AWNtotalsLOS AWNtotalsRICE] = ... 
    condAWN(tkappa, fadingblocksN, frac1 / fadingblocksN, ...
	    signallength, symbolN, K, fadingmean);
  toc

  ## Scale the envelopes appropriately  
  AEAsLOS = AEAsLOS / sqrt( 1 + 4 * tkappa * 1 ^ 2 );
  AEAsRICE = AEAsRICE / sqrt( 1 + 4 * tkappa * fadingmean ^ 2);
  
  ## Numerically estimate the AEA autocovariance functions
  maxlag = 10 * sqrt((log(2) / (2 * dd))); # Maximum lag 
  [autocovarianceLoS, lag] = ... # The AEA autocovariance in LoS channel
    xcov(AEAsLOS, [ceil((maxlag / orbitalspeed) ...
			* signallength / cohtime)], "unbiased");
  [autocovarianceRICE, lag] = ... # The AEA autocovariance in Rice channel
    xcov(AEAsRICE, [ceil((maxlag / orbitalspeed) ...
			 * signallength /  cohtime)], "unbiased");

  ## Modulate the samples into an analog signal
  frac2 = 1 / 10; # Limit the length of the modulated signal
  T = 0 : 0.01 : frac2 * fadingblocksN * signallength;
  xbRICE = digitaltoanalog(T, AWNtotalsRICE(1:frac2 * fadingblocksN ...
					      * signallength)');
  xbLOS = digitaltoanalog(T, AWNtotalsLOS(1:frac2 * fadingblocksN ...
					    * signallength)');  
  T1 = linspace(0, fadingblocksN * cohtime, length(AEAsLOS));
  T2 = linspace(0, frac2 * fadingblocksN * cohtime, length(T));
  figure;
  hold on;
  ## Plot the Rice channel AEA
  plot(T1, AEAsRICE', 'linewidth', 3);
  ## Plot the LoS channel AEA
  plot(T1, AEAsLOS', 'linewidth', 3);
  title('Intracell plus intercell interference AEA',...
	'fontname', 'DejaVu Serif');
  legend("Rice-1 channel  AEA_1", ...
	 "LoS channel channel AEA_1",...
	 'fontname', 'DejaVu Serif');
  xlabel("time (s)", 'fontname', 'DejaVu Serif');
  ylabel("Average envelope amplitude (non-dimensional)",...
	 'fontname', 'DejaVu Serif');
  set(gca, 'fontsize', 8);
  grid on;
  hold off;

  figure;
  ## plot(T2,xbLOS, 'linewidth',2) 
  hold on;
  ## Plot the analog Rician response waveform
  plot(T2, xbRICE, 'linewidth', 2) 
    title('Intracell plus intercell interference',...
	'fontname', 'DejaVu Serif');
  legend('AWN total interference waveform',...
	 'fontname','DejaVu Serif')
  xlabel('time (s)', 'fontname', 'DejaVu Serif');
  ylabel("Amplitude (non-dimensional)",...
	   'fontname', 'DejaVu Serif');
  set(gca, 'fontsize', 8);
  grid on;
  hold off;

  ## Plot the theoretical and simulated autocovariance functions 
  ## for the Rice channel
  figure;
  hold on;
  scaledtimelags = lag ./ max(lag) * (maxlag / orbitalspeed);
  actualcohtime = cohtime;
  plot(scaledtimelags, autocovarianceRICE, '-s', 'linewidth',3);
  ## plot(scaledtimelags, autocovarianceLoS, '-s','linewidth',3);
  axis([-maxlag / orbitalspeed maxlag / orbitalspeed 0 0.3]);
  x = linspace(-maxlag / orbitalspeed, maxlag / orbitalspeed, 10000);

  ## The theoretical autocorrelation function
  y1 =  (1 - fadingvar) *  tkappa ... # The LoS covariance component
   	* exp(-dd * (orbitalspeed * x) .^ 2 / 4) ...
	/ ( 1 + 4 * tkappa * fadingmean ^ 2 );
  ## Define the triangular fading autocovariance function
  triangr = tkappa * fadingvar * (actualcohtime - x(find(x > 0))) ...
   	    / actualcohtime .* (abs(x(find(x > 0))) < actualcohtime) ...
	    / ( 1 + 4 * tkappa * fadingmean ^ 2 );
  ## Alternatively, plot the autocovariance of the intercell interference
  ## The footprint restriction of the interferer domain   
  ## must be set in the function GP(.) for the simulation:
  ## y1 =  (1 - fadingvar) * 1 / 2 * tkappa * ... 
  ##       exp(-dd * (orbitalspeed * abs(x)) * 1 / 50 * h) ...
  ## 	/ ( 1 + 4 * tkappa * fadingmean ^ 2 ); 
  ## triangr = 1 / 2 * tkappa * fadingvar ...
  ## 	    * (actualcohtime - x(find(x > 0))) / actualcohtime ...
  ## 	    .* (abs(x(find(x > 0))) < actualcohtime) ...
  ## 	    / ( 1 + 4 * tkappa * fadingmean ^ 2 );

  triangl = flip(triangr);
  ## Plot the autocovariance function
  plot(x(find(x <= 0)), (triangl ... # Plot the left side of the triangle
			 + y1(find(x <= 0))), 'linewidth', 3, 'color', 'b');
  plot(x(find(x > 0)), (triangr ...  # Plot the right side of the triangle
			+ y1(find(x > 0))), 'linewidth', 3, 'color', 'b');
  xlabel('time lag (s)', 'fontname', 'DejaVu Serif');
  string = 'Autocovariance of the intercell interference AEA_2';
  title(string, 'fontname', 'DejaVu Serif');
  set(gca, 'fontsize', 8);
  legend('Simulated autocovariance', 'Theory', 'fontname', 'DejaVu Serif');
  grid on;
  hold off;

endfunction

function [AEAsLOS AEAsRICE AWNtotalsLOS AWNtotalsRICE] ...
	 = condAWN(tkappa, fadingblocksN, cohtime, signallength, ...
		   symbolN, K, fadingmean)
  a = 10 ^ 5; # Scaling for the PPP sampling region [-0.505, 0.505] ^ 2
  ## Initialize the memory for the waveform vector
  AWNtotalsLOS = zeros(fadingblocksN, signallength / symbolN); 
  AWNtotalsRICE = zeros(fadingblocksN, signallength / symbolN);
  ## Initialize the fading gains for the LoS and Rician channels
  AEAsLOS = zeros(fadingblocksN, 1);
  AEAsRICE = zeros(fadingblocksN, 1);  
  refs = planarHPPrefs(tkappa * a / pi); # The Earth transmitter locations
  refs = [refs(1, :) - 0.5; refs(2, :) + 0.5];
  for iii = 1:fadingblocksN
    if(mod(iii, 250) == 0) # Observe the progression
      disp([num2str(iii), "/", num2str(fadingblocksN)]);
    end
    [ithAEAsLOS ithAEAsRICE ithAWNsLOS ithAWNsRice] = ...
      nthAEAsandAWNs(iii, refs, cohtime, signallength / symbolN, K, ...
		     tkappa, fadingmean);
    AEAsRICE(iii) = ithAEAsRICE; # Rice AEA
    AEAsLOS(iii) = ithAEAsLOS; # LoS AEA
    AWNtotalsLOS(iii, 1:signallength / symbolN) = ithAWNsLOS;
    AWNtotalsRICE(iii, 1:signallength / symbolN) = ithAWNsRice;
  end
  ## AEA gains 
  AEAsLOS = AEAsLOS'(:);
  AEAsRICE = AEAsRICE'(:);
  AEAsLOS = repmat(AEAsLOS', signallength, 1)(:);
  AEAsRICE = repmat(AEAsRICE', signallength, 1)(:);
  ## The aggregate AWN signals in a vector
  AWNtotalsLOS = AWNtotalsLOS'(:);
  AWNtotalsRICE = AWNtotalsRICE'(:);
  ## The signals are constant within each symbol, so duplicate the vectors
  AWNtotalsLOS = repmat(AWNtotalsLOS', symbolN, 1)(:); 
  AWNtotalsRICE = repmat(AWNtotalsRICE', symbolN, 1)(:);
endfunction

## Modulate the given digital samples to an analog signal 
function xb = digitaltoanalog(T, digitalsignal)
  xb = zeros(0, length(T)); 
  iii = 1;
  for t = T
    xb(iii) = ... # Modulated analog signal value at t
	sum(digitalsignal .* sinc(t - (0:length(digitalsignal) - 1)));
    iii = iii + 1;
  end
endfunction

## Locations of the interferers
function refs = planarHPPrefs(density)
  yMin = -0.505; yMax = 0.505; # Scaled dimensions of the rectangle
  xMin = -0.505; xMax = 0.505;
  xDelta = xMax - xMin; yDelta = yMax - yMin; # Rectangle side length
  ## Number of points in the area is a Poisson variable of the given
  ## density
  numbPoints = poissrnd(density + density * (1.01 * 1.01 - 1));    
  ## Pick points from uniform distribution
  x = xDelta * (rand(numbPoints, 1)) + xMin;    
  ## Map referencepoints to geographical coordinates
  y = yDelta * (rand(numbPoints, 1)) + yMin; 
  refs = [x'; y'];  
endfunction

## Derive the AEAs and the interference waveforms at nth fading block
function [nthAEAsLOS nthAEAsRICE nthAWNsLOS nthAWNsRICE] ...
	 = nthAEAsandAWNs(n,refs, cohtime, signallength, K, tkappa, ...
			  fadingmean)
  refs = [refs(1, :) + n * cohtime; refs(2, :)]; # Move refs by n steps
  GPrefs = GP(refs, tkappa, fadingmean); # Interferer locations
  ## Determine the amplitude fading parameters for normalized power fading
  s = sqrt(K / (1 + K)); # Noncentrality parameter
  sigma = sqrt(1 / (2 * (1 + K))); # Scale parameter
  if s == 0 # Rayleigh fading case
    fadings = raylrnd(sigma, 1, length(GPrefs));
  else # Else, Rician samples (for some reason, ricernd can not handle s==0)
    fadings = ricernd(s, sigma, 1, length(GPrefs));
  endif
  nthAEAsRICE = sum(sqrt(GPrefs) .* fadings); # Interference AEA
  nthAEAsLOS = sum(sqrt(GPrefs));

  ## Alternatively, return the interference powers:
  ## nthAEAsRICE = sum(GPrefs .* fadings .^ 2); 
  ## nthAEAsLOS = sum(GPrefs);
  
  ## Simulate the AWGN waveform samples in each block
  signals = normrnd(0, 1, signallength, length(GPrefs));
  ## Alternatively, simulate binary white noise samples:
  ##  signals = ...
  ##  discrete_rnd([-1, 1], [0.5,0.5], signallength, length(GPrefs));
  nthAWNsRICE = ...
    sum((sqrt(GPrefs) .* signals .* fadings)')'; # Rice-K interference
  nthAWNsLOS = sum((sqrt(GPrefs) .* signals)')'; # The LoS interference
endfunction

## Produces a realization of the GP and sum the signals given refs
function GPrefs = GP(refs, tkappa, fadingmean)
  a = 10 ^ 5; # A scaling factor corresponding to refs in [-0.5 0.5] ^ 2
  GPrefs = exp(-a * norm(refs, 2, "cols") .^ 2);
  ## Here put a footprint restriction (for intercell, GPrefs < 0.5)
  GPrefs = GPrefs(find(GPrefs < 0.5));  
  ## ... or take only the nearest transmitter
  ## (here, the scaling corresponds to the Rice fading):
  ##  GPrefs = [0 max(GPrefs) * (1 + 4 * tkappa * fadingmean ^ 2)]; 
endfunction
