 % Open a progress bar   
cpb = ConsoleProgressBar();
    
    
% Set progress bar parameters
cpb.setLeftMargin(4);   % progress bar left margin
cpb.setTopMargin(1);    % rows margin

cpb.setLength(40);      % progress bar length: [.....]
cpb.setMinimum(0);      % minimum value of progress range [min max]
cpb.setMaximum(N);      % maximum value of progress range [min max]

cpb.setElapsedTimeVisible(1);
cpb.setRemainedTimeVisible(1);

cpb.setElapsedTimePosition('left');
cpb.setRemainedTimePosition('right');