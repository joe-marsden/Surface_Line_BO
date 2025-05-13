%% Third Optimisation - simulation
% Joe Marsden last updated 9/5/24

% This script will attempt to optimise the discovery of a line of solutions
% on a surface

% Here different grid sizes will be compared to see how long they take to
% converge

% Stopping criteria - a change in distance of less than 40 for 3
% consecutive iterations

%% Define the bounds

flowbounds = [0.5, 2.5];
filtbounds = [1, 60];

aa = 0.3403; % surface fitting values - defining aim
bb = 13.92;
cc = 2.41;
dd = 3.056;

%% Iterate through different grid sizes

% Data to be stored
requirediterations = [];

maxgridsize = 10;

for gridsize = 3  % Change to 1:gridsize to loop

    %% Define grid

    if gridsize == 1 % one is a special case
        flowincs = mean(flowbounds);
        filtincs = mean(filtbounds);
    else
        gridincs = linspace(0,1,gridsize);
        flowincs = flowbounds(1) + (flowbounds(2)- flowbounds(1)) * gridincs;
        filtincs = filtbounds(1) + (filtbounds(2)- filtbounds(1)) * gridincs;
    end

    flowlist = zeros(gridsize^2,1); 
    filtlist = zeros(gridsize^2,1); 

    % Generate list of conditions
    for i = 1:gridsize
        for j = 1:gridsize
            flowlist((i-1)*gridsize + j) = flowincs(i);
            filtlist((i-1)*gridsize + j) = filtincs(j);
        end
    end

    % Sample surface with initial data
    initdiss = sampleSurface(flowlist,filtlist);
    initdissadj = sampleSurface(flowlist+0.01,filtlist+0.01); % Adjacent points for case that algorithm is given single initial point

    %% Set up the optimiser

    nquant = 2; % number of continuous variables (res time, temp, equiv)
    nqual = 1; % number of discrete variables
    bounds = [flowbounds',filtbounds',[1; 1]]; % columns are lower and upper bounds
    
    levels = [1]; % number of levels of the qualitative variable
    dim_qual = [3]; % position of first qualitative variable
    
    if gridsize == 1 % Algorithm doesn't like single point input
        optimiser = LVBayesianOptimiser('AEI', bounds, [[flowlist+0.01,filtlist+0.01,ones(length(flowlist),1)];[flowlist,filtlist,ones(length(flowlist),1)]], [-initdissadj;-initdiss], dim_qual, levels); % creating the optimiser object
    else
        optimiser = LVBayesianOptimiser('AEI', bounds, [flowlist,filtlist,ones(length(flowlist),1)], -initdiss, dim_qual, levels); % creating the optimiser object
    end
    %% Create iteration while loop
    
    iteration = 0
    stopcondition = "not met";
    flowvalues = linspace(flowbounds(1),flowbounds(2),100)';      

    % Known curve
    filts2 = (bb.*exp(-cc .* flowvalues) + dd) .^ (1/(1-aa));
    filts3 = ((bb*exp(-cc .* flowvalues) + dd) .* (2 * flowvalues)).^(1/(1-aa));

    known1 = []; % initialising variables
    known2 = [];
    unknown1 = [];
    unknown2 = [];
    dist1 = [];
    dist2 = [];        
    change1 = [];
    change2 = [];
    stopcheck1 = [];
    stopcheck2 = [];

    while stopcondition == "not met"
        iteration = iteration + 1
        
        filtsol = []; % initialising variables to be reset each iteration
        filtsol2 = [];
    
        % Initial values
        x0 = 20;
        x1 = 20;

        for i = 1:length(flowvalues)
            % for each flow, find filt such that gpPredict - filt = 0
            fun =  @(x)abs(-optimiser.revertY(gpPredict(optimiser.mdl,[flowvalues(i),x,1]))-x);
            % also find filt such that gpPredict*flow - 0.5*filt = 0
            fun2 = @(x)abs(-optimiser.revertY(gpPredict(optimiser.mdl,[flowvalues(i),x,1]))*flowvalues(i)-0.5*x);
            % Solve functions
            filtsol(i,1) = fmincon(fun,x0,[],[],[],[],filtbounds(1),filtbounds(2));
            filtsol2(i,1) = fmincon(fun2,x1,[],[],[],[],filtbounds(1),filtbounds(2));
            % start next iteration from last point to encourage smooth
            % line
            x0 = filtsol(i,1);
            x1 = filtsol2(i,1);
        end

        % Prepare for distance measurements
        known1 = [known1, filts2];
        known2 = [known2, filts3];
        unknown1 = [unknown1, filtsol];
        unknown2 = [unknown2, filtsol2];
   
        if iteration >= 2

            % Distance between known and solved point
            dist1(1,iteration) = computedistance(unknown1(:,iteration),unknown1(:,iteration -1));
            dist2(1,iteration) = computedistance(unknown2(:,iteration),unknown2(:,iteration -1));

            % Change in distance
            change1(iteration-1) = dist1(iteration) - dist1(iteration - 1);
            change2(iteration-1) = dist2(iteration) - dist2(iteration - 1);

            % Check if distance below threshold value
            if abs(change1(iteration-1)) < 40
                stopcheck1(iteration - 1) = 1
            else
                stopcheck1(iteration - 1) = 0
            end            
            if abs(change2(iteration-1)) < 40
                stopcheck2(iteration - 1) = 1
            else
                stopcheck2(iteration - 1) = 0
            end
        end

        % Find mean and uncertainty (a,b) along the two lines
        [a,b] = gpPredict(optimiser.mdl,[flowvalues,filtsol,ones(size(flowvalues))]);
        [c,d] = gpPredict(optimiser.mdl,[flowvalues,filtsol2,ones(size(flowvalues))]);

        % Assess the lines
        
        % Find the point on the line with the greatest uncertainty
        [mx, idx] = max(diag(b));
        [mx2, idx2] = max(diag(d));
    
        % Check which line has greater uncertainty
        if mx2 > mx
            nextflow = flowvalues(idx2);
            nextfilt = filtsol2(idx2);
        else
            nextflow = flowvalues(idx);
            nextfilt = filtsol(idx);
        end

        % Sample the surface
        nextdiss = sampleSurface(nextflow,nextfilt);
    
        % Add to the optimiser 
        optimiser = optimiser.addData([nextflow, nextfilt, 1],-nextdiss);

        % Check for stop condition
        % Last three points were below required distance
        if (length(stopcheck1) >= 3) && (isequal(stopcheck1(end-2:end),[1,1,1])) && (isequal(stopcheck2(end-2:end),[1,1,1]))
            stopcondition = "met";
        end

        % Store number of iterations
        requirediterations(gridsize) = iteration;

        % Make a tiled figure plotting the uncertainty
        h(iteration) = figure;
        t = tiledlayout(2,4);
        t.Padding = 'compact';
        t.TileSpacing = "compact";
        set(gcf,'color','w');

        nexttile([2,2]) % Square figure

        set(0, 'defaultAxesFontSize', 16);
        set(0, 'DefaultLineLineWidth', 2);
        set(0, 'DefaultAxesLineWidth', 2);
        set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 34, 16])
        set(gcf,'color','w');
        set(gca, 'Layer', 'top')
        box on

        hold on

        % Define lines
        filts2 = (bb.*exp(-cc .* flowvalues) + dd) .^ (1/(1-aa));
        filts3 = ((bb*exp(-cc .* flowvalues) + dd) .* (2 * flowvalues)).^(1/(1-aa));

        % Plot known lines
        plot(flowvalues,filts3,Color='#7500D6',LineWidth=0.5)
        plot(flowvalues,filts2,Color='#7500D6',LineWidth=0.5)

        % Plot solved lines
        plot(flowvalues,filtsol,'-k',LineWidth=4)
        plot(flowvalues,filtsol2,'-k',LineWidth=4)

        % Plot initial and iterated data from optimiser
        if gridsize>1
            plot(optimiser.X(1:gridsize^2,1),optimiser.X(1:gridsize^2,2),'o',Color='#D600BD',markersize=8)
            plot(optimiser.X(gridsize^2+1:end,1),optimiser.X(gridsize^2+1:end,2),'x',Color='#D600BD',markersize=15)
        else
            plot(optimiser.X(1:gridsize^2,1),optimiser.X(1:gridsize^2,2),'o',Color='#D600BD',markersize=8)
            plot(optimiser.X(gridsize^2+2:end,1),optimiser.X(gridsize^2+2:end,2),'x',Color='#D600BD',markersize=15)
        end

        strgrid = sprintf('%d',gridsize^2+iteration); % add experiment label 
        annotation('textbox',[0.1 0.9 1 0],'String',strgrid,'FontSize',20,EdgeColor='none')

        ylabel('Filtration Time / min',FontWeight='bold');
        ylim([0 60])
        xlabel(t,'Solvent Flow / mL min^{-1}',FontWeight='bold',fontsize=16);
        xlim([0.5 2.5])

        hold off
        
        nexttile([1 2])
        % Top uncertainty
        hold on
        box on
        set(gca, 'YScale', 'log')

        ylabel('Line 1 Uncertainty',FontWeight='bold');
        ylim([0.0001 1])

        plot(flowvalues,diag(d),color='#7500D6')

        if mx2 > mx % if most recent sample point was top
            xline(nextflow,'r',LineWidth=2) % show this
        end

        hold off
        nexttile([1 2])
        % Bottom uncertainty

        hold on
        box on
        set(gca, 'YScale', 'log')

        ylabel('Line 2 Uncertainty',FontWeight='bold');
        ylim([0.0001 1])

        if mx2 < mx % if most recent was bottom
            xline(nextflow,'r',LineWidth=2)
        end

        plot(flowvalues,diag(b),color='#7500D6')
        
        hold off

        filename = sprintf('Gridsize %d result %d',gridsize,iteration);
        saveas(h(iteration),filename,'jpeg')
    end

    % Stop condition met - spit out a graph

    f(gridsize) = figure;
    
    set(0, 'defaultAxesFontSize', 16);
    set(0, 'DefaultLineLineWidth', 2);
    set(0, 'DefaultAxesLineWidth', 2);
    set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 17, 16])
    set(gcf,'color','w');
    set(gca, 'Layer', 'top')
    box on

    hold on
    
    % Known lines
    filts2 = (bb.*exp(-cc .* flowvalues) + dd) .^ (1/(1-aa));
    filts3 = ((bb*exp(-cc .* flowvalues) + dd) .* (2 * flowvalues)).^(1/(1-aa));
    plot(flowvalues,filts3,Color='#7500D6',LineWidth=0.5)
    plot(flowvalues,filts2,Color='#7500D6',LineWidth=0.5)

    % solved lines
    plot(flowvalues,filtsol,'-k',LineWidth=4)
    plot(flowvalues,filtsol2,'-k',LineWidth=4)

    % algorithm data
    if gridsize>1
        plot(optimiser.X(1:gridsize^2,1),optimiser.X(1:gridsize^2,2),'o',Color='#D600BD',markersize=8)
        plot(optimiser.X(gridsize^2+1:end,1),optimiser.X(gridsize^2+1:end,2),'x',Color='#D600BD',markersize=15)
    else
        plot(optimiser.X(1:gridsize^2,1),optimiser.X(1:gridsize^2,2),'o',Color='#D600BD',markersize=8)
        plot(optimiser.X(gridsize^2+2:end,1),optimiser.X(gridsize^2+2:end,2),'x',Color='#D600BD',markersize=15)
    end

    % annotate with grid size
    strgrid = sprintf('%d x %d',gridsize,gridsize);
    annotation('textbox',[0.2 0.9 1 0],'String',strgrid,'FontSize',20,EdgeColor='none')

    ylabel('Filtration Time / min',FontWeight='bold');
    ylim([0 60])
    xlabel('Solvent Flow / mL min^{-1}',FontWeight='bold');
    xlim([0.5 2.5])
    hold off
    
    filename = sprintf('Gridsize %d result',gridsize);
    saveas(f(gridsize),filename,'jpeg')

    % Save optimisers 
    alloptimisers(gridsize) = optimiser;
end