%%Alice Henman%%
%%Evolutionary Branching in Models for Symmetric and Asymmetric
%%Competition%%

%%Default graphing layout from Denis
set(0,'defaultTextFontName', 'Arial');
set(0,'defaultaxesfontsize', 25); % 25 for 1X3, 20 for 1X2 figures
%set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','none');
set(0,'defaulttextinterpreter','none');
set(0,'defaultAxesXGrid','off');
set(0,'defaultAxesYGrid','off');
set(0,'defaultAxesTickDir','out');
set(0,'defaultAxesLineWidth',1.5);

C = orderedcolors("gem");

%% functions

%(Mexican hat) Wavelet Competition function
function alpha=alpha(x,y,sigma,beta)
    alpha = exp(sigma.^2.*beta.^2/2).* ...
(1.-((x-y+sigma.^2.*beta)./sigma).^2).*exp(-(x-y+sigma.^2.*beta).^2./(2*sigma^2));
end

%dalpha(x-y)/dx, dalpha(y,x,...) = dalpha(y-x)/dy = -dalpha(y-x)/dx
function dalpha=dalpha(x,y,sigma,beta)
    dalpha = exp(sigma.^2.*beta.^2/2).* ...
(-3+((x-y+sigma.^2.*beta)./sigma).^2).*exp(-(x-y+sigma.^2.*beta).^2./(2*sigma^2)).*(x-y+sigma.^2.*beta)./sigma.^2;
end
%d^2alpha(x-y)/dx^2, ddalpha(y,x,...) = d^2alpha(y-x)/dy^2
function ddalpha=ddalpha(x,y,sigma,beta)
A = (x-y+sigma^2*beta);
    ddalpha = exp(sigma.^2.*beta.^2/2-A.^2./(2*sigma^2)).* ...
((-3+(A./sigma).^2)+ ...
(A./sigma^2).^2.*(5-(A./sigma).^2));
end

%Gaussian centred on x0
function Gauss=G(sigma,x0,x)
    Gauss = exp(-(x-x0).^2./(2*sigma^2));
end
%derivative of gaussian
function dGauss=dG(sigma,x0,x)
    dGauss = -(x-x0)./(sigma^2).*exp(-(x-x0).^2./(2*sigma^2));
end
%second derivative of gaussian
function ddGauss=ddG(sigma,x0,x)
    ddGauss = 1/(sigma^2)*(-1+(x-x0).^2/(sigma^2)).*exp(-(x-x0).^2./(2*sigma^2));
end
%K is sum of 2 gaussians
function K=K(C1,C2,sigma1,sigma2,x0,x)
    K = C1*G(sigma1,x0,x)-C2*G(sigma2,x0,x);
end
%derivative of K (when alpha peaks at y=x, the solutions to this are x*)
function dK=dK(C1,C2,sigma1,sigma2,x0,x)
    dK = C1*dG(sigma1,x0,x)-C2*dG(sigma2,x0,x);
end

%second derivative of K
function ddK=ddK(C1,C2,sigma1,sigma2,x0,x)
    ddK = C1*ddG(sigma1,x0,x)-C2*ddG(sigma2,x0,x);
end

%% Parameters
C1 = 1;%2.5;%5*rand;
C2 = 0;%2;%C1*rand;
sigma = 0.5;%1;
sigma1 = 0.5;%0.5*rand;
sigma2 = 0.35;%sigma1*sqrt(C2/C1)*rand; % condition for there to be 2 humps
x0 = 1.5;%10;%1+5*rand;
beta=0;
r=rand+0.001;
m=0.01; %max size of mutation
epsM = 1e-6; %fraction of population that mutates
xmin=0;
xmax=3;

%% Equilibria
N1 = @(x) K(C1,C2,sigma1,sigma2,x0,x)/alpha(0,0,sigma,beta); %1 resident equilibrium
%equilibrium for 2 residents 
%The equilibrium pop density of resident with trait x1, swap x1 and x2 to
%get pop dens of resident with trait x2
function N = N(x1,C1,C2,sigma1,sigma2,x0,sigma,beta,x2)
    N = (K(C1,C2,sigma1,sigma2,x0,x1).*alpha(0,0,sigma,beta)- ...
        K(C1,C2,sigma1,sigma2,x0,x2).*alpha(x1,x2,sigma,beta))./( ...
        alpha(0,0,sigma,beta).^2-alpha(x1,x2,sigma,beta).*alpha(x2,x1,sigma,beta));
end
Nfun = @(x,xb) N(x,C1,C2,sigma1,sigma2,x0,sigma,beta,xb);

%popdens vector for n traits, traits must be a column vector [x1;x2;...;xn]
function Nn = Nn(traits,C1,C2,sigma1,sigma2,x0,sigma,beta)

    A = alpha(traits,traits.',sigma,beta);%vector Aij=alpha(xi-xj)
    Kvec = K(C1,C2,sigma1,sigma2,x0,traits);
    Kvec = Kvec(:);%turns Kvec into column vector

    Nn = A\Kvec; %finds vector N=[N1;N2;...Nn] by solving K=AN
end
%% Invasion fitness
f = @(x,y) r*(1-alpha(y,x,sigma,beta).*K(C1,C2,sigma1,sigma2,x0,x) ...
    ./(alpha(0,0,sigma,beta).*K(C1,C2,sigma1,sigma2,x0,y)));

%fitness for 2 residents
function fit12 =fit12(r,sigma,beta,C1,C2,sigma1,sigma2,x0,y,x1,x2)
    fit12 = r*(1-(alpha(y,x1,sigma,beta).*N(x1,C1,C2,sigma1,sigma2,x0,sigma,beta,x2)+ ...
        alpha(y,x2,sigma,beta).*N(x2,C1,C2,sigma1,sigma2,x0,sigma,beta,x1))./K(C1,C2,sigma1,sigma2,x0,y));
end
f2 = @(x1,x2,y) fit12(r,sigma,beta,C1,C2,sigma1,sigma2,x0,y,x1,x2);

% fitness for n residents
function fitn=fitn(r,sigma,beta,C1,C2,sigma1,sigma2,x0,Nvec,y,traits)
    A = alpha(y,traits,sigma,beta);
    Ky = K(C1,C2,sigma1,sigma2,x0,y);
    fitn = r*(1-sum(A.*Nvec)./Ky);
end
fn = @(Nvec,y,traits) fitn(r,sigma,beta,C1,C2,sigma1,sigma2,x0,Nvec,y,traits);

%% Selection gradient
%selection gradient for branch with resident x1
function secgrad1 = secgrad1(r,sigma,beta,C1,C2,sigma1,sigma2,x0,x1,x2)
    secgrad1 = r*(-(dalpha(0,0,sigma,beta).*N(x1,C1,C2,sigma1,sigma2,x0,sigma,beta,x2)+ ...
        dalpha(x1,x2,sigma,beta).*N(x2,C1,C2,sigma1,sigma2,x0,sigma,beta,x1))./K(C1,C2,sigma1,sigma2,x0,x1)+ ...
        dK(C1,C2,sigma1,sigma2,x0,x1).*(alpha(0,0,sigma,beta).*N(x1,C1,C2,sigma1,sigma2,x0,sigma,beta,x2)+ ...
        alpha(x1,x2,sigma,beta).*N(x2,C1,C2,sigma1,sigma2,x0,sigma,beta,x1))./K(C1,C2,sigma1,sigma2,x0,x1).^2);
end

%selection gradient for branch with resident x2
function secgrad2 = secgrad2(r,sigma,beta,C1,C2,sigma1,sigma2,x0,x1,x2)
    secgrad2 = r*(-(dalpha(x2,x1,sigma,beta).*N(x1,C1,C2,sigma1,sigma2,x0,sigma,beta,x2)+ ...
        dalpha(0,0,sigma,beta).*N(x2,C1,C2,sigma1,sigma2,x0,sigma,beta,x1))./K(C1,C2,sigma1,sigma2,x0,x2)+ ...
        dK(C1,C2,sigma1,sigma2,x0,x2).*(alpha(x2,x1,sigma,beta).*N(x1,C1,C2,sigma1,sigma2,x0,sigma,beta,x2)+ ...
        alpha(0,0,sigma,beta).*N(x2,C1,C2,sigma1,sigma2,x0,sigma,beta,x1))./K(C1,C2,sigma1,sigma2,x0,x2).^2);
end

sx1 = @(x1,x2) secgrad1(r,sigma,beta,C1,C2,sigma1,sigma2,x0,x1,x2);
sx2 = @(x1,x2) secgrad2(r,sigma,beta,C1,C2,sigma1,sigma2,x0,x1,x2);
%Evolutionary singular strategies
x1=x0;
%x2=x0+sqrt(2*sigma1^2*sigma2^2*log(C1*sigma2^2/(C2*sigma1^2))/(sigma2^2-sigma1^2));
%x3=x0-sqrt(2*sigma1^2*sigma2^2*log(C1*sigma2^2/(C2*sigma1^2))/(sigma2^2-sigma1^2));

%% Numerical Simulation (up to 2 branches)

%{
numberofrealisations = 1; % number of sample paths
T = 2000; % end time for simulation
xinitial = 1; %initial resident trait value
% matrix to hold trait value at each timestep, for each realisation
xplot = zeros(T+1,numberofrealisations);
Nplot = zeros(T+1,numberofrealisations);
Nbplot = zeros(T+1,numberofrealisations);
xbplot = zeros(T+1,numberofrealisations);
timeplot = zeros(T+1,numberofrealisations);
for i = 1:numberofrealisations %Run numberofrealisations simulations
    x=xinitial;
    xb=xinitial; %resident of the branch
    Nr=N1(x);
    Nb=Nr;        %initially there is no branch population
    time=0;
    j=1;
    xplot(j,i)=x;       %To keep track of trait values for simulation i at "time" j
    xbplot(j,i)=xb;     %To keep track of branched trait values for simulation i at "time" j
    Nplot(j,i)=Nr;       %To keep track of population density for simulation i at "time" j
    Nbplot(j,i)=Nb;
    timeplot(j,i)=0;    %To keep track of time that trait x is obtained for simulation i
    branch = false;     %means the trait hasn't branched yet

    while time<T        %goes until specified time T
        %only runs for the trait value sufficiently far from x*
        if branch==false
            y = x+2*m*rand-m; %mutant value is resident +-random number <= m
            if f(x,y) > 0
                if f(y,x)>0
                    xb=y; %update branch value if mutally invasive
                    Nr=Nfun(x,xb);  %branch point affects both equilibrium densitities
                    Nb=Nfun(xb,x);
                    branch=true;
                else
                    x=y; %update trait value if y has positive invasion fitness
                    xb=x;
                    Nr=N1(x); %update new population
                    Nb=Nr;
                end
            end
        %This runs once we have got sufficiently close to "branching" point
        else
            if rand<=Nr/(Nr+Nb) %pick which subpopulation is mutating weighted by its frequency in the population
                y = x+2*m*rand-m;
                if f2(x,xb,y)>0
                    x=y;
                    Nr=Nfun(x,xb);
                    Nb=Nfun(xb,x);
                    %need to check mutations don't push traits outside of
                    %mutualinvasibility zone
                    if Nr<=0 %branch with trait x dies, only xb remains
                        x=xb;
                        Nr=N1(x);
                        Nb=Nr;
                        branch = false;
                    elseif Nb<=0
                        xb=x;
                        Nr=N1(x);
                        Nb=Nr;
                        branch = false;
                    end
                end
            else
                yb = xb+2*m*rand-m;
                if f2(x,xb,yb)>0
                    xb=yb;
                    Nr=Nfun(x,xb);
                    Nb=Nfun(xb,x);
                    if Nr<=0 %branch with trait x dies, only xb remains
                        x=xb;
                        Nr=N1(x);
                        Nb=Nr;
                        branch = false;
                    elseif Nb<=0
                        xb=x;
                        Nr=N1(x);
                        Nb=Nr;
                        branch = false;
                    end
                end
            end
        end
        time=time+1;%update timestep
        %update matrices
        j=j+1;
        xplot(j,i)=x;
        xbplot(j,i)=xb;
        Nplot(j,i)=Nr;
        Nbplot(j,i)=Nb;
        timeplot(j,i)=time;
    end
end

%Plotting
figure;
C = orderedcolors("gem");
%{
%green colour theme for poster
C=[0.757 0.820 0.122
    0.596 0.788 0.169
    0.431 0.753 0.027
    0.318 0.694 0.020
    0.204 0.635 0.012
    0.224 0.580 0.027
    0.239 0.525 0.043
    0.122 0.455 0.051
    0.000 0.380 0.055];
%}
hold on;
%plot trait value by time for each realisation
for i=1:numberofrealisations
    h=stairs(timeplot(:,i),xplot(:,i),'Color',C(mod(i,7)+1,:));
    k=stairs(timeplot(:,i),xbplot(:,i),'Color',C(mod(i,7)+1,:));
    set(h,'Linewidth',2);
    set(k,'Linewidth',2);
end
mean(mean(abs(xplot(end-200:end,:)-xbplot(end-200:end,:)))) %finds mean separation between branches
hold on;
%plots singular strategies
hx1=yline(x1,'--',LineWidth=2);
%hx2=plot([0,T],x2*ones(1,2),'r-.',LineWidth=2);
%hx3=plot([0,T],x3*ones(1,2),'b--',LineWidth=2);
legend(Location="southeast")
%legend([hx2 hx1 hx3], {sprintf('x_0+x_1 = %.2f', x2), sprintf('x_0 = %.2f', x1), sprintf('x_0-x_1 = %.2f', x3)});
%legend([hx1 hx3], {sprintf('x_0   = Repeller'), sprintf('x_0-x_1 = Branch point')});
legend(hx1,sprintf('x* = %.2f', x1))
xlabel('Time','FontSize',25);
ylabel('Trait value','FontSize',25);
axis([0 T 0 max(x0+1,xinitial)]);
set(gca,'linewidth',1.5);
set(gca,'FontSize',20);
grid on;
axis square;

%pop dens
figure;
Cb = orderedcolors("glow");
hold on;
%plots population density of each branch over time
for i=1:numberofrealisations
    h=stairs(timeplot(:,i),Nplot(:,i),'Color',C(mod(i,7)+1,:));
    k=stairs(timeplot(:,i),Nbplot(:,i),'Color',Cb(mod(i,7)+1,:));
    set(h,'Linewidth',2);
    set(k,'Linewidth',2);
end
hold on;
xlabel('Time','FontSize',25);
ylabel('Population Density','FontSize',25);
set(gca,'linewidth',1.5);
set(gca,'FontSize',20);
grid on;
axis square;

%total pop dens
figure;
hold on;
%plots population density of each branch over time
for i=1:numberofrealisations
    h=stairs(timeplot(:,i),Nplot(:,i)+Nbplot(:,i),'Color',C(mod(i,7)+1,:));
    set(h,'Linewidth',2);
end
hold on;
xlabel('Time','FontSize',25);
ylabel('Population Density','FontSize',25);
set(gca,'linewidth',1.5);
set(gca,'FontSize',20);
grid on;
axis square;

%trait evolution on a TEP

z=linspace(xmin,xmax,2000);
[X1,X2] = meshgrid(z, z);

F1 = f(X1,X2); %x1 resident, x2 mutant
F2 = f(X2,X1); %x2 resident, x1 mutant
S1 = sx1(X1,X2);
S2 = sx2(X1,X2);
%selection gradient
sx = @(x) r*(-dalpha(0,0,sigma,beta0)./alpha(0,0,sigma,beta0)+dK(C1,C2,sigma1,sigma2,x0,x) ...
    ./K(C1,C2,sigma1,sigma2,x0,x));
xstar=roots(chebfun(sx,[1e-6 5])); %singular strategies

F = zeros(size(X1));
F(F1>0)=1;   %fx(y)>0
F(F2>0)=2;   %fy(x)>0
F(F1>0&F2>0)=3;   %mutual invasibility

%don't care about direction outside of coexistence region
%S1(f(X1,X2)<0|f(X2,X1)<0)=0;
%S2(f(X1,X2)<0|f(X2,X1)<0)=0;

figure;
hold on
% TEP
contourf(X1,X2,F,'LineStyle','none');
contour(X1,X2,F,[0.5 1.5 2.5],'k','LineWidth',2); %black line for f=0
colormap([1 1 1; 0.8 0.9 0.9; 1 0.9 0.7; 0.6 0.9 0.3]); %blue for fx(y)>0, yellow for fy(x)>0, green for both
xlabel('$\textsf{Resident trait value}$ ($x_1$)','Interpreter','latex','FontSize',20);
ylabel('$\textsf{Resident trait value}$ ($x_2$)','Interpreter','latex','FontSize',20);
axis equal
%selection gradient nullclines
contour(X1,X2,S1,[0 0],'r:','LineWidth',2,'DisplayName','s_{x_1,x_2}(x_1)=0');
contour(X1,X2,S2,[0 0],'b:','LineWidth',2,'DisplayName','s_{x_1,x_2}(x_2)=0');
%singular strategies
xline(xstar(1),'LineStyle','--','LineWidth',1.5,'Color',C(3,:),'DisplayName',['B1=',num2str(xstar(1),3)])
%xline(xstar(2),'LineStyle','--','LineWidth',1.5,'Color',C(2,:),'DisplayName',['D=',num2str(xstar(2),3)])
xline(xstar(end),'LineStyle','--','LineWidth',1.5,'Color',C(1,:),'DisplayName',['B2=',num2str(xstar(end),3)])

hold on;
%plots population density of each branch over time
for i=1:numberofrealisations
    h=stairs(xplot(:,i),xbplot(:,i),'Color',C(mod(i,7)+1,:));
    set(h,'Linewidth',2);
end
hold on;

legend("off");
axis equal
box on
%}
%% Numerical Simulation (up to n branches)

numberofrealisations = 5; % number of sample paths
T = 100000; % end time for simulation
xinitial = 1; %initial resident trait value
% matrix to hold trait values at each timestep, for each realisation for
% each branch
traitHistories = cell(numberofrealisations,1);%starts  with 1 for 1 branch
densityHistories = cell(numberofrealisations,1);
timeHistories = cell(numberofrealisations,1);
for i = 1:numberofrealisations %Run numberofrealisations simulations
    traits=xinitial; %vector column of all current traits
    Nlist=N1(xinitial); %vector column of all current densities
    time=0;
    j=1;
    branch = 1; %initially 1 branch only
    
    %list of all traits and densities for all time
    xplot = nan(T+1,1); %xplot(1,1)=trait value of branch 1 at time 0
    Nplot = nan(T+1,1); %Nplot(T,5)=trait value of branch 5 at time T-1
    %keeps track of time each trait/density held at
    timeplot = nan(T+1,1);

    xplot(j,1)=traits;       %To keep track of trait values at "time" j
    Nplot(j,1)=Nlist;       %To keep track of population densities at "time" j
    timeplot(j)=time;    %To keep track of time that trait x is obtained
    
    while time<T        %goes until specified time T
        %only runs for the trait value sufficiently far from x*
        k = find(rand<=cumsum(Nlist./sum(Nlist)),1);%adds each density to fraction until it passes rand, this is the branch that mutates
        x=traits(k);
        y = x+m*(2*rand-1);
        if fn(Nlist,y,traits)>0 %if y has a positive invasion fitness in environment of residents
            swappedtrait = traits;
            swappedtrait(k) = y; %replaces x with y
            Nswapped = Nn(swappedtrait,C1,C2,sigma1,sigma2,x0,sigma,beta);
            if fn(Nswapped,x,swappedtrait)>0 %if x has +ve f when mutant is resident
                %we have mutual invasibility, add y to list of branches
                traits = [traits;y];
                branch = numel(traits);
                if branch > size(xplot,2) %only add a column if needed, don't wipe an existing one
                    xplot(:,branch)=nan;
                    Nplot(:,branch)=nan;
                end
            else %not mutual invasibility
                traits = swappedtrait;%y replaces x
            end
        %else x doesn't change
        end
        Nnew = Nn(traits,C1,C2,sigma1,sigma2,x0,sigma,beta);%find new densities
        extinct = Nnew<=0;%checks if any branches died
        while any(extinct)
            traits(extinct) = [];%remove extinct branches
            Nnew = Nn(traits,C1,C2,sigma1,sigma2,x0,sigma,beta);%find new densities
            extinct = Nnew<=0;%check again
        end
        Nlist = Nnew;
        time=time+1%update timestep
        %update matrices
        j=j+1;
        xplot(j,:)=nan;
        Nplot(j,:)=nan;
        branch=numel(traits);
        for b=1:branch
            xplot(j,b)=traits(b);
            Nplot(j,b) = Nlist(b);
        end
        timeplot(j)=time;
    end
    traitHistories{i} = xplot;
    densityHistories{i} = Nplot;
    timeHistories{i} = timeplot;
end

%Plotting
figure;
hold on;
%plot trait value by time for each branch
for i=1:numberofrealisations
    traitHistory = traitHistories{i};
    timeHistory = timeHistories{i};
    colour = C(mod(i,size(C,1))+1,:);
    for b = 1:size(traitHistory,2)
        plot(timeHistory,traitHistory(:,b),'.','Color',colour,'LineWidth',1);
    end
end
hold on;
%plots singular strategies
hx1=yline(x0,'--',LineWidth=2);
legend(Location="southeast")
legend(hx1,sprintf('x* = %.2f', x1))
xlabel('Time','FontSize',25);
ylabel('Trait value','FontSize',25);
%axis([0 T 0 max(x0+1,xinitial)]);
set(gca,'linewidth',1.5);
set(gca,'FontSize',20);
grid on;
axis square;
%}

%% ChatGPT sim

%{
% PARAMETERS FOR THE SIMULATION

numberofrealisations = 3;
m = 0.05;                 % maximum mutation step
T = 400;                % total evolutionary time
epsM = 1e-6;              % initial mutant density relative to parent
TecoStart = 100;          % initial ecological integration time
tolEco = 1e-10;           % tolerance for ecological equilibrium
maxEcoExtensions = 10;   % maximum times to extend ecological integration
extinctionTol = 1e-8;     % density below this is treated as extinct

% ODE OPTIONS
% AbsTol is made much smaller than the initial mutant density
odeOptionsBase = odeset( ...
    'RelTol',1e-9, ...
    'AbsTol',1e-12);
% STORAGE FOR ALL REALISATIONS
traitHistories = cell(numberofrealisations,1);
densityHistories = cell(numberofrealisations,1);
timeHistories = cell(numberofrealisations,1);
parentHistories = cell(numberofrealisations,1);
birthTimeHistories = cell(numberofrealisations,1);
finalTraitValues = cell(numberofrealisations,1);
finalDensityValues = cell(numberofrealisations,1);

% RUN ALL REALISATIONS
for realisation = 1:numberofrealisations
    % INITIAL POPULATION
    % Start with one resident
    traits = xinitial;

    % Start the ecological dynamics from a positive density.
    % It does not have to be the exact equilibrium because ode45 will
    % relax the population to equilibrium.
    initialDensity = K(C1,C2,sigma1,sigma2,x0,xinitial);
    initialDensities = initialDensity;

    % Find the actual ecological equilibrium using ode45
    densities = getEquilibriumDensities( ...
        traits,...
        initialDensities,...
        C1,C2,sigma1,sigma2,x0,sigma,beta,r,...
        TecoStart,tolEco,maxEcoExtensions,odeOptionsBase);

    % BRANCH IDs
    % The original resident is branch 1
    branchID = 1;
    % The next new branch will receive this ID
    nextID = 2;
    % Parent of every branch
    % parentOf(1) = 0 because branch 1 is the original resident
    parentOf = 0;
    % Evolutionary time at which each branch first appeared
    birthTime = 0;

    % HISTORY ARRAYS
    % Maximum possible number of evolutionary time steps
    maxSteps = ceil(T) + 1;
    % Start with one column because initially there is one branch.
    % New columns are automatically added when new branches appear.
    traitHistory = nan(maxSteps,1);
    densityHistory = nan(maxSteps,1);
    timeHistory = nan(maxSteps,1);

    % RECORD INITIAL STATE
    step = 1;
    time = 0;
    traitHistory(step,branchID) = traits;
    densityHistory(step,branchID) = densities;
    timeHistory(step) = time;

    % EVOLUTIONARY SIMULATION
    while time < T
        % CHOOSE WHICH EXISTING BRANCH MUTATES
        % Mutation is more likely to come from a branch with a larger
        % population.
        mutationProbability = densities ./ sum(densities);
        cumulativeProbability = cumsum(mutationProbability);
        parentIndex = find(rand <= cumulativeProbability,1);

        % Trait of the branch that produced the mutant
        parentTrait = traits(parentIndex);
        parentBranchID = branchID(parentIndex);

        % GENERATE MUTANT
        % Mutant lies uniformly within +/- m of parent trait
        mutantTrait = parentTrait + 2*m*rand - m;

        % CALCULATE INVASION FITNESS
        invasionFitnessValue = calcInvasionFitness( ...
            mutantTrait,...
            traits,...
            densities,...
            C1,C2,sigma1,sigma2,x0,sigma,beta,r);

        % MUTANT CAN ONLY ESTABLISH IF ITS INVASION FITNESS IS POSITIVE
        if invasionFitnessValue > 0
            % TEMPORARILY ADD MUTANT TO THE COMMUNITY
            % Give the mutant a new permanent ID
            newBranchID = nextID;
            % Add mutant to the list of resident traits
            trialTraits = [traits; mutantTrait];
            % Existing branches retain their current densities.
            %
            % Mutant starts at a very small density.
            mutantInitialDensity = epsM * densities(parentIndex);
            trialDensities0 = [densities; mutantInitialDensity];
            % IDs corresponding to trialTraits
            trialBranchID = [branchID; newBranchID];

            % RUN ECOLOGICAL DYNAMICS WITH MUTANT INCLUDED
            [trialTraits,trialDensities,trialBranchID] = ...
                getEquilibriumWithExtinction( ...
                trialTraits,...
                trialDensities0,...
                trialBranchID,...
                C1,C2,sigma1,sigma2,x0,sigma,beta,r,...
                TecoStart,tolEco,maxEcoExtensions,...
                extinctionTol,odeOptionsBase);

            % ACCEPT THE NEW ECOLOGICAL STATE
            traits = trialTraits;
            densities = trialDensities;
            branchID = trialBranchID;

            % CHECK WHETHER THE NEW MUTANT SURVIVED
            newBranchSurvives = any(branchID == newBranchID);
            if newBranchSurvives
                % Create a new history column filled with NaNs
                traitHistory(:,newBranchID) = NaN;
                densityHistory(:,newBranchID) = NaN;
                % The new branch now permanently exists
                nextID = nextID + 1;
                % Record which branch produced it
                parentOf(newBranchID) = parentBranchID;
                % Record when it appeared
                birthTime(newBranchID) = time + 1;
            end
        end

        % ADVANCE EVOLUTIONARY TIME
        time = time + 1

        % RECORD THE CURRENT STATE
        % Clear this row first.
        % This means extinct branches become NaN from this point onward.
        traitHistory(step+1,:) = nan;
        densityHistory(step+1,:) = nan;

        % Put every currently existing branch in its permanent column
        traitHistory(step+1,branchID) = traits;
        densityHistory(step+1,branchID) = densities;
        % Record time
        timeHistory(step+1) = time;
        % Move to next row
        step = step + 1;
    end


    % REMOVE UNUSED ROWS
    traitHistory = traitHistory(1:step,:);
    densityHistory = densityHistory(1:step,:);
    timeHistory = timeHistory(1:step);

    % SAVE THIS REALISATION
    traitHistories{realisation} = traitHistory;
    densityHistories{realisation} = densityHistory;
    timeHistories{realisation} = timeHistory;
    parentHistories{realisation} = parentOf;
    birthTimeHistories{realisation} = birthTime;

    % SAVE FINAL STATE
    finalTraitValues{realisation} = traits;
    finalDensityValues{realisation} = densities;
end



% PLOT TRAIT VALUES USING STAIRS
figure;
hold on;
for realisation = 1:numberofrealisations
    traitHistory = traitHistories{realisation};
    timeHistory = timeHistories{realisation};
    % Use the same colour for all branches belonging to the same
    % realisation
    colour = C(mod(realisation-1,size(C,1))+1,:);
    % Every column corresponds to one permanent branch ID
    for branch = 1:size(traitHistory,2)
        stairs( ...
            timeHistory,...
            traitHistory(:,branch),...
            'Color',colour,...
            'LineWidth',1);
    end
end
xlabel('Time');
ylabel('Trait');
box on;
hold off;

% PLOT POPULATION DENSITIES USING STAIRS
figure;
hold on;
for realisation = 1:numberofrealisations
    densityHistory = densityHistories{realisation};
    timeHistory = timeHistories{realisation};
    colour = C(mod(realisation-1,size(C,1))+1,:);
    for branch = 1:size(densityHistory,2)
        stairs( ...
            timeHistory,...
            densityHistory(:,branch),...
            'Color',colour,...
            'LineWidth',1);
    end
end
xlabel('Time');
ylabel('Population density');
box on;
hold off;

%% ========================================================================
%                       LOCAL FUNCTIONS
%  ========================================================================
% INVASION FITNESS

function s = calcInvasionFitness( ...
    mutantTrait,...
    traits,...
    densities,...
    C1,C2,sigma1,sigma2,x0,sigma,beta,r)
    % Carrying capacity of the mutant
    mutantK = K( ...
        C1,C2,sigma1,sigma2,x0,mutantTrait);
    % Competition exerted on the mutant by EVERY resident
    competitionFromResidents = ...
        alpha(mutantTrait,traits,sigma,beta) .* densities;
    % Integral approximated by a sum
    totalCompetition = sum(competitionFromResidents);
    % Invasion fitness
    s = r * ...
        (1 - totalCompetition ./ mutantK);
end



% ECOLOGICAL ODE

function dDdt = ecologicalODEGeneral(t,densities,traits,C1,C2,sigma1,sigma2,x0,sigma,beta,r)
    % Competition matrix
    %
    % A(i,j) = alpha(trait_i,trait_j)
    competitionMatrix = ...
        alpha(traits,traits.',sigma,beta);
    % Carrying capacity of every resident
    Kvec = K( ...
        C1,C2,sigma1,sigma2,x0,traits);
    % Ecological dynamics
    %
    % dN_i/dt =
    % r*N_i*(1 - sum_j alpha(x_i,x_j)*N_j/K_i)
    dDdt = r .* densities .* ...
        (1 - (competitionMatrix * densities) ./ Kvec);
end

% FIND ECOLOGICAL EQUILIBRIUM USING ODE45
function equilibriumDensities = getEquilibriumDensities( ...
    traits,...
    initialDensities,...
    C1,C2,sigma1,sigma2,x0,sigma,beta,r,...
    TecoStart,tolEco,maxEcoExtensions,odeOptionsBase)
    % Make sure both are columns
    traits = traits(:);
    initialDensities = initialDensities(:);
    % Start with the initial ecological time
    Teco = TecoStart;
    % Repeatedly integrate until equilibrium is reached
    for extension = 1:maxEcoExtensions
        % Number of resident branches
        numberOfBranches = numel(traits);
        % Prevent ode45 from producing negative densities
        odeOptions = odeset( ...
            odeOptionsBase,...
            'NonNegative',1:numberOfBranches);
        % Integrate ecological dynamics
        [~,densityTrajectory] = ode45( ...
            @(t,densities) ecologicalODEGeneral( ...
                t,...
                densities,...
                traits,...
                C1,C2,sigma1,sigma2,x0,sigma,beta,r),...
            [0 Teco],...
            initialDensities,...
            odeOptions);
        % Take the final density as the candidate equilibrium
        equilibriumDensities = ...
            densityTrajectory(end,:).';
        % Evaluate the ecological derivative at this point
        dDdt = ecologicalODEGeneral( ...
            0,...
            equilibriumDensities,...
            traits,...
            C1,C2,sigma1,sigma2,x0,sigma,beta,r);
        % Measure how far we are from equilibrium
        equilibriumError = max(abs(dDdt));
        % Stop if sufficiently close to equilibrium
        if equilibriumError < tolEco
            return
        end
        % If equilibrium has not been reached,
        % use the current final state as the new initial state
        initialDensities = equilibriumDensities;
        % Double the integration time
        Teco = 2*Teco;
    end
    % Warn if equilibrium still has not been reached
    warning( ...
        'Ecological equilibrium was not reached within the allowed integration time.');
end

% ECOLOGICAL EQUILIBRIUM WHILE REMOVING EXTINCT BRANCHES

function [remainingTraits,remainingDensities,remainingBranchID] = ...
    getEquilibriumWithExtinction( ...
    traits,...
    initialDensities,...
    branchID,...
    C1,C2,sigma1,sigma2,x0,sigma,beta,r,...
    TecoStart,tolEco,maxEcoExtensions,...
    extinctionTol,odeOptionsBase)
    % Force vectors to be columns
    traits = traits(:);
    initialDensities = initialDensities(:);
    branchID = branchID(:);
    while true

        % FIND ECOLOGICAL EQUILIBRIUM

        equilibriumDensities = getEquilibriumDensities( ...
            traits,...
            initialDensities,...
            C1,C2,sigma1,sigma2,x0,sigma,beta,r,...
            TecoStart,tolEco,maxEcoExtensions,...
            odeOptionsBase);


        % IDENTIFY EXTINCT BRANCHES
        extinct = equilibriumDensities <= extinctionTol;

        % IF NO BRANCHES ARE EXTINCT, WE ARE DONE
        if ~any(extinct)
            remainingTraits = traits;
            remainingDensities = equilibriumDensities;
            remainingBranchID = branchID;
            return
        end

        % REMOVE EXTINCT BRANCHES
        traits(extinct) = [];
        branchID(extinct) = [];
        initialDensities = ...
            equilibriumDensities(~extinct);
        % CHECK WHETHER ANY BRANCHES REMAIN
        if isempty(traits)
            remainingTraits = [];
            remainingDensities = [];
            remainingBranchID = [];
            return
        end
    end
end
%}

%% Branching video
%{
 vidfile = VideoWriter('testmovie.mp4(1)','MPEG-4');
 vidfile.FrameRate=40;
 open(vidfile);
 f=figure;
 for j = 1:T-3 %because I plot 3 steps at a time
    f.Name = ['Simulation time: t = ', num2str((j-1))];
    %Plotting
    hold on;
    for i=1:numberofrealisations
        h=stairs(timeplot([j:j+3],i),xplot([j:j+3],i),'Color',C(mod(i,7)+1,:));
        k=stairs(timeplot([j:j+3],i),xbplot([j:j+3],i),'Color',C(mod(i,7)+1,:));
        set(h,'Linewidth',2);
        set(k,'Linewidth',2);
    end
    hold on;
    hx1=plot([0,T],x1*ones(1,2),'k-',LineWidth=2);
    hx2=plot([0,T],x2*ones(1,2),'r-.',LineWidth=2);
    hx3=plot([0,T],x3*ones(1,2),'b--',LineWidth=2);
    legend(Location="southeast")
    %legend([hx2 hx1 hx3], {sprintf('x_0+x_1 = %.2f', x2), sprintf('x_0 = %.2f', x1), sprintf('x_0-x_1 = %.2f', x3)});
    legend([hx1 hx3], {sprintf('x_0   = repeller'), sprintf('x_0-x_1 = branch point')});
    xlabel('time','FontSize',25);
    ylabel('Trait value','FontSize',25);
    axis([0 T 0 x0+1]);
    set(gca,'linewidth',1.5);
    set(gca,'FontSize',20);
    grid on;
    axis square;
    writeVideo(vidfile, getframe(gcf));%takes screenshot of the 3 new steps
    each time
 end
close(vidfile)
%}