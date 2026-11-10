%%Alice Henman%%
%%Evolutionary Branching Paper - Figure 6 i-v%%

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
function alpha=alpha(x,y,sigma,gamma)
    alpha = (1-gamma*((x-y)./sigma).^2).*exp(-(x-y).^2./(2*sigma.^2));
end

%dalpha(x-y)/dx, dalpha(y,x,...) = dalpha(y-x)/dy = -dalpha(y-x)/dx
function dalpha=dalpha(x,y,sigma,gamma)
    dalpha = (-1-2*gamma+gamma*((x-y)./sigma).^2).*exp(-(x-y).^2./(2*sigma.^2)).*(x-y)./sigma.^2;
end
%d^2alpha(x-y)/dx^2, ddalpha(y,x,...) = d^2alpha(y-x)/dy^2
function ddalpha=ddalpha(x,y,sigma,gamma)
z = (x-y);
    ddalpha = exp(-z.^2./(2*sigma.^2)).* ...
(-1-2*gamma+(z./sigma).^2.*(1+5*gamma-gamma*(z./sigma).^2))./sigma.^2;
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
rho=0.5;
%i: rho=1.5
%ii: rho=1.75
%iii: rho=2
%iv: rho=0.7
%v: rho=0.5
sigma = 0.5*sqrt(rho);%sqrt(1.5/4.5);%A=(1+2gamma)/sigma^2=4.5>4
sigma1 = 0.5;%sqrt(0.27/(1+2gamma))%sqrt(0.5/(1+2gamma))
sigma2 = 0.35;%sigma1*sqrt(C2/C1)*rand; % condition for there to be 2 humps
x0 = 1.5;%10;%1+5*rand;
beta=0;
gamma=0.6;
%i: gamma=0.1
%ii: gamma=0.4
%iii: gamma=0.75
%iv: gamma=0.2
%v: gamma=0.6
r=0.5;
m=0.01; %max size of mutation
epsM = 1e-6; %fraction of population that mutates
xinitial=1; %initial resident trait value
xmin=0;
xmax=3;

%% Equilibria
N1 = @(x) K(C1,C2,sigma1,sigma2,x0,x)/alpha(0,0,sigma,gamma); %1 resident equilibrium
%equilibrium for 2 residents 
%The equilibrium pop density of resident with trait x1, swap x1 and x2 to
%get pop dens of resident with trait x2
function N = N(x1,C1,C2,sigma1,sigma2,x0,sigma,gamma,x2)
    N = (K(C1,C2,sigma1,sigma2,x0,x1).*alpha(0,0,sigma,gamma)- ...
        K(C1,C2,sigma1,sigma2,x0,x2).*alpha(x1,x2,sigma,gamma))./( ...
        alpha(0,0,sigma,gamma).^2-alpha(x1,x2,sigma,gamma).*alpha(x2,x1,sigma,gamma));
end
Nfun = @(x,xb) N(x,C1,C2,sigma1,sigma2,x0,sigma,gamma,xb);

%popdens vector for n traits, traits must be a column vector [x1;x2;...;xn]
function Nn = Nn(traits,C1,C2,sigma1,sigma2,x0,sigma,gamma)

    A = alpha(traits,traits.',sigma,gamma);%vector Aij=alpha(xi-xj)
    Kvec = K(C1,C2,sigma1,sigma2,x0,traits);
    Kvec = Kvec(:);%turns Kvec into column vector

    Nn = A\Kvec; %finds vector N=[N1;N2;...Nn] by solving K=AN
    if ~all(isfinite(Nn)) %A singular to working precision (e.g. a tight cluster of branches)
        Nn = lsqminnorm(A,Kvec); %smallest solution of A*N=K, shares density within the cluster
    end
end
%% Invasion fitness
f = @(x,y) r*(1-alpha(y,x,sigma,gamma).*K(C1,C2,sigma1,sigma2,x0,x) ...
    ./(alpha(0,0,sigma,gamma).*K(C1,C2,sigma1,sigma2,x0,y)));

%fitness for 2 residents
function fit12 =fit12(r,sigma,gamma,C1,C2,sigma1,sigma2,x0,y,x1,x2)
    fit12 = r*(1-(alpha(y,x1,sigma,gamma).*N(x1,C1,C2,sigma1,sigma2,x0,sigma,gamma,x2)+ ...
        alpha(y,x2,sigma,gamma).*N(x2,C1,C2,sigma1,sigma2,x0,sigma,gamma,x1))./K(C1,C2,sigma1,sigma2,x0,y));
end
f2 = @(x1,x2,y) fit12(r,sigma,gamma,C1,C2,sigma1,sigma2,x0,y,x1,x2);

% fitness for n residents
function fitn=fitn(r,sigma,gamma,C1,C2,sigma1,sigma2,x0,Nvec,y,traits)
    A = alpha(y,traits,sigma,gamma);
    Ky = K(C1,C2,sigma1,sigma2,x0,y);
    fitn = r*(1-sum(A.*Nvec)./Ky);
end
fn = @(Nvec,y,traits) fitn(r,sigma,gamma,C1,C2,sigma1,sigma2,x0,Nvec,y,traits);

%produces a vector of the values of fyy*K(x) at each branch value
function fityyn = fityyn(r,sigma,gamma,C1,C2,sigma1,sigma2,x0,Nvec,traits)
    ddA = ddalpha(traits,traits.',sigma,gamma);
    ddKx = ddK(C1,C2,sigma1,sigma2,x0,traits);
    Kx = K(C1,C2,sigma1,sigma2,x0,traits);
    fityyn = r*(ddKx-ddA*Nvec);
end

%% Selection gradient
%selection gradient for branch with resident x1
function secgrad1 = secgrad1(r,sigma,gamma,C1,C2,sigma1,sigma2,x0,x1,x2)
    secgrad1 = r*(-(dalpha(0,0,sigma,gamma).*N(x1,C1,C2,sigma1,sigma2,x0,sigma,gamma,x2)+ ...
        dalpha(x1,x2,sigma,gamma).*N(x2,C1,C2,sigma1,sigma2,x0,sigma,gamma,x1))./K(C1,C2,sigma1,sigma2,x0,x1)+ ...
        dK(C1,C2,sigma1,sigma2,x0,x1).*(alpha(0,0,sigma,gamma).*N(x1,C1,C2,sigma1,sigma2,x0,sigma,gamma,x2)+ ...
        alpha(x1,x2,sigma,gamma).*N(x2,C1,C2,sigma1,sigma2,x0,sigma,gamma,x1))./K(C1,C2,sigma1,sigma2,x0,x1).^2);
end

%selection gradient for branch with resident x2
function secgrad2 = secgrad2(r,sigma,gamma,C1,C2,sigma1,sigma2,x0,x1,x2)
    secgrad2 = r*(-(dalpha(x2,x1,sigma,gamma).*N(x1,C1,C2,sigma1,sigma2,x0,sigma,gamma,x2)+ ...
        dalpha(0,0,sigma,gamma).*N(x2,C1,C2,sigma1,sigma2,x0,sigma,gamma,x1))./K(C1,C2,sigma1,sigma2,x0,x2)+ ...
        dK(C1,C2,sigma1,sigma2,x0,x2).*(alpha(x2,x1,sigma,gamma).*N(x1,C1,C2,sigma1,sigma2,x0,sigma,gamma,x2)+ ...
        alpha(0,0,sigma,gamma).*N(x2,C1,C2,sigma1,sigma2,x0,sigma,gamma,x1))./K(C1,C2,sigma1,sigma2,x0,x2).^2);
end

sx1 = @(x1,x2) secgrad1(r,sigma,gamma,C1,C2,sigma1,sigma2,x0,x1,x2);
sx2 = @(x1,x2) secgrad2(r,sigma,gamma,C1,C2,sigma1,sigma2,x0,x1,x2);
%Evolutionary singular strategies
x1=x0;
%x2=x0+sqrt(2*sigma1^2*sigma2^2*log(C1*sigma2^2/(C2*sigma1^2))/(sigma2^2-sigma1^2));
%x3=x0-sqrt(2*sigma1^2*sigma2^2*log(C1*sigma2^2/(C2*sigma1^2))/(sigma2^2-sigma1^2));

%produces a vector of the values of s(x) at each branch value, using observed N
function secgradn = secgradn(r,sigma,gamma,C1,C2,sigma1,sigma2,x0,Nvec,traits)
    dKx = dK(C1,C2,sigma1,sigma2,x0,traits);
    dA = dalpha(traits,traits.',sigma,gamma);
    secgradn = r*(dKx - dA*Nvec);
end

%produces a vector of the values of s(x) at each branch value, calculating N
function sxn_fun = sxn_fun(X, r, sigma, gamma, C1, C2, sigma1,sigma2, x0)
    X  = X(:);
    A  = alpha(X,X.',sigma,gamma);
    dA = dalpha(X,X.',sigma,gamma);
    Kx = K(C1,C2,sigma1,sigma2,x0,X);
    dKx = dK(C1,C2,sigma1,sigma2,x0,X);
    Nn = A\Kx; %finds vector N=[N1;N2;...Nn] by solving K=AN
    if ~all(isfinite(Nn)) %A singular to working precision (e.g. a tight cluster of branches)
        Nn = lsqminnorm(A,Kx); %smallest solution of A*N=K, shares density within the cluster
    end
    sxn_fun  = r*(dKx - dA*Nn);        % = K .* s(x)
end

%% Numerical Simulation (up to n branches)
%
numberofrealisations = 1; % number of sample paths
T = 1e8; % end time for simulation
% matrix to hold trait values at each timestep, for each realisation for
% each branch
traitHistories = cell(numberofrealisations,1);%starts  with 1 for 1 branch
densityHistories = cell(numberofrealisations,1);
timeHistories = cell(numberofrealisations,1);
ws = warning('off','MATLAB:nearlySingularMatrix');%to stop flooding command window
numberofpoints = min(T,1e3);%5;%
stride = max(1,floor((T+1)/numberofpoints)); %plot at most ~50000 time points
numberofrows = floor(T/stride) + 1;
for i = 1:numberofrealisations %Run numberofrealisations simulations
    traits=xinitial; %vector column of all current traits
    Nlist=N1(xinitial); %vector column of all current densities
    Nnew=Nlist; %initial density
    time=0;
    j=1;
    branch = 1; %initially 1 branch only
    
    %list of all traits and densities for all time
    xplot = nan(numberofrows,1); %xplot(1,1)=trait value of branch 1 at time 0
    Nplot = nan(numberofrows,1); %Nplot(T,5)=trait value of branch 5 at time T-1
    %keeps track of time each trait/density held at
    timeplot = nan(numberofrows,1);

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
            Nswapped = Nn(swappedtrait,C1,C2,sigma1,sigma2,x0,sigma,gamma);
            if fn(Nswapped,x,swappedtrait)>0 %if x has +ve f when mutant is resident
                %we have mutual invasibility, add y to list of branches
                traits = [traits;y];
                branch = numel(traits);
                if branch > size(xplot,2) %only add a column if needed, don't wipe an existing one
                    xplot(:,branch)=nan;
                    Nplot(:,branch)=nan;
                end
                Nnew = Nn(traits,C1,C2,sigma1,sigma2,x0,sigma,gamma);%find new densities
                %this method ensures the branch removed was the correct choice
                %(one negative N means the eqm system is wrong not just that it should be removed)
                removed = []; %traits of branches removed this step
                for pass = 1:50 %cap on passes, in case near-twins flip back and forth
                    while any(Nnew<=0) %checks if any branches died
                        [~,w] = min(Nnew);              %remove only the most negative branch,
                        removed = [removed; traits(w)]; %others may survive once it's gone
                        traits(w) = [];
                        Nnew = Nn(traits,C1,C2,sigma1,sigma2,x0,sigma,gamma);%find new densities
                    end
                    if isempty(removed), break; end
                    fitBack = arrayfun(@(z) fn(Nnew,z,traits), removed); %could a removed branch grow back?
                    [fbest,w] = max(fitBack);
                    if fbest <= 1e-6, break; end        %no (ignoring rounding-level twin cases)
                    traits = [traits; removed(w)];      %yes: put it back and re-solve
                    removed(w) = [];
                    Nnew = Nn(traits,C1,C2,sigma1,sigma2,x0,sigma,gamma);
                end
                while any(Nnew<=0) %checks if any branches died
                    [~,w] = min(Nnew);              %remove only the most negative branch,
                    traits(w) = [];
                    Nnew = Nn(traits,C1,C2,sigma1,sigma2,x0,sigma,gamma);%find new densities
                end
            else %not mutual invasibility
                traits = swappedtrait;%y replaces x
                Nnew=Nswapped;
            end
        %else x doesn't change
        end
        Nlist = Nnew;
        time=time+1;%update timestep
        if mod(time,500000)==0, fprintf('\n %.3f',100*time/T); end %progress every 5000 steps
        %update matrices (reduces points being stored)
        if mod(time,stride)==0
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
        
    end
    traitHistories{i} = xplot(1:j,:);
    densityHistories{i} = Nplot(1:j,:);
    timeHistories{i} = timeplot(1:j);
end

%Plotting

tol=0.01;
%merges traits closer than tol into one: density = sum, trait = density-weighted mean
function [xc,Nc] = mergeClose(x,N,tol)
    [x,o] = sort(x(:)); N = N(o); N = N(:);
    edges = [0; find(diff(x) > tol); numel(x)]; %a gap bigger than tol starts a new cluster
    nc = numel(edges)-1;
    xc = zeros(nc,1); Nc = zeros(nc,1);
    for k = 1:nc
        idx = edges(k)+1:edges(k+1);
        Nc(k) = sum(N(idx));
        xc(k) = sum(x(idx).*N(idx))/Nc(k);
    end
end

figure;
hold on;
cmap=jet;
%plot trait value by time for each branch
%
for i=1:numberofrealisations
    X = traitHistories{i}; Nh = densityHistories{i}; t = timeHistories{i};
    tp = nan(numel(X),1); xp = tp; np = tp; c = 0; %preallocated, filled row by row
    for j = 1:size(X,1)
        keep = ~isnan(X(j,:));
        if isnan(t(j)) || ~any(keep), continue, end
        [xc,Nc] = mergeClose(X(j,keep),Nh(j,keep),tol);
        n = numel(xc);
        tp(c+1:c+n) = t(j); xp(c+1:c+n) = xc; np(c+1:c+n) = Nc; c = c+n;
    end
    scatter(tp(1:c),xp(1:c),25,log10(np(1:c)),"filled");
end
colormap("turbo");
cb=colorbar;
cb.Label.String = "log_{10}(Density)";
%plots singular strategies
hx1=yline(x0,'--',LineWidth=2);
%legend(Location="northwest")
%legend(hx1,sprintf('x* = x_0'))
xlabel('Time');
ylabel('Trait value');
xlim([0 T]);
grid on;
axis square;
box on;

%show value of s(x) and final derivative of each branch
keep=~isnan(X(end,:));
finalx=X(end,keep);
finalN=Nh(end,keep);
[finalx,o] = sort(finalx(:))
finalN = finalN(o);
finalN = finalN(:);
secgradn(r,sigma,gamma,C1,C2,sigma1,sigma2,x0,finalN,finalx)
fityyn(r,sigma,gamma,C1,C2,sigma1,sigma2,x0,finalN,finalx)

[finalx,finalN] = mergeClose(finalx,finalN,tol);
finalx
secgradn(r,sigma,gamma,C1,C2,sigma1,sigma2,x0,finalN,finalx)
fityyn(r,sigma,gamma,C1,C2,sigma1,sigma2,x0,finalN,finalx)

opts = optimset('Display','off', 'TolX',1e-12, 'TolFun',1e-12, 'MaxIter',20000, 'MaxFunEvals',100000);
sx = @(X) sxn_fun(X, r, sigma, gamma, C1, C2, sigma1,sigma2, x0);
S=fsolve(sx,finalx,opts)
secgradn(r,sigma,gamma,C1,C2,sigma1,sigma2,x0,finalN,finalx)
fityyn(r,sigma,gamma,C1,C2,sigma1,sigma2,x0,finalN,finalx)
%}