%%Alice Henman%%
%%Pigolotti but with facilitation

%%Default graphing layout from Denis
set(0,'defaultTextFontName', 'Arial');
set(0,'defaultaxesfontsize', 20); % 25 for 1X3, 20 for 1X2 figures
%set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','none');
set(0,'defaulttextinterpreter','none');
set(0,'defaultAxesXGrid','off');
set(0,'defaultAxesYGrid','off');
set(0,'defaultAxesTickDir','out');
set(0,'defaultAxesLineWidth',1.5);

C = orderedcolors("gem12");

%% parameters
r=1;
x0 = 1.5;
sigmak = 0.1/sqrt(2.5);%sigma/sqrt(rho)
sigma = 0.1;%sigmak*sqrt(1.5);
gamma = 1;
nT=1; %threshold population
ai=0.1;
C1 = nT/ai;
I=0.004;
xmax = 4*sigmak+x0;%38.6*sigmak+x0;%1;%
xmin = x0-4*sigmak;%x0-38.6*sigmak;%0;%
L=xmax-xmin; %length of domain
numspecies = 200;
T = 5e5;%1e4;%
dt = 0.1;
xinitial = xmin + (xmax - xmin)*rand(numspecies,1);
ninitial = nT + rand(numspecies,1);
%ninitial = K(C1,sigmak,x0,xinitial);

%model setup
periodic=false;
pig=false;
krondelt = false;


%% functions

%(Mexican hat) Wavelet Competition function
function alpha=alpha(x,y,sigma,gamma)
    alpha = (1-gamma*((x-y)./sigma).^2).*exp(-(x-y).^2./(2*sigma.^2));
end

function alpha=alpha1(z,R)%here z=the difference xi-xj
    alpha = exp(-abs(z)/R);
end

function alpha=alpha2(z,R)%here z=the difference xi-xj
    alpha = exp(-(abs(z)/R).^4);
end

%K is scarrying capacity
function K=K(C,sigma,x0,x)
    K = C*exp(-(x-x0).^2./(2*sigma^2));
end

function ddt = odefcn1(N,X,r,C1,x0,sigmak,sigma,gamma,L,periodic,pig,krondelt)
    X=X(:);
    N=N(:);
    D=X-X.';
    if periodic 
        D = mod(D + L/2, L) - L/2;
    end
    if pig
        %A = alpha1(D,sigma);
        A = alpha2(D,sigma);
        if krondelt
            A = A+eye(numel(X));
        end
        ddt = r*N.*(1-A*N./C1);
    else
        A = alpha(D, 0, sigma, gamma);
        ddt = r*N.*(1-A*N./K(C1,sigmak,x0,X));
    end
end

f = @(N,X) odefcn1(N,X,r,C1,x0,sigmak,sigma,gamma,L,periodic,pig,krondelt);

function n = rk2(f,N,X,dt)
    k1 = f(N,X);
    k2 = f(N+dt*k1,X);
    n = N+dt*(k1+k2)/2;
end

numberofpoints = min(T/dt,100);%5;%
stride = max(1,floor((T/dt+1)/numberofpoints)); %plot at most ~100 time points
recordedtimes = floor(T/(dt*stride)) + 1;

xvals = cell(recordedtimes,1);  
nvals = cell(recordedtimes,1);
tvals = cell(recordedtimes,1);
x=xinitial;
xvals{1} = x;
n=ninitial;
nvals{1} = n;
tvals{1} = zeros(numel(x),1);
j=1;
for step = 1:round(T/dt)
    t=step*dt;
    nold=n;
    n=rk2(f,nold,x,dt);
    while any(n<=nT) %checks if any branches died
        [~,w] = min(n);  %remove only the smallest branch,
        x(w) = [];
        nold(w) = [];
        n=rk2(f,nold,x,dt);
    end
    if t<T/2 && rand <I*dt %if we have immigration
        x = [x;xmin+rand*(xmax-xmin)]; %add new species at random point
        n = [n;nT+rand]; %with population a bit larger than threshold
    end
    if isempty(n)
        warning('Everything went extinct at t = %g', t);
        break;
    end
    if mod(step,stride)==0
        j=j+1;
        nvals{j} = n;
        xvals{j} = x;
        tvals{j} = t*ones(numel(x),1);
    end
end

%% plotting final trait distribution
figure;
stem(x(:), n(:), 'filled', 'LineWidth', 1.5);
xlabel('Niche position x'); ylabel('Population n');
xlim([xmin xmax]);

%% plotting traits over time
%
figure;
scatter(vertcat(tvals{:}), vertcat(xvals{:}), 4, log(vertcat(nvals{:})), 'filled');
colormap turbo; cb = colorbar; ylabel(cb, 'log(Density)');
xlabel('Time'); ylabel('Niche position x');
xline(T/2, 'k--', 'Immigration off', 'LineWidth', 1.5);
%}