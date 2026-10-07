%% 
%%Alice Henman%%
%%Evolutionary Branching in Models for Symmetric and Asymmetric
%%Competition%%

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
C1 = 1;%2.5;
C2 = 0;%2;
sigma1 = 0.5;%0.5;%
sigma2 = 0.15;%0.35;%
x0 = 1.5;
sigma = 0.5; %0.25;%width
beta0 = 0;   %symmetrical

%% (Mexican hat) Wavelet Competition function

% z represents x-y
z=linspace(-1.5,1.5,500);

%parameters
sigma0 = 0.5;
beta1 = 1.5; %asymmetric

function alpha=alpha(x,y,sigma,beta)
    alpha = exp(sigma.^2.*beta.^2/2).* ...
(1-((x-y+sigma.^2.*beta)./sigma).^2).*exp(-(x-y+sigma.^2.*beta).^2./(2*sigma.^2));
end
%The beta parameter can add asymmetry to the competition
%but I'm just keeping it as beta=0 for now

%dalpha(x-y)/dx, dalpha(y,x,...) = dalpha(y-x)/dy = -dalpha(y-x)/dx
function dalpha=dalpha(x,y,sigma,beta)
    dalpha = exp(sigma.^2.*beta.^2/2).* ...
(-3+((x-y+sigma.^2.*beta)./sigma).^2).*exp(-(x-y+sigma.^2.*beta).^2./(2*sigma.^2)).*(x-y+sigma.^2.*beta)./sigma.^2;
end
%d^2alpha(x-y)/dx^2, ddalpha(y,x,...) = d^2alpha(y-x)/dy^2
function ddalpha=ddalpha(x,y,sigma,beta)
A = (x-y+sigma^2*beta);
    ddalpha = exp(sigma.^2.*beta.^2/2-A.^2./(2*sigma.^2)).* ...
(-3./sigma.^2+(A./sigma.^2).^2.*(6-(A./sigma).^2));
end

%z=x-y
%{
figure;
%Plots of competition for 2 different widths
plot(z,alpha(z,0,sigma,beta0),'DisplayName', sprintf('σ_α = %.2f',sigma), ...
    "LineWidth",1.5);
hold on;
plot(z,alpha(z,0,sigma0,beta0),'r--','DisplayName', sprintf('σ_α = %.2f',sigma0), ...
    "LineWidth",1.5);
xlabel('Difference in trait value (x-y)');
ylabel('Strength of competition α(x − y)');
grid on;
ax=gca;
ax.FontSize = 20;
legend(FontSize=20);
legend('boxoff');
%}
%% Double Gaussian Carrying capacity K(x)
%Gaussian centred on x0
function Gauss=G(sigma,x0,x)
    Gauss = exp(-(x-x0).^2./(2*sigma.^2));
end
%derivative of gaussian
function dGauss=dG(sigma,x0,x)
    dGauss = -(x-x0)./(sigma.^2).*exp(-(x-x0).^2./(2*sigma.^2));
end
%second derivative of gaussian
function ddGauss=ddG(sigma,x0,x)
    ddGauss = 1./(sigma.^2).*(-1+(x-x0).^2/(sigma.^2)).*exp(-(x-x0).^2./(2*sigma.^2));
end
%K is sum of 2 gaussians
function K=K(C1,C2,sigma1,sigma2,x0,x)
    K = C1*G(sigma1,x0,x)-C2*G(sigma2,x0,x);
end
%derivative of K (when alpha peaks at y=x, the solutions to this are x*)
function dK=dK(C1,C2,sigma1,sigma2,x0,x)
    dK = C1*dG(sigma1,x0,x)-C2*dG(sigma2,x0,x);
end

dKx = @(x) dK(C1,C2,sigma1,sigma2,x0,x);
%second derivative of K
function ddK=ddK(C1,C2,sigma1,sigma2,x0,x)
    ddK = C1*ddG(sigma1,x0,x)-C2*ddG(sigma2,x0,x);
end

x=linspace(0,3,500);
%
figure;
plot(x,K(C1,C2,sigma1,sigma2,x0,x),"LineWidth",3);
hold on;
xlabel('Trait value x');
ylabel('Carrying Capacity, K(x)');
grid on;
ax=gca;
ax.FontSize = 20;
%}
%% singular stategies 

%Parameters
r=rand+0.001; %growth rate
xmin=0;
xmax=3;
%singular strategies
%{
x2=x0+sqrt(2*sigma1^2*sigma2^2*log(C1*sigma2^2/(C2*sigma1^2))/(sigma2^2-sigma1^2));
x3=x0-sqrt(2*sigma1^2*sigma2^2*log(C1*sigma2^2/(C2*sigma1^2))/(sigma2^2-sigma1^2));
%}

%selection gradient
sx = @(x) r*(-dalpha(0,0,sigma,beta0)./alpha(0,0,sigma,beta0)+dK(C1,C2,sigma1,sigma2,x0,x) ...
    ./K(C1,C2,sigma1,sigma2,x0,x));
%these are general but for symmetric alpha dalpha(x,x) is 0 so these can be
%simplified and alpha(x,x) can be taken to be 1
%in general we expect 1 or 3 singular strategies as a bimodal K has 1 or 3
%stat points and the alpha term just adds/subtracts some amount

%alternatively by finding roots of K
%{
%can also find selection gradient and substitute dKx by sx (selection
%gradient function)
%xstar=roots(chebfun(sx,[1e-6 5]));

inter = find(diff(sign(sx(z))));%finds location of sign changes in linspace
n=numel(inter);%number of sign changes
H(k)=n;
xstar = zeros(1,n);
for i =1:n
    xstar(i)=fzero(sx,z(inter(i)));%finds location of singular strategies by solving s(x)=0
end
%}

%showing they align with the stationary points of K(x)
%{
legend("boxoff");
h1=plot(x0*ones(1,2),[0,1.5],'k-');
h2=plot(x2*ones(1,2),[0,1.5],'r-.');
h3=plot(x3*ones(1,2),[0,1.5],'b--');
legend([h1,h2,h3], {sprintf('x_0 = %.2f', x0) sprintf('x_0+x_1 = %.2f', x2) sprintf('x_0-x_1 = %.2f', x3)});
%}
%% The Mathssss (The PIP)

%Invasion fitness
%for bimodal K
f = @(x,y) r*(1-alpha(y,x,sigma,beta0).*K(C1,C2,sigma1,sigma2,x0,x) ...
    ./(alpha(0,0,sigma,beta0).*K(C1,C2,sigma1,sigma2,x0,y)));

y = linspace(xmin,xmax,500);
[x,y] = meshgrid(y, y);
F = f(x,y);

figure
hold on
% PIP
contourf(x,y,F>0,1,'LineStyle','none')
contour(x,y,F,[0 0],'k','LineWidth',2) %black line for f=0
colormap([1 1 1; 0.6 0.9 0.3]); %green for f>0 white for f<0
xlabel('$\textsf{Resident trait value}$ ($x$)','Interpreter','latex','FontSize',17);
ylabel('$\textsf{Mutant trait value}$ ($y$)','Interpreter','latex','FontSize',17);
title(['PIP for σ_{α} = ',num2str(sigma)],'Interpreter','tex');
xline(x0,DisplayName=sprintf('x_0 = %.2f', x0));
%{
h1=plot(x0*ones(1,2),[y(1),y(end)],'k-');
h2=plot(x2*ones(1,2),[y(1),y(end)],'r-.');
h3=plot(x3*ones(1,2),[y(1),y(end)],'b--');

legend([h1,h2,h3], {sprintf('x_0 = %.2f', x0) sprintf('x_0+x_1 = %.2f', x2) sprintf('x_0-x_1 = %.2f', x3)});
%}
legend(Location="northwest");
axis equal
box on

%% 2 residents

%The equilibrium pop density of resident with trait x1, swap x1 and x2 to
%get pop dens of resident with trait x2
function N = N(x1,C1,C2,sigma1,sigma2,x0,sigma,beta,x2)
    N = (K(C1,C2,sigma1,sigma2,x0,x1).*alpha(0,0,sigma,beta)- ...
        K(C1,C2,sigma1,sigma2,x0,x2).*alpha(x1,x2,sigma,beta))./( ...
        alpha(0,0,sigma,beta).^2-alpha(x1,x2,sigma,beta).*alpha(x2,x1,sigma,beta));
end
%dN1/dx1 swapping x1 and x2 is dN2/dx2
function dN1x1 = dN1x1(x1,C1,C2,sigma1,sigma2,x0,sigma,beta,x2)
    K1 = K(C1,C2,sigma1,sigma2,x0,x1);
    K2 = K(C1,C2,sigma1,sigma2,x0,x2);
    dK1 = dK(C1,C2,sigma1,sigma2,x0,x1);
    a0 = alpha(0,0,sigma,beta);
    a12 = alpha(x1,x2,sigma,beta);
    a21 = alpha(x2,x1,sigma,beta);
    da12 = dalpha(x1,x2,sigma,beta);%dalpha(x1-x2)/dx1
    da21 = dalpha(x2,x1,sigma,beta);%dalpha(x2-x1)/dx2 = -dalpha(x2-x1)/dx1
    A = a0^2-a12.*a21;
    dA = -da12.*a21+a12.*da21;

    dN1x1 = (a0*dK1-K2.*da12)./A-(K1*a0-K2.*a12).*dA./A.^2;
end
%d^2N1/dx1^2 swapping x1 and x2 is d^2N2/dx2^2
function ddN1x1 = ddN1x1(x1,C1,C2,sigma1,sigma2,x0,sigma,beta,x2)
    K1 = K(C1,C2,sigma1,sigma2,x0,x1);
    K2 = K(C1,C2,sigma1,sigma2,x0,x2);
    dK1 = dK(C1,C2,sigma1,sigma2,x0,x1);
    ddK1 = ddK(C1,C2,sigma1,sigma2,x0,x1);
    a0 = alpha(0,0,sigma,beta);
    a12 = alpha(x1,x2,sigma,beta);
    a21 = alpha(x2,x1,sigma,beta);
    da12 = dalpha(x1,x2,sigma,beta);%dalpha(x1-x2)/dx1
    da21 = dalpha(x2,x1,sigma,beta);%dalpha(x2-x1)/dx2
    dda12 = ddalpha(x1,x2,sigma,beta);%d^2alpha(x1-x2)/dx1^2
    dda21 = ddalpha(x2,x1,sigma,beta);%d^2alpha(x2-x1)/dx2^2=d^2alpha(x2-x1)/dx1^2
    A = a0^2-a12.*a21;
    dA = -da12.*a21+a12.*da21;
    ddA = -dda12.*a21+2*da12.*da21-a12.*dda21;

    ddN1x1 = (ddK1*a0-K2.*dda12)./A-2*(dK1*a0-K2.*da12).*dA./A.^2-(K1*a0-K2.*a12).*ddA./A.^2+2*(K1*a0-K2.*a12).*dA.^2./A.^3;
end
%d^2N1/dx2^2 swapping x1 and x2 is d^2N2/dx1^2
function ddN1x2 = ddN1x2(x1,C1,C2,sigma1,sigma2,x0,sigma,beta,x2)
    K1 = K(C1,C2,sigma1,sigma2,x0,x1);
    K2 = K(C1,C2,sigma1,sigma2,x0,x2);
    dK2 = dK(C1,C2,sigma1,sigma2,x0,x2);
    ddK2 = ddK(C1,C2,sigma1,sigma2,x0,x2);
    a0 = alpha(0,0,sigma,beta);
    a12 = alpha(x1,x2,sigma,beta);
    a21 = alpha(x2,x1,sigma,beta);
    da12 = dalpha(x1,x2,sigma,beta);%dalpha(x1-x2)/dx1 = -dalpha(x1-x2)/dx2
    da21 = dalpha(x2,x1,sigma,beta);%dalpha(x2-x1)/dx2
    dda12 = ddalpha(x1,x2,sigma,beta);%d^2alpha(x1-x2)/dx1^2
    dda21 = ddalpha(x2,x1,sigma,beta);%d^2alpha(x2-x1)/dx2^2=d^2alpha(x2-x1)/dx1^2
    A = a0^2-a12.*a21;
    dA = da12.*a21-a12.*da21;
    ddA = -dda12.*a21+2*da12.*da21-a12.*dda21;
    k = K1*a0-K2.*a12;
    dk = -dK2.*a12+K2.*da12;
    ddk = -ddK2.*a12+2*dK2.*da12-K2.*dda12;

    ddN1x2 = ddk./A-2*dk.*dA./A.^2-k.*ddA./A.^2+2*k.*dA.^2./A.^3;
end

%fitness for 2 residents
function fit12 =fit12(r,sigma,beta,C1,C2,sigma1,sigma2,x0,y,x1,x2)
    fit12 = r*(1-(alpha(y,x1,sigma,beta).*N(x1,C1,C2,sigma1,sigma2,x0,sigma,beta,x2)+ ...
        alpha(y,x2,sigma,beta).*N(x2,C1,C2,sigma1,sigma2,x0,sigma,beta,x1))./K(C1,C2,sigma1,sigma2,x0,y));
end
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

sx1 = @(x1,x2) secgrad1(r,sigma,beta0,C1,C2,sigma1,sigma2,x0,x1,x2);
sx2 = @(x1,x2) secgrad2(r,sigma,beta0,C1,C2,sigma1,sigma2,x0,x1,x2);

%% TEP

z=linspace(xmin,xmax,2000);
[X1,X2] = meshgrid(z, z);

F1 = f(X1,X2); %x1 resident, x2 mutant
F2 = f(X2,X1); %x2 resident, x1 mutant
S1 = sx1(X1,X2);
S2 = sx2(X1,X2);
%xstar=roots(chebfun(sx,[1e-6 5])); %singular strategies

F = zeros(size(X1));
F(F1>0)=1;   %fx(y)>0
F(F2>0)=2;   %fy(x)>0
F(F1>0&F2>0)=3;   %mutual invasibility
F(X1>X2)=0;       %shows upper half only
S1(X1>=X2)=0;       %shows upper half only
S2(X1>=X2)=0;       %shows upper half only

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
%
%so we can have sparse arrows even if lots of data collected for accuracy
step = 80;%(stepsize of z)/25
X11= X1(1:step:end,1:step:end);
X22= X2(1:step:end,1:step:end);
S11=S1(1:step:end,1:step:end);
S22=S2(1:step:end,1:step:end);

%normalise arrow size
normy = sqrt(S11.^2+S22.^2);
normy(normy==0) = 1;
L = 0.08;
%vector field
quiver(X11,X22,S11./normy,S22./normy,0.6,"LineWidth",1,"Color","k","Alignment","center");
legend('','');
%selection gradient nullclines
contour(X1,X2,S1,[0 0],'r:','LineWidth',2,'DisplayName','s_{x_1,x_2}(x_1)=0');
contour(X1,X2,S2,[0 0],'b:','LineWidth',2,'DisplayName','s_{x_1,x_2}(x_2)=0');
%singular strategies
%xline(xstar(1),'LineStyle','--','LineWidth',1.5,'Color',C(3,:),'DisplayName',['B1=',num2str(xstar(1),3)])
%xline(xstar(2),'LineStyle','--','LineWidth',1.5,'Color',C(2,:),'DisplayName',['D=',num2str(xstar(2),3)])
%xline(xstar(end),'LineStyle','--','LineWidth',1.5,'Color',C(1,:),'DisplayName',['B2=',num2str(xstar(end),3)])

legend(Location="southeast");
axis equal
box on

%% Colour plot to show that the conditions for branching are satisfied
%{
figure;
z=linspace(0.001,1,1000); %range of sigmak1 and k2
[sk1,sk2] = meshgrid(z,z);
C=C1/C2; %ratio of C1/C2>1
h = @(x,y) 3*(y.^2-x.^2)./(2.*log(C*y.^2./x.^2));
H=h(sk1,sk2); %values of the RHS of branching condition
[h,w]=size(H);
for i = 2:h
    for j = 1:i
        if i>=j
            H(i,j)=-3;%makes sure all values with sigmak1<k2 are undefined
        end
    end
end
hlims=[-3 3]; %stops large values drowning out small variations close to 0
imagesc(z,z,H,hlims);
hold on;
cmap = [abyss;autumn];%makes difference between positive and negative obvious
colormap(cmap); 
colorbar;
%Marks all negatives as undefined, makes it obvious that values of "3"
%could be bigger
colorbar('Ticks',[-1.5,0,0.5,1,1.5,2,2.5,3],...
         'TickLabels',{'undefined','0','0.5','1','1.5','2','2.5','3+'})
xlabel("\sigma_{k1}",'Interpreter','tex');
ylabel("\sigma_{k2}",'Interpreter','tex');
set(gca,'YDir','normal'); %makes axis look normal
% title('C* colour plot on \sigma_{k1},\sigma_{k2} axis','Interpreter','tex');
%}
%% 2nd deriv
%second derivative of the 2 resident fitness evaluated at y=x1*
function fit1yy=fit1yy(r,sigma,beta,C1,C2,sigma1,sigma2,x0,x1,x2)
    A11 = alpha(0,0,sigma,beta);
    A12 = alpha(x1,x2,sigma,beta);
    dA11 = dalpha(0,0,sigma,beta);
    dA12 = dalpha(x1,x2,sigma,beta);
    ddA11 = ddalpha(0,0,sigma,beta);
    ddA12 = ddalpha(x1,x2,sigma,beta);
    N1 = N(x1,C1,C2,sigma1,sigma2,x0,sigma,beta,x2);
    N2 = N(x2,C1,C2,sigma1,sigma2,x0,sigma,beta,x1);
    K1 = K(C1,C2,sigma1,sigma2,x0,x1);
    dK1 = dK(C1,C2,sigma1,sigma2,x0,x1);
    ddK1 = ddK(C1,C2,sigma1,sigma2,x0,x1);

    fit1yy = r*(-(ddA11.*N1+ddA12.*N2)./K1+2*dK1.*(dA11.*N1+dA12.*N2)./K1.^2 + ...
        (ddK1.*K1-2*dK1.^2).*(A11.*N1+A12.*N2)./K1.^3);
end
%second derivative of the 2 resident fitness evaluated at y=x2*
function fit2yy=fit2yy(r,sigma,beta,C1,C2,sigma1,sigma2,x0,x1,x2)
    A22 = alpha(0,0,sigma,beta);
    A21 = alpha(x2,x1,sigma,beta);
    dA22 = dalpha(0,0,sigma,beta);
    dA21 = dalpha(x2,x1,sigma,beta);
    ddA22 = ddalpha(0,0,sigma,beta);
    ddA21 = ddalpha(x2,x1,sigma,beta);
    N1 = N(x1,C1,C2,sigma1,sigma2,x0,sigma,beta,x2);
    N2 = N(x2,C1,C2,sigma1,sigma2,x0,sigma,beta,x1);
    K2 = K(C1,C2,sigma1,sigma2,x0,x2);
    dK2 = dK(C1,C2,sigma1,sigma2,x0,x2);
    ddK2 = ddK(C1,C2,sigma1,sigma2,x0,x2);

    fit2yy = r*(-(ddA21.*N1+ddA22.*N2)./K2+2*dK2.*(dA21.*N1+dA22.*N2)./K2.^2 + ...
        (ddK2.*K2-2*dK2.^2).*(A21.*N1+A22.*N2)./K2.^3);
end
%d^2f12/dx1^2 evaluated at y=x1, swap x1 and x2 to get d^2f12/dx2^2 at y=x2
function fitxx=fitxx(r,sigma,beta,C1,C2,sigma1,sigma2,x0,x1,x2)
    a0 = alpha(0,0,sigma,beta);
    a12 = alpha(x1,x2,sigma,beta);
    da0 = dalpha(0,0,sigma,beta);
    dda0 = ddalpha(0,0,sigma,beta);
    N1 = N(x1,C1,C2,sigma1,sigma2,x0,sigma,beta,x2);
    dN1 = dN1x1(x1,C1,C2,sigma1,sigma2,x0,sigma,beta,x2);%dN1/dx1
    ddN1 = ddN1x1(x1,C1,C2,sigma1,sigma2,x0,sigma,beta,x2);%d^2N1/dx1^2
    ddN2 = ddN1x2(x2,C1,C2,sigma1,sigma2,x0,sigma,beta,x1);%d^2N2/dx1^2
    K1 = K(C1,C2,sigma1,sigma2,x0,x1);

    fitxx = -r*(dda0*N1-2*da0*dN1+a0*ddN1+a12.*ddN2)./K1;
end
function S = sx12_fun(X,r,sigma,beta,C1,C2,sigma1,sigma2,x0)
    x1 = X(1);
    x2 = X(2);

    S = [ ...
        secgrad1(r,sigma,beta,C1,C2,sigma1,sigma2,x0,x1,x2); ...
        secgrad2(r,sigma,beta,C1,C2,sigma1,sigma2,x0,x1,x2) ...
    ];
end
%% Tests location and type of singular strategies
%
%guesses
%for C1=2.5, C2=2, x0=1.5
%for sigma1=0.5, sigma2=0.35, sigm=0.5
x1=1;%0.481741;%0.590045;%1.05228;%1.75563
x2=2;%1.50225;%1.24537;%1.94872;%2.41096

%solving selection gradients=0 near these guesses
opts = optimset('Display','off', 'TolX',1e-12, 'TolFun',1e-12, 'MaxIter',20000, 'MaxFunEvals',100000);
sx12 = @(X) sx12_fun(X,r,sigma,beta0,C1,C2,sigma1,sigma2,x0);
S=fsolve(sx12,[x1,x2],opts);
fprintf('x1=%.11f',S(1))
fprintf(', x2=%.11f',S(2))

%Printing second derivatives at the singular strategy
fprintf('\n fyy(x1)=%f',fit1yy(r,sigma,beta0,C1,C2,sigma1,sigma2,x0,S(1),S(2)))
fprintf('\n fyy(x2)=%f',fit2yy(r,sigma,beta0,C1,C2,sigma1,sigma2,x0,S(1),S(2)))
fprintf('\n fx1x1(x1)=%f',fitxx(r,sigma,beta0,C1,C2,sigma1,sigma2,x0,S(1),S(2)))
fprintf('\n fx2x2(x2)=%f',fitxx(r,sigma,beta0,C1,C2,sigma1,sigma2,x0,S(2),S(1)))


%guesses
x1=0.481741;%1.1362;%0.590045;%1.05228;%1.75563
x2=1.50225;%1.8638;%1.24537;%1.94872;%2.41096

%solving selection gradients=0 near these guesses
S=fsolve(sx12,[x1,x2],opts);
fprintf('\n \nx1=%.11f',S(1))
fprintf(', x2=%.11f',S(2))

%Printing second derivatives at the singular strategy
fprintf('\n fyy(x1)=%f',fit1yy(r,sigma,beta0,C1,C2,sigma1,sigma2,x0,S(1),S(2)))
fprintf('\n fyy(x2)=%f',fit2yy(r,sigma,beta0,C1,C2,sigma1,sigma2,x0,S(1),S(2)))
fprintf('\n fx1x1(x1)=%f',fitxx(r,sigma,beta0,C1,C2,sigma1,sigma2,x0,S(1),S(2)))
fprintf('\n fx2x2(x2)=%f',fitxx(r,sigma,beta0,C1,C2,sigma1,sigma2,x0,S(2),S(1)))

%guesses
%{
zmin = sqrt(3)*sigma;
dels=delstar(sigma,beta0,C1,C2,sigma1,sigma2,x0,zmin);
x1=x0-dels/2;
x2=x0+dels/2;

%solving selection gradients=0 near these guesses
S=fsolve(sx12,[x1,x2],opts);
fprintf('\n \nx1=%.11f',S(1))
fprintf(', x2=%.11f',S(2))

%Printing second derivatives at the singular strategy
fprintf('\n fyy(x1)=%f',fit1yy(r,sigma,beta0,C1,C2,sigma1,sigma2,x0,S(1),S(2)))
fprintf('\n fyy(x2)=%f',fit2yy(r,sigma,beta0,C1,C2,sigma1,sigma2,x0,S(1),S(2)))
fprintf('\n fx1x1(x1)=%f',fitxx(r,sigma,beta0,C1,C2,sigma1,sigma2,x0,S(1),S(2)))
fprintf('\n fx2x2(x2)=%f',fitxx(r,sigma,beta0,C1,C2,sigma1,sigma2,x0,S(2),S(1)))
%}
%}
%% Example PIP to show branching conditions
%{
Y = linspace(1,1.48,500);
X = linspace(1,2,1000);
%Y = linspace(x3-0.1,x3+0.1,500);
[x,y] = meshgrid(X, X);
F = f(x,y);

%neighbours x' and y'
a=x0-0.15;
b=x0+0.15;

figure
hold on

%PIP
contourf(x,y,F>0,1,'LineStyle','none')
contour(x,y,F,[0 0],'k','LineWidth',2) %black line for f=0
colormap([1 1 1; 0.6 0.9 0.3]); %green for f>0 white for f<0
xlabel('$\textsf{Resident trait value}$ ($x$)','Interpreter','latex','FontSize',20);
ylabel('$\textsf{Mutant trait value}$ ($y$)','Interpreter','latex','FontSize',20);

%lines through a,b (x',y) and the singular strategy
h1=xline(x0,'--','f_x(y)<0','FontSize',20,'LabelOrientation','horizontal','LabelHorizontalAlignment','left','LineWidth',2);
h2=plot((b)*ones(1,2),[x(1),a],':',LineWidth=2,Color=C(1,:));
h3=plot((a)*ones(1,2),[x(1),b],'r-.',LineWidth=2,Color=C(2,:));
h4=plot([x(1),a],(b)*ones(1,2),':',LineWidth=2,Color=C(1,:));
h5=plot([x(1),b],(a)*ones(1,2),'r-.',LineWidth=2,Color=C(2,:));
h6=plot(X,2*x0-X,'--',LineWidth=2,Color=C(6,:)); %mutual invasibility

plot(b,a,'k.','MarkerSize',20);
plot(a,b,'k.','MarkerSize',20);
text(b,a,"f_{b}(a)>0",'Interpreter','tex','FontSize',20,'VerticalAlignment','bottom');
text(a,b,"f_{a}(b)>0",'Interpreter','tex','FontSize',20,'HorizontalAlignment','right','VerticalAlignment','top');

%legend(Location="northwest");
legend([h1,h6], {sprintf('x^*') sprintf("y=2x^*-x")},'Interpreter','tex');
axis equal
set(gca,'XTick',[a,x0,b], 'YTick', [a,b]);
set(gca,'XTickLabel',["a","x*","b"], 'YTickLabel', ["a","b"]);
box on

%% PIP with trait substitution sequence added

% Parameters
n_steps = 5; % number of mutations
step_size = 0.15;
initial_trait = 1.1;%x;
initial_trait2 = 1.9;%x;

% Trait values and meshgrid
traits = linspace(1,2,1000);
[X,Y] = meshgrid(traits, traits);

% Simulate TSS
[residents, mutants] = simulate_tss_single_mutant(initial_trait, n_steps, step_size, traits, f);
[residents2, mutants2] = simulate_tss_single_mutant(initial_trait2, n_steps, step_size, traits, f);

% Compute PIP
F = f(X,Y);

% Plot PIP and TSS
figure;
hold on;
contourf(X,Y,F>0,1,'LineStyle','none');
contour(X,Y,F,[0 0],'k','LineWidth',2); %black line for f=0
colormap([1 1 1; 0.6 0.9 0.3]); %green for f>0 white for f<0

scatter(residents(1), residents(1), 250, 'k','filled',...
    'MarkerEdgeColor','k'); %plots initial value in black
scatter(residents(1:end-1), mutants, 110,'filled'); % Mutant steps
scatter(residents(2:end), residents(2:end), 110,'filled'); % Resident steps
%scatter(residents(end), residents(end), 180, 'red','filled',...
%    'MarkerEdgeColor','k');

scatter(residents2(1), residents2(1), 250, 'k','filled',...
    'MarkerEdgeColor','k'); %plots initial value in black
scatter(residents2(1:end-1), mutants2, 110,C(2,:),'filled'); % Mutant steps
scatter(residents2(2:end), residents2(2:end),110, C(3,:),'filled'); % Resident steps
%scatter(residents2(end), residents2(end), 180, 'red','filled',...
%    'MarkerEdgeColor','k');
xlabel('$\textsf{Resident trait value}$ ($x$)','Interpreter','latex');
ylabel('$\textsf{Mutant trait value}$ ($y$)','Interpreter','latex');
%title('Trait Substitution Sequence on Pairwise Invasibility Plot');
legend('','','Starting trait value','Mutant trait','Resident trait', ...
    '','','','Location','northwest');
axis equal;
hold off;

% Compute the PIP
function [residents, mutants] = simulate_tss_single_mutant(initial_trait, n_steps, step_size, traits, f)
    residents = zeros(n_steps+1, 1);
    mutants = zeros(n_steps, 1);
    residents(1) = initial_trait;
    for i = 1:n_steps
        res = residents(i);
        % Generate a single candidate mutant nearby
        mut = res + (2*rand-1) * step_size;
        mut = min(max(mut, traits(1)), traits(end));%makes sure mutant value is in range of the trait values
        mutants(i) = mut;
        % If mutant is successful, update resident; else, resident doesn't change
        if f(res, mut) > 0
            residents(i+1) = mut;
        else
            residents(i+1) = res;
        end
    end
end

%}

%% Colour plot to show difference in delta* and delta_obs
%{
% parameters
C1 = 1;%2.5;
C2 = 0;      %symmetrical
x0 = 1.5;
beta0 = 0;   %symmetrical

%initial guesses that are close to where the x0+-del/2 ss will be
x1=1;
x2=2;
opts = optimset('Display','off', 'TolX',1e-13, 'TolFun',1e-13, 'MaxIter',200000, 'MaxFunEvals',1000000);

%delta star for a general symmetric system, zmin is the minimum of the competition
%function for mexican hat that is sqrt(3)*sigma
function delstar = delstar(sigma,beta,C1,C2,sigma1,sigma2,x0,zmin)
    k = K(C1,C2,sigma1,sigma2,x0,x0+zmin/2);
    dk = dK(C1,C2,sigma1,sigma2,x0,x0+zmin/2);
    ddk = ddK(C1,C2,sigma1,sigma2,x0,x0+zmin/2);
    a0 = alpha(0,0,sigma,beta);
    az = alpha(zmin,0,sigma,beta);
    dda = ddalpha(zmin,0,sigma,beta);
    A = dk./k;%lnK'
    B = ddk./k-(dk./k).^2;%lnK''
    C = dda./(az+a0);
    delstar = zmin - A./(B/2-C);%zmin-epsilon
end
z=linspace(0.02,1,1000); %range of sigmak and sigmaa
[ska,sk] = meshgrid(z,z);
H = nan(size(sk));
for i = 1:length(z) %y-axis
    sigma1 = z(i)
    for j = 1:length(z) %z-axis
        sigma = z(j);
        %A>kappa condition for branching turns into 3sigmak^2>sigmaa^2
        if 3*sigma1^2>sigma^2 %if branching happens, then we find delta
            zmin = sqrt(3)*sigma;
            dels=delstar(sigma,beta0,C1,C2,sigma1,sigma2,x0,zmin);
            %finds deterministic values of x1 and x2
            sx12 = @(X) sx12_fun(X,r,sigma,beta0,C1,C2,sigma1,sigma2,x0);
            S=fsolve(sx12,[x0-dels/2,x0+dels/2],opts);
            delobs = abs(S(2)-S(1)); %finds delta_obs
            if norm(sx12(S)) > 1e-10 || delobs <1e-6
                S = fsolve(sx12,[x1,x2],opts);
                delobs = abs(S(2)-S(1)); %finds delta_obs
            end
            H(i,j) = abs(delobs-dels);
        else
            H(i,j) = -0.01;
        end
    end
end
figure;
imagesc(z,z,H);
hold on;
set(gca,'YDir','normal'); %makes axis look normal
cmap = [0.7 0.7 0.7; parula(255)];
colormap(cmap);
clim([-0.01 max(H(:))]);
plot(z,z/sqrt(3),'w-','LineWidth',1.5);
colorbar;
ylabel("\sigma_{k}",'Interpreter','tex');
xlabel("\sigma_{\alpha}",'Interpreter','tex');
%}

%% Deltag vs deltam

function deltag=deltag(p)
    deltag=sqrt(2*log(2./p-1));
end

function deltam=deltam(p)
    deltam=sqrt(3)-(sqrt(3)*p*(1-2*exp(-3/2)))./(p*(1-2*exp(-3/2))+12*exp(-3/2));
end

function deltam_exact=deltam_exact(pspan)
    deltam_exact = zeros(size(pspan));
    for i = 1:numel(pspan)
        p=pspan(i);
        if p==0
            deltam_exact(i)=sqrt(3);
        elseif p==3
            deltam_exact(i)=0;
        else
            fun = @(x) p-(6-2*x.^2)./(1-x.^2+exp(x.^2/2));
            deltam_exact(i) = fzero(fun,[0,sqrt(3)]);
        end
    end
end

pmspan = linspace(0,3,20000);
pgspan = linspace(0,1,20000);

figure;
plot(pgspan,deltag(pgspan),'LineWidth',3,'DisplayName',"$\Delta_G^*$");
hold on;
plot(pmspan,deltam(pmspan),'LineWidth',3,'DisplayName',"$\Delta^*$ approximation");
plot(pmspan,deltam_exact(pmspan),':','LineWidth',3,'DisplayName',"$\Delta^*$ exact");
yline(sqrt(3),'--','LineWidth',1.5,'DisplayName',"$\Delta^*=z_{min}$");
xlabel("$\rho = \frac{\sigma_\alpha^2}{\sigma_k^2}$",Interpreter="latex");
ylabel("$\frac{\Delta^*}{\sigma_\alpha}$",Interpreter="latex");
legend('Interpreter','latex');
%{
pmspaninv = linspace(1/3,10000,2000000);
pgspaninv = linspace(1,10000,2000000);

figure;
plot(sqrt(pgspaninv),deltag(1./pgspaninv),'LineWidth',3);
hold on;
plot(sqrt(pmspaninv),deltam(1./pmspaninv),'LineWidth',3);
plot(sqrt(pmspaninv),deltam_exact(1./pmspaninv),':','LineWidth',3);
yline(sqrt(3),'--','LineWidth',1.5);
xlabel("$\sqrt{\frac{1}{\rho}} = \frac{\sigma_k}{\sigma_\alpha}$",Interpreter="latex");
ylabel("$\frac{\Delta^*}{\sigma_\alpha}$",Interpreter="latex");
%}