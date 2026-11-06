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

%z1=(x1-x0)/sigma
%z2=(x2-x0)/sigma
%assuming alpha and K symmetric

%% alpha function

function alpha=alpha(z1,z2,gamma)
    alpha = (1-gamma*(z1-z2).^2).*exp(-(z1-z2).^2./2);
end

%dalpha(x-y)/dx, dalpha(y,x,...) = dalpha(y-x)/dy = -dalpha(y-x)/dx
function dalpha=dalpha(z1,z2,gamma)
    dalpha = (z1-z2).*(-1-2*gamma+gamma*(z1-z2).^2).*exp(-(z1-z2).^2./2);
end
%d^2alpha(x-y)/dx^2, ddalpha(y,x,...) = d^2alpha(y-x)/dy^2
function ddalpha=ddalpha(z1,z2,gamma)
z = (z1-z2);
    ddalpha = exp(-z.^2./2).* ...
(-1-2*gamma+z.^2.*(1+5*gamma-gamma*z.^2));
end

%% Double Gaussian Carrying capacity K(x)
%K is sum of 2 gaussians
function K=K(C,rho,z)
    K = C*exp(-z.^2.*rho/2);
end
%derivative of K (when alpha peaks at y=x, the solutions to this are x*)
function dK=dK(C,rho,z)
    dK = -C*z.*rho.*exp(-z.^2.*rho/2);
end

%second derivative of K
function ddK=ddK(C,rho,z)
    ddK = C*rho.*(z.^2.*rho-1).*exp(-z.^2.*rho/2);
end
%% N
%N2=N1(C,rho,gamma,z2,z1)
function N1=N1(C,rho,gamma,z1,z2)
    a=alpha(z1,z2,gamma);
    denom = 1-a.^2;
    N1 = (K(C,rho,z1)-K(C,rho,z2).*a)./denom;
end
%dN1/dz1, dN2/dz2 = dN1z1(C,rho,gamma,z2,z1)
function dN1z1=dN1z1(C,rho,gamma,z1,z2)
    a=alpha(z1,z2,gamma);
    da=dalpha(z1,z2,gamma);
    denom = 1-a.^2;
    dN1z1 = (dK(C,rho,z1)-K(C,rho,z2).*da)./denom + 2*(K(C,rho,z1)-K(C,rho,z2).*a).*a.*da./denom.^2;
end
%dN1/dz2, dN2/dz1 = dN1z2(C,rho,gamma,z2,z1)
function dN1z2=dN1z2(C,rho,gamma,z1,z2)
    a=alpha(z1,z2,gamma);
    da=dalpha(z1,z2,gamma);
    denom = 1-a.^2;
    dN1z2 = (K(C,rho,z2).*da-dK(C,rho,z2).*a)./denom - 2*(K(C,rho,z1)-K(C,rho,z2).*a).*a.*da./denom.^2;
end

%% selection gradients and deltas

%simplified reduction of selection gradients=0 eqn
function S = sx12_fun(Z,rho,gamma,C)
    z1 = Z(1);
    z2 = Z(2);

    denom = 1-alpha(z1,z2,gamma)^2;

    if abs(denom) < 1e-10
        S = [1e6; 1e6];
        return
    end
    N1 = (K(C,rho,z1)-K(C,rho,z2).*alpha(z1,z2,gamma))./denom;
    N2 = (K(C,rho,z2)-K(C,rho,z1).*alpha(z1,z2,gamma))./denom;
    
    S = [ ...
        dK(C,rho,z1)-dalpha(z1,z2,gamma).*N2; ...
        dK(C,rho,z2)+dalpha(z1,z2,gamma).*N1 ...
    ];
end
%K1/r d^2f12/dy^2(y=x1), K2/r d^2f12/dy^2(y=x2)=ddfdy([z2,z1],rho,gamma,C)
function ddfdy= ddfdy(Z,rho,gamma,C)
    z1=Z(1);
    z2=Z(2);
    ddfdy = -ddalpha(0,0,gamma).*N1(C,rho,gamma,z1,z2)-ddalpha(z1,z2,gamma).*N1(C,rho,gamma,z2,z1)+ddK(C,rho,z1);
end
%1/r dS1/dz1, 1/r dS2/dz2 = dS1dz1([z2,z1],rho,gamma,C)
function dS1dz1 = dS1dz1(Z,rho,gamma,C)
    z1=Z(1);
    z2=Z(2);
    K1=K(C,rho,z1);
    dK1=dK(C,rho,z1);
    ddK1=ddK(C,rho,z1);
    da=dalpha(z1,z2,gamma);
    dda=ddalpha(z1,z2,gamma);
    N2 = N1(C,rho,gamma,z2,z1);
    dN2dz1 = dN1z2(C,rho,gamma,z2,z1);

    dS1dz1 = -dda.*N2./K1-da.*(dN2dz1./K1-N2.*dK1./K1.^2)+ddK1./K1-(dK1./K1).^2;
end
%1/r dS1/dz2, 1/r dS2/dz1 = dS1dz2([z2,z1],rho,gamma,C)
function dS1dz2 = dS1dz2(Z,rho,gamma,C)
    z1=Z(1);
    z2=Z(2);
    K1=K(C,rho,z1);
    da=dalpha(z1,z2,gamma);
    dda=ddalpha(z1,z2,gamma);
    N2 = N1(C,rho,gamma,z2,z1);
    dN2dz2 = dN1z1(C,rho,gamma,z2,z1);

    dS1dz2 = dda.*N2./K1-da.*dN2dz2./K1;
end
%Delta/sigma. for symmetric branching z1=-delta/2sigma, z2=delta/2sigma
function deltam_exact=deltam_exact(rho,gamma)
    if gamma==0
        deltam_exact=sqrt(2*log(2./rho-1));
    elseif rho==0
        deltam_exact=sqrt(2+1./gamma);
    elseif rho==1+2*gamma
        deltam_exact=0;
    else
        fun = @(x) rho-2*(1+2*gamma-gamma.*x.^2)./(1-gamma.*x.^2+exp(x.^2/2));
        deltam_exact = fzero(fun,[0,sqrt(2+1./gamma)]);
    end
end
%% parameters
C1 = 1;%2.5;
C2 = 0;%2;
sigma1 = sqrt(0.5);%0.5;%sqrt(0.5/(1+2gamma))
sigma2 = 0.15;%0.35;%
x0 = 1.5;
sigma = 0.5; %0.25;%width
beta0 = 0;   %symmetrical
gamma = 0;

%% Colour plot to show how many singular strategies there are
%
rhospan=linspace(0,3,5000); 
gammaspan = linspace(0,1,5000);
H = zeros(length(rhospan),length(gammaspan));
[height,width]=size(H);
% finding if symmetric sc is convergent or not, then if convergent sc is stable or not
%
opts = optimoptions('fsolve','Display','off');
for i = 1:height
    fprintf("i=%d\n",i);
    for j = 1:width
        solfound=false;
        rho = rhospan(i);
        gamma = gammaspan(j);
        if rho>=1+2*gamma %if no branching
            H(i,j)=0;
        else
            sx12 = @(Z) sx12_fun(Z,rho,gamma,C1);
            delta = deltam_exact(rho,gamma);
            %symmetric solutions for z1=a, z2=b
            a=-delta/2;
            b=delta/2;
            H(i,j)=1;%always 1 for gamma==0, always at least 1 in general
            if gamma>0 && rho>0 %rho=0: K is flat, so shifted copies of the symmetric solution all count as one
                Z1=[a,b];
                Z2=[b,a];
                detJ = dS1dz1(Z1,rho,gamma,C1).*dS1dz1(Z2,rho,gamma,C1)-dS1dz2(Z1,rho,gamma,C1).*dS1dz2(Z2,rho,gamma,C1);
                tr = dS1dz1(Z1,rho,gamma,C1)+dS1dz1(Z2,rho,gamma,C1);
                if detJ>0 && tr<0 %if convergent, z1=a,z2=b
                    F1yy = ddfdy(Z1,rho,gamma,C1);
                    F2yy = ddfdy(Z2,rho,gamma,C1);
                    if F1yy<0 && F2yy<0
                        H(i,j)=2;
                    end
                else %if not convergent (det>0&tr<0)
                    H(i,j)=3;
                    z1grid = linspace(0,a,10);
                    for k=1:length(z1grid)
                        z1=z1grid(k);
                        z2grid = linspace(z1+sqrt(2+1./gamma),b,10);
                        for l=1:length(z2grid)
                            z2=z2grid(l);
                            [Zsol,~,exitflag]=fsolve(sx12,[z1,z2],opts); %Zsol=[z1;z2] of a candidate singular coalition
                            res = norm(sx12(Zsol)./[K(C1,rho,Zsol(1)); K(C1,rho,Zsol(2))]); %is it really a root? (this checks actaul s, since the functions only check the numerator=0)
                            al12 = alpha(Zsol(1),Zsol(2),gamma);
                            N1s = (K(C1,rho,Zsol(1))-K(C1,rho,Zsol(2))*al12)/(1-al12^2); %densities of the coalition
                            N2s = (K(C1,rho,Zsol(2))-K(C1,rho,Zsol(1))*al12)/(1-al12^2);
                            zs = sort(Zsol); %order doesn't matter: (b,a) is the same coalition as (a,b)
                            isNew = abs(zs(1)-a)>1e-4 || abs(zs(2)-b)>1e-4; %differs from the symmetric solution
                            if exitflag>0 && res<1e-8 && N1s>0 && N2s>0 && abs(Zsol(2)-Zsol(1))>1e-3 && isNew %if (z1,z2) is indeed a new singular coalition
                                Z1=zs;
                                c=zs(1);
                                d=zs(2);
                                Z2=[d,c];
                                F1yy = ddfdy(Z1,rho,gamma,C1);
                                F2yy = ddfdy(Z2,rho,gamma,C1);
                                if F1yy<0 && F2yy<0
                                    H(i,j)=4;
                                end
                                solfound=true;
                                break;
                            end
                        end
                        if solfound
                            break;
                        end
                    end
                end
            end
        end
    end
end
%}

%by finding if symmetric sc is convergent
%{
for i = 1:height
    fprintf("i=%d\n",i);
    for j = 1:width
        rho = rhospan(i);
        gamma = gammaspan(j);
        if rho>=1+2*gamma %if no branching
            H(i,j)=0;
        else
            delta = deltam_exact(rho,gamma);
            %symmetric solutions for z1=a, z2=b
            a=-delta/2;
            b=delta/2;
            H(i,j)=1;%always 1 for gamma==0, always at least 1 in general
            if gamma>0 && rho>0 %rho=0: K is flat, so shifted copies of the symmetric solution all count as one
                Z1=[a,b];
                Z2=[b,a];
                detJ = dS1dz1(Z1,rho,gamma,C1).*dS1dz1(Z2,rho,gamma,C1)-dS1dz2(Z1,rho,gamma,C1).*dS1dz2(Z2,rho,gamma,C1);
                tr = dS1dz1(Z1,rho,gamma,C1)+dS1dz1(Z2,rho,gamma,C1);
                if detJ<0 || tr>0 %if not convergent (det>0&tr<0)
                    H(i,j)=3;
                end
            end
        end
    end
end
%}
%plots number of sign changes (number of singular strategies)
figure;
imagesc(gammaspan,rhospan,H);
hold on;
%colormap sky
colormap([0.9,0.9,0.9;0.6 0.9 0.4; 1 0.9 0.5;0.4 0.8 0.1; 0.7 0.8 0.8])%grey no branch;light green sym branch;yellow stable sym;dark green asym;blue stable asym
colorbar;
plot(0,25/27,'k.',"MarkerSize",20);
text(0,25/27,"A","FontSize",20);
plot(0.25,25/18,'k.',"MarkerSize",20);
text(0.25,25/18,"B","FontSize",20,"HorizontalAlignment","right");
plot(0.5,50/27,'k.',"MarkerSize",20);
text(0.5,50/27,"C","FontSize",20,"HorizontalAlignment","right");
plot(0.75,125/54,'k.',"MarkerSize",20);
text(0.75,125/54,"D","FontSize",20,"HorizontalAlignment","right");
plot(1,25/9,'k.',"MarkerSize",20);
text(1,25/9,"E","FontSize",20,"HorizontalAlignment","right");
xlabel("\textsf{Strength of facilitation } ($\gamma$)",'Interpreter','latex');
ylabel("$\rho=\frac{\sigma_\alpha^2}{\sigma_k^2}$",'Interpreter','latex');
set(gca,'YDir','normal'); %makes axis look normal 