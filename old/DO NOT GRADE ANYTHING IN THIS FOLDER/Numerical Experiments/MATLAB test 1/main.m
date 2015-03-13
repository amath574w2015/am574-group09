%Solves the shallow water equations with an augmented solver
%Note: this process will be modularized later once tests are complete
%Units:
%x: meters
%t: seconds
function main
clc;
g=9.81;

%setup domain
xmin=-100;
xmax=100;
t0=0;
tfinal=120;

nx=200; %number of cells
cfl=0.9; %Courant number

%very crude upper bound on largest wave speed
c=8;

%set up grids
dx=(xmax-xmin)/nx;
xcenters=linspace(xmin+dx/2,xmax-dx/2,nx)';
nt=ceil(c*(tfinal-t0)/(((xmax-xmin)/nx)*cfl));
dt=(tfinal-t0)/nt;
tgrid=linspace(t0,tfinal,nt+1);


%make topography vector
b=0*(0.02*xcenters).^2;
b=b+0*0.5*cos(0.5*xcenters).^2.*exp(-(0.1*xcenters).^2);

%make initial condition
h_0=max(0,1-b)+0.2*exp(-(0.1*xcenters).^2);
hu_0=0*0.1*h_0;


%plot initial conditions
figure(1);
hold on;
eta=h_0+b;
plot(xcenters,eta,'b');
plot(xcenters,b,'r');
hold off;

h_n=[0;0;h_0;0;0];
hu_n=[0;0;hu_0;0;0];
b=[b(end-1);b(end);b;b(1);b(2)];
%time steps
for n=1:150
    %update ghost cells (periodic)
    h_n(1)=h_n(end-3);
    h_n(2)=h_n(end-2);
    h_n(end-1)=h_n(3);
    h_n(end)=h_n(4);
    hu_n(1)=hu_n(end-3);
    hu_n(2)=hu_n(end-2);
    hu_n(end-1)=hu_n(3);
    hu_n(end)=hu_n(4);
    
    %compute velocity and momentum flux
    u_n=velocity(h_n,hu_n);
    phi_n=momentum_flux(h_n,hu_n);
    
    h_np=h_n;
    hu_np=hu_n;
    
    %for each cell interface
    for i=3:nx+2
        %at the i-1st interface:
        %skip iteration if left and right states are both dry
        if(h_n(i)<1e-16&&h_n(i-1)<1e-16)
            continue;
        end
        
        %make jump vector
        jump=[h_n(i)-h_n(i-1);hu_n(i)-hu_n(i-1);phi_n(i)-phi_n(i-1);b(i)-b(i-1)];
        
        %define wave basis
        
        %stationary wave
        Hbar=0.5*(h_n(i-1)+h_n(i));
        lambar=0.25*(u_n(i-1)+u_n(i)).^2-g*Hbar;
        lamtil=max(0,u_n(i-1)*u_n(i))-g*Hbar;
        Htil=Hbar*lamtil/lambar;
        w0=[g*Hbar/lambar;0;-g*Htil;1];
        
        %other waves
        lambdam=u_n(i-1)-sqrt(g*h_n(i-1));
        lambdap=u_n(i)+sqrt(g*h_n(i));
        %roe linearization values
        uhat=(u_n(i-1)*sqrt(h_n(i-1))+u_n(i)*sqrt(h_n(i)))/(sqrt(h_n(i-1))+sqrt(h_n(i)));
        hbar=0.5*(h_n(i-1)+h_n(i));
        chat=sqrt(g*hbar);
        lambdahatm=uhat-chat;
        lambdahatp=uhat+chat;
        %Einfeldt speeds
        sm=min(lambdam,lambdahatm);
        sp=max(lambdap,lambdahatp);
        %finally define the waves
        ss=zeros(3,1);
        w1=[1;sm;sm^2;0];
        ss(1)=sm;
        w3=[1;sp;sp^2;0];
        ss(3)=sp;
        
        %momentum flux wave
        w2=[0;0;1;0];
        ss(2)=0.5*(sm+sp);
        
        %solve linear system to decompose jump into waves
        R=[w0,w1,w2,w3];
        alpha=R\jump;  %find a better way to do this
        alpha=max(0,alpha)+min(0,alpha); %hack to get rid of NaN
        
        
        %fluctuations
        M=[zeros(2,1),eye(2),zeros(2,1)];
        W=[alpha(1)*w0,alpha(2)*w1,alpha(3)*w2,alpha(4)*w3];
        fwaves=M*W;
        ApQ=[0,0]';
        AmQ=[0,0]';
        for j=1:3
            if ss(j)>0
                ApQ=ApQ+fwaves(:,j);
            elseif ss(j)<0 
                AmQ=AmQ+fwaves(:,j);
            end
        end
        
        ApQ
        AmQ
        
        
        %let's try the first order method first
        
        %add waves to the left cell
        h_np(i-1)=h_np(i-1)-dt/dx*AmQ(1);
        hu_np(i-1)=hu_np(i-1)-dt/dx*AmQ(2);
        
        %add waves to the right cell
        h_np(i)=h_np(i)-dt/dx*ApQ(1);
        hu_np(i)=hu_np(i)-dt/dx*ApQ(2);
        
    end
    
    h_n=max(h_np,0);
    hu_n=hu_np;
    
    
    %close all;
    %continue;
    
    %animated plot
    figure(2);
    clf(2);
    hold on;
    eta=h_n+b;
    plot(xcenters,eta(3:end-2),'b');
    plot(xcenters,b(3:end-2),'r');
    hold off;
end

end


%obtain velocity from height and momentum
function u=velocity(h,hu)

m=numel(h);
u=zeros(m,1);

for i=1:m
    if(h(i)~=0&&abs(hu(i)/h(i))<1e10)
        u(i)=hu(i)/h(i);
    end
end

end

%function for calculating the momentum flux from the shallow water
%equations
function phi=momentum_flux(h,hu)
if h<=0
    phi=0*hu;
else
    g=9.81;
    phi=hu.^2./h+0.5*g*h.^2;
end
end