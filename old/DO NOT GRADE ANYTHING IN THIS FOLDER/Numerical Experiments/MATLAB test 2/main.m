function main
%physical parameters
g=9.81;

%setup domain
xmin=-50;
xmax=150;
t0=0;
tfinal=300;

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

%make topography grid
b=0.5*exp(-0.02*(xcenters-50).^2);

%initial conditions
h_0=max(0,1-b)+0.8*exp(-(0.1*xcenters).^2);
hu_0=0*0.1*h_0;

%plot initial conditions
figure(1);
hold on;
eta=h_0+b;
plot(xcenters,eta,'b');
plot(xcenters,b,'r');
hold off;


%add ghost cells
h_n=[0;0;h_0;0;0];
hu_n=[0;0;hu_0;0;0];
b=[b(end-1);b(end);b;b(1);b(2)];

%time steps
for n=1:nt
    
    
    %update ghost cells (periodic)
    h_n(1)=h_n(end-3);
    h_n(2)=h_n(end-2);
    h_n(end-1)=h_n(3);
    h_n(end)=h_n(4);
    hu_n(1)=hu_n(end-3);
    hu_n(2)=hu_n(end-2);
    hu_n(end-1)=hu_n(3);
    hu_n(end)=hu_n(4);
    
    
    %compute velocity:
    %PROBLEM IF DRY
    u_n=velocity(h_n,hu_n);
    %compute momentum flux:
    %PROBLEM IF DRY
    phi_n=momentum_flux(h_n,hu_n);
    
    
    h_np=h_n;
    hu_np=hu_n;
    
    %loop over cell interfaces, updating left and right cells
    for i=3:nx+3
        
        %Roe wave speeds:
        hbar=0.5*(h_n(i-1)+h_n(i));
        %PROBLEM IF DRY:
        uhat=(sqrt(h_n(i-1))*u_n(i-1)+sqrt(h_n(i))*u_n(i))/(sqrt(h_n(i-1))+sqrt(h_n(i)));
        chat=sqrt(g*hbar);
        
        %eigenvalues and eigenvectors
        lamm_left=u_n(i-1)-sqrt(g*h_n(i-1));%outflow at left
        lamp_right=u_n(i)+sqrt(g*h_n(i));%outflow at right
        roe_lamm=uhat-chat;
        roe_lamp=uhat+chat;
        s1=min(lamm_left,roe_lamm);%Einfeldt speed
        s3=max(lamp_right,roe_lamp);%Einfeldt speed
        r1=[1;s1;s1^2];
        r3=[1;s3;s3^2];
        
        %momentum flux wave
        r2=[0;0;1];
        s2=0.5*(s1+s3);
        
        
        %the jump:
        delta=[h_n(i)-h_n(i-1);hu_n(i)-hu_n(i-1);phi_n(i)-phi_n(i-1);b(i)-b(i-1)];
        
        %first, compute/subtract off the stationary jump
        lamlambar=0.25*(u_n(i-1)+u_n(i))^2-0.5*g*(h_n(i-1)+h_n(i));
        lamlamtil=max(0,u_n(i)*u_n(i-1))-0.5*g*(h_n(i-1)+h_n(i));
        %PROBLEM IF DRY:
        htil=hbar*lamlamtil/lamlambar;
        w0=[g*hbar/lamlambar;0;-g*htil;1];
        %apply bounds to the first term of w0 to prevent negative depth
        if b(i)>b(i-1)
            %PROBLEM IF DRY:
            lb=(hu_n(i-1)-hu_n(i)+s3*h_n(i)-s1*h_n(i-1))/s1;
        else
            %PROBLEM IF DRY:
            lb=(hu_n(i-1)-hu_n(i)+s3*h_n(i)-s1*h_n(i-1))/s3;
        end
        w0(1)=max(w0(1),lb);
        w0(1)=min(w0(1),-1);
        %get wave
        W0=(b(i)-b(i-1))*w0;
        
        %subtract off this wave and only decompose what remains
        delta=delta-W0;
        
        
        
        
        %decompose the jump into waves
        %PROBLEMS (2) IF DRY:
        alpha1=((uhat+chat)*delta(1)-delta(2))/(2*chat);
        alpha3=(-(uhat-chat)*delta(1)+delta(2))/(2*chat);
        alpha2=delta(3)-alpha1*s1^2-alpha3*s3^2;
        
        %define waves
        W1=alpha1*r1;
        W2=alpha2*r2;
        W3=alpha3*r3;
        
        %define fluctuations
        Z1=W1(2:3);
        Z2=W2(2:3);
        Z3=W3(2:3);
        
        
        AmQ=[0;0];
        ApQ=[0;0];
        
        %AmQ=AmQ+min(s1,0)*W1+min(s3,0)*W3;
        %ApQ=ApQ+max(s1,0)*W1+max(s3,0)*W3;
        AmQ=AmQ+(s1<0)*Z1+(s2<0)*Z2+(s3<0)*Z3;
        ApQ=ApQ+(s1>0)*Z1+(s2>0)*Z2+(s3>0)*Z3;
        
        
        %update left and right cells
        h_np(i-1)=h_np(i-1)-dt/dx*AmQ(1);
        hu_np(i-1)=hu_np(i-1)-dt/dx*AmQ(2);
        h_np(i)=h_np(i)-dt/dx*ApQ(1);
        hu_np(i)=hu_np(i)-dt/dx*ApQ(2);
        
    end
    
    %update values
    h_n=h_np;
    hu_n=hu_np;
    
    
    %plot current frame
    %animated plot
    figure(2);
    clf(2);
    title(strcat('t=',num2str(tgrid(n))))
    hold on;
    eta=h_n+b;
    plot(xcenters,eta(3:end-2),'b');
    plot(xcenters,b(3:end-2),'r');
    hold off;
    
end

end




%obtain velocity from height and momentum
function u=velocity(h,hu)

u=hu./h;

end

%calculate the momentum flux
function phi=momentum_flux(h,hu)
g=9.81;
phi=hu.^2./h+0.5*g*h.^2;
end