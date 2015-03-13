function main
global drytol g

%physical parameters
g=9.81;

%dry cell tolerance
drytol=1e-3;

%setup domain
xmin=-50;
xmax=50;
t0=0;
tfinal=180;

%set boundary condition types
%types: wall, extrapolate, periodic
lbc='extrapolate';
rbc='extrapolate';

nx=100; %number of cells
cfl_target=0.9; %target courant Courant number

%set up grids
dx=(xmax-xmin)/nx;
xcenters=linspace(xmin+dx/2,xmax-dx/2,nx)';

%very crude upper bound on largest wave speed (for the first step)
c=5;
dt=cfl_target*dx/c;
%make empty array of times and
tgrid=[0];%linspace(t0,tfinal,nt+1);
current_t=0;
nmax=3000; %maximum number of time steps before stopping

%make topography grid
b=1.3/(xmax-xmin)*(xcenters-xmin);
%b=0.6*exp(-0.005*(xcenters-50).^2).*cos(0.5*xcenters).^2;

%initial conditions
h_0=max(0,1-b)+0.2*exp(-(0.1*xcenters).^2);
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
b=[0;b(1);b;b(end);0];

%time steps
for n=1:nmax
    
    
    %update ghost cells based on boundary condition
    switch lbc
        case 'wall'
            %update ghost cells (solid wall)
            h_n(1)=h_n(3);
            h_n(2)=h_n(3);
            hu_n(1)=-hu_n(3);
            hu_n(2)=-hu_n(3);
        case 'extrapolate'
            %update ghost cells (outflow)
            h_n(1)=h_n(3);
            h_n(2)=h_n(3);
            hu_n(1)=hu_n(3);
            hu_n(2)=hu_n(3);
        otherwise %just use periodic
            h_n(1)=h_n(end-3);
            h_n(2)=h_n(end-2);
            hu_n(1)=hu_n(end-3);
            hu_n(2)=hu_n(end-2);
    end
    
    switch rbc
        case 'wall'
            h_n(end-1)=h_n(end-2);
            h_n(end)=h_n(end-2);
            hu_n(end-1)=-hu_n(end-2);
            hu_n(end)=-hu_n(end-2);
        case 'extrapolate'
            h_n(end-1)=h_n(end-2);
            h_n(end)=h_n(end-2);
            hu_n(end-1)=hu_n(end-2);
            hu_n(end)=hu_n(end-2);
        otherwise %just use periodic
            h_n(end-1)=h_n(3);
            h_n(end)=h_n(4);
            hu_n(end-1)=hu_n(3);
            hu_n(end)=hu_n(4);
    end
    
    
    h_np=h_n;
    hu_np=hu_n;
    
    %make to store wave speeds/fwaves, etc. (note: the array is larger than
    %it needs to be to avoid index confusion.)
    s=zeros(nx+4,3);  %s(grid point to right, wave number)
    fwave=zeros(nx+4,3,2); %fwave(grid point to right, wave number, component)
    wallL=ones(nx+4);
    wallR=ones(nx+4);
    
    %loop over cell interfaces, updating left and right cells
    %note that we compute fwaves even in boundary ghost cells so that we
    %can use limiters.
    for i=2:nx+4
        
        hL=h_n(i-1);
        hR=h_n(i);
        huL=hu_n(i-1);
        huR=hu_n(i);
        bL=b(i-1);
        bR=b(i);
        
        %compute velocities/momentum fluxes
        if hL<drytol
            hL=0;
            huL=0;
            uL=0;
            phiL=0;
        else
            uL=huL/hL;
            phiL=huL^2/hL+0.5*g*hL^2;
        end
        if hR<drytol
            hR=0;
            huR=0;
            uR=0;
            phiR=0;
        else
            uR=huR/hR;
            phiR=huR^2/hR+0.5*g*hR^2;
        end
        
        %If both states are dry, skip this interface
        if hR<drytol&&hL<drytol
            continue;
        end
        
        %wall indicator: 1 indicates to update the cell; 2 indicates to
        %leave it as a wall, so no update
        
        %If one side is dry, do something special
        
        if hL<drytol
            
            %just do a solid wall problem
            uhat=0;
            chat=sqrt(g*hR);
            lamm_left=-uR-sqrt(g*hR);
            lamp_right=uR+sqrt(g*hR);
            roe_lamm=uhat-chat;
            roe_lamp=uhat+chat;
            
            s1=min(lamm_left,roe_lamm);%Einfeldt speed
            s3=max(lamp_right,roe_lamp);%Einfeldt speed
            
            
            %estimate the middle state for the dry wall problem (HLLE
            %estimate; maybe not the best choice?)
            Hstar=(-2*huR+(s3-s1)*hR)/(s3-s1);
            %if the middle state does not exceed the jump in bathymetry,
            %continue the solid wall problem and use this to update the
            %left cell
            
            if Hstar<bL-bR
                %repose the problem, then continue on
                hL=hR;
                huL=-huR;
                uL=-uR;
                phiL=phiR;
                bL=bR;
                
                %set wall flag not to update the left cell
                wallL(i)=0;
                
                %Now, run the normal method with this ghost cell problem.
            end
            %otherwise, run the normal method
            
        elseif hR<drytol
            
            %just do a solid wall problem
            uhat=0;
            chat=sqrt(g*hL);
            lamm_left=uL-sqrt(g*hL);
            lamp_right=-uL+sqrt(g*hL);
            roe_lamm=uhat-chat;
            roe_lamp=uhat+chat;
            
            s1=min(lamm_left,roe_lamm);%Einfeldt speed
            s3=max(lamp_right,roe_lamp);%Einfeldt speed
            
            
            %estimate the middle state for the dry wall problem
            Hstar=(2*huL+(s3-s1)*hL)/(s3-s1);
            %if the middle state does not exceed the jump in bathymetry,
            %continue the solid wall problem and use this to update the
            %left cell
            
            if Hstar<bR-bL
                %repose the problem
                hR=hL;
                huR=-huL;
                uR=-uL;
                phiR=phiL;
                bR=bL;
                
                %set wall flag not to update the right cell
                wallR(i)=0;
                
                %Now, run the normal method with this ghost cell problem
            end
            %otherwise, run the normal method
        end
        
        %Beyond this point, just carry out the Roe solver.
        %can assume at least one cell is not dry
        
        %Roe wave speeds:
        hbar=0.5*(hL+hR);
        %PROBLEM IF DRY:
        uhat=(sqrt(hL)*uL+sqrt(hR)*uR)/(sqrt(hL)+sqrt(hR));
        chat=sqrt(g*hbar);
        
        %eigenvalues and eigenvectors
        lamm_left=uL-sqrt(g*hL);%outflow at left
        lamp_right=uR+sqrt(g*hR);%outflow at right
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
        delta=[hR-hL;huR-huL;phiR-phiL;bR-bL];
        
        %first, compute/subtract off the stationary jump
        lamlambar=0.25*(uL+uR)^2-0.5*g*(hL+hR);
        lamlamtil=max(0,uR*uL)-0.5*g*(hL+hR);
        %PROBLEM IF DRY:
        htil=hbar*lamlamtil/lamlambar;
        w0=[g*hbar/lamlambar;0;-g*htil;1];
        %apply bounds to the first term of w0 to prevent negative depth
        if bR>bL
            %PROBLEM IF DRY:
            lb=(huL-huR+s3*hR-s1*hL)/s1;
        else
            %PROBLEM IF DRY:
            lb=(huL-huR+s3*hR-s1*hL)/s3;
        end
        w0(1)=max(w0(1),lb);
        w0(1)=min(w0(1),-1);
        %get wave
        W0=(bR-bL)*w0;
        
        %subtract off this wave and only decompose what remains
        delta=delta-W0;
        
        
        %decompose the jump into waves
        alpha1=-(s3*delta(1)-delta(2))/(s1-s3);
        alpha3=-(delta(2)-s1*delta(1))/(s1-s3);
        alpha2=delta(3)-(s1+s3)*delta(2)+s1*s3*delta(1);
        
        
        %define waves
        W1=alpha1*r1;
        W2=alpha2*r2;
        W3=alpha3*r3;
        
        
        %define flux waves
        Z1=W1(2:3);
        Z2=W2(2:3);
        Z3=W3(2:3);
        
        %add flux waves and wave speeds to the relevant arrays
        s(i,1)=s1; s(i,2)=s2; s(i,3)=s3;
        fwave(i,1,:)=Z1; fwave(i,2,:)=Z2; fwave(i,3,:)=Z3;
        
        %Compute speeds and unlimited flux waves for all other interfaces
        %before updating cell values/applying limiters
    end
    
    
    %compute Courant number and new time step based on it
    largest_s=max(max(abs(s)));
    cfl=largest_s*dt/dx;
    dtn=cfl_target*dx/largest_s;
    
    %check to see if the CFL condition was held
    if cfl>=1
        strcat('Repeating step with smaller time step.  CFL=',num2str(cfl),' at t=',num2str(current_t))
        %if the cfl condition is not held, then repeat the step
        dt=dtn;
        continue;
    end
    
    
    %Now that we have determined that we will not have to retake the step
    %and we have fwaves at each interface, loop over them to determine
    %limited values and update appropriately
    
    for i=3:nx+3
        %compute fluctuations
        AmQ=[0;0];
        ApQ=[0;0];
        for p=1:3
            if s(i,p)>0
                ApQ=ApQ+squeeze(fwave(i,p,:));
            else
                AmQ=AmQ+squeeze(fwave(i,p,:));
            end
        end
        
        %apply limiters to waves
        Z1=squeeze(fwave(i,1,:));
        Z2=squeeze(fwave(i,2,:));
        Z3=squeeze(fwave(i,3,:));
        
        I1=i-sign(s(i,1));
        I2=i-sign(s(i,2));
        I3=i-sign(s(i,3));
        theta1=(dot(squeeze(fwave(I1,1,:)),Z1))/norm(Z1,2);
        theta2=(dot(squeeze(fwave(I2,2,:)),Z2))/norm(Z2,2);
        theta3=(dot(squeeze(fwave(I3,3,:)),Z3))/norm(Z3,2);
        
        lim1=limiter(theta1);
        lim2=limiter(theta2);
        lim3=limiter(theta3);
        
        Z1lim=lim1.*Z1; Z2lim=lim2.*Z2; Z3lim=lim3.*Z3; 
        %use limited waves for second order corrections
        F_corrector=[0;0];
        F_corrector=F_corrector+0.5*(sign(s1)*(1-abs(s1)*dt/dx)*Z1lim);
        F_corrector=F_corrector+0.5*(sign(s2)*(1-abs(s2)*dt/dx)*Z2lim);
        F_corrector=F_corrector+0.5*(sign(s3)*(1-abs(s3)*dt/dx)*Z3lim);
        
        %update left and right cells
        h_np(i-1)=h_np(i-1)-wallL(i)*dt/dx*(AmQ(1)+F_corrector(1));
        hu_np(i-1)=hu_np(i-1)-wallL(i)*dt/dx*(AmQ(2)+F_corrector(2));
        h_np(i)=h_np(i)-wallR(i)*dt/dx*(ApQ(1)-F_corrector(1));
        hu_np(i)=hu_np(i)-wallR(i)*dt/dx*(ApQ(2)-F_corrector(2));
    end
    
    
    %update values
    h_n=h_np;
    hu_n=hu_np;
    
    
    %reduce all cells below the dry tolerance to zero
    for i=3:numel(h_n-2)
        if h_n(i)<drytol
            h_n(i)=0;
            hu_n(i)=0;
        end
    end
    
    %update time
    current_t=current_t+dt;
    dt=dtn;
    tgrid(numel(tgrid)+1)=current_t;
    
    
    %plot current frame
    %animated plot
    figure(2);
    clf(2);
    title(strcat('t=',num2str(current_t)))
    hold on;
    eta=h_n+b;
    plot(xcenters,eta(3:end-2),'b');
    plot(xcenters,b(3:end-2),'r');
    ylim([0,2]);
    hold off;
    if n==2
        pause;
    end
    
    
    if current_t>tfinal
        return;
    end
    continue;
    
    figure(3);
    clf(3);
    hold on;
    plotyy(xcenters,hu_n(3:end-2),xcenters,h_n(3:end-2));
    hold off;
    
end

end

%function that applies the MC limiter
function phi=limiter(theta)
%phi=[0;0];
phi=max(min(min(2*theta,0.5*(1+theta)),2),0);
end