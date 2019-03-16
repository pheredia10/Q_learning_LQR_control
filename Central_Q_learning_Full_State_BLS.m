close all
clc
clear
% Initial Structures
T=200; % number of time steps
dt=1e-2; % time step size
time=T*dt; %elapsed time
t_vector=0:dt:time;
m=1; % number of agents
n=2; % number of states
n_bar=1; % number of control commands
xx=n*(n+1)/2; %combination of states with replacement
uu=n_bar*(n_bar+1)/2;  % n_bar*(n_bar+1)/2 is combinations of control commands with replacement
tau=(n+n_bar)*(n+n_bar+1)/2; % number of weights
w=zeros(tau,T,m); % weights , Note: weights that correspnd to x*x elements in phi(x) go with S_xx and so on
x=zeros(n,T,m); %states
psi=zeros(tau,T+1,m); % combined weights
mu=.0000001*ones(1,T); % optimization algorithm step-size
epsilon=.01; % convergance tolerance  .0000000001;
h_coef=zeros(n_bar,n,m); % control policy coefficient

% Initial Conditions
x(:,1,:)=[10 -10 ]';
k=ones(m,1);
zbar(1:tau,1:21)=0;
d_target(1:21,1)=0;
H=[1 0 0;0 1 0;1 0 1];
Hyy=H(3,3);Hyx=H(3,1:2);
h_coef(:,:,1)=-inv(Hyy)*Hyx;
L=-inv(Hyy)*Hyx;


R=1;P=eye(2,2);G=[P [0 ; 0];[0 0] R];

for t=1:T
    
    for i =1:m
        
        if t<50
            a1=0;
            a2=0.97;
        else
            a1=0;
            a2=0;
        end
        
        
        x_prev=x(:,t,i);
        
        
        u_prev=h_coef(:,:,i)*x(:,t,i) + a2*(0.5*sin(2.0*t)^2*cos(10.1*t)+0.9*sin(1.102*t)^2*cos(4.001*t)+0.3*sin(1.99*t)^2*cos(7*t)+0.3*sin(10.0*t)^3+0.7*sin(3.0*t)^2*cos(4.0*t)+0.3*sin(3.00*t)*1*cos(1.2*t)^2+0.400*sin(1.12*t)^2+0.5*cos(2.4*t)*sin(8*t)^2+0.3*sin(1.000*t)^1*cos(0.799999*t)^2+0.3*sin(4*t)^3+0.4*cos(2*t)*1*sin(5*t)^4+0.3*sin(10.00*t)^3);
        
        
        A=[1 1;0.16 1]; % Unstable
        B=[0;1];
        
        
        % Discrete Simulation
        x(:,t+1,i)=A*x_prev+B*u_prev;
        
        %                 % Continous Simulation
        %         dx=A*x_prev+B*u_prev;
        %         x(:,t+1,i)=x(:,t,i)+dx*dt;
        
        x_now=x(:,t+1,i); % updated states
        
        u_now=h_coef(:,:,i)*x_now;
        
        
        if t==1
            d0= [x_prev; u_prev]'*G*[x_prev; u_prev]+[x_now; u_now]'*H(:,:,i)*[x_now; u_now];   %d_target(8,1);   r(x_prev,u_prev)
            Jnew(i)=(d0-[x_prev; u_prev]'*H(:,:,i)*[x_prev; u_prev])^2;
        end
        
        if t<50
            
            d_target(1,1)=d_target(2,1);
            d_target(2,1)=d_target(3,1);
            d_target(3,1)=d_target(4,1);
            d_target(4,1)=d_target(5,1);
            d_target(5,1)=d_target(6,1);
            d_target(6,1)=d_target(7,1);
            d_target(7,1)= [x_prev; u_prev]'*G*[x_prev; u_prev]+[x_now; u_now]'*H*[x_now; u_now];   %d_target(8,1);   r(x_prev,u_prev)
            
            
            zbar(:,1)=zbar(:,2);
            zbar(:,2)=zbar(:,3);
            zbar(:,3)=zbar(:,4);
            zbar(:,4)=zbar(:,5);
            zbar(:,5)=zbar(:,6);
            zbar(:,6)=zbar(:,7);
            zbar(:,7)=[x_prev(1)^2;x_prev(1)*x_prev(2);x_prev(1)*u_prev;x_prev(2)^2;x_prev(2)*u_prev;u_prev^2];%phi(x_prev,u_prev)  %  zbar(:,8);
            
            
            if mod(t,7)==0
                zbar2=zbar*zbar';
                q=zbar*d_target;
                w(:,t+1,i)=inv(zbar2)*q;
                
                H=[w(1,t+1,i) w(2,t+1,i)/2 w(3,t+1,i)/2; w(2,t+1,i)/2 w(4,t+1,i) w(5,t+1,i)/2;w(3,t+1,i)/2 w(5,t+1,i)/2 w(6,t+1,i)];
                Hyy=H(3,3);Hyx=H(3,1:2);
                h_coef(:,:,i)=-inv(Hyy)*Hyx;
                Jnew(i)=(d_target(7,1,i)-[x_prev; u_prev]'*H(:,:)*[x_prev; u_prev])^2;
                k(i)=k(i)+1;
                
                
            end
            
        end
        J(t,i)=Jnew(i);
    end
end

plot((1:T+1)-1,x(:,1:T+1,1))
xlabel('Time Step')
ylabel('System States')

figure
loglog((1:T)-1,(J(1:T,1)))
xlabel('Time Step')
ylabel('J(w)')



