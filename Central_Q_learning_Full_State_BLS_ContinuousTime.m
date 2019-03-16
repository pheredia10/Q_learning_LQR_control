close all
clc
clear
% Initial Structures
T=5000; % number of time steps
dt=1e-3; % time step size
time=T*dt; %elapsed time
t_vector=0:dt:time;
m=1; % number of agents
n=2; % number of states
n_bar=1; % number of control commands
xx=n*(n+1)/2; %combination of states with replacement
uu=n_bar*(n_bar+1)/2;  % n_bar*(n_bar+1)/2 is combinations of control commands with replacement
tau=(n+n_bar)*(n+n_bar+1)/2; % number of weights
w=zeros(tau,T,m); % weights , Note: weights that correspnd to x*x elements in phi(x) go with S_xx and so on
%C=doubly_stochastic(m); % communication matrix
x=zeros(n,T,m); %states
psi=zeros(tau,T+1,m); % combined weights
mu=.0000001*ones(1,T); % optimization algorithm step-size
epsilon=1; % convergance tolerance  .0000000001;
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
sample=10;


t=2;
count=0;
for tc=1:T
    
    for i =1:m
        
        if count<200
            a1=0;
            a2=0.97;
        else
            a1=0;
            a2=0;
        end
        
        
        x_prev=x(:,tc,i);
        
        u_prev=h_coef(:,:,i)*x(:,tc,i) + a2*(0.5*sin(2.0*tc)^2*cos(10.1*tc)+0.9*sin(1.102*tc)^2*cos(4.001*tc)+0.3*sin(1.99*tc)^2*cos(7*tc)+0.3*sin(10.0*tc)^3+0.7*sin(3.0*tc)^2*cos(4.0*tc)+0.3*sin(3.00*tc)*1*cos(1.2*tc)^2+0.400*sin(1.12*tc)^2+0.5*cos(2.4*tc)*sin(8*tc)^2+0.3*sin(1.000*tc)^1*cos(0.799999*tc)^2+0.3*sin(4*tc)^3+0.4*cos(2*tc)*1*sin(5*tc)^4+0.3*sin(10.00*tc)^3);
        u(tc)=u_prev;
        A=[1 1;0.16 1]; % Unstable
        B=[0;1];
        
        
        %         % Discrete Simulation
        %         x(:,t+1,i)=A*x_prev+B*u_prev;
        
        % Continous Simulation
        dx=A*x_prev+B*u_prev;
        x(:,tc+1,i)=x(:,tc,i)+dx*dt;
        
        
        
        if mod(tc,sample)==0
            t=tc-sample+1;
            r_integral=0;
            for ti=t:t+sample-1
                x_prev=x(:,ti,i);
                u_prev=u(ti);
                r_integral= r_integral+([x_prev; u_prev]'*G*[x_prev; u_prev])*dt;
            end
            r_integral_vec(count+1)=r_integral;
            x_prev=x(:,t,i);
            u_prev=u(t);%h_coef(:,:,i)*x(:,t,i) + a2*(0.5*sin(2.0*t)^2*cos(10.1*t)+0.9*sin(1.102*t)^2*cos(4.001*t)+0.3*sin(1.99*t)^2*cos(7*t)+0.3*sin(10.0*t)^3+0.7*sin(3.0*t)^2*cos(4.0*t)+0.3*sin(3.00*t)*1*cos(1.2*t)^2+0.400*sin(1.12*t)^2+0.5*cos(2.4*t)*sin(8*t)^2+0.3*sin(1.000*t)^1*cos(0.799999*t)^2+0.3*sin(4*t)^3+0.4*cos(2*t)*1*sin(5*t)^4+0.3*sin(10.00*t)^3);
            x_now=x(:,t+sample,i); % updated states
            u_now=h_coef(:,:,i)*x_now;
            if count<200
                
                d_target(1,1)=d_target(2,1);
                d_target(2,1)=d_target(3,1);
                d_target(3,1)=d_target(4,1);
                d_target(4,1)=d_target(5,1);
                d_target(5,1)=d_target(6,1);
                d_target(6,1)=d_target(7,1);
                d_target(7,1)= r_integral +[x_now; u_now]'*H*[x_now; u_now];   %d_target(8,1);   r(x_prev,u_prev)
                
                
                zbar(:,1)=zbar(:,2);
                zbar(:,2)=zbar(:,3);
                zbar(:,3)=zbar(:,4);
                zbar(:,4)=zbar(:,5);
                zbar(:,5)=zbar(:,6);
                zbar(:,6)=zbar(:,7);
                zbar(:,7)=[x_prev(1)^2;x_prev(1)*x_prev(2);x_prev(1)*u_prev;x_prev(2)^2;x_prev(2)*u_prev;u_prev^2];%phi(x_prev,u_prev)  %  zbar(:,8);
                
                
                if mod(count,7)==0
                    
                    zbar2=zbar*zbar';
                    q=zbar*d_target;
                    w(:,t+1,i)=inv(zbar2)*q;
                    if abs(w(:,t+1,i)-w(:,t,i))< epsilon*ones(tau,1)
                        H=[w(1,t+1,i) w(2,t+1,i)/2 w(3,t+1,i)/2; w(2,t+1,i)/2 w(4,t+1,i) w(5,t+1,i)/2;w(3,t+1,i)/2 w(5,t+1,i)/2 w(6,t+1,i)];
                        Hyy=H(3,3);Hyx=H(3,1:2);
                        h_coef(:,:,i)=-inv(Hyy)*Hyx;
                        k(i)=k(i)+1;
                    else
                    end
                    
                    
                end
                
            end
            count=count+1;
        end
        
    end
end


plot(dt*((1:T+1)-1),x(:,1:T+1,1))
xlabel('Time(s)')
ylabel('System States')



