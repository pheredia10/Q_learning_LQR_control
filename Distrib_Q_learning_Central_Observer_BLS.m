close all
clc
clear
% Initial Structures
T=200; % number of time steps
dt=1e-3; % time step size
time=T*dt; %elapsed time
t_vector=0:dt:time;
m=3; % number of agents
n=2; % number of states
n_bar=1; % number of control commands
xx=n*(n+1)/2; %combination of states with replacement
uu=n_bar*(n_bar+1)/2;  % n_bar*(n_bar+1)/2 is combinations of control commands with replacement
tau=(n+n_bar)*(n+n_bar+1)/2; % number of weights
w=zeros(tau,T,m); % weights , Note: weights that correspnd to x*x elements in phi(x) go with S_xx and so on
G_c=[1,1,0;1,1,1;0,1,1];% communication matrix  doubly_stochastic(m);
x=zeros(n,T,m); %states
psi=zeros(tau,T+1,m); % combined weights
mu=.0000001*ones(1,T); % optimization algorithm step-size
epsilon=.01; % convergance tolerance  .0000000001;
h_coef=zeros(n_bar,n,m); % control policy coefficient

% Initial Conditions
x(:,1,1)=[-.1;.1];
x_hat(:,1,1)=[.01;0];
x_hat(:,1,2)=[0;.01];
x_hat(:,1,3)=[0;0];

k=ones(m,1);
zbar(1:tau,1:7,m)=0;
d_target(1:7,1,m)=0;
H=zeros(3,3,m);
H(:,:,1)=[1 0 0;0 1 0;1 0 1];
H(:,:,2)=[1 0 0;0 1 0;1 0 1];
H(:,:,3)=[1 0 0;0 1 0;1 0 1];
Hyy=H(3,3,1);Hyx=H(3,1:2,1);
h_coef(:,:,1)=-inv(Hyy)*Hyx;
h_coef(:,:,2)=-inv(Hyy)*Hyx;
h_coef(:,:,3)=-inv(Hyy)*Hyx;
L=-inv(Hyy)*Hyx;
R=1;P=eye(2,2);G=[P [0 ; 0];[0 0] R];

% System Description
A=[1 1;0.16 1]; % Unstable
B=[0;1];
C=[1,0];
D=0;

% Observer Design
poles=[.1,.2];
L_o=-place(A',C',poles)';
wnew(:,1,:)=zeros(tau,1,m);

for t=1:T
    
    if t<50
        a1=0;
        a2=0.97;
    else
        a1=0;
        a2=0;
    end
    noise=a2*(0.5*sin(2.0*t)^2*cos(10.1*t)+0.9*sin(1.102*t)^2*cos(4.001*t)+0.3*sin(1.99*t)^2*cos(7*t)+0.3*sin(10.0*t)^3+0.7*sin(3.0*t)^2*cos(4.0*t)+0.3*sin(3.00*t)*1*cos(1.2*t)^2+0.400*sin(1.12*t)^2+0.5*cos(2.4*t)*sin(8*t)^2+0.3*sin(1.000*t)^1*cos(0.799999*t)^2+0.3*sin(4*t)^3+0.4*cos(2*t)*1*sin(5*t)^4+0.3*sin(10.00*t)^3);
    
    h_coef_bar=zeros(n_bar,n,1);
    
    u_prev=zeros(n_bar,1);
    % u_bar Calculation
    for ind=1:m
        x_prev=x_hat(:,t,ind);
        %h_coef_bar=((h_coef(:,:,ind)./m)+h_coef_bar);
        u_prev=((h_coef(:,:,ind)*x_prev +noise)./m +u_prev);
    end
    
    % Discrete Simulation
    x(:,t+1,1)=A*x(:,t,1)+B*u_prev;
    y(t,1)=C*x(:,t,1)+D*u_prev;
    
    for i =1:m
        
        x_prev=x_hat(:,t,i);
        %Central Observer
        y_hat(t,i)=C*x_hat(:,t,i);
        x_hat(:,t+1,i)=A*x_hat(:,t,i)+ B*u_prev + L_o*(y_hat(t,i)-y(t,1));
        x_tilde(:,t+1)=(x_hat(:,t+1,1)-x(:,t+1,1));
        % u_bar_now Calculation
        u_now=zeros(n_bar,1);
        for ind=1:m
            x_now=x_hat(:,t+1,ind);
            %h_coef_bar=((h_coef(:,:,ind)./m)+h_coef_bar);
            u_now=((h_coef(:,:,ind)*x_now)./m +u_now);
        end
        x_now=x_hat(:,t+1,i); % updated states
        
        if t==1
            d0= [x_prev; u_prev]'*G*[x_prev; u_prev]+[x_now; u_now]'*H(:,:,i)*[x_now; u_now];   %d_target(8,1);   r(x_prev,u_prev)
            Jnew(i)=(d0-[x_prev; u_prev]'*H(:,:,i)*[x_prev; u_prev])^2;
        end
        
        if t<50
            
            d_target(1,1,i)=d_target(2,1,i);
            d_target(2,1,i)=d_target(3,1,i);
            d_target(3,1,i)=d_target(4,1,i);
            d_target(4,1,i)=d_target(5,1,i);
            d_target(5,1,i)=d_target(6,1,i);
            d_target(6,1,i)=d_target(7,1,i);
            d_target(7,1,i)= [x_prev; u_prev]'*G*[x_prev; u_prev]+[x_now; u_now]'*H(:,:,i)*[x_now; u_now];   %d_target(8,1);   r(x_prev,u_prev)
            
            
            zbar(:,1,i)=zbar(:,2,i);
            zbar(:,2,i)=zbar(:,3,i);
            zbar(:,3,i)=zbar(:,4,i);
            zbar(:,4,i)=zbar(:,5,i);
            zbar(:,5,i)=zbar(:,6,i);
            zbar(:,6,i)=zbar(:,7,i);
            zbar(:,7,i)=[x_prev(1)^2;x_prev(1)*x_prev(2);x_prev(1)*u_prev;x_prev(2)^2;x_prev(2)*u_prev;u_prev^2];%phi(x_prev,u_prev)  %  zbar(:,8);
            
            
            
            if mod(t,7)==0
                zbar2=zbar(:,:,i)*zbar(:,:,i)';
                q=zbar(:,:,i)*d_target(:,:,i);
                wnew(:,k(i)+1,i)=inv(zbar2)*q;
                w(:,k(i)+1,i)=inv(zbar2)*q;
                H_i(:,:)=[w(1,k(i)+1,i) w(2,k(i)+1,i)/2 w(3,k(i)+1,i)/2; w(2,k(i)+1,i)/2 w(4,k(i)+1,i) w(5,k(i)+1,i)/2;w(3,k(i)+1,i)/2 w(5,k(i)+1,i)/2 w(6,k(i)+1,i)];
                Jnew(i)=(d_target(7,1,i)-[x_prev; u_prev]'*H_i(:,:)*[x_prev; u_prev])^2;
                k(i)=k(i)+1;
                
                
            end
            
        end
        
        J(t,i)=Jnew(i);
        wt(:,t,i)=wnew(:,k(i),i);
    end
    
    if mod(t,7)==0
        for i=1:m
            w_bar=zeros(tau,T,1);
            for j=1:m
                w_bar=(G_c(i,j)*(w(:,k(j),j))+w_bar);
            end
            H(:,:,i)=[w_bar(1,k(1)) w_bar(2,k(1))/2 w_bar(3,k(1))/2; w_bar(2,k(1))/2 w_bar(4,k(1)) w_bar(5,k(1))/2;w_bar(3,k(1))/2 w(5,k(1))/2 w_bar(6,k(1))];
            Hyy=H(3,3,i);Hyx=H(3,1:2,i);
            h_coef(:,:,i)=-inv(Hyy)*Hyx;
        end
    end
    hvec_coef1(:,t)=h_coef(:,:,1);
    hvec_coef2(:,t)=h_coef(:,:,2);
    h_vec(:,t)=h_coef_bar;
    
end


plot((1:T+1)-1,x(:,1:T+1,1))
hold
plot((1:T+1)-1,x(:,1:T+1,2))
xlabel('Time Step')
ylabel('System States')
figure
loglog(1:T,J)
xlabel('Time Step')
ylabel('J(w)')


