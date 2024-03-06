%This script generates Fig. 1a and Fig. 1b from
%Random Evolutionary Dynamics in Predator-Prey Systems Yields Large, Clustered Ecosystems
% by Hamster, Schaap, Van Heijster and Dijksman.
%Parameters are as described in Table 1 and the equations as outlined in
%section 2.
%
%%Comments and questions can be send to c.h.s.hamster@uva.nl


% time settings of simulation
Tend = 8e4;  % Length of simulation

% New species spawning parameters
exist_thresh = 0.001;            % threshold below which species is considered extinct
eta = 0.02;                      % evolutionary noise parameter when creating new species

NspecP = 1;   % Number of starting prey species - default 1
NspecR = 1;   % Number of starting predator species - default 1
K=10;         %carrying capacity  prey
expecnumspec = 200; % max number of species, i.e. 100 prey and 100 predator species are allowed

%parameters
r0=1.1;
r = zeros(expecnumspec,1); r(1:NspecP) = r0;
delta=0.4;
Delta=zeros(expecnumspec,1); Delta(expecnumspec/2+1:end)=delta;
beta=1/4;

Tvector=0;  %vector that will store all separate t values

%if species j is alive, we set Alive(j)=j, otherwise 0.
Alive=zeros(1,expecnumspec); Alive(1:NspecP)=1; 
Alive(expecnumspec/2+1:expecnumspec/2+NspecR)=expecnumspec/2+1; 

%generate interaction matrix A
A=zeros(expecnumspec);
A(1:NspecP,expecnumspec/2+1:expecnumspec/2+NspecR)=-0.4; %embed intial matrix in larger matrix
A(expecnumspec/2+1:expecnumspec/2+NspecR,1:NspecP)=transpose(A(1:NspecP,expecnumspec/2+1:expecnumspec/2+NspecR))*-beta;
A=sparse(A);

xold=zeros(1,expecnumspec); %xold is the startingpoint for ode45
xold(1:NspecP)=4; xold(expecnumspec/2+1:expecnumspec/2+NspecR)=2; %set IC
x=zeros(expecnumspec*50,expecnumspec); %preallocating vector x that stores all solutions of each speciation step;
x(1,:)=xold; %add IC to x

B0=sum(xold); %initial biomass

NtotalP=NspecP;  %total number of prey species, start with NspecP
NtotalR=NspecR;  %total number of predator species, start with NspecR

pm=0.001;  %set exponential parameter of speciation probability
told=0; %start time
tnew=exprnd(1/pm); %time until new speciation event
i=2; %set counter
options = odeset('RelTol',1e-6,'AbsTol',1e-6); %set ode45 options
while told<Tend && max(NtotalP,NtotalR)<expecnumspec/2
    [t,X]=ode45(@(t,x) LV(t,x,r,K,A,Delta),[told,tnew],xold,options);  %ode45 until next event
    tl=length(t);
    x(i:i+tl-2,:)=single(X(2:tl,:)); %concatenate solutions, and ignore the first value, bc its the last of the previous
    i=i+tl-1; %update counter
    Tvector=cat(1,Tvector,t(2:end)); %concatenate time
    xnew=X(end,:); %endpoint of simulation, i.e. initial condition of next step
    xnew=xnew.*(xnew>exist_thresh); %set small species to 0
    Alive=Alive.*(xnew>0);         %set extinct species in Alive to 0
    Al=nonzeros(Alive)';           %get vector of alive species
    A(:,Alive==0)=0;               %delete old species from matrix to keep it sparse
    A(Alive==0,:)=0;               %delete old species from matrix to keep it sparse
    j=randsample(Al,1,true,xnew(Al)); %choose random ancestor
    if j<expecnumspec/2+1 %update A for new species if ancestor is prey
        xnew(NtotalP+1)=0.05*xnew(j); xnew(j)=0.95*xnew(j);
        A(NtotalP+1,nonzeros(Alive(expecnumspec/2+1:end)))=min(A(j,nonzeros(Alive(expecnumspec/2+1:end)))+eta*randn(1,length(nonzeros(Alive(expecnumspec/2+1:end)))),0);
        A(expecnumspec/2+1:end,1:expecnumspec/2)=transpose(A(1:expecnumspec/2,expecnumspec/2+1:end))*-beta;
        r(NtotalP+1)=r(j);
        Alive(NtotalP+1)=NtotalP+1;
        NtotalP=NtotalP+1;
    else %update A for new species if ancestor is predator
        xnew(expecnumspec/2+NtotalR+1)=0.05*xnew(j); xnew(j)=0.95*xnew(j);
        A(expecnumspec/2+NtotalR+1,nonzeros(Alive(1:expecnumspec/2)))=max(A(j,nonzeros(Alive(1:expecnumspec/2)))+beta*eta*randn(1,length(nonzeros(Alive(1:expecnumspec/2)))),0);
        A(1:expecnumspec/2,expecnumspec/2+1:end)=transpose(A(expecnumspec/2+1:end,1:expecnumspec/2))/(-beta);
        r(NtotalR+1)=r(j);
        Alive(expecnumspec/2+NtotalR+1)=expecnumspec/2+NtotalR+1;
        NtotalR=NtotalR+1;
    end

    pmb=pm*sum(xnew)/sum(x(1,:)); %create new exponential parameter
    told=tnew;
    tnew=told+exprnd(1/pmb); %create new random time
    xold=xnew;
end

S=uint16((x(1:i-1,:)>0));  %number of species, as integer
S=uint16(sum(S,2));
Sp=uint16((x(1:i-1,1:expecnumspec/2)>0));  %number of prey species, as integer
Sp=uint16(sum(Sp,2));
Sr=uint16((x(1:i-1,expecnumspec/2+1:end)>0));  %number of predator species, as integer
Sr=uint16(sum(Sr,2));

B=sum(x(1:i-1,:),2); %Total Biomass
Bp=sum(x(1:i-1,1:expecnumspec/2),2); %Biomass prey
Br=sum(x(1:i-1,expecnumspec/2+1:end),2); %Biomass predator

figure(1)
hold on
xnan=x;
xnan(xnan==0)=nan; %zeros are set to nan so the inital zero part of each species does not show.
ylim([0 4.75])
plot(Tvector,xnan(1:i-1,:),'linewidth',3)
ylabel('$x_i,y_j$','Interpreter','latex','Fontsize',45)
xlabel('$t$','Interpreter','latex','Fontsize',45)
xlim([0 Tvector(i-1)])
hold off


figure(2)
hold on
xnan=x;
xnan(xnan==0)=nan;
ylim([0 4.75])
plot(Tvector,xnan(1:i-1,1:NtotalP),'linewidth',3,'color','b')
plot(Tvector,xnan(1:i-1,expecnumspec/2+1:expecnumspec/2+NtotalR),'linewidth',3,'color','r')
xlabel('$t$','Interpreter','latex','Fontsize',45)
ylabel('$x_i,y_j$','Interpreter','latex','Fontsize',45)
xlim([0 Tvector(i-1)])
hold off
