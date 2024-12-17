%This script generates Figs. 3-14 from
%Random Evolutionary Dynamics in Predator-Prey Systems Yields Large, Clustered Ecosystems
% by Hamster, Schaap, Van Heijster and Dijksman.
%Parameters are as described in Table 1 and the equations as outlined in
%section 2. To computes averages such as Fig.4, this script must be run
%multiple times. For Fig. 6, the script must be run multiple times for
%different eta values. For the figures with variable r_i, the noise in line
%80 must be uncommented. 
%
%Comments and questions can be sent to c.h.s.hamster@uva.nl

% time settings of simulation
Tend = 3e6;                 % Length of simulation

% New species spawning parameters
exist_thresh = 0.001;            % threshold below which species is considered extinct
eta = 0.02;                      % evolutionary noise parameter when creating new species
frac = 0.05;                     % fraction of species biomass that is transferred to daughter species

NspecP = 1;   % Number of starting prey species - default 1
NspecR = 1;   % Number of starting predator species - default 1
K=10;         %carrying capacity
expecnumspec = 8000; % max number of species, i.e. 4000 prey and 4000 predator species are allowed

%parameters
r0=1.1; 
r = zeros(expecnumspec,1); r(1:NspecP) = r0;
delta=0.4; 
Delta=zeros(expecnumspec,1); Delta(expecnumspec/2+1:end)=delta;
beta=1/4;

% Preallocate arrays
Tvector=zeros(expecnumspec,1);  %vector that will store all speciation times
TvectorP=zeros(expecnumspec/2,1);  %vector that will store all speciation times for prey
TvectorR=zeros(expecnumspec/2,1);  %vector that will store all speciation times forpredators

%if species j is alive, we set Alive(j)=j, otherwise 0.
Alive=zeros(1,expecnumspec); Alive(1:NspecP)=1; 
Alive(expecnumspec/2+1:expecnumspec/2+NspecR)=expecnumspec/2+1; 

%vector that stores for each species the index of the ancestor
Ancestor=zeros(expecnumspec,1);
Ancestor(1)=1; Ancestor(expecnumspec/2+1)=expecnumspec/2+1;

%generate interaction matrix A
A=zeros(expecnumspec); A(1:NspecP,expecnumspec/2+1:expecnumspec/2+NspecR)=-0.4; %embed intial matrix in larger matrix
A(expecnumspec/2+1:expecnumspec/2+NspecR,1:NspecP)=transpose(A(1:NspecP,expecnumspec/2+1:expecnumspec/2+NspecR))*-beta;
A=sparse(A);

xold=zeros(1,expecnumspec); %xold is the startingpoint for ode45
xold(1:NspecP)=4; xold(expecnumspec/2+1:expecnumspec/2+NspecR)=2; %set IC
x=zeros(expecnumspec,expecnumspec); %preallocating x that stores the solution at speciation events;
x(1,:)=xold; %store IC

B0=sum(xold); %initial biomass

NtotalP=NspecP;  %total number of prey species, start with NspecP
NtotalR=NspecR;  %total number of predator species, start with NspecR

pm=0.001; %set exponential parameter of speciation probability
told=0;   %start time
tnew=exprnd(1/pm); %time until new speciation event
Tvector(2)=tnew;
i=2; %set counter
options = odeset('RelTol',1e-6,'AbsTol',1e-6); %set options for ode45

while told<Tend && max(NtotalP,NtotalR)<expecnumspec/2
    [~,X]=ode45(@(t,x) LV(t,x,r,K,A,Delta),[told,tnew],xold,options);  %ode45 until next event
    xnew=X(end,:);                   %endpoint of simulation
    x(i,:)=xnew;                     %store endpoint of simulation
    xnew=xnew.*(xnew>exist_thresh);  %set small species to 0
    Alive=Alive.*(xnew>0);           %set extinct species in Alive to 0
    Al=nonzeros(Alive)';             %get vector of alive species
    A(:,Alive==0)=0;                 %delete old species from matrix to keep it sparse
    A(Alive==0,:)=0;                 %delete old species from matrix to keep it sparse
    j=randsample(Al,1,true,xnew(Al));%choose random ancestor
    if j<expecnumspec/2+1  %update A for new species if ancestor is prey
        xnew(NtotalP+1)=0.05*xnew(j); xnew(j)=0.95*xnew(j);
        A(NtotalP+1,nonzeros(Alive(expecnumspec/2+1:end)))=min(A(j,nonzeros(Alive(expecnumspec/2+1:end)))+eta*randn(1,length(nonzeros(Alive(expecnumspec/2+1:end)))),0);
        A(expecnumspec/2+1:end,1:expecnumspec/2)=transpose(A(1:expecnumspec/2,expecnumspec/2+1:end))*-beta;
        r(NtotalP+1)=r(j);%+eta*randn(1); 
        NtotalP=NtotalP+1;
        TvectorP(NtotalP)=tnew;
        Ancestor(NtotalP)=j;
        Alive(NtotalP)=NtotalP;
    else  %update A for new species if ancestor is predator
        xnew(expecnumspec/2+NtotalR+1)=0.05*xnew(j); xnew(j)=0.95*xnew(j);
        A(expecnumspec/2+NtotalR+1,nonzeros(Alive(1:expecnumspec/2)))=max(A(j,nonzeros(Alive(1:expecnumspec/2)))+beta*eta*randn(1,length(nonzeros(Alive(1:expecnumspec/2)))),0);
        A(1:expecnumspec/2,expecnumspec/2+1:end)=transpose(A(expecnumspec/2+1:end,1:expecnumspec/2))/(-beta);
        Alive(expecnumspec/2+NtotalR+1)=expecnumspec/2+NtotalR+1;
        NtotalR=NtotalR+1;
        TvectorR(NtotalR)=tnew;
        Ancestor(expecnumspec/2+NtotalR)=j;
    end
    i=i+1
    pmb=pm*sum(xnew)/B0; %create new exponential parameter
    told=tnew;
    tnew=told+exprnd(1/pmb); %create new random time
    Tvector(i)=tnew; %concatenate time
    xold=xnew;
end

S=uint16((x(1:i-1,:)>0));  %Total number of species, as integer
S=uint16(sum(S,2));
Sp=uint16((x(1:i-1,1:expecnumspec/2)>0));  %number of prey species, as integer
Sp=uint16(sum(Sp,2));
Sr=uint16((x(1:i-1,expecnumspec/2+1:end)>0));  %number of predator species, as integer
Sr=uint16(sum(Sr,2));

B=sum(x(1:i-1,:),2); %Total Biomass
Bp=sum(x(1:i-1,1:expecnumspec/2),2); %Biomass Prey
Br=sum(x(1:i-1,expecnumspec/2+1:end),2); %Biomass Preydator

figure(1)
hold on
plot(Tvector(1:i-1),Sr,'linewidth',3,'color','r')
plot(Tvector(1:i-1),Sp,'linewidth',3,'color','b')
ylabel('$S$','Interpreter','latex','Fontsize',45)
xlabel('$t$','Interpreter','latex','Fontsize',45)
xlim([0 Tvector(i-1)])
hold off

figure(2)
hold on
plot(Tvector(1:i-1),Br,'linewidth',3,'color','r')
plot(Tvector(1:i-1),Bp,'linewidth',3,'color','b')
ylabel('$B$','Interpreter','latex','Fontsize',45)
xlabel('$t$','Interpreter','latex','Fontsize',45)
xlim([0 Tvector(i-1)])
hold off


%compute final run to equilibrium for Fig. 4a
Aend=A;
Aend(:,Alive==0)=[];    %delete old species from matrix to keep it sparse
Aend(Alive==0,:)=[];    %delete old species from matrix to keep it sparse
xend=nonzeros(xnew);
rend=r;
rend(Alive==0)=[];
Deltaend=Delta;
Deltaend(Alive==0)=[];

%run final simulation 
[t,X]=ode45(@(t,x) LV(t,x,rend,K,Aend,Deltaend),[0,10^5],xold(Alive>0),options); 
figure(3)
hold on
plot(t,X,'linewidth',2)
xlabel('$t$','Interpreter','latex','Fontsize',45)
ylabel('Biomass','Interpreter','latex','Fontsize',45)
hold off

%set new initial condition for Perturbation experiment in Fig 4b
xnew=X(end,:);
xnew=xnew.*(xnew>exist_thresh); %set small species to 0
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[tpert,Xpert]=ode45(@(t,x) LV(t,x,rend,K,Aend,Deltaend),[0,2000],2*xnew,options);  
figure(4)
hold on
plot(tpert,Xpert,'linewidth',2)
ylabel('Biomass','Interpreter','latex','Fontsize',45)
xlabel('$t$','Interpreter','latex','fontsize',45)
hold off


Clusters=GCModulMax1(Aend); %this computes the clustering for a modularity algorithm as in Fig. 9
RemSpec=length(Aend); %Number of remaining species
NumbAl=1:RemSpec;
Sort=cat(1,NumbAl(Clusters==1)',NumbAl(Clusters==2)',NumbAl(Clusters==3)');
Asortend=Aend(Sort,Sort);

figure(5)
mycolormap=[parula];
mycolormap(52,:)=[0.6350 0.0780 0.1840]; %52 is here the index of zero, might change from iteration to interation
hold on
colormap(mycolormap);
colorbar
subplot(1,2,1)
pcolor(-flip(Aend,1)); shading flat;
xlim([1 RemSpec]); ylim([1 RemSpec])
axis off
subplot(1,2,2)
pcolor(-flip(Asortend,1)); shading flat;
xlim([1 RemSpec]); ylim([1 RemSpec])
axis off
hold off

%make graph of ancestoral distances for prey species
adjP=zeros(NtotalP);
for k=1:NtotalP
    adjP(k,Ancestor(k))=abs(TvectorP(k)-TvectorP(Ancestor(k))); adjP(Ancestor(k),k)=adjP(k,Ancestor(k));
end
adjP=sparse(adjP);
GP=graph(adjP);

NtotalPA=length(nonzeros(rend)); %number of remaining prey species
DistPend=zeros(NtotalPA);       %matrix that stores genealogical distances between remaining prey species
IndexAlP=nonzeros(Alive(1:expecnumspec/2));
for i=1:NtotalPA
    for j=i+1:NtotalPA
        ShortDist=cell2mat(allpaths(GP,IndexAlP(i),IndexAlP(j)));
        [~,I]=min(ShortDist);
        if I==1
            DistPend(i,j)=tnew-TvectorP(ShortDist(2));
            DistPend(j,i)=DistPend(i,j);
        else
            DistPend(i,j)=tnew-min(TvectorP(ShortDist(I-1)),TvectorP(ShortDist(I+1)));
            DistPend(j,i)=DistPend(i,j);
        end
    end
end

%make graph of ancestoral distances for predartor species
adjR=zeros(NtotalR);
for k=1:NtotalR
    ancR=Ancestor(expecnumspec/2+k)-expecnumspec/2;
    adjR(k,ancR)=abs(TvectorR(k)-TvectorR(ancR));
    adjR(ancR,k)=adjR(k,ancR);
end
adjR=sparse(adjR);
GR=graph(adjR);

IndexAlR=nonzeros(Alive(expecnumspec/2+1:end));  
NtotalRA=length(IndexAlR); %number of remaining predator species
DistRend=zeros(NtotalRA);  %matrix that stores genealogical distances between remaining predator species

for i=1:NtotalRA
    for j=i+1:NtotalRA
        ShortDist=cell2mat(allpaths(GR,IndexAlR(i)-expecnumspec/2,IndexAlR(j)-expecnumspec/2));
        [~,I]=min(ShortDist);
        if I==1
            DistRend(i,j)=tnew-TvectorR(ShortDist(2));
            DistRend(j,i)=DistRend(i,j);
        else
            DistRend(i,j)=tnew-min(TvectorR(ShortDist(I-1)),TvectorR(ShortDist(I+1)));
            DistRend(j,i)=DistRend(i,j);
        end
    end
end

%turn distance matrix into an linkage object for matlab's dendogram 
ZP=linkage(squareform(DistPend));
figure(6)
%H1=dendrogram(ZP,0,'ColorThreshold',1e6);
H1=dendrogram(ZP,0);
set(H1,'LineWidth',3);
set(gca,'XTick',[])
ylabel('Distance to Ancestor','Interpreter','Latex','Fontsize',45)
xlabel('Prey','Interpreter','Latex','Fontsize',45)

ZR = linkage(squareform(DistRend));
figure(7)
H2=dendrogram(ZR,0,'ColorThreshold',1e6);
set(H2,'LineWidth',3);
set(gca,'XTick',[])
ylabel('Distance to Ancestor','Interpreter','Latex','Fontsize',45)
xlabel('Predator','Interpreter','Latex','Fontsize',45)

P=Aend(1:NtotalPA,NtotalPA+1:end);
%kmeans clustering on functional differences Prey and evaluation using
%silhoutette. 5 is the number of cluster it is forced to find
ClusterDistPphen=kmeans(full(P),5,'Distance','sqeuclidean','Replicates',500);
[silh] = silhouette(full(P),ClusterDistPphen,'sqeuclidean');
mean(silh) %the higher the value, the better

%kmeans clustering on functional differences Predator and evaluation using
%silhoutette. 5 is the number of cluster it is forced to find
ClusterDistRphen=kmeans(beta*full(P'),5,'Distance','sqeuclidean','Replicates',500);
[silh] = silhouette(beta*full(P'),ClusterDistRphen,'sqeuclidean');
mean(silh)


