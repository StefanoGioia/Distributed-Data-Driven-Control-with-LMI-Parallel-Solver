%STEFANO GIOIA, DISTRIBUTED DATA-DRIVEN CONTROL AND OPTIMIZATION
%DISTRIBUTED IDENTIFICATION - INCREASED T RUNS

yalmip('clear')
clear all
close all

tic

%% Matrices to identify

%Agents dynamics in DT 
A_agent = [1.0008 1e-3; 3e-4 1.0001];
Acoupl = [1e-4 1.5e-4; 2e-5 1e-4];
B_agent = [5e-4 1e-3; 8e-4 1.1e-3]; %2 inputs and 2 states
% 
% Bcoupl

N_agents= 3;

%Block diagonal structure
IdAg= eye(N_agents);
%Adjaceny matrix - states (state coupling)
AD_A = [0 1 0 ; 0 0 1; 1 0 0]; %graph

%Full A
% kron (structure, matrix at)
A= kron(IdAg, A_agent) + kron(AD_A, Acoupl) ;

%Full B, each agent depends on own input(that can depend on coupled states)
B = kron(IdAg,B_agent) ;%+ kron(Coupling, Bcoupl) ;%

n_ag=2;
m_ag=2;
n= N_agents*n_ag;
m= N_agents*m_ag;

%% Generate collected data PE
iterTOT=15; 
incrT=50;
%Margin to satisfy T inequality
marg = 2 +incrT*(iterTOT-1);

T= (m+1)*n + m+ marg;


    
% Generate input for the global system

U0Tfix= rand(m,T);  % m: #inputs of the global system 
% U0Tfix= (-0.5+rand(m,T)).*10.^(-5+ 10*rand(m,T));  %wide


%% Build Hankel

order= n+1; % n: #states of the global system

for iter=0:iterTOT-1
    
    HankelU =[];
    consideredT= 50 + incrT*iter;
    U0T = U0Tfix(:,1:consideredT);
    

%Timestep indices
first_col_Hankel =1:1:order;
last_row_Hankel = order:1:consideredT;
hankel_idx = hankel(first_col_Hankel,last_row_Hankel);

%Place input samples
%Fill row-wise
for i = 1:order
    HankelRow_idx = hankel_idx(i,:);
    HankelU((i-1)*m+1 : i*m, :) = U0T(:,HankelRow_idx); %TO CLEAR IF TEST DIFFERENT T
end

%PE Test
if rank(HankelU)== m*(n+1) 
    disp('Persistently exciting input')
else
    error('Input not P.E.')
end

end

%% Generate collected data

x0 = (10.^(-0.5 + rand(n,1))).*(-0.5 + rand(n,1)) ;%10*(-0.5 + rand(n,1)) ;

%Computation of X0T and X1T
Xd(:,1) = x0;

for i = 1:T
    
    Xd(:,i+1) = A*Xd(:,i)+ B*U0Tfix(:,i);
    
end

X0T = Xd(:,1:T);
X1T = Xd(:,2:T+1);

X1Tfix= X1T;
X0Tfix= X0T;

%% Distributed identification

% %Centralized Identification with pseudoinverse
BA_ID_inv = X1T*pinv([U0T;X0T]);
diff_inv= BA_ID_inv - [B A];

for iter=1:iterTOT

yalmip('clear')

consideredT = 50 + (iter-1)*incrT;

X1T= X1Tfix(:,1:consideredT);
X0T= X0Tfix(:,1:consideredT);
U0T= U0Tfix(:,1:consideredT);

%Centralized identification with minimization
%selection: structure of BA, with inf in place of ones
selection=double([B A]~=0)./double([B A]==0);

BA = sdpvar(n,n+m);
cost = sum(sum( (X1T-BA*([U0T;X0T])).^2 ))
Constraints = []; % [-selection <= BA <= selection]; 
sol = optimize(Constraints, cost);
BA_ND = value(BA);    

maxabsdiff_ND = max(max(abs(BA_ND-[B A])));

%modified for zeros
BA_ND_noZeros = BA_ND.*double([B A]~=0) + double([B A]==0);
BA_real_noZeros = [B A].*double([B A]~=0) + double([B A]==0);

maxErrRel_ND(iter)= max(max(abs((BA_ND_noZeros-BA_real_noZeros)./BA_real_noZeros)));


%Adiacency matrices if B A not known a priori
self_dep = eye(N_agents);

A_dep    = AD_A + self_dep;

AD_B      = zeros(N_agents);
B_dep    = AD_B + self_dep;
 
%States and inputs interdependencies
for j = 1: N_agents
    states_dep (:,n_ag*(j-1)+1:n_ag*j) = repmat(A_dep(:,j),1,2);
    inputs_dep (:,m_ag*(j-1)+1:m_ag*j) = repmat(B_dep(:,j),1,2);
end


clear BA cost Constraints sol nonzeromin
yalmip('clear')
close all

BAsave = [];
for j=1:N_agents % = number of agents

    %Data needed by agent j, where n_ag: #states of each agent
    agent_state_rows = (j-1)*n_ag+1 : 1 : j*n_ag; 
    X1T_agent = X1T(agent_state_rows,:);
    
    dep_states_rows = find(states_dep(j,:));
    dep_input_rows  = find(inputs_dep(j,:));
    
    X0T_Dagent = X0T(dep_states_rows,:);
    U0T_Dagent = U0T(dep_input_rows,:);
    
    %YALMIP toolbox settings and solve
    
    %Number of states and inputs on which agent j depends on
    N_dep = sum(states_dep(j,:))+sum(inputs_dep(j,:));
    %Variable declaration
    BAagent= sdpvar(n_ag, N_dep);
    
    %Frobenius norm to be minimized (5.1)
    cost = sum(sum((X1T_agent-BAagent*([U0T_Dagent;X0T_Dagent])).^2 ))
    
    %Solve
    sol(j,:)=optimize([],cost);

    %Save result
    BA.(["agent"]+num2str(j)) = value(BAagent); 
    
    BAsave(agent_state_rows,:)= value(BAagent);    
    
    
    
end



%Put identified matrices on orde
BA_ID= zeros(n,n+m);
for j=1:N_agents
    
    dep_states_rows = find(states_dep(j,:))+m;
    dep_input_rows  = find(inputs_dep(j,:));

    BA_ID([1:1:n_ag]+n_ag*(j-1),[dep_input_rows,dep_states_rows]) = BAsave([1:1:n_ag]+n_ag*(j-1),:);
    
end

maxabsdiff_ID = max(max(abs(BA_ID-[B A])));

%modified for zeros
BA_ID_noZeros = BA_ID.*double([B A]~=0) + double([B A]==0);

maxErrRel_ID(iter)= max(max(abs((BA_ID_noZeros-BA_real_noZeros)./BA_real_noZeros)));
ccc(iter)=consideredT;
maxrapp(iter)=maxErrRel_ND(iter)/maxErrRel_ID(iter); %max?


end

toc


figure()
plot(50:incrT:(iterTOT-1)*incrT+50,maxErrRel_ID)
hold on
xlabel('T')
hold on
ylabel('maxErrD')
hold off

figure()
plot(50:incrT:(iterTOT-1)*incrT+50,log10(maxErrRel_ID))
hold on
xlabel('T')
hold on
ylabel('log_{10} maxErrD')
hold off


