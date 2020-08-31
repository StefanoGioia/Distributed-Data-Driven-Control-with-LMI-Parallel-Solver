%STEFANO GIOIA, DISTRIBUTED DATA-DRIVEN CONTROL AND OPTIMIZATION
%DISTRIBUTED STABILIZATION

clear all
close all

dt = 1e-3;

%  cvx_precision default

%% Dynamics
%Agents dynamics in DT 
A_agent = [1.0008 1e-3; 3e-4 1.0001];
Acoupl = [1e-4 1.5e-4; 2e-5 1e-4];
B_agent = [5e-4 1e-3; 8e-4 1.1e-3]; %2 inputs and 2 states
% 
% Bcoupl

N_agents= 3;%9;

%Block diagonal structure
IdAg= eye(N_agents);
%Adjaceny matrix - states (state coupling)
%  load('AD_A.mat'); %9 agents
AD_A = [0 1 0 ; 0 0 1; 1 0 0]; %graph

%%%%%%%%%%%%%%%%%%%%%%%%%

% h= WattsStrogatz(9,2,0.2)
% plot(h,'NodeColor','k','Layout','circle');
% title('Watts-Strogatz Graph with $N = 9$ nodes, $K = 2$, and $\beta = 0.2$', ...
%     'Interpreter','latex')
% load('edges.mat');
% edges=table2array(edges);
% %build adjacency
% AD_A = zeros(N_agents);
% %convert to linear idx
% linidxUP = N_agents*(edges(:,2)-1)+edges(:,1);
% linidxLOW = N_agents*(edges(:,1)-1)+edges(:,2);
% AD_A([linidxUP;linidxLOW])=1;


%%%%%%%%%%%%%%%%%%%%%%%%%


%Full A
% kron (structure, matrix agent)
A= kron(IdAg, A_agent) + kron(AD_A, Acoupl) ;

%Full B, each agent depends on own input(that can depend on coupled states)
B = kron(IdAg,B_agent) ;%+ kron(Coupling, Bcoupl) ;%

n_ag=2;
m_ag=2;
n= N_agents*n_ag;
m= N_agents*m_ag;

%% F function set up 

max_abs_eig_A = max(abs(eig(A))) ;

desired_LHS = 0.5; % Value for the LHS of (5.9) 
alpha = desired_LHS/max_abs_eig_A ; 

A1 = alpha*A; 

A2 = (1-alpha)*A;    

%Y initialization

%Y size = n x m (n,m are the global system's states and inputs)
rowY = m_ag*N_agents ; %m_ag: each agent inputs, N_agents: number of agents
colY = n_ag*N_agents ; %n_ag: each agent states
Y = -0.5+rand(rowY,colY);  

% %Gain simple structure(no n_ag,m_ag)
% K_struct_simple = AD_A + IdAg; 
% Y = kron(K_struct, - 0.5+rand(m_ag,n_ag));

%Load a saved initial Y
% initialY=Y;
% load('initialY3ag.mat'); %3ag or 9ag
% Y=initialY;

%For P and inv(P) computation 
%P size = n x n (n:global system's states)
dimP = n_ag*N_agents ; %n_ag: each agent states, N_agents: number of agents

%CVX can't solve strict inequalities
tol = 1e-8;  %tolY?
I_LMI_F = eye(2*dimP); % >= tol*I_LMI_F instead of > 0
I_LMI_P = eye(dimP);

I_P = eye(dimP);  %For the inequalities, dimP = length of P

P= I_P; %P initialization

DPinv_each= []; %Variable collecting every agent variation

%% Number of iterations and other(gamma) parameters %%
%Homotopy loops
M=1

cSAME=1.2e-8;
DSAME=2;

D_P = DSAME; 
c_P= cSAME; 

c_Y=  cSAME;
D_Y = DSAME;

c_inv = 0.8;
D_inv =60; %100

%% Other counters and optionz

% put 1 To test with "centralized" (no distributed alg)
PNotDistr=0;
YNotDistr=0;
InvNotDistr=0;

No_Previous_Pinv =1;

%General scheme counter
i= 1;

%global counter for k storage
l=0;

% counters for when no agent finds solution ("ag"[] for inverse)
totbadY=0;
totbadP=0;

%Agents results %
agCVXres = '';
agcomb   ='';
%counter for ag variables, global
invtest=1;


%counter for storage in while loop %
cc=1;


%symmetric assignments of P for odd N agents 
symmEqualPassign= diag(1:1:N_agents);
symmEqualPassign=symmEqualPassign+ tril(ones(N_agents),-1);
[freerow,freecol]= find(symmEqualPassign==0);
freeidx = [freerow,freecol];

ASSIGMORE = (N_agents-1)/2;
Assignment=0;
for aa=1:N_agents
    
    for as=1:ASSIGMORE
           
      Assignment=Assignment+1;
      %linear indices in columns
        symmEqualPassign((freeidx(Assignment,2)-1)*N_agents+freeidx(Assignment,1))= aa;
        symmEqualPassign((freeidx(Assignment,1)-1)*N_agents+freeidx(Assignment,2))= aa;
    end
 
end

%

assignP= kron(symmEqualPassign,ones(n_ag));

assignPinv_simple= zeros(N_agents);
assignY_simple=zeros(N_agents);
for ii=1:N_agents

   assignPinv_simple(:,ii) = ii;
   assignY_simple(ii,:) = ii;

end
assignPinv =  kron(assignPinv_simple,ones(n_ag)); %column wise
assignY = kron(assignY_simple,ones(m_ag,n_ag)); %row wise

%graph
K_struct = kron(AD_A + IdAg, ones(m_ag,n_ag));

%% Computation

%Store largest eigenvalues in magnitude, P, u, Y
eigHomotopystory=zeros(n_ag*N_agents,M+1);
eigHomotopystoryu=zeros(n_ag*N_agents,M+1);
ustory=zeros(1,M+1);
Pstory=[];
Ystory=[];


while i<= M
    
    u=(i-1)/M;
    ustory(i)=u;

disp('computation of P');

    if PNotDistr==1
        cvx_begin sdp
        variable P(dimP,dimP) symmetric %diagonal
        
        [           P                 ( A1*P  + u*(A2*P + B*Y) )';...
           (A1*P  + u*(A2*P + B*Y))               P          ] >= tol*I_LMI_F


        P >= tol*I_LMI_P 

        cvx_end

    else
       tic
     %Distributed P computation  %%%%%%%%%%%%%%%%%%%%%%
        for d=0:D_P

          g = c_P*(1 - d/D_P);  %gamma (5.10)
          DP=zeros(dimP);   %initialize computed global variation 
          PFIX=P;           %P from the previous iteration
          baditP = 0;
          
          for agent=1:N_agents
          
          %  AGENT P computation %
   
          indices=(assignP~=agent);
            cvx_begin sdp
                variable P(dimP,dimP) symmetric 
                [           P                 ( A1*P  + u*(A2*P + B*Y) )';...
                   (A1*P  + u*(A2*P + B*Y))               P          ] >= -g*I_LMI_F + tol*I_LMI_F

                P >= tol*I_LMI_P
                %Assignment : fix other agents entries
                
                P(indices) == PFIX(indices)
         
            cvx_end
            %Agent contribute to global P variation DP
            if sum(sum(isnan(P)))>0%cvx_status == "Infeasible"||cvx_status == "Failed"||cvx_status=="Inaccurate/Infeasible"
                DP = DP+0;
                baditP=baditP+1;
            else
                DP= DP+full(P)-PFIX;
            end
         end

            if baditP==N_agents
        %       pause(0.7);
        %       error('P');
                totbadP=totbadP+1;
               warning('No step for P');
            end

            %New P 
            P = PFIX + DP;  %full() ?
            Pstory(:,:,d+1)=P;
        end

    end
timeExecP=toc;
    
disp('computation of the inverse of P');

    if InvNotDistr==1
        Pinv = inv(P); %Is positive definite when P is
        % Pinv  =zeros(dimP);
        % for iii=1:length(P)
        %     Pinv(iii,iii) = 1/P(iii,iii);
        % end

    else
         % Distributed inverse computation %%%%%%%%%%
         
invstory=[];
tic
        for d=0:D_inv

            if No_Previous_Pinv ==1 % No previously computed inverse available
                Pinv = I_P; 
                No_Previous_Pinv = 0;
            end

            gamma=c_inv*(1-d/D_inv) ; 
            PinvFIX=Pinv;

            made_variations=0;  %Agents variations counter
            DPINV= zeros(dimP); %Will contain the total variation at the end of the loop

          for agent=1:N_agents           
            
            %AGENT 
             indices=(assignPinv~=agent);
            
            cvx_begin sdp
                variable DPinv(dimP,dimP) 

                P*(DPinv+PinvFIX) >=  (1-gamma)*I_P 
                P*(DPinv+PinvFIX) <=  (1+gamma)*I_P  

                DPinv(indices) == 0 %Assignment, here column wise

            cvx_end

            if sum(sum(isnan(DPinv)))>0%cvx_status == "Infeasible"||cvx_status == "Failed"||cvx_status=="Inaccurate/Infeasible" 
               DPinv_each(:,:,agent) = zeros(dimP);
            else
               DPinv_each(:,:,agent) = full(DPinv);
               made_variations = made_variations + 1;
            end
            agCVXres{invtest,agent}=cvx_status;
          end
    
            %Total variation made by agents   
            for agent = 1 : N_agents
               DPINV = DPinv_each(:,:,agent) + DPINV; 
            end

            %Average (4.36) , avoid 0/0 division
            if made_variations ~= 0 
               DPINV =DPINV/made_variations;
            end

            Pinv = PinvFIX + DPINV ;  

            %test of resulting inverse
             Pinvtest=Pinv;
     
                 cvx_begin sdp
                    variable Pinvt(dimP,dimP) 
                   
                    P*Pinvt >=  (1-gamma)*I_P  
                    P*Pinvt <=  (1+gamma)*I_P  
                    Pinvt==Pinvtest

                 cvx_end
               invstory(:,:,d+1)=Pinv;
                 agcomb{invtest}=cvx_status;
                 ACVAR(invtest) =made_variations;
                 invtest=invtest+1;
        end

    end
    timeExecInv=toc;
    eigHomotopystory(:,i)= eig(A+B*Y*Pinv);
    eigHomotopystoryu(:,i)= eig((A1  + u*(A2 + B*Y*Pinv)));
    
    
    i=i+1;
    u= (i-1)/M;
    ustory(i)=u;

disp('computation of Y');


 if YNotDistr==1
     
    cvx_begin sdp
        variable Y(rowY,colY)
      
        [           P                 ( A1*P  + u*(A2*P + B*Y) )';...
           (A1*P  + u*(A2*P + B*Y))               P          ] >= tol*I_LMI_F

       %3 agent system case %
        Y(1:2,:)  *Pinv(:,5:end) == 0
        Y(3:4,:)  *Pinv(:,1:2)   == 0 %%%%%%%%%%
        Y(5:end,:)*Pinv(:,3:4)   == 0


    cvx_end
    
 else
      % Distributed Y computation %%%%%%%%%%% 
 Ystory=[];
 tic
    for d=0:D_Y
      g = c_Y*(1 - d/D_Y);      %gamma (5.12)
      DY=zeros(rowY,colY);   %initialize computed collective variation of Y
      YFIX=Y;           %Y from the previous iteration

      %Store K
        l=l+1;
        store_K(:,:,l)= Y*Pinv;

        baditY=0;
        
      for agent=1:N_agents
      %  AGENT   Y computation %
              
              indices=(assignY~=agent);  
              rowsag= (1:1:n_ag)+(agent-1)*n_ag;

              [ridx,cidx]=find( K_struct(rowsag,:)==0);
              cidx = unique(cidx);

        cvx_begin sdp
            variable Y(rowY,colY)  
            [           P                 ( A1*P  + u*(A2*P + B*Y) )';...
               (A1*P  + u*(A2*P + B*Y))               P          ] >= -g*I_LMI_F + tol*I_LMI_F

            %Assignment : fix other agents rows 
            Y(indices)==YFIX(indices)

            %Constraint on agent 1 gain (same rows of Y)
            %K(1:2,5:end) == 0, agent 1 input doesn't depend on agent 3 states

            Y(rowsag,:)*Pinv(:,cidx) == 0 %Pinv is the inverse of P

        cvx_end
        %Agent contribute to the collective variation DY
        if sum(sum(isnan(Y)))>0%cvx_status == "Infeasible"||cvx_status == "Failed"||cvx_status=="Inaccurate/Infeasible"
            DY = DY+0;
            baditY=baditY+1;
        else
            DY= DY+full(Y)-YFIX;
        end

        %New Y 
        Y = YFIX + DY;
         Ystory(:,:,d+1) =Y;
        if baditY==N_agents

            totbadY=totbadY+1;
            %       pause(0.7);
            %       error('Y');
            warning('No step for Y');
        end


     end

    end
    timeExecY=toc;

 end


pp(:,:,cc)= full(P);
yy(:,:,cc)= full(Y);

cc=cc+1;
eigHomotopystory(:,i)= eig(A+B*Y*Pinv);
eigHomotopystoryu(:,i)= eig((A1  + u*(A2 + B*Y*Pinv)));
i=i+1;
      
    
end

ABsvaleigs= abs(eigHomotopystory); %sqrt(eigHomotopystory.*conj(eigHomotopystory)); %


end_=1

%% STEP4 Gain extraction

K = (full(Y)*Pinv).*K_struct;

%test
TESTDYN= (A+B*K)*P*(A+B*K)'-P;
TESTDYN=-TESTDYN ;%pos def test

try chol(TESTDYN)
    disp('Matrix TESTDYN is symmetric positive definite.')
catch ME
    disp('Matrix TESTDYN is not symmetric positive definite')
    error('NO')
end



