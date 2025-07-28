%
% Source code for "Pushing the Boundary of Quantum Advantage in Hard Combinatorial Optimization with Probabilistic Computer"
% by Shuvro Chowdhury et al.
%
% Date: July 2025
%
%
% This code performs Simulated Quantum Annealing on 3D spin glass instances. 

clc; clearvars; close all;

%*********************************** Inputs ******************************% 
% % loading single replica J, h
logical = 1;    % 1 for logical, 0 for embedded
instanceSize = 15;
instanceID = 55;
run_per_instance = 50;
total_real_num_sweeps = 2000; % how many total sweeps in MCMC chain
betaAll = 0.5;
base_seed  = 223;
GammaX = 3.0;
num_replicas = 2850;

if(logical)
	pathName = ['./logicalMatInstances/size' num2str(instanceSize,'%02d') '/'];
	colorMapConstant = 3;
else
	pathName = ['./matInstances/size' num2str(instanceSize,'%02d') '/'];
	colorMapConstant = 5;
end

% load saved J matrices
JFileName = ['JOrig_' num2str(instanceID,'%04d') '.mat'];
load(JFileName,'J');
cMapFileName = [pathName 'colorMap.csv'];

%******************************* Pre-processing **************************% 
addpath(pathName);

W_org = sparse(J); %original values from 1 replica
num_pbits_per_replica = length(J);

% now work on replicas
JT = sparse(eye(num_pbits_per_replica));
JJ1 = kron(sparse(eye(num_replicas)),W_org/num_replicas);
JJ2O = kron(diag(sparse(ones(1,num_replicas-1)),+1),JT)+kron(diag(sparse(ones(1,num_replicas-1)),-1),JT);

% PBC along replica direction
if(num_replicas>1)
	JJ2O(1:num_pbits_per_replica,(num_replicas-1)*num_pbits_per_replica+1:num_pbits_per_replica*num_replicas) = JT;
	JJ2O((num_replicas-1)*num_pbits_per_replica+1:num_pbits_per_replica*num_replicas,1:num_pbits_per_replica) = JT;
end

clear J0 JT J;

%%   Work on the biases now
% for the DWAVE problem, h is all 0
GammaZ = 0;
h_org = sparse(zeros(num_pbits_per_replica,1));
hhr = zeros(1,num_pbits_per_replica*num_replicas);
hir = GammaZ*(ones(1,num_pbits_per_replica))/num_replicas;
hir =repmat(hir,1,num_replicas);
h = (hhr+hir)';
clear hhr hir;

num_pbits = num_pbits_per_replica*num_replicas; % # of magnets
seq_time = total_real_num_sweeps;

% placeholders
m_last3 = zeros(num_pbits/num_replicas,length(seq_time),run_per_instance);
E_last3 = zeros(1,length(seq_time),run_per_instance);

BestState = zeros(length(h_org),1);
BestStateAll = zeros(length(h_org),run_per_instance);
BestEnAll = zeros(1,run_per_instance,length(seq_time));

colorMap_org = readmatrix(cMapFileName);

colorMap = zeros(1,num_pbits_per_replica*num_replicas);
colorMap(1:num_pbits_per_replica) = colorMap_org;
for ii = 2:2:num_replicas
colorMap((ii-1)*num_pbits_per_replica+1:ii*num_pbits_per_replica) = colorMapConstant - colorMap_org;
end
for ii = 3:2:num_replicas
colorMap((ii-1)*num_pbits_per_replica+1:ii*num_pbits_per_replica) = colorMap_org;
end

JJ1_sparsed = sparse(JJ1); % using sparse J for fast performance
JJ2O_sparsed = sparse(JJ2O);

if size(h,1)==1, h =h';end % handling dimension error
required_colors = length(unique(colorMap)); % num of colors required
Groups = cell(1,required_colors); % groups to be updated in parallel
for k = 1:required_colors
    Groups{k} = find(colorMap==k);
    JJJ1{k} = JJ1_sparsed(Groups{k},:);
	JJJ2{k} = JJ2O_sparsed(Groups{k},:);
    hhh{k} = h(Groups{k});
end

clear JJ1_sparsed JJ2O_sparsed;

numWorkers = 6;
% If a parallel pool exists, shut it down
poolobj = gcp('nocreate');  % Get the current pool if it exists
if ~isempty(poolobj)
    delete(poolobj);
end

% Start a new parallel pool with the desired number of workers
parpool('local', numWorkers);
	
fprintf('Simulation started: %s\n',datetime("now"));


% ******************** Start of the simulation ***************************%
parfor runID = 1:run_per_instance
    %tic
    runID
    seedNum = base_seed+runID;
	sstream = RandStream('mt19937ar','Seed',seedNum);
    RandStream.setGlobalStream(sstream);
    timeIndx = 1;
    BestEn = 1000;
	x = zeros(num_pbits_per_replica*num_replicas,1);
	
	m_last_temp2 = zeros(num_pbits/num_replicas,length(seq_time));
	E_last_temp2 = zeros(1,length(seq_time));
    BestState = zeros(length(h_org),1);

    %% RANDOMIZE THE EXPERIMENT INITIALS
    m_start = sign(2*rand(num_pbits_per_replica*num_replicas,1)-1);


    %tic		
    for jj=1:total_real_num_sweeps % time step loop
        %% GC for loop
		GammaXt = GammaX*(1-jj/total_real_num_sweeps)+0.00001;
		KN = -1*log(tanh(GammaXt*betaAll));  
		for ijk = 1:1:required_colors  % color group loop
			x(Groups{ijk}) = (betaAll*JJJ1{ijk}*num_replicas+KN*JJJ2{ijk})*m_start + betaAll*hhh{ijk};
			m_start(Groups{ijk}) = sign(tanh(x(Groups{ijk}))-2*rand(length(Groups{ijk}),1)+1);
		end
        
        if(ismember(jj,seq_time))
			%tic
			% compute energies of all the replicas
			m_reshaped = reshape(m_start,num_pbits_per_replica,num_replicas);
			EnAll = - 0.5*dot(m_reshaped,W_org*m_reshaped) - h_org'*m_reshaped;

			[minE,minE_idx] = min(EnAll);
			m_last_temp2(:,timeIndx) = m_reshaped(:,minE_idx);
		    E_last_temp2(:,timeIndx) = minE; 
            timeIndx = timeIndx + 1;
        end
    end
	m_last3(:,:,runID) = m_last_temp2;
	E_last3(:,:,runID) = E_last_temp2;
    BestStateAll(:,runID) = BestState;
    %toc
end

% post-processing
m_last = permute(m_last3,[1,3,2]);
E_last = permute(E_last3,[1,3,2]);

% save results
if(logical)
	fname = ['./logicalResultsSQA/size' num2str(instanceSize,'%02d') '/m_instance_' num2str(instanceID,'%04d') '_' num2str(total_real_num_sweeps) '_' num2str(base_seed) '_2_ICM.mat'];
else
	fname = ['./ResultsSQA/size' num2str(instanceSize,'%02d') '/m_instance_' num2str(instanceID,'%04d') '_' num2str(total_real_num_sweeps) '_' num2str(base_seed) '_2_ICM.mat'];
end
save(fname,'m_last','E_last','BestEnAll','BestState','-v7.3');
fprintf('Simulation ended: %s\n',datetime("now"));

