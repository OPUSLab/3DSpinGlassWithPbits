%
% Source code for "Pushing the Boundary of Quantum Advantage in Hard Combinatorial Optimization with Probabilistic Computer"
% by Shuvro Chowdhury et al.
%
% Date: July 2025
%
%
% This code performs Adaptive Parallel Tempering on 3D spin glass instances. 


clc; clearvars; close all;

%*********************************** Inputs ******************************% 
% % loading single replica J, h
logical = 1;   % 1 for logical, 0 for embedded
instanceSize = 16;  
instanceID = 55;
run_per_instance = 50;
total_real_num_sweeps = 10000; % how many total sweeps in MCMC chain
num_internal_replicas = 4; % number of replicas that runs at same temperature per temperature
num_sweeps_per_swap = 1; %how many sweeps in MCMC chain before a swap
num_sweeps_read_per_swap  = 1; %how many sweeps you will read per swap
base_seed  = 223;  % seed to random number generator (arbitrarily chosen)

if(logical)
	pathName = ['./logicalMatInstances/size' num2str(instanceSize,'%02d') '/'];
    betaFilename = ['./logicalSchedules/size' num2str(instanceSize,'%02d')...
        '/instance_' num2str(instanceID,'%04d')...
        '/beta_schedule_alpha_1.25_parallelComputations_100_MCMCSweeps_10000_EnergySweeps_1000.mat'];
else
	pathName = ['./matInstances/size' num2str(instanceSize,'%02d') '/'];
    betaFilename = ['./Schedules/size' num2str(instanceSize,'%02d')...
        '/instance_' num2str(instanceID,'%04d')...
        '/beta_schedule_alpha_0.62_parallelComputations_100_MCMCSweeps_10000_EnergySweeps_1000.mat'];
end

% load saved J matrices
JFileName = ['JOrig_' num2str(instanceID,'%04d') '.mat'];
load (JFileName,'J');

cMapFileName = [pathName 'colorMap.csv'];

% MCS steps where we will computeresidual energy
seq_time = [2,4,10,12,20,30,40,50,60,70,80,90,100,128,150,176,200,250,300,...
    400, 500,600,700,800, 900, 1000, 1280, 1500, 1750, 2000, 2500, 3000, ...
    4000, 5000, 6000, 7000,8000,9000,10000];

%******************************* Pre-processing **************************% 
addpath(pathName);

% for the DWAVE problem, h is all 0
h = zeros(length(J),1);

norm_factor = full(max(max(abs((J)))));
W_org = sparse(J); h_org = sparse(h); % original values from 1 replica
clear J h;

load(betaFilename);
num_replicas = length(saveSigma);
startingBeta = 1;
selectedBeta = startingBeta: startingBeta+num_replicas-1; % choose beta that you select

I0 = beta(selectedBeta);
I0 =I0/norm_factor;

betaAll = kron(I0,ones(1,num_internal_replicas));
num_betas = num_replicas*num_internal_replicas;

num_pbits_per_replica = length(W_org);

for ii = 1:1:num_pbits_per_replica
    [~,nbr{ii},~] = find(W_org(ii,:));
end

num_swap_attempts= round(total_real_num_sweeps/num_sweeps_per_swap);
num_sweeps_read = num_sweeps_read_per_swap*num_swap_attempts; %how many sweeps you will read in total

num_pbits= length(W_org)*num_betas;

save seq_time_size_20_516_ICM.mat seq_time

% placeholders
m_last3 = zeros(num_pbits/num_replicas/num_internal_replicas,length(seq_time),run_per_instance);
E_last3 = zeros(1,length(seq_time),run_per_instance);

BestState = zeros(length(h_org),1);
BestStateAll = zeros(length(h_org),run_per_instance);
BestEnAll = zeros(1,run_per_instance,length(seq_time));

colorMap = readmatrix(cMapFileName);

JJ2 = sparse(W_org); % using sparse J for fast performance
h = h_org;
if size(h,1)==1, h =h';end % handling dimension error
required_colors = length(unique(colorMap)); % number of colors required
Groups = cell(1,required_colors); % groups to be updated in parallel
for k = 1:required_colors
    Groups{k} = find(colorMap==k);
    JJJ{k} = JJ2(Groups{k},:);
    hhh{k} = h(Groups{k});
end

clear JJ2 h

numWorkers = 10;
% If a parallel pool exists, shut it down
poolobj = gcp('nocreate');  % Get the current pool if it exists
if ~isempty(poolobj)
    delete(poolobj);
end

% Start a new parallel pool with the desired number of workers
parpool('local', numWorkers);
	
fprintf('Simulation started: %s\n',datetime("now"));

% ******************** Start of the simulation ***************************%
parfor (runID = 1:run_per_instance, numWorkers)
    %tic
    runID
    seedNum = base_seed+runID;
	sstream = RandStream('mt19937ar','Seed',seedNum);
    RandStream.setGlobalStream(sstream);
    timeIndx = 1;
    BestEn = 1000;
	x_reshaped = zeros(num_pbits_per_replica,num_replicas*num_internal_replicas);
	
	m_last_temp2 = zeros(num_pbits/num_replicas/num_internal_replicas,length(seq_time));
	E_last_temp2 = zeros(1,length(seq_time));
    BestState = zeros(length(h_org),1);

    % generate all possible consecutive pairs of replicas
    all_pairs = reshape(1:num_replicas-1, [], 1) + [0 1];

    % preallocate cell array to store selected pairs for all attempts
    all_selected_pairs = cell(num_swap_attempts, 1);

    % loop over swap attempts
    for attempt = 1:num_swap_attempts
        if(mod(attempt,2))
            selected_pairs = all_pairs(1:2:end,:);
        else
            selected_pairs = all_pairs(2:2:end,:);
        end
        all_selected_pairs{attempt} = selected_pairs;
    end

    %% RANDOMIZE THE EXPERIMENT INITIALS
    m_reshaped = sign(2*rand(num_pbits_per_replica,num_replicas*num_internal_replicas)-1);

    for ii = 1:1:num_swap_attempts
	%ii
        %tic		
        for jj=1:num_sweeps_per_swap % time step loop
            %% GC for loop
			for bij = 1:1:num_betas
				for ijk = 1:1:required_colors  % color group loop
					x_reshaped(Groups{ijk},bij) = betaAll(bij)*(JJJ{ijk}*m_reshaped(:,bij) + hhh{ijk});
					m_reshaped(Groups{ijk},bij) = sign(tanh(x_reshaped(Groups{ijk},bij))-2*rand(length(Groups{ijk}),1)+1);
				end
			end
        end
        %t2imes(ii) = toc;
        
		%tic
        % now we have to perform the Houdayer move
        for t_reps = 1:1:num_replicas   % loop over temperature replicas
            permuted_i_reps = randperm(num_internal_replicas);  % randomly permute internal replicas
            for i_reps = 1:2:num_internal_replicas
				start_1 = (t_reps-1)*num_internal_replicas+(permuted_i_reps(i_reps));
                start_2 = (t_reps-1)*num_internal_replicas+(permuted_i_reps(i_reps+1));
                [replicas_1,replicas_2,~] = houdayer_cluster_move(nbr,m_reshaped(:,start_1),m_reshaped(:,start_2));
                m_reshaped(:,start_1) = replicas_1;
                m_reshaped(:,start_2) = replicas_2;
				
            end
        end
		%t3imes(ii) = toc;
		
		%tic
		% compute energies of all the replicas
		EnAll = - 0.5*dot(m_reshaped,W_org*m_reshaped) - h_org'*m_reshaped;
        %t5imes(ii) = toc;
		
		%tic
        for jj = 1:1:size(all_selected_pairs{ii},1)
            % Select a random replica
            sel = all_selected_pairs{ii}(jj,1);
            next = all_selected_pairs{ii}(jj,2);
			
			start_1 = (sel-1)*num_internal_replicas;
			start_2 = (next-1)*num_internal_replicas;

            for i_reps = 1:1:num_internal_replicas

				E_sel = EnAll(start_1+i_reps);
				E_next = EnAll(start_2+i_reps);
				
                if(E_sel<BestEn)
                    BestEn = E_sel;
                    BestState = m_reshaped(:,start_1+i_reps);
                end
                if(E_next<BestEn)
                    BestEn = E_next;
                    BestState = m_reshaped(:,start_2+i_reps);
                end
				
                beta_sel = I0(sel);
                beta_next = I0(next);

                DeltaE = E_next-E_sel;
                DeltaB = beta_next - beta_sel;

                if rand<min(1,exp(DeltaB*DeltaE))
					m_reshaped(:,[start_1+i_reps,start_2+i_reps]) = m_reshaped(:,[start_2+i_reps,start_1+i_reps]);
					
					EnAll([start_1+i_reps,start_2+i_reps]) = [E_next E_sel];
                end
            end

        end
        %t4imes(ii) = toc;

        if(ismember(ii*num_sweeps_per_swap,seq_time))
            min_indices = zeros(num_replicas,1);
            Energy = zeros(num_replicas,1);
            for look_replica = 1:num_replicas
                [Energy(look_replica),min_indices(look_replica)] = min(EnAll((look_replica-1)*num_internal_replicas+1: look_replica*num_internal_replicas));
            end

			[minE,minE_idx] = min(Energy);
			m_last_temp2(:,timeIndx) = m_reshaped(:,(minE_idx-1)*num_internal_replicas+(min_indices(minE_idx)));
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
	fname = ['./logicalResultsAPT/size' num2str(instanceSize,'%02d') '/m_instance_' num2str(instanceID,'%04d') '_' num2str(total_real_num_sweeps) '_' num2str(base_seed) '_2_ICM.mat'];
else
	fname = ['./ResultsAPT/size' num2str(instanceSize,'%02d') '/m_instance_' num2str(instanceID,'%04d') '_' num2str(total_real_num_sweeps) '_' num2str(base_seed) '_2_ICM.mat'];
end
save(fname,'m_last','E_last','BestEnAll','BestState','-v7.3');

fprintf('Simulation ended: %s\n',datetime("now"));
