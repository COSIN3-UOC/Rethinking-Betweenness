% clc;
addpath('./mcm');
clear variables;
% %close all;
seed = ceil(rand()*100);

% if isempty(gcp('nocreate'))
%     parpool(6);
% end

numLayers = 1;
timeSteps = 10000;
% max number of packets thacet can jump from each node in 1 time step
tau = 1;


%tol_ls = [30];
tol_ls = [0, 10, 20, 30];
    
%possGenRates = logspace(-2,log10(3),100);
possGenRates = linspace(0.0001,0.008,10);
numRhos = numel(possGenRates);

nstat        = zeros(numRhos,1);
nstat_std    = zeros(numRhos,1);
noutstat     = zeros(numRhos,1);
noutstat_std = zeros(numRhos,1);

str_mat = sprintf('./mcm/other_graphs/bcn_amb_full/pen_path_ls_t_03_k_20_p_01.h5');

PathLisArr = h5paths_to_cell(str_mat, '/paths');
numNodes = length(PathLisArr)

for tol_i = 1:length(tol_ls)

    tol = tol_ls(tol_i);

    
    % counter for completed iterations

    %h = waitbar(0, 'Initializing waitbar...');  

    DelivExp   = zeros(numRhos,1);
    DeltaNExp  = zeros(numNodes,numRhos);

    % parfor rhosIndex = 1:numRhos
    for rhosIndex = 1:numRhos
        rhosIndex
        %waitbar_label = sprintf('Processing rho %g for tol %d.', ...
           % possGenRates(rhosIndex) ,tol);
        
        %waitbar(rhosIndex / numRhos, ...
        %        h, waitbar_label)
        % setting node generation rates
        genRatePerMinuteGates = zeros(numNodes,1);
        genRatePerMinuteGates(:) = possGenRates(rhosIndex);       

        %define the "OD" matrix somehow, this is all to all but with the
        %given generationrate

        [possDest,possDestCumProb,possDestProb] =...
            matrixOfAccesibleDestinations_AllToAll(...
                find(genRatePerMinuteGates>0),...
                find(genRatePerMinuteGates==0),...
                genRatePerMinuteGates,...
                numNodes);   
        genRatePerMinute = genRatePerMinuteGates;        
        processingRatePerMinute = ones(numNodes,1);      

        % monte carlo simulations adapted for tolerance bwss
        [numPackets,numPacketOutUnitTime,B_i_ext,e_i_ext,B_i,e_i,s_i,PExp,sigma,DExp,...
        extractedPackets,packetsThatContinue]...
            = cSPCongestion_statMem_dir_weighted_local_list_paths_no_pad(...
                genRatePerMinute,...
                processingRatePerMinute,...
                tol,...%SPModelData.DistRG,...
                [],...%SPModelData.NumPathsRG,...
                PathLisArr,...%path list cell array
                possDest,...
                possDestCumProb,...
                numNodes,1,timeSteps,seed); 
          
        % --- statistics for this rho ---

        % Add stationary number of packets and out rate
        tot_n = sum(numPackets, 1, "default");

        nstat(rhosIndex) = mean(tot_n(end-999:end));
        nstat_std(rhosIndex) = std(tot_n(end-999:end));

        noutstat(rhosIndex) = mean(numPacketOutUnitTime(end-999:end));
        noutstat_std(rhosIndex) = std(numPacketOutUnitTime(end-999:end));
        % compute microscopic vars obtained by the montecarlo
    
        normFactor = 1./sum(PExp);
        normFactor(normFactor == Inf) = 0;
        % Order parameter
        PExp = PExp*diag(normFactor);    
        pExp = packetsThatContinue./extractedPackets;
        DExp = mean(DExp,2);

        % compute queue increments obtained by the montecarlo
        for j=1:numNodes
            p = polyfit(1:size(numPackets(j,:),2),numPackets(j,:),1);
            DeltaNExp(j,rhosIndex) = (p(1)); 
        end      

        DelivExp(rhosIndex) = sum(numPacketOutUnitTime);%/sum(s_i, "all");
    end

    % Write file for each graph

    fname = sprintf('./sim_results/all_graphs/penal/bcn_full/tol%d_vx2_hires.mat', tol);
    save(fname, 'DeltaNExp','DelivExp', ...
    'nstat','nstat_std','noutstat', ...
    'noutstat_std')

    %close(h);

end

