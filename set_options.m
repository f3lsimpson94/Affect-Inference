function options = set_options(task)

%% set_options.m specifies all the options that are needed to run the simulation 
%               (e.g., paths, simulation settings, etc.)
%   SYNTAX:     options = set_options
%
%   OUT:        options: struct containing settings and information
%                        needed to run these simulations
%
%   AUTHOR:     Coded by: Katharina V. Wellstein, XX.2024
%                         katharina.wellstein@newcastle.edu.au 
%                         (based on tutorial code by Smith et al., 2022,
%                          https://www.sciencedirect.com/science/article/pii/S0022249621000973)
% -------------------------------------------------------------------------
% This file is released under the terms of the GNU General Public Licence
% (GPL), version 3. You can redistribute it and/or modify it under the
% terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any
% later version.
%
% This file is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
% more details.
% 
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <https://www.gnu.org/licenses/>.
% _______________________________________________________________________________%
 
%% SETTING INDIVIDUAL ENVIRONMENT OPTIONS

% add toolbox path
options.codeDir = pwd;
options.figDir  = [pwd,'/figures'];
addpath(genpath('/Users/kwellste/projects/Toolboxes/spm/spm12')) % SPM12 for model inversion, specifically VBA

if ~exist(options.figDir,'dir')
    mkdir(options.figDir);
end
%% SET SIMULATION TYPE

switch task
    case 'AAA'
        options.sim_agents = [1 1 1 1 1 1;    % Agent1: "agnostic"
                             2 2 2 2 0 0;     % Agent2: "adaptive"
                             0 2 0 2 1 1;     % Agent3: "rejectBias"
                             2 0 2 0 0 0;     % Agent4: "acceptBias"
                             0.5 2 0.5 2 1 1; % Agent5: "high blunted"
                             1 1 1 1 2 2];    % Agent6: "patient blunted"
        options.sim_nTrials = [140, 280, 420];
        options.states    = {'f1Welc', 'f1Rej', 'f1Unkown',...
                             'f2Welc', 'f2Rej', 'f2Unkown'};
        options.outcomes  = {'f1Welc_pAppr', 'f1Rej_pAppr', 'f1Unkown_pAvoid',...
                               'f2Welc_pAppr', 'f2Rej_pAppr', 'f2Unkown_pAvoid'};
        options.policies  = {'appr_welcomed', 'appr_rejected','avoid'};
        % options.outcomeMods = {'face','response'};
        % options.outcome(1).type = {'face1','face2'};
        % options.outcome(2).type = {'welcome','reject','unknown'};


    case 'SAP'
        options.sim_agents = [1 1 1 1 1 1;     % Agent1: "agnostic"
                             2 2 2 2 0 0;     % Agent2: "adaptive"
                             0 2 0 2 1 1;     % Agent3: "rejectBias"
                             2 0 2 0 0 0;     % Agent4: "acceptBias"
                             0.5 2 0.5 2 1 1; % Agent5: "high blunted"
                             1 1 1 1 2 2];    % Agent6: "patient blunted"
        options.sim_nTrials = [120, 240, 360];
end



% If options.doSimTyp = 1, simulate single trial. This will reproduce fig. 8. (Although
% note that, for this and the following simulations, results
% will vary each time due to random sampling)

% If options.doSimTyp = 2, simulate multiple trials where the left context is active
% (D{1} = [1 0]'). This will reproduce fig. 10.

% If options.doSimTyp = 3, simulate reversal learning, where the left context is active
% (D{1} = [1 0]') in early trials and then reverses in later
% trials (D{1} = [0 1]'). This will reproduce fig. 11.

% If options.doSimTyp = 4, run parameter estimation on simulated data with reversal
% learning. This will reproduce the top panel of fig. 17.

% If options.doSimTyp = 5, run parameter estimation on simulated data with reversal
% learning from multiple participants under different models
% (i.e., different parameter values) and perform model comparison.
% This will reproduce the bottom panel of fig. 17. This option
% will also save two structures that include results of model
% comparison, model fitting, parameter recoverability analyses,
% and inputs needed for group (PEB) analyses.

options.doSimType = 3;

%% SET SIMULATION PARAMETERS
options.rs1 = 4; % Risk-seeking parameter (set to the variable rs below)

% To reproduce fig. 8, use values of 4 or 8 (with Sim = 1)
% To reproduce fig. 10, use values of 3 or 4 (with Sim = 2)
% To reproduce fig. 11, use values of 3 or 4 (with Sim = 3)
% This will have no effect on Sim = 4 or Sim = 5

% When Sim = 5, if PEB = 1 the script will run simulated group-level
% (Parametric Empirical Bayes) analyses.

options.doPEB = 0; 
% Note: GCM_2 and GCM_3 (the inputs to PEB; see below) are saved
% after running Sim = 5 to avoid needing to re-run it each time
% you want to use PEB (i.e., because Sim = 5 takes a long time).
% After running Sim = 5 once, you can simply load GCM_2 and GCM_3 and
% run the PEB section separately if you want to come back to it later.

% You can also run the sections separately after building the model by
% simply clicking into that section and clicking 'Run Section' above
end
