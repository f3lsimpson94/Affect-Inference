%% AI_Simulate_ApproachAvoidTask.m contains the main simulation
%                                             script to simulate the SAP task
%   SYNTAX:     AI_Simulate_ApproachAvoidTask
%
%   AUTHOR:     Coded by: Ryan Smith, Karl J. Friston, Christopher J. Whyte, 2022
%               Originally supplementary Code for: 
%                          "A Step-by-Step Tutorial on Active Inference Modelling and its 
%                          Application to Empirical Data"
%                          https://www.sciencedirect.com/science/article/pii/S0022249621000973)
%
%               Amended by: Katharina V. Wellstein, XX.2024
%                           katharina.wellstein@newcastle.edu.au
%                           Felicity M. Simpson 15.2025
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First, you need to add SPM12, the DEM toolbox of SPM12 and the
% folder with the example scripts to your path in Matlab.

clear all
close all      % These commands clear the workspace and close any figures

rng('shuffle') % This sets the random number generator to produce a different 
               % random sequence each time, which leads to variability in 
               % repeated simulation results (you can alse set to 'default'
               % to produce the same random sequence each time)

options = set_options('AAA');

%% 1. Set up model structure

% Number of time points or 'epochs' within a trial: T
% =========================================================================

% Trials take from HGF optimized trace, based on current simulation, T=160. 

for t = 1:numel(options.sim_nTrials)

    T = options.sim_nTrials(t);
    if t==1
        u = readmatrix('/Users/kwellste/projects/SEPAB/tasks/approach_avoid_task/task/+eventCreator/input_sequence.csv');

        f1_idx      = find(u(:,1)==1);
        f2_idx      = find(u(:,1)==2);
        f1_vec      = zeros(options.sim_nTrials(t),1);
        f1_vec(f1_idx,1) = 1;
        f2_vec      = zeros(options.sim_nTrials(t),1);
        f2_vec(f2_idx,1) = 1;
        welcome_idx = find(u(:,2)==1);
        reject_idx  = find(u(:,2)==0);
        welcome_vec = zeros(options.sim_nTrials(t),1);
        welcome_vec(welcome_idx,1) = 1;
        reject_vec  = zeros(options.sim_nTrials(t),1);
        reject_vec(reject_idx,1) = 1;
    else
        % code for all trial structures
    end

% Priors about initial states: D and d
% =========================================================================

%--------------------------------------------------------------------------
% Specify prior probabilities about initial states in the generative 
% process (D)
% Note: By default, these will also be the priors for the generative model
%--------------------------------------------------------------------------

for a = 1:size(options.sim_agents,1)
    D{a} = options.sim_agents(a,:);
end

d = D;

% Here we seperate the generative process (the capital D)
% from the generative model (the lower case d) allowing learning to occur
% (i.e. to acccumulate concentration paramaters) in the generative model, 
% independent of the generative process.

% probabilistic (likelihood) mapping from hidden states to outcomes: A
%--------------------------------------------------------------------------

%  initialize A matrix
A{1} = zeros(6,6,options.sim_nTrials(t));
for i = 1:options.sim_nTrials(t)
% If no f1 or not the defined response in this trial this will be 0 
    s1 = f1_vec(i)*welcome_vec(i); 
    s2 = f1_vec(i)*reject_vec(i);  
    s3 = f1_vec(i); % unknown is always possible
    s4 = f2_vec(i)*welcome_vec(i); 
    s5 = f2_vec(i)*reject_vec(i);  
    s6 = f2_vec(i); % unknown is always possible
                % f1Welc|f1Rej|f1Unkn|f2Welc|f2Rej|f2Unkn
    A{1}(:,:,i) =   [s1    0     0      0      0     0;   % f1Welc+pAppr
                   0     s2    0      0      0     0;   % f1Rej_pAppr
                   0     0     s3     0      0     0;   % f1Unkown_pAvoid
                   0     0     0      s4     0     0;   % f2Welc_pAppr
                   0     0     0      0      s5    0;   % f2Rej_pAppr
                   0     0     0      0      0     s6]; % f2Unkown_pAvoid
end
   % options.outcomes    = {'f1Welc_pAppr', 'f1Rej_pAppr', 'f1Unkown_pAvoid',...
   %                             'f2Welc_pAppr', 'f2Rej_pAppr', 'f2Unkown_pAvoid'};
        % options.outcomeMods = {'face','response'};
% seperate generative model from generative process
a = A;

% reduce precision
pr1 = 2; % precision (inverse termperature) parameter (lower = less precise)
a{1} = spm_softmax(pr1*log(A{1}+exp(-4)));

a = a{1}*100;

% By passing the a matrix  through a softmax function with a precision paramater of 2
% we slightly reduce the precision of the generative model, analagous to introducing 
% a degree of noise into our model of tone perception. We then multiply it
% by 100 so that the level of noise stays constant across trials.

% Transitions between states: B
%--------------------------------------------------------------------------

% Controlled transitions and transition beliefs : B{:,:,u} and b(:,:,u)
%==========================================================================

%--------------------------------------------------------------------------
% Next, we have to specify the probabilistic transitions between hidden states
% under each action (sometimes called 'control states'). 
% Note: By default, these will also be the transitions beliefs 
% for the generative model
%--------------------------------------------------------------------------

% Columns are states at time t. Rows are states at t+1.

% The agent cannot control the context state, so there is only 1 'action',
% indicating that contexts remain stable within a trial:

B{1}(:,:,1) = [1 0;  % See face1 context
               0 1]; % See face2 Context
           
% Move to approach_welcomed state from any other state
B{2}(:,:,2) = [0.5 0.5 1;  % appr_welcomed
               0 0 0;  % appr_rejected
               0 0 0]; % avoid

% Move to approach_welcomed state from any other state
B{3}(:,:,3) = [0 0 0;  % appr_welcomed
               0.5 0.5 1;  % appr_rejected
               0 0 0]; % avoid

% Move to approach_welcomed state from any other state
B{4}(:,:,4) = [0 0 0;  % appr_welcomed
               0 0 0;  % appr_rejected
               0.5 0.5 1]; % avoid

end