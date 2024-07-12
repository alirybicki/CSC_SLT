function [traj, infStates] = tapas_rw_social_shift_vol(r, p, varargin)
% Calculates the trajectories of v under the Rescorla-Wagner learning model
%
% This function can be called in two ways:
%
% (1) tapas_rw_social_reward(r, p)
%
%     where r is the structure generated by tapas_fitModel and p is the parameter vector in native space;
%
% (2) tapas_rw_social_reward(r, ptrans, 'trans')
%
%     where r is the structure generated by tapas_fitModel, ptrans is the parameter vector in
%     transformed space, and 'trans' is a flag indicating this.
%
% --------------------------------------------------------------------------------------------------
% Copyright (C) 2012-2013 Christoph Mathys, Andreea Diaconescu TNU, UZH & ETHZ
%
% This file is part of the HGF toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.

% Transform paramaters back to their native space if needed
if ~isempty(varargin) && strcmp(varargin{1},'trans');
    p = tapas_rw_social_reward_unrew_vol_transp(r, p);
end

r_vol = r.vol{1};
a_vol = r.vol{2};
choice = [0; r.y]; % add dummy "zeroth" trial

% Unpack parameters
vr_0 = p(1);
al_shift_v_r  = p(2);
al_stay_v_r  = p(3);
al_shift_s_r  = p(4);
al_stay_s_r  = p(5);
va_0 = p(6);
al_shift_v_a  = p(7);
al_stay_v_a  = p(8);
al_shift_s_a  = p(9);
al_stay_s_a  = p(10);

% Add dummy "zeroth" trial
u_r = [0; r.u(:,1)]; % JC 1/12/15 changed from u to u_r
u_a = [0; r.u(:,2)]; % JC 1/12/15 added u_a
n = length(u_r);
r_vol = [0 r_vol];
a_vol = [0 a_vol];

% Initialize updated quantity: value
v_r  = NaN(n,1);
da_r = NaN(n,1);

v_a  = NaN(n,1);
da_a = NaN(n,1);

% Priors
v_r(1) = vr_0;
v_a(1) = va_0;

% Pass through value update loop
for k = 2:1:n
    if not(ismember(k, r.ign))
        
        if choice(k) == choice(k-1) % if this is a stay trial, then use that learning rate
            if r_vol(k) % if this is a volatile trial, then use that learning rate
                al_r = al_stay_v_r;
            else
                al_r = al_stay_s_r;
            end
            
            if a_vol(k) % if this is a volatile trial, then use that learning rate
                al_a = al_stay_v_a;
            else
                al_a = al_stay_s_a;
            end
            
        else % if this is a shift trial, then use that learning rate
            if r_vol(k) % if this is a volatile trial, then use that learning rate
                al_r = al_shift_v_r;
            else
                al_r = al_shift_s_r;
            end
            
            if a_vol(k) % if this is a volatile trial, then use that learning rate
                al_a = al_shift_v_a;
            else
                al_a = al_shift_s_a;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%
        % Effect of input u(k)
        %%%%%%%%%%%%%%%%%%%%%%
        
        % Prediction error
        da_r(k) = u_r(k)-v_r(k-1);
        da_a(k) = u_a(k)-v_a(k-1); % JC 1/12/15 added u_a
        
        % Value
        v_r(k) = v_r(k-1)+al_r*da_r(k);
        v_a(k) = v_a(k-1)+al_a*da_a(k);
    else
        da_r(k) = 0;
        da_a(k) = 0;
        v_r(k)  = v_r(k-1);
        v_a(k)  = v_a(k-1);
    end
end

% Predicted value
vhat_r = v_r;
vhat_r(end) = [];
vhat_a = v_a;
vhat_a(end) = [];

% Remove representation priors
v_r(1)  = [];
da_r(1) = [];
v_a(1)  = [];
da_a(1) = [];

% Create result data structure
traj = struct;

traj.v_r     = v_r;
traj.vhat_r  = vhat_r;
traj.da_r    = da_r;

traj.v_a     = v_a;
traj.vhat_a  = vhat_a;
traj.da_a    = da_a;

% Create matrix (in this case: vector) needed by observation model
infStates = NaN(n-1,1,4); % trials, levels, trajectories
infStates(:,:,1) = traj.vhat_r;
infStates(:,:,2) = traj.vhat_a;

return;
