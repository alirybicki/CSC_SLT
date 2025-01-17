function [y, prob] = rw_softmax_social_reward_sim(r, infStates, p)
% Simulates observations from a Bernoulli distribution
%
% --------------------------------------------------------------------------------------------------
% Copyright (C) 2012-2013 Christoph Mathys, Andreea Diaconescu TNU, UZH & ETHZ
%
% This file is part of the HGF toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.

x_r = infStates(:,1,1);
x_a = infStates(:,1,2);


ze1 = log(p(1));
beta = log(p(2));

% Precision (i.e., Fisher information) vectors
px = 1./(x_a.*(1-x_a));
pc = 1./(x_r.*(1-x_r));

% Weight vectors
wx = ze1.*px./(ze1.*px + pc); % precision first level
wc = pc./(ze1.*px + pc);

% Belief vector
b = wx.*x_a + wc.*x_r;


% Apply the unit-square sigmoid to the inferred states
prob = b.^(exp(beta))./(b.^(exp(beta))+(1-b).^(exp(beta)));

y = binornd(1, prob);

return;
