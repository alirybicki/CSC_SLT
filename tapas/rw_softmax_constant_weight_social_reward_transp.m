function [pvec, pstruct] = rw_softmax_constant_weight_social_reward_transp(r, ptrans)
% --------------------------------------------------------------------------------------------------
% Copyright (C) 2012 Christoph Mathys, TNU, UZH & ETHZ
%
% This file is part of the HGF toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.

pvec    = NaN(1,length(ptrans));
pstruct = struct;

pvec(1)     = sgm(ptrans(1),1);       % ze1
%pvec(1)     = ptrans(1); %when ze1 is fixed
pstruct.ze1 = pvec(1);

pvec(2)     = exp(ptrans(2));         % ze2
%pvec(2)     = ptrans(2);  %when ze2 is fixed 
pstruct.ze2 = pvec(2);

return;