function pstruct = tapas_rw_reward_unrew_vol_namep(pvec)
% --------------------------------------------------------------------------------------------------
% Copyright (C) 2012-2013 Christoph Mathys, TNU, UZH & ETHZ
%
% This file is part of the HGF toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.

pstruct = struct;


pstruct.v_0   = pvec(1);
pstruct.al_unrew_v    = pvec(2);
pstruct.al_rew_v    = pvec(3);
pstruct.al_unrew_s    = pvec(4);
pstruct.al_rew_s   = pvec(5);

return;
