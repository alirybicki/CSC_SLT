function [pvec, pstruct] = tapas_rw_social_reward_unrew_vol_transp(r, ptrans)
% --------------------------------------------------------------------------------------------------
% Copyright (C) 2012-2013 Christoph Mathys, Andreea Diaconescu TNU, UZH & ETHZ
%
% This file is part of the HGF toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.

pvec    = NaN(1,length(ptrans));
pstruct = struct;

pvec(1)       = tapas_sgm(ptrans(1),1); % vr_0
pstruct.vr_0  = pvec(1);
pvec(2)       = tapas_sgm(ptrans(2),1); % alpha_r_unrew_v
pstruct.al_unrew_v_r  = pvec(2);
pvec(3)       = tapas_sgm(ptrans(3),1); % alpha_r_volatile
pstruct.al_rew_v_r  = pvec(3);
pvec(4)       = tapas_sgm(ptrans(4),1); % alpha_r_s
pstruct.al_unrew_s_r  = pvec(4);
pvec(5)       = tapas_sgm(ptrans(5),1); % alpha_r_s
pstruct.al_rew_s_r  = pvec(5);

pvec(6)       = tapas_sgm(ptrans(6),1); % va_0
pstruct.va_0  = pvec(6);
pvec(7)       = tapas_sgm(ptrans(7),1); % alpha_a_stable
pstruct.al_unrew_v_a  = pvec(7);
pvec(8)       = tapas_sgm(ptrans(8),1); % alpha_a_volatile
pstruct.al_rew_v_a  = pvec(8);
pvec(9)       = tapas_sgm(ptrans(9),1); % alpha_a_volatile
pstruct.al_unrew_s_a  = pvec(9);
pvec(10)       = tapas_sgm(ptrans(10),1); % alpha_a_volatile
pstruct.al_rew_s_a  = pvec(10);

return;
