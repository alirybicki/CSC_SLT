function [pvec, pstruct] = tapas_rw_social_shift_transp(r, ptrans)
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
pvec(2)       = tapas_sgm(ptrans(2),1); % alpha_r_shift
pstruct.al_shift_r  = pvec(2);
pvec(3)       = tapas_sgm(ptrans(3),1); % alpha_r_stay
pstruct.al_stay_r  = pvec(3);
pvec(4)       = tapas_sgm(ptrans(4),1); % va_0
pstruct.va_0  = pvec(4);
pvec(5)       = tapas_sgm(ptrans(5),1); % alpha_a_shift
pstruct.al_shift_a  = pvec(5);
pvec(6)       = tapas_sgm(ptrans(6),1); % alpha_a_stay
pstruct.al_stay_a  = pvec(6);

return;