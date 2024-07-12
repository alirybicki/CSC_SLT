function pstruct = tapas_rw_social_shift_vol_namep(pvec)
% --------------------------------------------------------------------------------------------------
% Copyright (C) 2012-2013 Christoph Mathys, TNU, UZH & ETHZ
%
% This file is part of the HGF toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.

pstruct = struct;

pstruct.vr_0  = pvec(1);
pstruct.al_shift_v_r  = pvec(2);
pstruct.al_stay_v_r  = pvec(3);
pstruct.al_shift_s_r  = pvec(4);
pstruct.al_stay_s_r  = pvec(5);

pstruct.va_0  = pvec(6);
pstruct.al_shift_v_a  = pvec(7);
pstruct.al_stay_v_a  = pvec(8);
pstruct.al_shift_s_a  = pvec(9);
pstruct.al_stay_s_a  = pvec(10);

return;

