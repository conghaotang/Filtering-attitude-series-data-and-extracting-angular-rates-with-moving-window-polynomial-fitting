function qnb = qnormlz(qnb)
% Quaternion normalization, so ||qnb||=1.
%
% Prototype: qnb = qnormlz(qnb)
% Input: qnb - input quaternion whose norm may not be 1
% Output: qnb - input quaternion whose norm equals 1
%
% See also  vnormlz, mnormlz.
 
    norq = norm(qnb);
    qnb = qnb/norq;