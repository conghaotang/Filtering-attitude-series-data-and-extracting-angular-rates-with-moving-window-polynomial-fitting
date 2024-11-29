function [Qdiag] = diagq(QZ)

% Qdiag=zeros(15);

 for i=1:length(QZ)/3
     Qdiag((3*i-2):(3*i),(3*i-2):(3*i))=QZ((3*i-2):(3*i),:);
 end