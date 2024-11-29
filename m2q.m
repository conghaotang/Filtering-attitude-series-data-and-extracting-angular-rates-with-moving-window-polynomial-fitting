function [q,max_num]=m2q(A)
%需要再检查 一下
% q4 为 标量
A11=A(1,1);  A22=A(2,2);
A33=A(3,3);  trA = trace(A);
A12=A(1,2);  A21=A(2,1); 
A13=A(1,3);  A31=A(3,1);
A32=A(3,2);  A23=A(2,3);
AA=[A11 A22 A33 trA; 1:4 ];
[AAmax,max_num,max_value]=shengxu(AA,4);

if  max_num==1
    q1=1+2*A11-trA;
    q2=A12+A21;
    q3=A13+A31;
    q4=A23-A32;
else if max_num==2
    q1=A21+A12;
    q2=1+2*A22-trA;
    q3=A23+A32;
    q4=A31-A13;
    else if max_num==3
          q1=A31+A13;
          q2=A32+A23;
          q3=1+2*A33-trA;
          q4=A12-A21;
        else % 即 max_num==4
            q1=A23-A32;
            q2=A31-A13;
            q3=A12-A21;
            q4=1+trA;
        end
    end
end

Q=[q1  q2  q3 q4]';  %q4 为 正 的情况下
q=Q./norm(Q);

%   s = sign([ A(2,3)-A(3,2); 
%              A(3,1)-A(1,3); 
%              A(1,2)-A(2,1);
%                1]); %求数字的符号函数。 正 或 负
%   qz = s.*q;



