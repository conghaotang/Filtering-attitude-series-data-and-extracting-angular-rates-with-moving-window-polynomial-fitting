function [g]=q2g(q)

  qv=q(1:3,1);
  g=qv./q(4,1);
  
end