function [g]=qua2g(q)

  qv=q(2:4,1);
  g=qv/q(1);
  
end