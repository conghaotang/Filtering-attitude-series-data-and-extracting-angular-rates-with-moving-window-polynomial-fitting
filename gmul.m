function g = gmul(g1,g2)
 % Gibbs  ∏¡ø œ‡≥À
 
  gfz = g2 + g1 + askew(g2) * g1;
  gfm = 1 - g2' * g1;
  
  g = gfz / gfm;
end