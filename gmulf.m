function gf = gmulf(g1,g2,g1f)
 % Gibbs  ∏¡ø œ‡≥À
 
  gfz = g2 + g1 + askew(g2) * g1;
  gfm = 1 - g2' * g1;
  
  gfzf =  g1f + askew(g2) * g1f ;
  gfmf = - g2' * g1f;
  
  gf = 1/gfm^2 * (gfzf * gfm - gfz *  gfmf);
  
end