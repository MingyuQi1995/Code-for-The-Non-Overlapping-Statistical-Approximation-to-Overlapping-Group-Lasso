% function to get the group information to apply SLEP package
% Input: G (group structure)
% Output: opt_ind (see SLEP: Sparse Learning with Efficient Projections for details)





function [opt_ind] = cvt_G_ind(G)
  B = mypar(G);
  p = size(B,2);
  g = size(B,1);
  opt_ind =  zeros(3, g);
  l0 = 1;
  H = sum(G);
  
 for i = 1:g
     tmp_idx = find(B(i,:) == 1);
     l1 = l0 + length(tmp_idx) ;
    opt_ind(1,i) = l0;
    opt_ind(2,i) = l1-1;
    opt_ind(3,i) = mean(H(tmp_idx));
    l0 = l1;

 end
  
  
end