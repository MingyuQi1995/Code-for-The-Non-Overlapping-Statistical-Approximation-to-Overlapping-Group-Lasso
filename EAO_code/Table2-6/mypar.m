%Function to parition groups
%Input: overlapping group structure G
$Output: non-overlapping group sturcture G1

function [G1] = mypar(G)
 % g = size(G,1);
  p = size(G,2);
  H = sum(G);
  
  G1 = zeros(p, p);
  k = 1;
  l0 = 1;
  l = 1;
  flag = 1;
  while(flag)
      if(l == p)
        G1(k, l0:l) = 1;
      flag = 0;
   else
       if(H(l + 1) == H(l))
      l = l + 1;
       else
      G1(k, l0:l) = 1;
      k = k + 1;
      l = l + 1;
      l0 = l; 
       end
     end
  end

 
  G1 = G1(1:k,:);
  
 
  
end