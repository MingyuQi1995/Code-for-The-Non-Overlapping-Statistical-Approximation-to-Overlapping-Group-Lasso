% function to get the group information to apply SPAM package
% Input: G (group structure)
% Output:ov,nov, eta,treegroups (see SPAM for details)


function [ov,nov,eta, tgroup] = treespamcvt_OG_ind(G)
  p =  size(G,2);
  g =  size(G,1);

 ov = int32(zeros(1, g));

  for i = 1:g
     j = g - i + 1;
     tmpd = find(G(j,:) == 1);
      ov(i) = p - max(tmpd );
  end  


 nov = int32(zeros(1, g));

 for i = 1:(g-1)
     j = g - i + 1;
     tmpd = G(j,:) - G((j-1),:);
     nov(i) = sum(tmpd );
 end  
  
nov(g) = sum(G(1,:) );
 
eta = flip(sqrt(sum(G ~= 0, 2)));
eta = eta.^4;
eta = 1./eta;
%eta = 1;

tgroupint = zeros(g, g);
tgroupsub = diag(ones((g-1), 1));
tgroupint(2:g,1:(g-1)) = tgroupsub;
tgroup=sparse(tgroupint);
end