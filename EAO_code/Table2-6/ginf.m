% Function to get the group information to apply SLEP package
% Input: G (group structure)
% Output: opt_G and opt_ind (see SLEP: Sparse Learning with Efficient Projections for details)

function [opt_G,opt_ind] = ginf(G,type)

  if strcmp(type, '0')
  p =  size(G,2);
  g =  size(G,1);
  opt_G = [];
  
 opt_ind = zeros(3, g);

  l0 = 1;
  
 for i = 1:g
     tmp_idx = find(G(i,:) == 1);
  
    l1 = l0 + length(tmp_idx) ;
    opt_G = [opt_G, tmp_idx];
    opt_ind(1,i) = l0;
    opt_ind(2,i) = l1-1;
    opt_ind(3,i) = sqrt(length(tmp_idx));
    l0 = l1;
 end

  elseif strcmp(type, '2')
     orw = sqrt(length( find(G(1,:) == 1)));
     G1 = mypar(G);
     H = sum(G);
     k = size(G1,1);
     col =  size(G1,2);
     H1 = [];
  
 for i = 1:k
     tmp_idx = find(G1(i,:) == 1);
     H1(i)  = mean(H(tmp_idx))*orw;

 end
     opt_ind= H1;
  
     opt_G = [];
  
  for i = 1:col
     idx = find(G1(:,i) == 1);
     opt_G(i) = idx;

  end

    elseif strcmp(type, '3')
     G1 = mypar(G);
     H = sum(G);
     k = size(G1,1);
     col =  size(G1,2);
     opt_ind= ones(1, k);
  
     opt_G = [];
  
  for i = 1:col
     idx = find(G1(:,i) == 1);
     opt_G(i) = idx;

  end

   elseif strcmp(type, '4')
     G1 = mypar(G);
     H = sum(G);
     k = size(G1,1);
     col =  size(G1,2);
     H1 = [];
  
 for i = 1:k
     tmp_idx = find(G1(i,:) == 1);
     H1(i)  = sqrt(length(tmp_idx));
 end
     opt_ind= H1;
  
     opt_G = [];
  
  for i = 1:col
     idx = find(G1(:,i) == 1);
     opt_G(i) = idx;

 end
  
  else
        disp('Invalid type. No action performed.');
    end
end