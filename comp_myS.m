function mS = comp_myS(x)

[ n, ~ ] = size(x); 
mS = eye(n-1,n);
mS = mS - [ zeros(n-1,1), mS(:,1:n-1) ];

end