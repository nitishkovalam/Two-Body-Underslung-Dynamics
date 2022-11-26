function [A,B,H]=abh(T,K_N,Ra)
    % A -> Transformation from generalized to configurational coordinates
    % B -> Inverse of A i.e Transformation from configurational to generalized coordinates
    % H -> Relates tension and cable force 
    
    A_23=-T.NB{1}*vp(Ra(:,1));
    A_24=T.NB{2}*vp(Ra(:,2));
    
    B_23=-T.BN{3}*A_23;
    B_24=-T.BN{3}*A_24;
    
    
    A=[eye(3)	 zeros(3)	zeros(3)	zeros(3);
       eye(3)    T.NB{3}     A_23        A_24;
       zeros(3)  zeros(3)    eye(3)      zeros(3);
       zeros(3)  zeros(3)    zeros(3)    eye(3)];
   
    B=[eye(3)    zeros(3)	zeros(3)	zeros(3);
      -T.BN{3}   T.BN{3}     B_23        B_24;
       zeros(3)  zeros(3)    eye(3)      zeros(3);
       zeros(3)  zeros(3)    zeros(3)    eye(3)];
   
    H=[K_N;
      -K_N;
       vp(Ra(:,1))*T.BN{1}*K_N;
      -vp(Ra(:,2))*T.BN{2}*K_N];   % Relates tension and cable force 
   
end