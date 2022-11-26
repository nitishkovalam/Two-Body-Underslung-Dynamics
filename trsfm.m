function [T,WI,K_N,phc] = trsfm(r,Ra,t)
% trsfm -> contains all the transformation between reference frames
% T.BN(i) -> Transforms from inertial frame to i-th body frame.
% T.NB(i) -> inverse of T.BN(i) i.e. transformation from i-th body frame to inertial frame
% i = 1 -> Helicopter
%   = 2 -> RLV
%   = 3 -> Cable (well, cable is not a body but transformation to that frame is required)
%
% WI(i) -> Transformation from body rates and Euler angle rates for i-th body
% W(i)  -> Transformation from Euler angle rates to body rates for i-th body
%
% K_N -> Cable direction vector in inertial frame
%        (along the cable from helicopter attachement point to RLV attachment point)
%
%

global n PLT;
    for ii=1:n
        ph=r(3*n+(ii-1)*3+1);   % phi
        th=r(3*n+(ii-1)*3+2);   % theta
        ps=r(3*n+(ii-1)*3+3);   % psi
        
        T.BN{ii}=[ cos(ps)*cos(th)                                sin(ps)*cos(th)                               -sin(th);
                   sin(ph)*cos(ps)*sin(th)-cos(ph)*sin(ps)     sin(ph)*sin(ps)*sin(th)+cos(ph)*cos(ps)    sin(ph)*cos(th);
                   cos(ph)*cos(ps)*sin(th)+sin(ph)*sin(ps)     cos(ph)*sin(ps)*sin(th)-sin(ph)*cos(ps)    cos(ph)*cos(th)];
        
        WI{ii}=[ 1      sin(ph)*tan(th)    cos(ph)*tan(th);
                 0      cos(ph)             -sin(ph) ;
                 0      sin(ph)/cos(th)    cos(ph)/cos(th) ]; 
        
        T.NB{ii}=T.BN{ii}';
    end
        
    %%

        K_N=1*(r(4:6)+T.NB{2}*Ra(:,2)-r(1:3)-T.NB{1}*Ra(:,1));
        leng=sqrt(sum((K_N).^2));
        K_N=K_N/leng;
        phc=asin(-K_N(2));
        thc=asin(K_N(1)/cos(phc));
        
        addpoints(PLT.c,t,thc*180/pi);
        addpoints(PLT.ph,t,phc*180/pi);
%         drawnow limitrate
        
        T.BN{n+1}=[ cos(thc)               0            -sin(thc);
                    sin(phc)*sin(thc)     cos(phc)    sin(phc)*cos(thc);
                    cos(phc)*sin(thc)    -sin(phc)    cos(phc)*cos(thc)];
        
%         W{n+1}=[ 1    0            -sin(thc);
%                  0    cos(phc)    sin(phc)*cos(thc) ;
%                  0    -sin(phc)   cos(phc)*cos(thc) ]; 
        
        
        T.NB{n+1}=T.BN{n+1}';
end