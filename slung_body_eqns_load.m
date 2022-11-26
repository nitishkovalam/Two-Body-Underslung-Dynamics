function dzdt = slung_body_eqns_load (t,z,Ra,L,D,DI,g,F,FA,vw)

global PLT

% % uncomment one of the below
% drawnow
% drawnow limitrate
% % ^ 'drawnow limitrate' gives faster performance but is only available in 2016 i think


% Writing the state variables in two body form ( Helicopter and RLV)
r=[F.r{1}(t);F.r{2}(t);F.r{3}(t);z(7:9);F.r{4}(t);F.r{5}(t);F.r{6}(t);z(10:12)];
v=[F.v{1}(t);F.v{2}(t);F.v{3}(t);z(1:3);F.v{4}(t);F.v{5}(t);F.v{6}(t);z(4:6)];
a_hel=[F.a{1}(t);
         F.a{2}(t);
         F.a{3}(t);
         F.a{4}(t);
         F.a{5}(t);
         F.a{6}(t);
        ];

 % trsfm -> contains all the transformation between reference frames
[T,WI,K_N,phc] = trsfm(r,Ra,t);

% abh -> contains the transformation matrix to the generalized co-ordinates
[A,B,H]=abh(T,K_N,Ra);


u0=B*v;             


Adot = dA_dt(u0,phc,L,T,Ra);

%% 
% fa ->  Aerodynamic Forces and Moments
% Aero DATA gives forces in body frame 
fa = force_aero_load(z,v,T,FA,vw,t);

% ficg -> Inertia Force + Coriolis Force + Gravity Force + Inertia Coupling
ficg = force_icg_load(g,D,Adot,u0);

f0=(fa)+ficg;
% Sum of External Forces 


C = invr_s_compute(H,D,A,DI);

    x = C*[ a_hel(1:3);
            f0(1:3);
            a_hel(4:6);
            f0(4:6)
            ];

dudt = [ a_hel(1:3);
         x(4:6);
         a_hel(4:6);
         x(10:12)
         ];
     
dbdt=Adot*u0+A*dudt;

dzdt=[dbdt([4:6 10:12]); 
      z(1:3);
      WI{2}*z(4:6)];

%   Ene=1/2*D(4,4)*z(1:3)'*z(1:3)+1/2*(D(10:12,10:12)*z(4:6))'*z(4:6)-g*D(4,4)*z(9);
%   addpoints(PLT.E,t,Ene);
  

        
    
end