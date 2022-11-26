%%GIVEN THE HELICOPTER TRAJECTORY, ACCELERATIONS OF THE LOAD AND CABLE
%%FORCES OF THE SYSTEM CAN BE OBTAINED
% Helicopter with constant forward velocity is taken
% Input parameters required 
% r0;v0;a0 of the helicopter
% D (Mass of helicopter and load)
% Fa (Aerodynamic Data)

% Transformation Matrics 
% Ra ,L -> cable attachment pt w.r.t body, cable length 
% 
% Euler angle triplet order = (phi,theta,psi)

%% Inputs


global n;    % no. of bodies (including helicopter)

n=2;     % no of bodies

L1=15;  % Helicopter to middle attachement point
L2=4;  % Middle attachement point to ris

g=9.81;
m1 = 9000;
m2 = 2.4742e+03;

J1=0.1*[36100   0   14800;
    0   191500   0;
    14800   0   179200];        % Inertia matrix of helicopter -> does not matter for straight flight

J2=[2572.71   53.15  161.08;
      53.15   7811.00  -55.06;
      161.08   -55.06   6596.88]; % Inertia matrix of RLV

% Ra=[0,0;
%     0,0;
%     1*5.7/2,-L2-1.037];    % Position of attachment points in respective body frames

Ra=[0,0;
    0,0;
    1*5.7/2,-1.037];    % Position of attachment points in respective body frames



%load heli_vel_30.mat;     % Input Helicopter Trajectory
%load heli_traj_with_turn.mat
load sortie_1.mat
load RLV_aero_coeff_CFD.mat;   % Input RLV aero coefficients

% Fa.s_ref=6;
% Fa.b_ref=3.6;
% Fa.rho=0.81;
% Fa.Cnr=-0.41;
% Fa.Clp=-0.104;
% Fa.l_ref=1.84;  

%vel=30;        % Initial forward velocity of RLV
vel = 0;
ph_c=-0*pi/180; % Euler angles of cables 
ph_l=0*pi/180;

r0=[Ra(1,1)+L1*sin(ph_c)+L2*sin(ph_l),0,Ra(3,1)+L1*cos(ph_c)+L2*cos(ph_l),0,ph_l,0];
v0=[vel,0,0,0,0,0];   

% r0 = [r_RLV_x,n , r_RLV_y,n , r_RLV_z,n , phi , theta , psi]
% v0 = [v_RLV_x,n , v_RLV_y,n , v_RLV_z,n , omga_RLV_x,rlv , omga_RLV_y,rlv , omga_RLV_z,rlv]
% ,n -> inertial frame
% ,RLV -> RLB body frame
% Euler angles are between body frame and inertial frame

%tspan=[0 F.r{1,1}.GridVectors{1}(end)];
tspan=[0 0.5];   % Time of simulation
vw = [0; 0 ;0];  % Wind velocity


% ~~~~~  Inputs ended ~~~~~~
% Everything else remains same except what needs to be plotted
% ^ Add RIS aero data through force_aero_load or directly to aero data of RLV
% ^^ currently RIS data has been added separately in force_aero_load
% %% Intermediate plotting for checking solution while solving
% global PLT;
% 
% % uncomment required lines in slung_body_eqns_load to display the plots
% 
% figure(5);
% title('\alpha \beta');
% % xlim([0,100])
% PLT.a=animatedline('LineStyle','-.','LineWidth',3);
% PLT.b=animatedline('LineStyle','-','LineWidth',3);
% grid on;
% %drawnow
% 
% figure(8);
% title('\theta_c \phi_c');
% PLT.c=animatedline('LineStyle','-.','LineWidth',3);
% PLT.ph=animatedline('LineStyle','-','LineWidth',3);
% grid on;
% % figure(9);
% % title('aero forces');
% % PLT.cs=animatedline('LineStyle','-.','LineWidth',3);
% % PLT.cr=animatedline('LineStyle','--','LineWidth',3);
% % PLT.cy=animatedline('LineStyle','-','LineWidth',3);

%% Solving

D = [ m1*eye(3)  zeros(3)    zeros(3)    zeros(3);
      zeros(3)   m2*eye(3)   zeros(3)    zeros(3);
      zeros(3)   zeros(3)    J1          zeros(3);
      zeros(3)   zeros(3)    zeros(3)    J2];     % Inertia matrix 
 
DI=inv(D);


z0=[v0,r0];     % Solution variable in configurational coordinates

 options = odeset('RelTol',1e-4,'AbsTol',1e-4,'MaxStep',1e-2);
%options = odeset('RelTol',1e-5,'AbsTol',1e-7);

sol=ode45(@(t,z)slung_body_eqns_load(t,z,Ra,L1,D,DI,g,F,Fa,vw), tspan, z0,options); % Ode45 

t=sol.x;
[z,dzdt]=deval(sol,t);   
% Returns solution variable and its derivative
% z = [ v r ] => dzdt = [a v] 

z=z';                    
dzdt=dzdt';
%% Generalized Co-ordinates
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
%% System Dynamics
function Adot = dA_dt(u,phc,L,T,Ra)
    phcdot=-u(5)/L;
    thcdot=u(4)/(L*cos(phc));
    
    Wc=[phcdot;
        thcdot*cos(phc);
        -thcdot*sin(phc)];
    
    Adot_22=T.NB{3}*vp(Wc);
    Adot_23=-T.NB{1}*vp(u(7:9))*vp(Ra(:,1));
    Adot_24=T.NB{2}*vp(u(10:12))*vp(Ra(:,2));
    
    Adot=[zeros(3,3*4);
         zeros(3)   Adot_22     Adot_23     Adot_24;
         zeros(3*2,3*4)];

end
%% Aerodynamics Subroutine
function Fa = force_aero_load (z,v,T,Fa,vw,t)


 global PLT;
   
    % Body experiences two "winds" (as seen by the body)
    % 1. due to its motion in a stationary atmosphere
    % 2. due to crosswinds
    % So the net wind in the body frame is vw_b
    
    
    va_b=T.BN{2}*(v(4:6));
    vw_b=va_b-0*([0;1;0]*(-5*1)*(heaviside(t+0.0001)-heaviside(t-1)));
    % ^ this is the wind experienced by the body
    % ^^ direction of this gives the x-direction of the wind-axis in the body frame 
    v_mag=sqrt(vw_b'*vw_b);
    

    % alp  -> angle of attack
    % beta -> angle of sideslip
    alp=atan2(vw_b(3),vw_b(1));
    beta=asin(vw_b(2)/sqrt(sum(vw_b.^2)));


    addpoints(PLT.a,t,alp*180/pi);
    addpoints(PLT.b,t,beta*180/pi);     
%     drawnow limitrate
    
% CA and CN are negated due to the difference in axis definition between aero data and this case

    force_b=0.5*Fa.rho*v_mag^2*Fa.s_ref*...
           [ (-Fa.ca(alp*180/pi,beta*180/pi))+(Fa.CLq*z(5)*Fa.l_ref/(2*v_mag)*sin(alp));
              Fa.cs(alp*180/pi,beta*180/pi);
             -Fa.cn(alp*180/pi,beta*180/pi)-(Fa.CLq*z(5)*Fa.l_ref/(2*v_mag)*cos(alp))
            ];
        
    force_iner = T.NB{2}*force_b;
    % Convert the forces into the inertial frame from body frame
    
    % Add the RIS component to the yawing moment co-efficient
    mom_b=0.5*Fa.rho*v_mag^2*Fa.s_ref* ...
           [ (Fa.crm(alp*180/pi,beta*180/pi)+1*Fa.Clp*z(4)*Fa.b_ref/(2*v_mag))*Fa.b_ref;
             (Fa.cpm(alp*180/pi,beta*180/pi)+1*Fa.Cmq*z(5)*Fa.l_ref/(2*v_mag))*Fa.l_ref;
             (Fa.cym(alp*180/pi,beta*180/pi)+1*Fa.Cnr*z(6)*Fa.b_ref/(2*v_mag) )*Fa.b_ref  
            ];
    

    
    Fa = [ force_iner ;
           mom_b
          ];

end
%% Inertial Forces and Moments Subroutine
function ficg = force_icg_load(g,D,Adot,u)

        % Force and Moments due to Gravity
        fg=[0;
            0;
            D(4,4)*g;
            zeros(3,1)];
        
        % Force and Moments due to Coriolis component
        fcori=-[zeros(3,1);
                vp(u(10:12))*D(10:12,10:12)*u(10:12)];
            
        fi=-D*Adot*u;
            
        ficg=fg+fi([4:6 10:12])+fcori;
            
end
%% Generalized Co-ordinates Transformation
function C = invr_s_compute(H,D,A,DI)

    P = eye(12) -H*((H'*DI*H)\(H'*DI)); 
    Q = D*A;
    R = [P(:,1:3) -Q(:,4:6) P(:,7:9) -Q(:,10:12) ] ;
    S = [-Q(:,1:3) P(:,4:6) -Q(:,7:9) P(:,10:12) ] ;

    C = -R\S ;
    
    
end
%% Differential Equation Subroutine
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
%% Reference Frames Subroutine
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
%% Cross Product of Angular Velocity
function S = vp(V)
    S=[0       -V(3)    V(2);
       V(3)     0      -V(1);
      -V(2)     V(1)    0];
end
