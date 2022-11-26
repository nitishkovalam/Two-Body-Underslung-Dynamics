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