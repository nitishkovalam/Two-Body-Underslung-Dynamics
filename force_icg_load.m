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