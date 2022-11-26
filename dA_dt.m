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