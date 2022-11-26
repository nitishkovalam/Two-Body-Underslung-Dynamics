function C = invr_s_compute(H,D,A,DI)

    P = eye(12) -H*((H'*DI*H)\(H'*DI)); 
    Q = D*A;
    R = [P(:,1:3) -Q(:,4:6) P(:,7:9) -Q(:,10:12) ] ;
    S = [-Q(:,1:3) P(:,4:6) -Q(:,7:9) P(:,10:12) ] ;

    C = -R\S ;
    
    
end