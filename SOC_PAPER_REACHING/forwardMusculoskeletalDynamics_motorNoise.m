function dX = forwardMusculoskeletalDynamics_motorNoise(X,u,T_EXT,wM,auxdata)
q = X(1:2);
qdot = X(3:4);

[Fa,Fp,~,~,~,~,~,~,~] = getMuscleForce(q,qdot,auxdata);

Fm = u.*Fa + Fp;
T = TorqueForceRelation(Fm,q,auxdata) + wM;

F_forceField = auxdata.forceField*(auxdata.l1*cos(q(1,:)) + auxdata.l2*cos(q(1,:)+q(2,:)));
T_forceField = -F_forceField*[auxdata.l2*sin(q(1,:)+q(2,:))+auxdata.l1*sin(q(1,:));auxdata.l2*sin(q(1,:)+q(2,:))];
              

ddtheta = armForwardDynamics(T,q(2),qdot,T_EXT+T_forceField,auxdata);



dX =  [qdot; ddtheta];