function functions = generateFunctions(auxdata)
import casadi.*

nDOF = auxdata.nDOF;
nMuscles = auxdata.nMuscles;
nMotorNoises = auxdata.nMotorNoises;
nStates = auxdata.nStates;
% dt = auxdata.dt;

% FORWARD DYNAMICS
u_MX = MX.sym('u_MX',nMuscles); % Controls (muscle excitations)
X_MX = MX.sym('X_MX',nStates); % States (joint kinematics)
dt_MX = MX.sym('dt_MX',1); % Time step
wM_MX = MX.sym('wM_MX',nMotorNoises); % Motor noise

% unperturbed stochastic forward dynamics
dX_MX = forwardMusculoskeletalDynamics_motorNoise(X_MX, u_MX, 0, wM_MX, auxdata); 
f_forwardMusculoskeletalDynamics = Function('f_forwardMusculoskeletalDynamics',{X_MX,u_MX,wM_MX},{dX_MX});  
functions.f_forwardMusculoskeletalDynamics = f_forwardMusculoskeletalDynamics;

% Jacobian of forward dynamics wrt state
DdX_DX_MX = jacobian(dX_MX, X_MX');
f_DdX_DX = Function('f_DdX_DX',{X_MX,u_MX,wM_MX},{DdX_DX_MX});
functions.f_DdX_DX = f_DdX_DX;

% Jacobian of forward dynamics wrt motor noise
DdX_DwM_MX = jacobian(dX_MX, wM_MX');
f_DdX_Dw = Function('f_DdX_Dw',{X_MX,u_MX,wM_MX},{DdX_DwM_MX});
functions.f_DdX_Dw = f_DdX_Dw;

% Trapezoidal integration scheme (implicit)
X_plus_MX = MX.sym('X_plus_MX',nStates); % States at mesh end
dX_MX = MX.sym('dX_MX',nStates); % State derivative
dX_plus_MX = MX.sym('dX_plus_MX',nStates); % State derivative at mesh end
f_G_Trapezoidal = Function('f_G_Trapezoidal',{X_MX,X_plus_MX,dX_MX,dX_plus_MX, dt_MX},{G_Trapezoidal(X_MX,X_plus_MX,dX_MX,dX_plus_MX,dt_MX)});
functions.f_G_Trapezoidal = f_G_Trapezoidal;

% Sensitivity of trapezoidal integration scheme to changes in initial state
DdX_DX_MX = MX.sym('DdX_DX_MX',nStates,nStates); 
DG_DX_MX = DG_DX_Trapezoidal(DdX_DX_MX,dt_MX);
f_DG_DX = Function('f_DG_DX',{DdX_DX_MX, dt_MX},{DG_DX_MX});
functions.f_DG_DX = f_DG_DX;

% Sensitivity of trapezoidal integration scheme to changes in final state
DG_DZ_MX = DG_DZ_Trapzoidal(DdX_DX_MX,dt_MX);
f_DG_DZ = Function('f_DG_DZ',{DdX_DX_MX, dt_MX},{DG_DZ_MX});
functions.f_DG_DZ = f_DG_DZ;

% % Sensitivity of trapezoidal integration scheme to changes in motor noise
% DdX_Dw_MX = MX.sym('DdX_Dw_MX',nStates,2); 
% DG_DW_MX = DG_DW_Trapezoidal(DdX_Dw_MX,dt);
% f_DG_DW = Function('f_DG_DW',{DdX_Dw_MX},{DG_DW_MX});
% functions.f_DG_DW = f_DG_DW;
% Sensitivity of trapezoidal integration scheme to changes in motor noise
DdX_Dw_MX = MX.sym('DdX_Dw_MX',nStates, nMotorNoises); 
DdX_plus_Dw_MX= MX.sym('DdX_plus_Dw_MX',nStates, nMotorNoises); 
DG_DW_MX = DG_DW_Trapezoidal(DdX_Dw_MX, DdX_plus_Dw_MX, dt_MX);
f_DG_DW = Function('f_DG_DW',{DdX_Dw_MX, DdX_plus_Dw_MX, dt_MX},{DG_DW_MX});
functions.f_DG_DW = f_DG_DW;

% Some other useful function definitions
% End effector position and variability
q_MX = MX.sym('q_MX',nDOF); 
P_q_MX = MX.sym('P_q_MX',nDOF,nDOF);
EEPos_MX = EndEffectorPos(q_MX,auxdata); f_EEPos = Function('f_EEPos',{q_MX},{EEPos_MX});
functions.f_EEPos = f_EEPos;
P_EEPos_MX = jacobian(EEPos_MX,q_MX)*P_q_MX*jacobian(EEPos_MX,q_MX)';
f_P_EEPos = Function('f_P_EEPos',{q_MX,P_q_MX},{P_EEPos_MX}); % Covariance of end effector position
functions.f_P_EEPos = f_P_EEPos;
% End effector velocity and variability
q_MX = MX.sym('q_MX',nDOF); qdot_MX = MX.sym('q_MX',nDOF); P_qdot_MX = MX.sym('P_qdot_MX',2*nDOF,2*nDOF);
EEVel_MX = EndEffectorVel(q_MX,qdot_MX,auxdata);
f_EEVel = Function('f_EEVel',{q_MX,qdot_MX},{EEVel_MX}); % End effector position
functions.f_EEVel = f_EEVel;
P_EEVel_MX = jacobian(EEVel_MX,[q_MX qdot_MX])*P_qdot_MX*jacobian(EEVel_MX,[q_MX qdot_MX])';
f_P_EEVel = Function('f_P_EEVel',{q_MX,qdot_MX,P_qdot_MX},{P_EEVel_MX}); % Covariance of end effector velocity
functions.f_P_EEVel = f_P_EEVel;
end