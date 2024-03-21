function functions = generateFunctions(auxdata)
import casadi.*

% FORWARD DYNAMICS
e_ff_MX = MX.sym('u_MX',6);       % Controls (muscle excitations - feedforward)
X_MX = MX.sym('X_MX',auxdata.nStates);    % States (muscle activations and joint kinematics)
wM_MX = MX.sym('wM_MX',2);        % Motor noise
wPq_MX = MX.sym('wPq_MX',2);      % Sensor position noise
wPqdot_MX = MX.sym('wPqdot_MX',2);% Sensor velocity noise
EE_ref_MX = MX.sym('EE_ref_MX', 4); % Reference end-effector kinematics

% unperturbed stochastic forward dynamics
dX_MX = forwardMusculoskeletalDynamics_motorNoise(X_MX, e_ff_MX, 0, wM_MX, auxdata); 
f_forwardMusculoskeletalDynamics = Function('f_forwardMusculoskeletalDynamics',{X_MX,e_ff_MX,wM_MX},{dX_MX});  
functions.f_forwardMusculoskeletalDynamics = f_forwardMusculoskeletalDynamics;

% Jacobian of forward dynamics wrt state
DdX_DX_MX = jacobian(dX_MX, X_MX');
f_DdX_DX = Function('f_DdX_DX',{X_MX,e_ff_MX,wM_MX},{DdX_DX_MX});
functions.f_DdX_DX = f_DdX_DX;

% Jacobian of forward dynamics wrt motor noise
DdX_DwM_MX = jacobian(dX_MX, wM_MX');
f_DdX_Dw = Function('f_DdX_Dw',{X_MX,e_ff_MX,wM_MX},{DdX_DwM_MX});
functions.f_DdX_Dw = f_DdX_Dw;

% Trapezoidal integration scheme (implicit)
X_plus_MX = MX.sym('X_plus_MX',auxdata.nStates); % States at mesh end
dX_MX = MX.sym('dX_MX',auxdata.nStates); % State derivative
dX_plus_MX = MX.sym('dX_plus_MX',auxdata.nStates); % State derivative at mesh end
f_G_Trapezoidal = Function('f_G_Trapezoidal',{X_MX,X_plus_MX,dX_MX,dX_plus_MX},{G_Trapezoidal(X_MX,X_plus_MX,dX_MX,dX_plus_MX,auxdata.dt)});
functions.f_G_Trapezoidal = f_G_Trapezoidal;

% Sensitivity of trapezoidal integration scheme to changes in initial state
DdX_DX_MX = MX.sym('DdX_DX_MX',auxdata.nStates,auxdata.nStates); 
DG_DX_MX = DG_DX_Trapezoidal(DdX_DX_MX,auxdata.dt);
f_DG_DX = Function('f_DG_DX',{DdX_DX_MX},{DG_DX_MX});
functions.f_DG_DX = f_DG_DX;

% Sensitivity of trapezoidal integration scheme to changes in final state
DG_DZ_MX = DG_DZ_Trapzoidal(DdX_DX_MX,auxdata.dt);
f_DG_DZ = Function('f_DG_DZ',{DdX_DX_MX},{DG_DZ_MX});
functions.f_DG_DZ = f_DG_DZ;

% Sensitivity of trapezoidal integration scheme to changes in motor noise
DdX_Dw_MX = MX.sym('DdX_Dw_MX',auxdata.nStates,2); 
DG_DW_MX = DG_DW_Trapezoidal(DdX_Dw_MX,auxdata.dt);
f_DG_DW = Function('f_DG_DW',{DdX_Dw_MX},{DG_DW_MX});
functions.f_DG_DW = f_DG_DW;

% Some other useful function definitions
% End effector position and variability
q_MX = MX.sym('q_MX',2); P_q_MX = MX.sym('P_q_MX',2,2);
EEPos_MX = EndEffectorPos(q_MX,auxdata); f_EEPos = Function('f_EEPos',{q_MX},{EEPos_MX});
functions.f_EEPos = f_EEPos;
P_EEPos_MX = jacobian(EEPos_MX,q_MX)*P_q_MX*jacobian(EEPos_MX,q_MX)';
f_P_EEPos = Function('f_P_EEPos',{q_MX,P_q_MX},{P_EEPos_MX}); % Covariance of end effector position
functions.f_P_EEPos = f_P_EEPos;
% End effector velocity and variability
q_MX = MX.sym('q_MX',2); qdot_MX = MX.sym('q_MX',2); P_qdot_MX = MX.sym('P_qdot_MX',4,4);
EEVel_MX = EndEffectorVel(q_MX,qdot_MX,auxdata);
f_EEVel = Function('f_EEVel',{q_MX,qdot_MX},{EEVel_MX}); % End effector position
functions.f_EEVel = f_EEVel;
P_EEVel_MX = jacobian(EEVel_MX,[q_MX qdot_MX])*P_qdot_MX*jacobian(EEVel_MX,[q_MX qdot_MX])';
f_P_EEVel = Function('f_P_EEVel',{q_MX,qdot_MX,P_qdot_MX},{P_EEVel_MX}); % Covariance of end effector velocity
functions.f_P_EEVel = f_P_EEVel;
end