function [X, M, EE_ref, Pmat] = forwardSim_no_fb(X_init,Pmat_init,e_ff, T, auxdata,functions)
import casadi.*;
Urf = MX.sym('Urf',auxdata.nStates*2);
X_i = X_init;
X = NaN(auxdata.nStates,auxdata.N+1); X(:,1) = X_init;
EE_ref = NaN(4,auxdata.N);
EE_ref(:,1) = [EndEffectorPos(X_i(1:2),auxdata);
                EndEffectorVel(X_i(1:2),X_i(3:4),auxdata)];

% f_forwardMusculoskeletalDynamics = functions.f_forwardMusculoskeletalDynamics;
f = @(x, u) forwardMusculoskeletalDynamics_motorNoise(x, u, 0, 0, auxdata);
nStates = auxdata.nStates;
dt = T/auxdata.N;
for i = 1:auxdata.N
    dX_i = f(X_i, e_ff(:,i));
    rf = rootfinder('rf','newton',struct('x',Urf,'g', ...
        [Urf(1:nStates) - (X_i + (dX_i + Urf(nStates+1:end))/2*dt);
        Urf(nStates+1:end) - f(Urf(1:nStates), e_ff(:,i+1))]),struct('abstol',1e-16));
    solution = rf([X_i;dX_i],[]);
    X_i = full(solution(1:nStates));
    X(:,i+1) = X_i;
    EE_ref(:,i+1) = [EndEffectorPos(X_i(1:2),auxdata); EndEffectorVel(X_i(1:2),X_i(3:4),auxdata)];
end

M = NaN(auxdata.nStates,auxdata.nStates*auxdata.N);
Pmat_i = Pmat_init;
Pmat = NaN(auxdata.nStates,auxdata.nStates,auxdata.N+1);
Pmat(:,:,1) = Pmat_i;

for i = 1:auxdata.N
    DdX_DX_i = functions.f_DdX_DX(X(:,i),e_ff(:,i),auxdata.wM);
    DdZ_DX_i = functions.f_DdX_DX(X(:,i+1),e_ff(:,i+1),auxdata.wM);
    DdX_Dw_i = functions.f_DdX_Dw(X(:,i),e_ff(:,i),auxdata.wM);
    DdZ_Dw_i = functions.f_DdX_Dw(X(:,i+1),e_ff(:,i+1),auxdata.wM);
    
    DG_DX_i = functions.f_DG_DX(DdX_DX_i, dt);
    DG_DZ_i = functions.f_DG_DZ(DdZ_DX_i, dt);
    DG_DW_i = functions.f_DG_DW(DdX_Dw_i, DdZ_Dw_i, dt);
    M_i = DG_DZ_i^(-1);
    M(:,(i-1)*auxdata.nStates + 1:i*auxdata.nStates) = full(M_i);
    Pmat_i = full(M_i*(DG_DX_i*Pmat_i*DG_DX_i' + DG_DW_i*auxdata.sigma_w*DG_DW_i')*M_i');
    Pmat(:,:,i+1) = Pmat_i;
end
