function DG_DW = DG_DW_Trapezoidal(DdX_DwM, DdX_plus_DwM,dt)

% DG_DW =  - (DdX_DwM*dt);
DG_DW = -(dt/2) * (DdX_DwM + DdX_plus_DwM);