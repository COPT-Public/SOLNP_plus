#pragma once
SOLNPWork* init_work_RESCUE
(
    SOLNPWork* w_old,
    SOLNPSettings* stgs
);
SOLNPWork* SOLNP_RESCUE_init(
    SOLNPWork* w_old,
    SOLNPSettings* stgs
    );
solnp_int update_work_rescue
(
    SOLNPWork* w,
    SOLNPSettings* stgs
);
solnp_int SOLNP_RESCUE_solve
(
    SOLNPWork* w,
    SOLNPSettings* stgs,
    SOLNPInfo* info
    );
SOLNPWork* SOLNP_RESCUE
(
    SOLNPWork* w_old,
    SOLNPSettings* stgs,
    SOLNPInfo* info // record running time
    );