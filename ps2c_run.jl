#Run this file for 2c

include("ps2b_run.jl")

using Pkg
Pkg.add("LinearAlgebra")
using LinearAlgebra

#Single-Value Decomposition will give:
    #A = U * Diagonal(S) * Vt

#Phase 1 Single-Value Decomposition
F_phase1 = svd(sensitivity_phase1)
U_phase1 = F_phase1.U
S_phase1 = F_phase1.S
Vt_phase1 = F_phase1.Vt
V_phase1 = F_phase1.V

#Phase 2early Single-Value Decomposition
F_phase2early = svd(sensitivity_phase2early)
U_phase2early = F_phase2early.U
S_phase2early = F_phase2early.S
Vt_phase2early = F_phase2early.Vt
V_phase2early = F_phase2early.V

#Phase 2late Single-Value Decomposition
F_phase2late = svd(sensitivity_phase2late)
U_phase2late = F_phase2late.U
S_phase2late = F_phase2late.S
Vt_phase2late = F_phase2late.Vt
V_phase2late = F_phase2late.V
