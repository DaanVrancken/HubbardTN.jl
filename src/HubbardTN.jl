module HubbardFunctions

# if occursin("amd",lowercase(Sys.cpu_info()[1].model))
#     using BLISBLAS
#     println("Using BLIS for BLAS operations.")
# else
#     using MKL
#     println("Using MKL for BLAS and LAPACK operations.")
# end

export OB_Sim, MB_Sim, OBC_Sim, MBC_Sim
export produce_groundstate, produce_excitations, produce_bandgap, produce_TruncState
export dim_state, density_spin, density_state, plot_excitations, plot_spin, extract_params, save_state, load_state

using DrWatson
using ThreadPinning
using Base.Threads
using LinearAlgebra
using MPSKit, MPSKitModels
using TensorKit
using KrylovKit
# using DataFrames
using Plots
using Plots.PlotMeasures
using TensorOperations
using JLD2

# function __init__()
#     LinearAlgebra.BLAS.set_num_threads(1)
#     if haskey(ENV, "SLURM_JOB_ID") || haskey(ENV, "JOB_ID") || haskey(ENV, "PBS_JOBID")
#         # Running on remote cluster
#         ThreadPinning.pinthreads(:affinitymask)
#     else
#         # Running locally
#         ThreadPinning.pinthreads(:cores)
#     end
#     MPSKit.Defaults.set_scheduler!(:dynamic)   # serial -> disable multithreading, greedy -> greedy load-balancing, dynamic -> moderate load-balancing
#     println("Running on $(Threads.nthreads()) threads.")
# end

include("simulations.jl")
include("hamiltonian.jl")
include("groundstate.jl")
include("excitations.jl")
include("tools.jl")
        
end