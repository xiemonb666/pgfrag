# TODO
function mpi_data_handler(ms::AbstractMultibodySetup, solver::AbstractTimeSolver)
    _reject_initial_policy_in_multibody(ms, "MPI")
    throw(ArgumentError("MPI multibody data handling is not implemented in this pinned fork"))
end
