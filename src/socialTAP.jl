using NPZ
using JuMP
using LinearAlgebra


function read_np_data()
	A = npzread("tmp_jl/A.npz")["arr_0"]
	B = npzread("tmp_jl/B.npz")["arr_0"]
	C = npzread("tmp_jl/C.npz")["arr_0"]
	h = npzread("tmp_jl/h.npz")["arr_0"]
	flow_k = npzread("tmp_jl/flow_k.npz")["arr_0"]
	data = npzread("tmp_jl/data.npz")["arr_0"]
	dxdb_k = npzread("tmp_jl/dxdb.npz")["arr_0"]
	dxdg_k = npzread("tmp_jl/dxdg.npz")["arr_0"]
	g_k = npzread("tmp_jl/gk.npz")["arr_0"]
	beta_k = npzread("tmp_jl/betak.npz")["arr_0"]
	pars = npzread("tmp_jl/par.npz")["arr_0"]
	return A, B, C, h, flow_k, data, dxdg_k, dxdb_k, g_k, beta_k, pars
end

