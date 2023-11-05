# This file contains the module dependencies
mod_bound.mod = $(./src/bound.f90)
mod_chkdiv.mod = $(./src/chkdiv.f90)
mod_chkdt.mod = $(./src/chkdt.f90)
mod_common_cudecomp.mod = $(./src/common_cudecomp.f90)
mod_common_mpi.mod = $(./src/common_mpi.f90)
mod_const.mod = $(./src/const.f90)
mod_correc.mod = $(./src/correc.f90)
mod_debug.mod = $(./src/debug.f90)
mod_fft.mod = $(./src/fft.f90)
mod_fftw_param.mod = $(./src/fftw.f90)
mod_fillps.mod = $(./src/fillps.f90)
mod_gpu_utils.mod = $(./src/gpu_utils.f90)
procedure.mod = $(./src/gpu_utils.f90)
mod_initflow.mod = $(./src/initflow.f90)
mod_initgrid.mod = $(./src/initgrid.f90)
mod_initmpi.mod = $(./src/initmpi.f90)
mod_initsolver.mod = $(./src/initsolver.f90)
mod_load.mod = $(./src/load.f90)
mod_mom.mod = $(./src/mom.f90)
mod_nvtx.mod = $(./src/nvtx.f90)
mod_output.mod = $(./src/output.f90)
mod_param.mod = $(./src/param.f90)
mod_post.mod = $(./src/post.f90)
mod_rk.mod = $(./src/rk.f90)
mod_sanity.mod = $(./src/sanity.f90)
mod_scal.mod = $(./src/scal.f90)
mod_solver.mod = $(./src/solver.f90)
mod_solver_gpu.mod = $(./src/solver_gpu.f90)
mod_timer.mod = $(./src/timer.f90)
mod_typedef.mod = $(./src/typedef.f90)
mod_updatep.mod = $(./src/updatep.f90)
mod_utils.mod = $(./src/utils.f90)
mod_wm.mod = $(./src/wallmodel.f90)
mod_workspaces.mod = $(./src/workspaces.f90)
$(./src/bound.f90) += $(mpi.mod)
$(./src/bound.f90) += $(mod_common_mpi.mod)
$(./src/bound.f90) += $(mod_const.mod)
$(./src/bound.f90) += $(mod_typedef.mod)
$(./src/bound.f90) += $(mod_wm.mod)
$(./src/bound.f90) += $(cudecomp.mod)
$(./src/bound.f90) += $(mod_common_cudecomp.mod)
$(./src/chkdiv.f90) += $(mpi.mod)
$(./src/chkdiv.f90) += $(mod_common_mpi.mod)
$(./src/chkdiv.f90) += $(mod_const.mod)
$(./src/chkdt.f90) += $(mpi.mod)
$(./src/chkdt.f90) += $(mod_common_mpi.mod)
$(./src/chkdt.f90) += $(mod_const.mod)
$(./src/common_cudecomp.f90) += $(mod_const.mod)
$(./src/common_cudecomp.f90) += $(cudecomp.mod)
$(./src/common_cudecomp.f90) += $(openacc.mod)
$(./src/common_cudecomp.f90) += $(mod_param.mod)
$(./src/const.f90) += $(mpi.mod)
$(./src/correc.f90) += $(mod_const.mod)
$(./src/debug.f90) += $(mpi.mod)
$(./src/debug.f90) += $(mod_common_mpi.mod)
$(./src/debug.f90) += $(mod_param.mod)
$(./src/debug.f90) += $(mod_const.mod)
$(./src/fft.f90) += $(mod_common_mpi.mod)
$(./src/fft.f90) += $(mod_fftw_param.mod)
$(./src/fft.f90) += $(mod_const.mod)
$(./src/fft.f90) += $(mod_utils.mod)
$(./src/fft.f90) += $(mod_param.mod)
$(./src/fft.f90) += $(mod_common_cudecomp.mod)
$(./src/fillps.f90) += $(mod_const.mod)
$(./src/gpu_utils.f90) += $(mod_const.mod)
$(./src/initflow.f90) += $(mpi.mod)
$(./src/initflow.f90) += $(mod_common_mpi.mod)
$(./src/initflow.f90) += $(mod_param.mod)
$(./src/initflow.f90) += $(mod_const.mod)
$(./src/initgrid.f90) += $(mod_param.mod)
$(./src/initgrid.f90) += $(mod_const.mod)
$(./src/initmpi.f90) += $(mpi.mod)
$(./src/initmpi.f90) += $(decomp_2d.mod)
$(./src/initmpi.f90) += $(mod_common_mpi.mod)
$(./src/initmpi.f90) += $(mod_const.mod)
$(./src/initmpi.f90) += $(mod_common_cudecomp.mod)
$(./src/initmpi.f90) += $(mod_param.mod)
$(./src/initsolver.f90) += $(mod_fft.mod)
$(./src/initsolver.f90) += $(mod_const.mod)
$(./src/initsolver.f90) += $(mod_param.mod)
$(./src/load.f90) += $(mpi.mod)
$(./src/load.f90) += $(mod_common_mpi.mod)
$(./src/load.f90) += $(mod_const.mod)
$(./src/load.f90) += $(mod_utils.mod)
$(./src/load.f90) += $(decomp_2d.mod)
$(./src/load.f90) += $(cudecomp.mod)
$(./src/load.f90) += $(mod_common_cudecomp.mod)
$(./src/main.f90) += $(mpi.mod)
$(./src/main.f90) += $(decomp_2d.mod)
$(./src/main.f90) += $(mod_bound.mod)
$(./src/main.f90) += $(mod_chkdiv.mod)
$(./src/main.f90) += $(mod_chkdt.mod)
$(./src/main.f90) += $(mod_common_mpi.mod)
$(./src/main.f90) += $(mod_correc.mod)
$(./src/main.f90) += $(mod_fft.mod)
$(./src/main.f90) += $(mod_fillps.mod)
$(./src/main.f90) += $(mod_initflow.mod)
$(./src/main.f90) += $(mod_initgrid.mod)
$(./src/main.f90) += $(mod_initmpi.mod)
$(./src/main.f90) += $(mod_initsolver.mod)
$(./src/main.f90) += $(mod_load.mod)
$(./src/main.f90) += $(mod_mom.mod)
$(./src/main.f90) += $(mod_rk.mod)
$(./src/main.f90) += $(mod_output.mod)
$(./src/main.f90) += $(mod_param.mod)
$(./src/main.f90) += $(mod_sanity.mod)
$(./src/main.f90) += $(mod_solver.mod)
$(./src/main.f90) += $(mod_solver_gpu.mod)
$(./src/main.f90) += $(mod_workspaces.mod)
$(./src/main.f90) += $(mod_common_cudecomp.mod)
$(./src/main.f90) += $(mod_timer.mod)
$(./src/main.f90) += $(mod_updatep.mod)
$(./src/main.f90) += $(mod_utils.mod)
$(./src/main.f90) += $(mod_const.mod)
$(./src/main.f90) += $(omp_lib.mod)
$(./src/mom.f90) += $(mpi.mod)
$(./src/mom.f90) += $(mod_common_mpi.mod)
$(./src/mom.f90) += $(mod_const.mod)
$(./src/mom.f90) += $(mod_param.mod)
$(./src/output.f90) += $(mpi.mod)
$(./src/output.f90) += $(decomp_2d_io.mod)
$(./src/output.f90) += $(mod_common_mpi.mod)
$(./src/output.f90) += $(mod_const.mod)
$(./src/output.f90) += $(decomp_2d.mod)
$(./src/output.f90) += $(mod_load.mod)
$(./src/param.f90) += $(mod_const.mod)
$(./src/param.f90) += $(mpi.mod)
$(./src/post.f90) += $(mod_const.mod)
$(./src/rk.f90) += $(mod_mom.mod)
$(./src/rk.f90) += $(mod_scal.mod)
$(./src/rk.f90) += $(mod_utils.mod)
$(./src/rk.f90) += $(mod_const.mod)
$(./src/sanity.f90) += $(mpi.mod)
$(./src/sanity.f90) += $(decomp_2d.mod)
$(./src/sanity.f90) += $(mod_bound.mod)
$(./src/sanity.f90) += $(mod_chkdiv.mod)
$(./src/sanity.f90) += $(mod_common_mpi.mod)
$(./src/sanity.f90) += $(mod_correc.mod)
$(./src/sanity.f90) += $(mod_debug.mod)
$(./src/sanity.f90) += $(mod_fft.mod)
$(./src/sanity.f90) += $(mod_fillps.mod)
$(./src/sanity.f90) += $(mod_initflow.mod)
$(./src/sanity.f90) += $(mod_initmpi.mod)
$(./src/sanity.f90) += $(mod_initsolver.mod)
$(./src/sanity.f90) += $(mod_param.mod)
$(./src/sanity.f90) += $(mod_solver.mod)
$(./src/sanity.f90) += $(mod_solver_gpu.mod)
$(./src/sanity.f90) += $(mod_const.mod)
$(./src/sanity.f90) += $(mod_workspaces.mod)
$(./src/sanity.f90) += $(openacc.mod)
$(./src/scal.f90) += $(mod_const.mod)
$(./src/scal.f90) += $(mpi.mod)
$(./src/scal.f90) += $(mod_param.mod)
$(./src/solver.f90) += $(decomp_2d.mod)
$(./src/solver.f90) += $(mod_fft.mod)
$(./src/solver.f90) += $(mod_const.mod)
$(./src/solver.f90) += $(mod_param.mod)
$(./src/solver_gpu.f90) += $(cudecomp.mod)
$(./src/solver_gpu.f90) += $(mod_fft.mod)
$(./src/solver_gpu.f90) += $(mod_const.mod)
$(./src/solver_gpu.f90) += $(mod_common_mpi.mod)
$(./src/solver_gpu.f90) += $(mod_common_cudecomp.mod)
$(./src/solver_gpu.f90) += $(mod_param.mod)
$(./src/timer.f90) += $(mpi.mod)
$(./src/timer.f90) += $(mod_nvtx.mod)
$(./src/typedef.f90) += $(mod_const.mod)
$(./src/updatep.f90) += $(mod_const.mod)
$(./src/utils.f90) += $(mpi.mod)
$(./src/utils.f90) += $(mod_const.mod)
$(./src/utils.f90) += $(mod_common_cudecomp.mod)
$(./src/wallmodel.f90) += $(mod_const.mod)
$(./src/wallmodel.f90) += $(mod_typedef.mod)
$(./src/workspaces.f90) += $(mod_const.mod)
$(./src/workspaces.f90) += $(mod_common_cudecomp.mod)
$(./src/workspaces.f90) += $(mod_utils.mod)
$(./src/workspaces.f90) += $(mod_common_mpi.mod)
$(./src/workspaces.f90) += $(mod_fft.mod)
$(./src/workspaces.f90) += $(mod_param.mod)
$(./src/workspaces.f90) += $(cudecomp.mod)
$(./src/workspaces.f90) += $(openacc.mod)
$(./src/workspaces.f90) += $(cufft.mod)
$(./src/main.f90) += out1d.h90
$(./src/main.f90) += out2d.h90
$(./src/main.f90) += out3d.h90