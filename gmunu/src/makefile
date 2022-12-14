################# Definitions for source files #############################
INCLUDES = mod_global_parameters.t mod_usr_methods.t\
mod_physicaldata.t mod_connectivity.t mod_constants.t mod_variables.t	\
mod_basic_types.t

FOBJECTS = mod_advance.t setdt.t mod_finite_volume.t mod_tvd.t mod_finite_difference.t \
connectivity.t boundary_conditions.t comm_lib.t mod_geometry.t	        \
mod_ghostcells_update.t mod_small_values.t mod_fix_conserve.t \
mod_source.t mod_initialize.t mod_gmunu.t mod_bc_data.t mod_lookup_table.t \
mod_constrained_transport.t mod_multigrid_coupling.t 

############################################################################
#  Dependencies of INCLUDES
############################################################################

mod_global_parameters.o: mod_connectivity.mod
mod_global_parameters.o: mod_physicaldata.mod
mod_global_parameters.o: mod_constants.mod
mod_global_parameters.o: mod_variables.mod
mod_global_parameters.o: mod_basic_types.mod
mod_variables.o: mod_basic_types.mod
mod_usr_methods.o: mod_global_parameters.mod

# Add multigrid files and dependencies. Due to the m_ prefix, the
# update_dependencies.sh script ignores these dependencies.
FOBJECTS += m_octree_mg_$(NDIM)d.t
mod_multigrid_coupling.o: m_octree_mg_$(NDIM)d.mod

# Update this file by running ./update_dependencies.sh in the src folder
# Don't change the text below or the script will break
#
# Dependencies of FOBJECTS (generated by list_module_deps.sh)

amr_grid.o: mod_advance.mod
amr_grid.o: mod_amr_fct.mod
amr_grid.o: mod_coarsen_refine.mod
amr_grid.o: mod_fix_conserve.mod
amr_grid.o: mod_ghostcells_update.mod
amr_ini.o: mod_constrained_transport.mod
amr_ini.o: mod_ghostcells_update.mod
amr_ini.o: mod_multigrid_coupling.mod
amr_neighbors.o: mod_geometry.mod
amr_solution_node.o: mod_geometry.mod
boundary_conditions.o: mod_bc_data.mod
boundary_conditions.o: mod_physics.mod
coarsen.o: mod_physics.mod
connectivity.o: mod_ghostcells_update.mod
convert.o: mod_calculate_xw.mod
convert.o: mod_ghostcells_update.mod
convert.o: mod_physics.mod
forest.o: mod_space_filling_curve.mod
gmunu.o: mod_advance.mod
gmunu.o: mod_fix_conserve.mod
gmunu.o: mod_ghostcells_update.mod
gmunu.o: mod_initialize.mod
gmunu.o: mod_input_output.mod
gmunu.o: mod_multigrid_coupling.mod
gmunu.o: mod_physics.mod
gmunu.o: mod_timing.mod
gmunu.o: mod_usr.mod
mod_advance.o: mod_cfc.mod
mod_advance.o: mod_finite_difference.mod
mod_advance.o: mod_finite_volume.mod
mod_advance.o: mod_fix_conserve.mod
mod_advance.o: mod_ghostcells_update.mod
mod_advance.o: mod_physics.mod
mod_advance.o: mod_small_values.mod
mod_advance.o: mod_source.mod
mod_advance.o: mod_tvd.mod
mod_amr_fct.o: mod_physics.mod
mod_bc_data.o: mod_lookup_table.mod
mod_calculate_xw.o: mod_geometry.mod
mod_calculate_xw.o: mod_limiter.mod
mod_cfc.o: mod_cfc_alp.mod
mod_cfc.o: mod_cfc_beta.mod
mod_cfc.o: mod_cfc_parameters.mod
mod_cfc.o: mod_cfc_psi.mod
mod_cfc.o: mod_cfc_weak_field.mod
mod_cfc.o: mod_geometry.mod
mod_cfc.o: mod_ghostcells_update.mod
mod_cfc.o: mod_interpolation.mod
mod_cfc.o: mod_multigrid_coupling.mod
mod_cfc.o: mod_physics.mod
mod_cfc_alp.o: mod_cfc_parameters.mod
mod_cfc_alp.o: mod_eos.mod
mod_cfc_alp.o: mod_multigrid_coupling.mod
mod_cfc_alp.o: mod_physics.mod
mod_cfc_beta.o: mod_cfc_parameters.mod
mod_cfc_beta.o: mod_geometry.mod
mod_cfc_beta.o: mod_multigrid_coupling.mod
mod_cfc_beta.o: mod_physics.mod
mod_cfc_psi.o: mod_cfc_parameters.mod
mod_cfc_psi.o: mod_geometry.mod
mod_cfc_psi.o: mod_multigrid_coupling.mod
mod_cfc_psi.o: mod_physics.mod
mod_cfc_weak_field.o: mod_cfc_parameters.mod
mod_cfc_weak_field.o: mod_geometry.mod
mod_cfc_weak_field.o: mod_multigrid_coupling.mod
mod_cfc_weak_field.o: mod_physics.mod
mod_coarsen_refine.o: mod_amr_fct.mod
mod_coarsen_refine.o: mod_ghostcells_update.mod
mod_coarsen_refine.o: mod_load_balance.mod
mod_coarsen_refine.o: mod_multigrid_coupling.mod
mod_coarsen_refine.o: mod_physics.mod
mod_coarsen_refine.o: mod_small_values.mod
mod_coarsen_refine.o: mod_space_filling_curve.mod
mod_collapse.o: mod_calculate_xw.mod
mod_collapse.o: mod_physics.mod
mod_collapse.o: mod_slice.mod
mod_constrained_transport.o: mod_fix_conserve.mod
mod_constrained_transport.o: mod_ghostcells_update.mod
mod_constrained_transport.o: mod_limiter.mod
mod_constrained_transport.o: mod_physics.mod
mod_eos.o: mod_gmunu.mod
mod_eos_hqs.o: mod_eos.mod
mod_eos_hqs.o: mod_gmunu.mod
mod_eos_hybrid.o: mod_eos.mod
mod_eos_hybrid.o: mod_gmunu.mod
mod_eos_idealgas.o: mod_eos.mod
mod_eos_idealgas.o: mod_gmunu.mod
mod_eos_piecewise.o: mod_eos.mod
mod_eos_piecewise.o: mod_gmunu.mod
mod_eos_polytrope.o: mod_eos.mod
mod_eos_polytrope.o: mod_gmunu.mod
mod_eos_readtable.o: mod_eos.mod
mod_eos_readtable.o: mod_eos_tabulated_parameters.mod
mod_eos_tabulated.o: mod_eos.mod
mod_eos_tabulated.o: mod_eos_readtable.mod
mod_eos_tabulated.o: mod_eos_tabulated_parameters.mod
mod_eos_tabulated.o: mod_interpolation.mod
mod_eos_tabulated.o: mod_rootfinding.mod
mod_finite_difference.o: mod_finite_volume.mod
mod_finite_difference.o: mod_limiter.mod
mod_finite_difference.o: mod_physics.mod
mod_finite_difference.o: mod_source.mod
mod_finite_volume.o: mod_cfc.mod
mod_finite_volume.o: mod_eos.mod
mod_finite_volume.o: mod_limiter.mod
mod_finite_volume.o: mod_physics.mod
mod_finite_volume.o: mod_source.mod
mod_finite_volume.o: mod_tvd.mod
mod_geometry.o: mod_limiter.mod
mod_geometry.o: mod_ppm.mod
mod_ghostcells_update.o: mod_amr_fct.mod
mod_ghostcells_update.o: mod_physics.mod
mod_global_parameters.o: mod_basic_types.mod
mod_gmunu.o: mod_geometry.mod
mod_gremhd.o: mod_gmunu.mod
mod_gremhd.o: mod_gremhd_phys.mod
mod_gremhd.o: mod_gremhd_ppm.mod
mod_gremhd_phys.o: mod_cfc.mod
mod_gremhd_phys.o: mod_eos.mod
mod_gremhd_phys.o: mod_geometry.mod
mod_gremhd_phys.o: mod_gremhd_phys_add_source.mod
mod_gremhd_phys.o: mod_gremhd_phys_convert.mod
mod_gremhd_phys.o: mod_gremhd_phys_divb_ct.mod
mod_gremhd_phys.o: mod_gremhd_phys_divb_mg.mod
mod_gremhd_phys.o: mod_gremhd_phys_flux.mod
mod_gremhd_phys.o: mod_gremhd_phys_implicit_update.mod
mod_gremhd_phys.o: mod_gremhd_phys_parameters.mod
mod_gremhd_phys.o: mod_physics.mod
mod_gremhd_phys_add_source.o: mod_geometry.mod
mod_gremhd_phys_add_source.o: mod_gremhd_phys_parameters.mod
mod_gremhd_phys_add_source.o: mod_physics.mod
mod_gremhd_phys_convert.o: mod_eos.mod
mod_gremhd_phys_convert.o: mod_geometry.mod
mod_gremhd_phys_convert.o: mod_gremhd_phys_parameters.mod
mod_gremhd_phys_convert.o: mod_physics.mod
mod_gremhd_phys_convert.o: mod_rootfinding.mod
mod_gremhd_phys_convert.o: mod_small_values.mod
mod_gremhd_phys_divb_ct.o: mod_constrained_transport.mod
mod_gremhd_phys_divb_ct.o: mod_geometry.mod
mod_gremhd_phys_divb_ct.o: mod_gremhd_phys_parameters.mod
mod_gremhd_phys_divb_ct.o: mod_physics.mod
mod_gremhd_phys_divb_mg.o: mod_geometry.mod
mod_gremhd_phys_divb_mg.o: mod_ghostcells_update.mod
mod_gremhd_phys_divb_mg.o: mod_gremhd_phys_parameters.mod
mod_gremhd_phys_divb_mg.o: mod_multigrid_coupling.mod
mod_gremhd_phys_divb_mg.o: mod_physics.mod
mod_gremhd_phys_flux.o: mod_geometry.mod
mod_gremhd_phys_flux.o: mod_gremhd_phys_parameters.mod
mod_gremhd_phys_flux.o: mod_physics.mod
mod_gremhd_phys_implicit_update.o: mod_eos.mod
mod_gremhd_phys_implicit_update.o: mod_eos_idealgas.mod
mod_gremhd_phys_implicit_update.o: mod_geometry.mod
mod_gremhd_phys_implicit_update.o: mod_gremhd_phys_parameters.mod
mod_gremhd_phys_implicit_update.o: mod_physics.mod
mod_gremhd_phys_implicit_update.o: mod_rootfinding.mod
mod_gremhd_phys_parameters.o: mod_eos.mod
mod_gremhd_phys_parameters.o: mod_geometry.mod
mod_gremhd_phys_parameters.o: mod_physics.mod
mod_gremhd_phys_parameters.o: mod_small_values.mod
mod_gremhd_ppm.o: mod_eos.mod
mod_gremhd_ppm.o: mod_physics.mod
mod_gremhd_ppm.o: mod_physics_ppm.mod
mod_grhd.o: mod_gmunu.mod
mod_grhd.o: mod_grhd_phys.mod
mod_grhd.o: mod_grhd_ppm.mod
mod_grhd_phys.o: mod_cfc.mod
mod_grhd_phys.o: mod_eos.mod
mod_grhd_phys.o: mod_geometry.mod
mod_grhd_phys.o: mod_grhd_phys_add_source.mod
mod_grhd_phys.o: mod_grhd_phys_convert.mod
mod_grhd_phys.o: mod_grhd_phys_flux.mod
mod_grhd_phys.o: mod_grhd_phys_one_d_core.mod
mod_grhd_phys.o: mod_grhd_phys_parameters.mod
mod_grhd_phys.o: mod_physics.mod
mod_grhd_phys_add_source.o: mod_geometry.mod
mod_grhd_phys_add_source.o: mod_grhd_phys_parameters.mod
mod_grhd_phys_add_source.o: mod_physics.mod
mod_grhd_phys_convert.o: mod_eos.mod
mod_grhd_phys_convert.o: mod_geometry.mod
mod_grhd_phys_convert.o: mod_grhd_phys_parameters.mod
mod_grhd_phys_convert.o: mod_physics.mod
mod_grhd_phys_convert.o: mod_rootfinding.mod
mod_grhd_phys_convert.o: mod_small_values.mod
mod_grhd_phys_flux.o: mod_geometry.mod
mod_grhd_phys_flux.o: mod_grhd_phys_parameters.mod
mod_grhd_phys_flux.o: mod_physics.mod
mod_grhd_phys_one_d_core.o: mod_grhd_phys_parameters.mod
mod_grhd_phys_one_d_core.o: mod_physics.mod
mod_grhd_phys_parameters.o: mod_eos.mod
mod_grhd_phys_parameters.o: mod_geometry.mod
mod_grhd_phys_parameters.o: mod_physics.mod
mod_grhd_phys_parameters.o: mod_small_values.mod
mod_grhd_ppm.o: mod_eos.mod
mod_grhd_ppm.o: mod_physics.mod
mod_grhd_ppm.o: mod_physics_ppm.mod
mod_grmhd.o: mod_gmunu.mod
mod_grmhd.o: mod_grmhd_phys.mod
mod_grmhd.o: mod_grmhd_ppm.mod
mod_grmhd_phys.o: mod_cfc.mod
mod_grmhd_phys.o: mod_eos.mod
mod_grmhd_phys.o: mod_geometry.mod
mod_grmhd_phys.o: mod_grmhd_phys_add_source.mod
mod_grmhd_phys.o: mod_grmhd_phys_convert.mod
mod_grmhd_phys.o: mod_grmhd_phys_divb_ct.mod
mod_grmhd_phys.o: mod_grmhd_phys_divb_mg.mod
mod_grmhd_phys.o: mod_grmhd_phys_flux.mod
mod_grmhd_phys.o: mod_grmhd_phys_parameters.mod
mod_grmhd_phys.o: mod_limiter.mod
mod_grmhd_phys.o: mod_physics.mod
mod_grmhd_phys_add_source.o: mod_geometry.mod
mod_grmhd_phys_add_source.o: mod_grmhd_phys_parameters.mod
mod_grmhd_phys_add_source.o: mod_physics.mod
mod_grmhd_phys_convert.o: mod_eos.mod
mod_grmhd_phys_convert.o: mod_geometry.mod
mod_grmhd_phys_convert.o: mod_grmhd_phys_parameters.mod
mod_grmhd_phys_convert.o: mod_physics.mod
mod_grmhd_phys_convert.o: mod_rootfinding.mod
mod_grmhd_phys_convert.o: mod_small_values.mod
mod_grmhd_phys_divb_ct.o: mod_constrained_transport.mod
mod_grmhd_phys_divb_ct.o: mod_geometry.mod
mod_grmhd_phys_divb_ct.o: mod_grmhd_phys_parameters.mod
mod_grmhd_phys_divb_ct.o: mod_physics.mod
mod_grmhd_phys_divb_mg.o: mod_geometry.mod
mod_grmhd_phys_divb_mg.o: mod_ghostcells_update.mod
mod_grmhd_phys_divb_mg.o: mod_grmhd_phys_parameters.mod
mod_grmhd_phys_divb_mg.o: mod_multigrid_coupling.mod
mod_grmhd_phys_divb_mg.o: mod_physics.mod
mod_grmhd_phys_flux.o: mod_grmhd_phys_parameters.mod
mod_grmhd_phys_flux.o: mod_physics.mod
mod_grmhd_phys_parameters.o: mod_eos.mod
mod_grmhd_phys_parameters.o: mod_geometry.mod
mod_grmhd_phys_parameters.o: mod_physics.mod
mod_grmhd_phys_parameters.o: mod_small_values.mod
mod_grmhd_ppm.o: mod_eos.mod
mod_grmhd_ppm.o: mod_physics.mod
mod_grmhd_ppm.o: mod_physics_ppm.mod
mod_grrhd.o: mod_gmunu.mod
mod_grrhd.o: mod_grrhd_phys.mod
mod_grrhd.o: mod_grrhd_ppm.mod
mod_grrhd_phys.o: mod_eos.mod
mod_grrhd_phys.o: mod_geometry.mod
mod_grrhd_phys.o: mod_grrhd_phys_add_source.mod
mod_grrhd_phys.o: mod_grrhd_phys_convert.mod
mod_grrhd_phys.o: mod_grrhd_phys_flux.mod
mod_grrhd_phys.o: mod_grrhd_phys_parameters.mod
mod_grrhd_phys.o: mod_physics.mod
mod_grrhd_phys_add_source.o: mod_geometry.mod
mod_grrhd_phys_add_source.o: mod_grrhd_phys_parameters.mod
mod_grrhd_phys_add_source.o: mod_physics.mod
mod_grrhd_phys_convert.o: mod_eos.mod
mod_grrhd_phys_convert.o: mod_geometry.mod
mod_grrhd_phys_convert.o: mod_grrhd_phys_parameters.mod
mod_grrhd_phys_convert.o: mod_physics.mod
mod_grrhd_phys_convert.o: mod_rootfinding.mod
mod_grrhd_phys_convert.o: mod_small_values.mod
mod_grrhd_phys_flux.o: mod_geometry.mod
mod_grrhd_phys_flux.o: mod_grrhd_phys_parameters.mod
mod_grrhd_phys_flux.o: mod_physics.mod
mod_grrhd_phys_parameters.o: mod_eos.mod
mod_grrhd_phys_parameters.o: mod_geometry.mod
mod_grrhd_phys_parameters.o: mod_physics.mod
mod_grrhd_phys_parameters.o: mod_rootfinding.mod
mod_grrhd_phys_parameters.o: mod_small_values.mod
mod_grrhd_ppm.o: mod_eos.mod
mod_grrhd_ppm.o: mod_physics.mod
mod_grrhd_ppm.o: mod_physics_ppm.mod
mod_initialize.o: mod_amr_fct.mod
mod_initialize.o: mod_bc_data.mod
mod_initialize.o: mod_fix_conserve.mod
mod_initialize.o: mod_geometry.mod
mod_initialize.o: mod_ghostcells_update.mod
mod_initialize.o: mod_input_output.mod
mod_initialize.o: mod_physics.mod
mod_input_output.o: mod_collapse.mod
mod_input_output.o: mod_geometry.mod
mod_input_output.o: mod_limiter.mod
mod_input_output.o: mod_physics.mod
mod_input_output.o: mod_slice.mod
mod_input_output.o: mod_small_values.mod
mod_input_output.o: mod_source.mod
mod_input_output.o: mod_timing.mod
mod_limiter.o: mod_mp5.mod
mod_limiter.o: mod_ppm.mod
mod_limiter.o: mod_venk.mod
mod_limiter.o: mod_weno.mod
mod_load_balance.o: mod_space_filling_curve.mod
mod_mp5.o: mod_physics.mod
mod_multigrid_coupling.o: mod_geometry.mod
mod_physics.o: mod_physics_hllc.mod
mod_physics.o: mod_physics_ppm.mod
mod_physics.o: mod_physics_roe.mod
mod_rootfinding.o: mod_lnsrch.mod
mod_rootfinding.o: mod_lu.mod
mod_slice.o: mod_basic_types.mod
mod_slice.o: mod_calculate_xw.mod
mod_source.o: mod_ghostcells_update.mod
mod_source.o: mod_physics.mod
mod_tvd.o: mod_limiter.mod
mod_tvd.o: mod_physics_roe.mod
mod_usr_template.o: mod_mhd.mod
mod_variables.o: mod_basic_types.mod
refine.o: mod_amr_fct.mod
refine.o: mod_physics.mod
setdt.o: mod_physics.mod
