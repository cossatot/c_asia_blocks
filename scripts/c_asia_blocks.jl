using Revise

using Oiler

using CSV
using JSON
using DataFrames:eachrow
using DataFrames, DataFramesMeta
using Setfield

using PyPlot

# options
fault_weight = 2.
save_results = true


# load data
cea_fault_file = "../block_data/c_asia_faults.geojson"
cea_block_file = "../block_data/c_asia_blocks.geojson"
cea_slip_rate_file = "../block_data/c_asia_geol_slip_rates.geojson"
tris_file = "../block_data/c_asia_sub_tris.geojson"

chn_block_file = "../../../fault_data/china/block_data/chn_blocks.geojson"
chn_fault_file = "../../../fault_data/china/block_data/chn_faults.geojson"
chn_slip_rate_file = "../../../fault_data/china/block_data/geol_slip_rate_pts.geojson"

ana_block_file = "../../../../geodesy/global_block_comps/anatolia/block_data/anatolia_blocks.geojson"
ana_fault_file = "../../../../geodesy/global_block_comps/anatolia/block_data/anatolia_faults.geojson"


gsrm_vels_file = "../gnss_data/gsrm_c_asia_vels.geojson"
comet_gnss_vels_file = "../gnss_data/c_asia_vels_rollins.geojson"
vel_field_file = "../../../fault_data/china/geod/tibet_vel_field_2021_12_06.geojson"

boundary_file = "../block_data/cea_gnss_block_domain.geojson"


@info "joining blocks"
cea_block_df = Oiler.IO.gis_vec_file_to_df(cea_block_file)
chn_block_df = Oiler.IO.gis_vec_file_to_df(chn_block_file)
ana_block_df = Oiler.IO.gis_vec_file_to_df(ana_block_file)
chn_block_df.fid = string.(chn_block_df.fid)
block_df = vcat(chn_block_df, 
                ana_block_df, 
                cea_block_df; 
                cols=:union)

@info "culling blocks"
println("n blocks before ", size(block_df, 1))
bound_df = Oiler.IO.gis_vec_file_to_df(boundary_file)
#block_df = Oiler.IO.get_blocks_in_bounds!(block_df, bound_df; epsg=2991)
println("n blocks after ", size(block_df, 1))

@info "doing faults"
fault_df, faults, fault_vels = Oiler.IO.process_faults_from_gis_files(
                                                        chn_fault_file,
                                                        ana_fault_file,
                                                        cea_fault_file;
                                                        block_df=block_df,
                                                        subset_in_bounds=true,
                                                        check_blocks=true,
                                                        )
fault_df[:,:fid] = string.(fault_df[:,:fid])
println("n faults: ", length(faults))
println("n fault vels: ", length(fault_vels))

@info "doing GNSS"
gsrm_vel_df = Oiler.IO.gis_vec_file_to_df(gsrm_vels_file)
comet_vel_df = Oiler.IO.gis_vec_file_to_df(comet_gnss_vels_file)

vel_field_df = Oiler.IO.gis_vec_file_to_df(vel_field_file)
vel_field_df[!,"station"] = string.(vel_field_df[!,:v_id])

comet_vel_df.e_err .* 2.
comet_vel_df.n_err .* 2.

@info "doing COMET GNSS vels"
@time comet_vels = Oiler.IO.make_vels_from_gnss_and_blocks(comet_vel_df, block_df;
    ve=:e_vel, vn=:n_vel, ee=:e_err, en=:n_err, name=:site,
    fix="1111"
)

@info "doing GSRM vels"
@time gsrm_vels = Oiler.IO.make_vels_from_gnss_and_blocks(gsrm_vel_df, block_df;
    ve=:e_vel, vn=:n_vel, ee=:e_err, en=:n_err, name=:station,
    fix="1111"
)

@info "doing COMET InSAR vels"
vel_field_vels = Oiler.IO.make_vels_from_gnss_and_blocks(vel_field_df, block_df;
    fix="1111")

#subsample
vel_field_vels = vel_field_vels[2:2:end]

gnss_vels = vcat(comet_vels, gsrm_vels)#, vel_field_vels)

println("n gnss vels: ", length(gnss_vels))
#println("n vel field vels: ", length(vel_field_vels))


@info "doing geol slip rates"
cea_slip_rate_df = Oiler.IO.gis_vec_file_to_df(cea_slip_rate_file)
chn_slip_rate_df = Oiler.IO.gis_vec_file_to_df(chn_slip_rate_file)
geol_slip_rate_df = vcat(cea_slip_rate_df, chn_slip_rate_df)

geol_slip_rate_df, geol_slip_rate_vels = Oiler.IO.make_geol_slip_rate_vels!(
                                                geol_slip_rate_df,
                                                fault_df)

println("n geol slip rates: ", length(geol_slip_rate_vels))

# tris
tri_json = JSON.parsefile(tris_file)
tris = Oiler.IO.tris_from_geojson(tri_json)

function set_tri_rates(tri; ds=5., de=1., ss=0., se=0.5)
    tri = @set tri.dip_slip_rate = ds
    tri = @set tri.dip_slip_err = de
    tri = @set tri.strike_slip_rate = ss
    tri = @set tri.strike_slip_err = se
    tri
end

tris = map(set_tri_rates, tris)


vels = vcat(fault_vels, gnss_vels, geol_slip_rate_vels, vel_field_vels)

vel_groups = Oiler.group_vels_by_fix_mov(vels);


# solve
results = Oiler.solve_block_invs_from_vel_groups(vel_groups; faults=faults,
                                                tris,
                                                sparse_lhs=true,
                                                weighted=true,
                                                elastic_floor=1e-2,
                                                regularize_tris=true,
                                                tri_priors=false,
                                                tri_distance_weight=20.,
                                                predict_vels=true,
                                                check_closures=true,
                                                pred_se=true,
                                                constraint_method="kkt_sym",
                                                factorization="lu")

Oiler.ResultsAnalysis.compare_data_results(results=results,
                                           vel_groups=vel_groups,
                                           geol_slip_rate_df=geol_slip_rate_df,
                                           geol_slip_rate_vels=geol_slip_rate_vels,
                                           fault_df=fault_df,
                                           )

block_bound_df = Oiler.IO.gis_vec_file_to_df("../block_data/c_asia_block_bounds.geojson")

function make_bound_fault(row)
    trace = Oiler.IO.get_coords_from_geom(row[:geometry])

    Oiler.Fault(trace=trace, dip_dir=row[:dip_dir], dip=89., hw=row[:hw],
                fw=row[:fw], fid=row[:fid])
end

#@info "getting block rates"
#bound_faults = []
#for i in 1:size(block_bound_df, 1)
#    push!(bound_faults, make_bound_fault(block_bound_df[i,:]))
#end
#
#block_bound_rates = Oiler.Utils.get_fault_slip_rates_from_poles(bound_faults,
#                                                                results["poles"],
#                                                                use_path=true)
if save_results == true
    Oiler.IO.write_tri_results_to_gj(tris, results, 
        "../results/makran_zagros_tri_results.geojson";
        name="Makran-Zagros tri results")
    
    Oiler.IO.write_fault_results_to_gj(results, 
        "../results/c_asia_fault_results.geojson",
        name="Central Asia fault results",
        calc_rake=true,
        calc_slip_rate=true)

    Oiler.IO.write_gnss_vel_results_to_csv(results, vel_groups;
        name="../results/c_asia_gnss_results.csv")

end

map_fig = Oiler.Plots.plot_results_map(results, vel_groups, faults, tris)

slip_rate_fig = Oiler.Plots.plot_slip_rate_fig(geol_slip_rate_df,
    geol_slip_rate_vels, fault_df, results)

Oiler.WebViewer.write_web_viewer(results=results, block_df=block_df,
                                 ref_pole="1111", directory="../web_viewer_insar")

show()


