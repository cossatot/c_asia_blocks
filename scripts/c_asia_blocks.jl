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

chn_block_file = "../../china/block_data/chn_blocks.geojson"
chn_fault_file = "../../china/block_data/chn_faults.geojson"
chn_slip_rate_file = "../../china/block_data/geol_slip_rate_pts.geojson"
zheng_gnss_file = "../../china/geod_data/gang_zheng_gnss.geojson"

ana_block_file = "../../anatolia/block_data/anatolia_blocks.geojson"
ana_fault_file = "../../anatolia/block_data/anatolia_faults.geojson"
weiss_vel_field_file = "../../anatolia/geod_data/weiss_et_al_2020_vels_down_100.geojson"

nea_block_file = "../../ne_asia_blocks/ne_asia_blocks.geojson"
nea_fault_file = "../../ne_asia_blocks/ne_asia_faults.geojson"
nea_slip_rate_file = "../../ne_asia_blocks/ne_asia_slip_rates.geojson"

gsrm_vels_file = "../gnss_data/gsrm_c_asia_vels.geojson"
comet_gnss_vels_file = "../gnss_data/c_asia_vels_rollins.geojson"
tibet_vel_field_file = "../../china/geod_data/tibet_insar_vels_2023_04.geojson"
zagros_vel_field_file = "../gnss_data/zagros_gacos_ml1_nonan_down_100.geojson"
iran_vel_field_file = "../gnss_data/iran_insar_spring_23_full.geojson"

#boundary_file = "../block_data/cea_gnss_block_domain.geojson"
boundary_file = "../block_data/cea_hazard_boundary.geojson"


@info "joining blocks"
cea_block_df = Oiler.IO.gis_vec_file_to_df(cea_block_file)
chn_block_df = Oiler.IO.gis_vec_file_to_df(chn_block_file)
ana_block_df = Oiler.IO.gis_vec_file_to_df(ana_block_file)
nea_block_df = Oiler.IO.gis_vec_file_to_df(nea_block_file)
chn_block_df.fid = string.(chn_block_df.fid)
block_df = vcat(chn_block_df, 
                ana_block_df,
                #nea_block_df,
                cea_block_df; 
                cols=:union)

@info "culling blocks"
println("n blocks before ", size(block_df, 1))
bound_df = Oiler.IO.gis_vec_file_to_df(boundary_file)
block_df = Oiler.IO.get_blocks_in_bounds!(block_df, bound_df)
println("n blocks after ", size(block_df, 1))

@info "doing faults"
fault_df, faults, fault_vels = Oiler.IO.process_faults_from_gis_files(
                                                        chn_fault_file,
                                                        ana_fault_file,
                                                        nea_fault_file,
                                                        cea_fault_file;
                                                        block_df=block_df,
                                                        subset_in_bounds=true,
                                                        check_blocks=true,
                                                        )
fault_df[:,:fid] = string.(fault_df[:,:fid])
println("n faults: ", length(faults))
println("n fault vels: ", length(fault_vels))

@info "doing non-fault block boundaries"
@time non_fault_bounds = Oiler.IO.get_non_fault_block_bounds(block_df, faults)
bound_vels = vcat(map(b->Oiler.Boundaries.boundary_to_vels(b, ee=2.0, en=2.0), 
                      non_fault_bounds)...)
println("n non-fault-bound vels: ", length(bound_vels))

@info "doing GNSS"
gsrm_vel_df = Oiler.IO.gis_vec_file_to_df(gsrm_vels_file)
comet_vel_df = Oiler.IO.gis_vec_file_to_df(comet_gnss_vels_file)


comet_vel_df.e_err .* 2.
comet_vel_df.n_err .* 2.

@info "doing COMET GNSS vels"
@time comet_vels = Oiler.IO.make_vels_from_gnss_and_blocks(comet_vel_df, block_df;
    ve=:e_vel, vn=:n_vel, ee=:e_err, en=:n_err, name=:site,
    fix="1111"
)

@info "doing Gang Zheng GNSS vels"
gang_df = Oiler.IO.gis_vec_file_to_df(zheng_gnss_file)
@time gang_vels = Oiler.IO.make_vels_from_gnss_and_blocks(gang_df, block_df;
    name=:name,
    fix="1111"
)


@info "doing GSRM vels"
@time gsrm_vels = Oiler.IO.make_vels_from_gnss_and_blocks(gsrm_vel_df, block_df;
    ve=:e_vel, vn=:n_vel, ee=:e_err, en=:n_err, name=:station,
    fix="1111"
)

@info "doing COMET InSAR vels"

tibet_vel_field_df = Oiler.IO.gis_vec_file_to_df(tibet_vel_field_file)
tibet_vel_field_df[!,"station"] = map(x->join(["tibet_insar_", x]),
                                      string.(tibet_vel_field_df[!,:fid]))

weiss_vel_field_df = Oiler.IO.gis_vec_file_to_df(weiss_vel_field_file)
weiss_vel_field_df[!,"station"] = map(x->join(["weiss_", x]), 
                                      string.(weiss_vel_field_df[!,:fid]))

@time tibet_vel_field_vels = Oiler.IO.make_vels_from_gnss_and_blocks(
    tibet_vel_field_df, block_df; 
    name=:station, fix="1111")

#subsample
tibet_vel_field_vels = tibet_vel_field_vels[1:10:end]
@info "$(length(tibet_vel_field_vels)) Tibet vels"

@info "doing weiss vels"
@time weiss_vel_field_vels = Oiler.IO.make_vels_from_gnss_and_blocks(
    weiss_vel_field_df, block_df;
    ve=:e_vel, vn=:n_vel, ee=:e_err, en=:n_err, name=:station,
    fix="1111"
)

#@info "doing zagros insar vels"
#zagros_vel_field_df = Oiler.IO.gis_vec_file_to_df(zagros_vel_field_file)
#zagros_vel_field_df[!,"station"] = map(x->join(["zagros_insar_", x]),
#                                       string.(zagros_vel_field_df[!,:fid]))
#zagros_vel_field_df[!,"e_err"] = ones(size(zagros_vel_field_df,1)) .* 2.0
#zagros_vel_field_df[!,"n_err"] = ones(size(zagros_vel_field_df,1)) .* 2.0
#
#zagros_vel_field_vels = Oiler.IO.make_vels_from_gnss_and_blocks(
#    zagros_vel_field_df, block_df; 
#    ve=:e_vel, vn=:n_vel, ee=:e_err, en=:n_err, name=:station,
#    fix="1111")

@info "doing Iran InSAR vels"
@info "  ...reading and processing"
iran_vel_field_df = Oiler.IO.gis_vec_file_to_df(iran_vel_field_file)
iran_vel_field_df[!, "station"] = map(x->join(["iran_insar_",x]),
                                      string.(iran_vel_field_df[!,:fid]))

iran_vel_field_df = iran_vel_field_df[1:10000:end,:]

@info "  ...making vels"
@time iran_vel_field_vels = Oiler.IO.make_vels_from_gnss_and_blocks(
    iran_vel_field_df, block_df;
    name=:station, fix="1111")

@info "$(length(iran_vel_field_vels)) Iran vels"
    
gnss_vels = vcat(comet_vels, 
                 gsrm_vels,
                 tibet_vel_field_vels,
                 weiss_vel_field_vels,
                 #zagros_vel_field_vels,
                 gang_vels,
                 iran_vel_field_vels,
                 )

println("n gnss vels: ", length(gnss_vels))
#println("n vel field vels: ", length(vel_field_vels))


@info "doing geol slip rates"
cea_slip_rate_df = Oiler.IO.gis_vec_file_to_df(cea_slip_rate_file)
chn_slip_rate_df = Oiler.IO.gis_vec_file_to_df(chn_slip_rate_file)
nea_slip_rate_df = Oiler.IO.gis_vec_file_to_df(nea_slip_rate_file)
geol_slip_rate_df = vcat(cea_slip_rate_df, 
                         chn_slip_rate_df, 
                         #nea_slip_rate_df
                         )

geol_slip_rate_df, geol_slip_rate_vels = Oiler.IO.make_geol_slip_rate_vels!(
                                                geol_slip_rate_df,
                                                fault_df)

println("n geol slip rates: ", length(geol_slip_rate_vels))

# tris
tri_json = JSON.parsefile(tris_file)
tris = Oiler.IO.tris_from_geojson(tri_json)

function set_tri_rates(tri; ds=15., de=1., ss=0., se=0.5)
    tri = @set tri.dip_slip_rate = ds
    tri = @set tri.dip_slip_err = de
    tri = @set tri.strike_slip_rate = ss
    tri = @set tri.strike_slip_err = se
    tri
end

tris = map(set_tri_rates, tris)


vels = vcat(fault_vels,
            bound_vels,
            gnss_vels, 
            geol_slip_rate_vels, 
            #vel_field_vels,
            )

vel_groups = Oiler.group_vels_by_fix_mov(vels);


# solve
results = Oiler.solve_block_invs_from_vel_groups(vel_groups; faults=faults,
                                                tris,
                                                sparse_lhs=true,
                                                weighted=true,
                                                elastic_floor=1e-3,
                                                regularize_tris=true,
                                                tri_priors=false,
                                                tri_distance_weight=20.,
                                                check_nans=true,
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
                                 ref_pole="1111", directory="../web_viewer")

show()


