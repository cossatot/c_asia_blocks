using DataFrames:eachrow
using Revise

using Oiler

using CSV
using JSON
using DataFrames, DataFramesMeta
using ArchGDAL
using Setfield
const AG = ArchGDAL

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

gsrm_vels_file = "../gnss_data/gsrm_c_asia_vels.geojson"
comet_gnss_vels_file = "../gnss_data/comet_c_asia_gnss_vels.geojson"
vel_field_file = "../../../fault_data/china/geod/tibet_vel_field.geojson"

boundary_file = "../block_data/cea_gnss_block_domain.geojson"


cea_fault_df = Oiler.IO.gis_vec_file_to_df(cea_fault_file)
cea_block_df = Oiler.IO.gis_vec_file_to_df(cea_block_file)
cea_slip_rate_df = Oiler.IO.gis_vec_file_to_df(cea_slip_rate_file)

chn_fault_df = Oiler.IO.gis_vec_file_to_df(chn_fault_file)
chn_block_df = Oiler.IO.gis_vec_file_to_df(chn_block_file)
chn_slip_rate_df = Oiler.IO.gis_vec_file_to_df(chn_slip_rate_file)

gsrm_vel_df = Oiler.IO.gis_vec_file_to_df(gsrm_vels_file)
comet_vel_df = Oiler.IO.gis_vec_file_to_df(comet_gnss_vels_file)

vel_field_df = Oiler.IO.gis_vec_file_to_df(vel_field_file)
vel_field_df[!,"station"] = string.(vel_field_df[!,:id])
bound_df = Oiler.IO.gis_vec_file_to_df(boundary_file)


## filter data by location

# filter CHN blocks and faults in model domain
function get_blocks_in_bounds(block_df, bdf)
    bound = bdf[1,:geometry]
    bgs = block_df[:,:geometry]
    idxs = falses(length(bgs))
    for (i, block_geom) in enumerate(bgs)
        if AG.intersects(bound, block_geom)
            idxs[i] = true
        end
    end
    block_df[idxs, :]
end

chn_block_df = get_blocks_in_bounds(chn_block_df, bound_df)
chn_fault_df = get_blocks_in_bounds(chn_fault_df, bound_df)

cea_block_df.color = ["" for i in 1:size(cea_block_df, 1)]

# join blocks and faults
chn_block_df.fid = string.(chn_block_df.fid)
chn_fault_df.fid = string.(chn_fault_df.fid)
block_df = vcat(chn_block_df, cea_block_df)
println("n blocks: ", size(block_df, 1))

chn_faults = map(Oiler.IO.row_to_fault, eachrow(chn_fault_df))
cea_faults = map(Oiler.IO.row_to_fault, eachrow(cea_fault_df))
faults = vcat(chn_faults, cea_faults)
fault_vels = reduce(vcat, map(Oiler.IO.fault_to_vels, faults))

println("n faults: ", length(faults))
println("n faults vels: ", length(fault_vels))

# get GNSS indices

comet_vel_df.sig_east .* 2.
comet_vel_df.sig_north .* 2.

comet_vels = Oiler.IO.make_vels_from_gnss_and_blocks(comet_vel_df, block_df;
    ve=:v_east, vn=:v_north, ee=:sig_east, en=:sig_north, name=:name,
    fix="1111"
)
gsrm_vels = Oiler.IO.make_vels_from_gnss_and_blocks(gsrm_vel_df, block_df;
    ve=:e_vel, vn=:n_vel, ee=:e_err, en=:n_err, name=:station,
    fix="1111"
)

vel_field_vels = Oiler.IO.make_vels_from_gnss_and_blocks(vel_field_df, block_df;
    fix="1111")
gnss_vels = vcat(comet_vels, gsrm_vels)

println("n gnss vels: ", length(gnss_vels))
println("n vel field vels: ", length(vel_field_vels))


# geol slip rates
chn_slip_rate_df = filter(x -> x.fault_seg in chn_fault_df.fid, 
    chn_slip_rate_df)

chn_slip_rates = Oiler.IO.make_geol_slip_rate_vel_vec(chn_slip_rate_df, 
   chn_fault_df; err_return_val=5.)
cea_slip_rates = Oiler.IO.make_geol_slip_rate_vel_vec(cea_slip_rate_df, 
   cea_fault_df; err_return_val=5.)
geol_slip_rate_vels = vcat(chn_slip_rates, cea_slip_rates)

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
                                                tris=tris,
                                                sparse_lhs=false,
                                                weighted=true,
                                                regularize_tris=true,
                                                tri_priors=true,
                                                tri_distance_weight=50.,
                                                predict_vels=true,
                                                check_closures=true,
                                                pred_se=true,
                                                factorization="lu")


if save_results == true
    Oiler.IO.write_tri_results_to_gj(tris, results, 
        "../results/makran_zagros_tri_results.geojson";
        name="Makran-Zagros tri results")
    
    Oiler.IO.write_fault_results_to_gj(results, 
        "../results/c_asia_fault_results.geojson",
        name="Central Asia fault results")
end

# plot

obs_vel_df = Oiler.Utils.make_df_from_vel_array(gnss_vels)
pred_vel_df = Oiler.Utils.get_gnss_results(results, vel_groups)

if save_results == true
    CSV.write("../results/pred_vels.csv", pred_vel_df)
end

figure(figsize=(14, 14))

cm = get_cmap(:viridis)

function get_tri_total_rate(tri)
    ds = results["tri_slip_rates"][tri.name]["dip_slip"]
    ss = results["tri_slip_rates"][tri.name]["strike_slip"]
    total_rate = sqrt(ds^2 + ss^2)
end

tri_rates = [get_tri_total_rate(tri) for tri in tris]
tri_rate_min = minimum(tri_rates)
tri_rate_max = maximum(tri_rates)

function plot_tri(tri; vmin=tri_rate_min, vmax=tri_rate_max)
    lons = [tri.p1[1], tri.p2[1], tri.p3[1], tri.p1[1]]
    lats = [tri.p1[2], tri.p2[2], tri.p3[2], tri.p1[2]]
    total_rate = get_tri_total_rate(tri)
    rate_frac = (total_rate - vmin) / (vmax - vmin)
    color = cm(rate_frac)
    # color = [0., 0., rate_frac]
    fill(lons, lats, color=color, alpha=0.25, zorder=0)
end

for tri in tris
    plot_tri(tri)
end


for fault in faults
    plot(fault.trace[:,1], fault.trace[:,2], "k-", lw=0.3)
end

quiver(obs_vel_df.lon, obs_vel_df.lat,
       obs_vel_df.ve, obs_vel_df.vn, color="b", scale=300)
quiver(pred_vel_df.lon, pred_vel_df.lat,
       pred_vel_df.ve, pred_vel_df.vn, color="r", scale=300)
axis("equal")
show()

