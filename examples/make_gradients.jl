
# This is a demonstration of how to create gradient maps that will be used for
# parcel creation, in the event that you need just a single gradient map
# (rather than an average of many gradient maps, as would be the case for
# doing a group parcellation).
#
# See config.json in this same directory for a structure for supplying input
# and output paths.
#
# You must have Connectome Workbench in your PATH, and the dconn input is expected
# to be a nvertices x nvertices matrix of correlations, Fisher z-transformed.
#
# With 64 CPU cores, the wb_command calls take about 30 minutes each.
#
# The dconn-related operations require high RAM (~48 GB).

using JSON
using CIFTI
using StatsBase: cor
using JSON

config = open("config.json", "r") do fid
	JSON.parse(fid)
end

dconn = CIFTI.load(config["dconn"])
dconn.data[.!isfinite.(dconn.data)] .= 0
corrofcorr = atanh.(cor(dconn.data))
dconn = nothing

# fill the diagonal with a reasonable high value instead of Inf
# (using 3.8002 to replicate what MATLAB does in atanh())
for v in 1:size(corrofcorr, 1)
	corrofcorr[v, v] = 3.8002
end

outname = config["outputs"]["simliarity map"]
CIFTI.save(outname, corrofcorr; template = config["dconn"])

# done with this; free memory
corrofcorr = zeros(Float32, 2, 2)
GC.gc()

run(`wb_command -cifti-gradient
	$(config["outputs"]["similarity map"])
	ROW
	$(config["outputs"]["gradient map"])
	-left-surface $(config["surfaces"]["for gradients"]["L"])
	-right-surface $(config["surfaces"]["for gradients"]["R"])
`)

run(`wb_command -cifti-smoothing
	$(config["outputs"]["gradient map"])
	$(config["kernel"])
	$(config["kernel"])
	ROW
	$(config["outputs"]["smoothed gradients"])
	-left-surface $(config["surfaces"]["for smoothing"]["L"])
	-right-surface $(config["surfaces"]["for smoothing"]["R"])
`)

