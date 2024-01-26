#!/bin/sh

export SIM_OUT_DIR=sim_output/cora_maps/

# Generate 21cm and foreground maps with cora
cora-makesky 21cm --nside 128 --freq 400. 500. 20 --freq-mode centre --pol full --filename ${SIM_OUT_DIR}in_map_21cm.h5
cora-makesky foreground --nside 128 --freq 400. 500. 20 --freq-mode centre --pol full --filename ${SIM_OUT_DIR}in_map_fg.h5

# Can also generate Gaussian foreground map
# cora-makesky gaussianfg --nside 128 --freq 400. 500. 20 --freq-mode centre --pol full  --filename ${SIM_OUT_DIR}in_map_gfg.h5
