# If you have not already done a general SL7 software setup:
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
export DUNELAR_VERSION=v10_00_04d00
export DUNELAR_QUALIFIER=e26:prof 
setup dunesw $DUNELAR_VERSION -q $DUNELAR_QUALIFIER
export METACAT_AUTH_SERVER_URL=https://metacat.fnal.gov:8143/auth/dune
export METACAT_SERVER_URL=https://metacat.fnal.gov:9443/dune_meta_prod/app 

# then you can set up metacat and rucio
setup metacat 
setup rucio

# metacat query "files from dune:all where core.file_type=detector and core.run_type=hd-protodune and core.data_tier=raw and core.data_stream=cosmics and core.runs[any]=28542" >> runs.txt
#rucio replica list file fardet-vd:prodmarley_nue_es_flat_radiological_decay0_dunevd10kt_1x8x14_3view_30deg_20250217T033222Z_gen_004122_supernova_g4stage1_g4stage2_detsim_reco.root --pfns --protocols=root

#metacat query "files from dune:all where core.file_type=mc and core.run_type=hd-protodune and core.data_tier=full-reconstructed"

#hd-protodune-det-reco:np04hd_raw_run029004_0066_dataflow3_datawriter_0_20240830T075834_reco_stage1_reco_stage2_20240830T130243_keepup.root
#hd-protodune-det-reco:np04hd_raw_run029004_0027_dataflow0_datawriter_0_20240830T074004_reco_stage1_reco_stage2_20240830T124435_keepup.root

#rucio replica list file hd-protodune-det-reco:np04hd_raw_run029004_0027_dataflow0_datawriter_0_20240830T074004_reco_stage1_reco_stage2_20240830T124435_keepup.root