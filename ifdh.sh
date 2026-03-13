source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
kx509

export ROLE=Analysis
voms-proxy-init -rfc -noregen -voms=dune:/dune/Role=$ROLE -valid 120:00
setup ifdhc
export IFDH_TOKEN_ENABLE=0
#ifdh cp root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/tape_backed/dunepro//hd-protodune/full-reconstructed/2024/detector/cosmics/hd-protodune-reco-keepup-v0/00/02/82/19/np04hd_raw_run028219_0000_dataflow3_datawriter_0_20240729T093114_reco_stage1_reco_stage2_20240729T113114_keepup.root /exp/dune/data/users/gabrielb

ifdh  cp root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/tape_backed/dunepro//hd-protodune/raw/2024/detector/cosmics/None/00/02/82/19/np04hd_raw_run028219_0000_dataflow5_datawriter_0_20240729T093114.hdf5 /exp/dune/data/users/gabrielb