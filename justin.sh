/cvmfs/dune.opensciencegrid.org/products/dune/justin/justin-sl7-setup
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup python v3_9_15
setup rucio
setup justin

justin get-token
export RUCIO_ACCOUNT=justinreadonly
rucio whoami
rucio list-scopes