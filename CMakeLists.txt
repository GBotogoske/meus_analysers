
add_subdirectory(utils)

cet_build_plugin(MyFlashMatching art::module
  LIBRARIES
    larcorealg::Geometry
    larcore::Geometry_Geometry_service
    lardataobj::RawData
    dunecore::DuneObj                       
    larsim::MCCheater_ParticleInventoryService_service  
    larsim::MCCheater_PhotonBackTrackerService_service
    larsim::PhotonPropagation
    larsim::PhotonPropagation_PhotonVisibilityService_service 
    lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
    art::Framework_Core
    art::Framework_Principal
    canvas::canvas
    duneopdet::OpticalDetector
    art_root_io::TFileService_service
    ROOT::Core
    ROOT::Tree
    duneopdet::OpticalDetector
    utils
)

install_headers()
install_fhicl()
install_source()
install_scripts()
