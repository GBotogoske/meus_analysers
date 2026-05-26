
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
    larevt::SpaceCharge
    art_root_io::TFileService_service
    ROOT::Core
    ROOT::Tree
    duneopdet::OpticalDetector
    utils
)

cet_build_plugin(MyFlashMatchingMC art::module
  LIBRARIES
    larcorealg::Geometry
    lardataobj::RawData
    dunecore::DuneObj                       
    larsim::MCCheater_ParticleInventoryService_service  
    larsim::MCCheater_PhotonBackTrackerService_service
    larsim::MCCheater_BackTrackerService_service
    lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
    lardataalg::DetectorInfo 
    art::Framework_Core
    art::Framework_Principal
    canvas::canvas
    duneopdet::OpticalDetector
    art_root_io::TFileService_service
    ROOT::Core
    ROOT::Tree
)

cet_build_plugin(getophitMCPE art::module
  LIBRARIES
    larcorealg::Geometry
    lardataobj::RawData
    dunecore::DuneObj                       
    larsim::MCCheater_ParticleInventoryService_service  
    larsim::MCCheater_PhotonBackTrackerService_service
    larsim::MCCheater_BackTrackerService_service
    lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
    lardataalg::DetectorInfo 
    art::Framework_Core
    art::Framework_Principal
    canvas::canvas
    duneopdet::OpticalDetector
    art_root_io::TFileService_service
    ROOT::Core
    ROOT::Tree
)

cet_build_plugin(GetMyWireData art::module
  LIBRARIES
    larcorealg::Geometry
    lardataobj::RawData
    dunecore::DuneObj                       
    larsim::MCCheater_ParticleInventoryService_service  
    larsim::MCCheater_PhotonBackTrackerService_service
    lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
    art::Framework_Core
    art::Framework_Principal
    canvas::canvas
    duneopdet::OpticalDetector
    art_root_io::TFileService_service
    ROOT::Core
    ROOT::Tree
)

cet_build_plugin(GetWfHitsAnalyzer art::module
  LIBRARIES
    larcorealg::Geometry
    lardataobj::RawData
    dunecore::DuneObj                       
    larsim::MCCheater_ParticleInventoryService_service  
    larsim::MCCheater_PhotonBackTrackerService_service
    lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
    art::Framework_Core
    art::Framework_Principal
    canvas::canvas
    duneopdet::OpticalDetector
    art_root_io::TFileService_service
    ROOT::Core
    ROOT::Tree
)

cet_build_plugin(getVis art::module
  LIBRARIES
    larcorealg::Geometry
    lardataobj::RawData
    dunecore::DuneObj                       
    larsim::MCCheater_ParticleInventoryService_service  
    larsim::MCCheater_PhotonBackTrackerService_service
    lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
    art::Framework_Core
    art::Framework_Principal
    canvas::canvas
    duneopdet::OpticalDetector
    art_root_io::TFileService_service
    ROOT::Core
    ROOT::Tree
    larcore::Geometry_Geometry_service
    larsim::PhotonPropagation
    larsim::PhotonPropagation_PhotonVisibilityService_service 
)

cet_build_plugin(CompareSimPhotonsOpHits art::module
  LIBRARIES
    larcorealg::Geometry
    lardataobj::RawData
    dunecore::DuneObj                       
    lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
    lardataalg::DetectorInfo 
    art::Framework_Core
    art::Framework_Principal
    canvas::canvas
    duneopdet::OpticalDetector
    art_root_io::TFileService_service
    ROOT::Core
    ROOT::Tree
)

cet_build_plugin(MyFlashMatchingGetClusters art::module
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
    larevt::SpaceCharge
    canvas::canvas
    duneopdet::OpticalDetector
    art_root_io::TFileService_service
    ROOT::Core
    ROOT::Tree
    duneopdet::OpticalDetector
    utils
)

cet_build_plugin(MyFlashMatchingGetClustersEfield art::module
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
    larevt::SpaceCharge
    canvas::canvas
    duneopdet::OpticalDetector
    art_root_io::TFileService_service
    ROOT::Core
    ROOT::Tree
    duneopdet::OpticalDetector
    utils
)

cet_build_plugin(GetHitFlashJoinAnalyzer art::module
  LIBRARIES
    lardataobj::RecoBase
    art::Framework_Core
    art::Framework_Principal
    art_root_io::TFileService_service
    canvas::canvas
    ROOT::Core
    ROOT::Tree
)

cet_build_plugin(PrintBeamInfo art::module
  LIBRARIES
    art::Framework_Core
    art::Framework_Principal
    art::Framework_Services_Registry
    canvas::canvas
    fhiclcpp::fhiclcpp
    messagefacility::MF_MessageLogger
    dunecore::DuneObj
    art_root_io::TFileService_service
    ROOT::Core
    ROOT::Tree
    ROOT::RIO
)

cet_build_plugin(MyWaveformAnalyzerDetSim art::module
  LIBRARIES
    larcorealg::Geometry
    lardataobj::RawData
    dunecore::DuneObj                       # << Novo correto para OpDetDivRec
    larsim::MCCheater_ParticleInventoryService_service  # Novo: ParticleInventoryService
    art::Framework_Core
    art::Framework_Principal
    canvas::canvas
    art_root_io::TFileService_service
    ROOT::Core
    ROOT::Tree
)

install_headers()
install_fhicl()
install_source()
install_scripts()
