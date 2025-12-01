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

cet_build_plugin(MyWaveformAnalyzerDetSimMC art::module
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

cet_build_plugin(GetMCInfoPDSAnalyzer art::module
  LIBRARIES
    larcorealg::Geometry
    lardataobj::RawData
    dunecore::DuneObj                       
    larsim::MCCheater_ParticleInventoryService_service  
    larsim::MCCheater_PhotonBackTrackerService_service
    art::Framework_Core
    art::Framework_Principal
    canvas::canvas
    duneopdet::OpticalDetector
    art_root_io::TFileService_service
    ROOT::Core
    ROOT::Tree
)

cet_build_plugin(GetAllTimesAnalyzer art::module
  LIBRARIES
    larcorealg::Geometry
    lardataobj::RawData
    dunecore::DuneObj                       
    larsim::MCCheater_ParticleInventoryService_service  
    larsim::MCCheater_PhotonBackTrackerService_service
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
    art::Framework_Core
    art::Framework_Principal
    canvas::canvas
    duneopdet::OpticalDetector
    art_root_io::TFileService_service
    ROOT::Core
    ROOT::Tree
)

cet_build_plugin(GetPandoraTracksAnalyzer art::module
  LIBRARIES
    larcorealg::Geometry
    lardataobj::RawData
    dunecore::DuneObj                       
    larsim::MCCheater_ParticleInventoryService_service  
    larsim::MCCheater_PhotonBackTrackerService_service
    art::Framework_Core
    art::Framework_Principal
    canvas::canvas
    duneopdet::OpticalDetector
    art_root_io::TFileService_service
    ROOT::Core
    ROOT::Tree
)

add_subdirectory(utils)

cet_build_plugin(MyFlashMatching art::module
  LIBRARIES
    larcorealg::Geometry
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
    utils
)

cet_build_plugin(GetFlashMatchingMC art::module
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

install_headers()
install_fhicl()
install_source()
install_scripts()
