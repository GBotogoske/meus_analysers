set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package(ROOT REQUIRED COMPONENTS Minuit)

# Em alguns setups isso já vem do projeto pai.
# Se ROOT_GENERATE_DICTIONARY não for reconhecido, descomente:
include(${ROOT_USE_FILE})


add_library(utils STATIC
    MyMatch.cc
    my_utils.cxx
)

# Caminho de includes para QUEM USA a lib
target_include_directories(utils
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include/my_dune/utils>
        ${ROOT_INCLUDE_DIRS}
)

target_link_libraries(utils
    PUBLIC
        ROOT::Minuit
        larsim::PhotonPropagation
        larsim::PhotonPropagation_PhotonVisibilityService_service
        larevt::SpaceCharge
        duneopdet::OpticalDetector 
)


ROOT_GENERATE_DICTIONARY(utilsDict
    my_utils.hh
    LINKDEF MyUtilsLinkDef.h
    MODULE utils
)

# Instala headers
install(FILES
    MyMatch.hh
    my_utils.hh
    MyUtilsLinkDef.h
    DESTINATION include/my_dune/utils
)

# Instalar a lib e exportar corretamente
install(TARGETS utils EXPORT my_dune)
