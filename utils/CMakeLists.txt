set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package(ROOT REQUIRED COMPONENTS Minuit)

add_library(utils STATIC
    MyMatch.cc
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
        duneopdet::OpticalDetector 
)

# Instala headers
install(FILES
    MyMatch.hh
    my_utils.hh
    DESTINATION include/my_dune/utils
)

# Instalar a lib e exportar corretamente
install(TARGETS utils EXPORT my_dune)
