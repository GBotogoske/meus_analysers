add_library(utils STATIC
    MyMatch.cc
)

# Instala headers
install(FILES
    MyMatch.hh
    my_utils.hh
    DESTINATION include/my_dune/utils
)

# Caminho de includes para QUEM USA a lib
target_include_directories(utils
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include/my_dune/utils>
)

# Instalar a lib e exportar corretamente
install(TARGETS utils EXPORT my_dune)
