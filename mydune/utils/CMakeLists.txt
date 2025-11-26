add_library(utils INTERFACE)

# Headers serão instalados via install_headers() no CMakeLists principal.
target_include_directories(utils
  INTERFACE
    $<INSTALL_INTERFACE:include>
)

install(TARGETS utils EXPORT my_dune)