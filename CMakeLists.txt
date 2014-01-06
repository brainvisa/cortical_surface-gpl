cmake_minimum_required( VERSION 2.6 )
find_package( brainvisa-cmake REQUIRED )
BRAINVISA_FIND_PACKAGE( brainvisa-share REQUIRED )
BRAINVISA_FIND_PACKAGE( axon REQUIRED )
BRAINVISA_FIND_PACKAGE( morphologist-gpl QUIET )
BRAINVISA_PROJECT()

BRAINVISA_DEPENDENCY( RUN DEPENDS "morphologist-gpl" RUN ">= ${morphologist-gpl_VERSION_MAJOR}.${morphologist-gpl_VERSION_MINOR}" )
BRAINVISA_DEPENDENCY( RUN DEPENDS "axon" RUN ">= ${axon_VERSION_MAJOR}.${axon_VERSION_MINOR}" )

BRAINVISA_COPY_PYTHON_DIRECTORY( "${CMAKE_CURRENT_SOURCE_DIR}/brainvisa"
                                 ${PROJECT_NAME} )

BRAINVISA_COPY_PYTHON_DIRECTORY( "${CMAKE_CURRENT_SOURCE_DIR}/python"
                                 ${PROJECT_NAME} )
                               
BRAINVISA_COPY_DIRECTORY( "${CMAKE_CURRENT_SOURCE_DIR}/bin"
                          bin
                          ${PROJECT_NAME} )
BRAINVISA_COPY_DIRECTORY( "${CMAKE_CURRENT_SOURCE_DIR}/share-data"
                          share/brainvisa-share-${brainvisa-share_VERSION_MAJOR}.${brainvisa-share_VERSION_MINOR}
                          ${PROJECT_NAME} )
BRAINVISA_COPY_DIRECTORY( "${CMAKE_CURRENT_SOURCE_DIR}/share/doc/images"
                          share/doc/axon-${axon_VERSION_MAJOR}.${axon_VERSION_MINOR}/images
                          ${PROJECT_NAME} )

BRAINVISA_CREATE_CMAKE_CONFIG_FILES()

