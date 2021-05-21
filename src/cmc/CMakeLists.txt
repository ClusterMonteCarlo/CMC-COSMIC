# add the executable
add_executable(cmc cmc.c)
# add library
add_library(cmc_library STATIC cmc_bhlosscone.c cmc_binbin.c cmc_binsingle.c cmc_core.c
              cmc_dynamics.c cmc_dynamics_helper.c cmc_bse_utils.c
              cmc_evolution_thr.c cmc_fits.c  
              cmc_io.c cmc_nr.c cmc_orbit.c
              cmc_remove_star.c cmc_search_grid.c cmc_sort.c cmc_sscollision.c
              cmc_stellar_evolution.c cmc_utils.c cmc_mpi.c)
# Include paths to headers
include_directories ("${PROJECT_SOURCE_DIR}/include/common")
include_directories ("${PROJECT_SOURCE_DIR}/include/cmc")
include_directories ("${PROJECT_SOURCE_DIR}/include/bse_wrap")
include_directories(SYSTEM ${MPI_INCLUDE_DIRS})
include_directories(SYSTEM ${GSL_INCLUDE_DIRS})
include_directories(SYSTEM ${CFITSIO_INCLUDE_DIR})
include_directories(SYSTEM ${ZLIB_INCLUDE_DIRS})
include_directories(SYSTEM ${HDF5_INCLUDE_DIRS})

# MPI compile options
target_compile_options(cmc PRIVATE ${MPI_COMPILE_FLAGS})

# link library to executable
target_link_libraries(cmc m)
target_link_libraries(cmc_library ${ZLIB_LIBRARIES})
target_link_libraries(cmc_library ${HDF5_LIBRARIES})
target_link_libraries(cmc_library ${HDF5_HL_LIBRARIES})
target_link_libraries (cmc cmc_library)
target_link_libraries (cmc fewbody)
target_link_libraries (cmc bsewrap)
target_link_libraries (cmc support)
target_link_libraries(cmc ${MPI_LIBRARIES} ${MPI_LINK_FLAGS})
target_link_libraries(cmc ${GSL_LIBRARIES})
target_link_libraries(cmc ${ZLIB_LIBRARIES})
target_link_libraries(cmc ${HDF5_LIBRARIES})
target_link_libraries(cmc ${HDF5_HL_LIBRARIES})

install(TARGETS cmc DESTINATION bin)
install(TARGETS cmc_library DESTINATION lib)
