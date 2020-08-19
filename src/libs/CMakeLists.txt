add_library(support STATIC fitslib.c taus113-v2.c)

# Include paths to headers
include_directories ("${PROJECT_SOURCE_DIR}/include/common")
include_directories(SYSTEM ${CFITSIO_INCLUDE_DIR})
include_directories(SYSTEM ${HDF5_INCLUDE_DIRS})

# link library to executable
target_link_libraries(support ${CFITSIO_LIBRARIES})
target_link_libraries(support ${HDF5_LIBRARIES})
