# add library
add_library(bsewrap STATIC bse_wrap.c)

# Include paths to headers
include_directories ("${PROJECT_SOURCE_DIR}/include/common")
include_directories ("${PROJECT_SOURCE_DIR}/include/bse_wrap")
# link library to executable
target_link_libraries(bsewrap bse)

add_subdirectory(bse)
