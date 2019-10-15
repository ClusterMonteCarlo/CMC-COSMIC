# add library
add_library(fewbody STATIC fewbody.c fewbody_classify.c fewbody_coll.c fewbody_hier.c fewbody_int.c fewbody_io.c fewbody_isolate.c fewbody_ks.c fewbody_nonks.c fewbody_scat.c fewbody_utils.c)

# Include paths to headers
include_directories ("${PROJECT_SOURCE_DIR}/include/fewbody-0.24")
include_directories(SYSTEM ${GSL_INCLUDE_DIRS})
# link library to executable
target_link_libraries(fewbody ${GSL_LIBRARIES})

install(TARGETS fewbody DESTINATION lib)
