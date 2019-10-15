# add the executable
add_executable(cmc_mkplummer mkplummer/mkplummer.c)
# add the path to the headers
include_directories ("${PROJECT_SOURCE_DIR}/include/common")
# link library to executable
target_link_libraries(cmc_mkplummer support)

# add the executable
add_executable(cmc_mkelson mkelson/mkelson.c)
# add the path to the headers
include_directories ("${PROJECT_SOURCE_DIR}/include/common")
include_directories(SYSTEM ${GSL_INCLUDE_DIRS})
# link library to executable
target_link_libraries(cmc_mkelson support)
target_link_libraries(cmc_mkelson ${GSL_LIBRARIES})

# add the executable
add_executable(cmc_mkking mkking/mkking.c)
# add the path to the headers
include_directories ("${PROJECT_SOURCE_DIR}/include/common")
include_directories(SYSTEM ${GSL_INCLUDE_DIRS})
# link library to executable
target_link_libraries(cmc_mkking support)
target_link_libraries(cmc_mkking ${GSL_LIBRARIES})

# add the executable
add_executable(cmc_setimf setimf/setimf.c setimf/genimf.c )
# add the path to the headers
include_directories ("${PROJECT_SOURCE_DIR}/include/common")
include_directories ("${PROJECT_SOURCE_DIR}/include/utils")
include_directories(SYSTEM ${GSL_INCLUDE_DIRS})
# link library to executable
target_link_libraries(cmc_setimf support)
target_link_libraries(cmc_setimf ${GSL_LIBRARIES})

# add the executable
add_executable(cmc_setunits setunits/setunits.c)
# add the path to the headers
include_directories ("${PROJECT_SOURCE_DIR}/include/common")
# link library to executable
target_link_libraries(cmc_setunits support)

# add the executable
add_executable(cmc_setstellar setstellar/setstellar.c)
# add the path to the headers
include_directories ("${PROJECT_SOURCE_DIR}/include/common")
include_directories ("${PROJECT_SOURCE_DIR}/include/cmc")
include_directories ("${PROJECT_SOURCE_DIR}/include/bse_wrap")
# link library to executable
target_link_libraries(cmc_setstellar support)
target_link_libraries(cmc_setstellar bsewrap)

# add the executable
add_executable(cmc_addbinaries addbinaries/addbinaries.c)
# add the path to the headers
include_directories ("${PROJECT_SOURCE_DIR}/include/common")
include_directories ("${PROJECT_SOURCE_DIR}/include/cmc")
# link library to executable
target_link_libraries(cmc_addbinaries support)
target_link_libraries(cmc_addbinaries bsewrap)

# add the executable
add_executable(cmc_ascii2fits ascii2fits/nbody2fits.c)
# add the path to the headers
include_directories ("${PROJECT_SOURCE_DIR}/include/common")
include_directories(SYSTEM ${GSL_INCLUDE_DIRS})
# link library to executable
target_link_libraries(cmc_ascii2fits support)
target_link_libraries(cmc_ascii2fits ${GSL_LIBRARIES})

install(TARGETS cmc_mkplummer DESTINATION bin)
install(TARGETS cmc_mkelson DESTINATION bin)
install(TARGETS cmc_mkking DESTINATION bin)
install(TARGETS cmc_setimf DESTINATION bin)
install(TARGETS cmc_setunits DESTINATION bin)
install(TARGETS cmc_setstellar DESTINATION bin)
install(TARGETS cmc_addbinaries DESTINATION bin)
install(TARGETS cmc_ascii2fits DESTINATION bin)
