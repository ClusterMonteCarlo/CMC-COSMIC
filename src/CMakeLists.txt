
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_subdirectory(bse_wrap)
add_subdirectory(fewbody-0.24)
add_subdirectory(libs)
add_subdirectory(cmc)
IF(FITS)
add_subdirectory(utils)
ENDIF(FITS)
