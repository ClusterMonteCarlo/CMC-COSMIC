# add library
add_library(bse STATIC ./COSMIC/cosmic/src/comenv.f ./COSMIC/cosmic/src/corerd.f ./COSMIC/cosmic/src/deltat.f ./COSMIC/cosmic/src/dgcore.f ./COSMIC/cosmic/src/evolv1.f ./COSMIC/cosmic/src/evolv2.f ./COSMIC/cosmic/src/gntage.f ./COSMIC/cosmic/src/hrdiag.f ./COSMIC/cosmic/src/instar.f ./COSMIC/cosmic/src/kick.f ./COSMIC/cosmic/src/mix.f ./COSMIC/cosmic/src/mlwind.f ./COSMIC/cosmic/src/mrenv.f ./COSMIC/cosmic/src/rl.f ./COSMIC/cosmic/src/star.f ./COSMIC/cosmic/src/zcnsts.f ./COSMIC/cosmic/src/zfuncs.f ./COSMIC/cosmic/src/concatkstars.f ./COSMIC/cosmic/src/bpp_array.f ./COSMIC/cosmic/src/checkstate.f ./COSMIC/cosmic/src/int64.f ./COSMIC/cosmic/src/tausworth.f ./COSMIC/cosmic/src/taus113-ran3.f ./COSMIC/cosmic/src/hrdiag_remnant.f ./COSMIC/cosmic/src/hrdiag.f ./COSMIC/cosmic/src/instar.f ./COSMIC/cosmic/src/kick.f ./COSMIC/cosmic/src/mix.f ./COSMIC/cosmic/src/mlwind.f ./COSMIC/cosmic/src/mrenv.f ./COSMIC/cosmic/src/rl.f ./COSMIC/cosmic/src/star.f ./COSMIC/cosmic/src/zcnsts.f ./COSMIC/cosmic/src/zfuncs.f ./COSMIC/cosmic/src/concatkstars.f ./COSMIC/cosmic/src/bpp_array.f ./COSMIC/cosmic/src/checkstate.f ./COSMIC/cosmic/src/int64.f ./COSMIC/cosmic/src/tausworth.f ./COSMIC/cosmic/src/taus113-ran3.f ./COSMIC/cosmic/src/assign_remnant.f  )

# Include paths to headers
include_directories (./COSMIC/cosmic/src/)

install(TARGETS bse DESTINATION lib)

IF(BUILD_COSMIC)
	install(CODE "execute_process(COMMAND ${Python3_EXECUTABLE} -m pip install numpy)")
	install(CODE "execute_process(COMMAND ${Python3_EXECUTABLE} -m pip install . WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/COSMIC/)")
ENDIF(BUILD_COSMIC)
