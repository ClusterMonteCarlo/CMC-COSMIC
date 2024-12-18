# add library
add_library(bse STATIC ./COSMIC/src/cosmic/src/comenv.f ./COSMIC/src/cosmic/src/corerd.f ./COSMIC/src/cosmic/src/deltat.f ./COSMIC/src/cosmic/src/dgcore.f ./COSMIC/src/cosmic/src/evolv1.f ./COSMIC/src/cosmic/src/evolv2.f ./COSMIC/src/cosmic/src/gntage.f ./COSMIC/src/cosmic/src/hrdiag.f ./COSMIC/src/cosmic/src/instar.f ./COSMIC/src/cosmic/src/kick.f ./COSMIC/src/cosmic/src/mix.f ./COSMIC/src/cosmic/src/mlwind.f ./COSMIC/src/cosmic/src/mrenv.f ./COSMIC/src/cosmic/src/rl.f ./COSMIC/src/cosmic/src/star.f ./COSMIC/src/cosmic/src/zcnsts.f ./COSMIC/src/cosmic/src/zfuncs.f ./COSMIC/src/cosmic/src/concatkstars.f ./COSMIC/src/cosmic/src/bpp_array.f ./COSMIC/src/cosmic/src/checkstate.f ./COSMIC/src/cosmic/src/int64.f ./COSMIC/src/cosmic/src/tausworth.f ./COSMIC/src/cosmic/src/taus113-ran3.f ./COSMIC/src/cosmic/src/hrdiag_remnant.f ./COSMIC/src/cosmic/src/hrdiag.f ./COSMIC/src/cosmic/src/instar.f ./COSMIC/src/cosmic/src/kick.f ./COSMIC/src/cosmic/src/mix.f ./COSMIC/src/cosmic/src/mlwind.f ./COSMIC/src/cosmic/src/mrenv.f ./COSMIC/src/cosmic/src/rl.f ./COSMIC/src/cosmic/src/star.f ./COSMIC/src/cosmic/src/zcnsts.f ./COSMIC/src/cosmic/src/zfuncs.f ./COSMIC/src/cosmic/src/concatkstars.f ./COSMIC/src/cosmic/src/bpp_array.f ./COSMIC/src/cosmic/src/checkstate.f ./COSMIC/src/cosmic/src/int64.f ./COSMIC/src/cosmic/src/tausworth.f ./COSMIC/src/cosmic/src/taus113-ran3.f ./COSMIC/src/cosmic/src/assign_remnant.f  )

# Include paths to headers
include_directories (./COSMIC/src/cosmic/src/)

install(TARGETS bse DESTINATION lib)

IF(BUILD_COSMIC)
	install(CODE "execute_process(COMMAND ${Python3_EXECUTABLE} -m pip install numpy)")
	install(CODE "execute_process(COMMAND ${Python3_EXECUTABLE} -m pip install . WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/COSMIC/)")
ENDIF(BUILD_COSMIC)
