# BSE/COSMIC + SSE + METISSE Fortran sources.
# Tracks COSMIC's src/cosmic/meson.build lib_source list. Keep in sync when bumping the submodule.

set(COSMIC_SRC ./COSMIC/src/cosmic/src)

set(COSMIC_CORE_SOURCES
    ${COSMIC_SRC}/hrdiag_remnant.f
    ${COSMIC_SRC}/assign_remnant.f
    ${COSMIC_SRC}/benchmarkevolv2.f
    ${COSMIC_SRC}/int64.f
    ${COSMIC_SRC}/corerd.f
    ${COSMIC_SRC}/comenv.f
    ${COSMIC_SRC}/dgcore.f
    ${COSMIC_SRC}/evolv2.f
    ${COSMIC_SRC}/gntage.f
    ${COSMIC_SRC}/instar.f
    ${COSMIC_SRC}/kick.f
    ${COSMIC_SRC}/mix.f
    ${COSMIC_SRC}/mrenv.f
    # ran3.f intentionally excluded — CMC's taus113-ran3.f provides
    # ran3 (Tausworth) plus the RandomNormal/RandomTruncatedNormal
    # helpers that the rest of COSMIC's Fortran calls.
    ${COSMIC_SRC}/rl.f
    ${COSMIC_SRC}/hrdiag.f
    ${COSMIC_SRC}/star.f
    ${COSMIC_SRC}/concatkstars.f
    ${COSMIC_SRC}/comprad.f
    ${COSMIC_SRC}/bpp_array.f
    ${COSMIC_SRC}/checkstate.f
    ${COSMIC_SRC}/deltat.f
    ${COSMIC_SRC}/mlwind.f
    ${COSMIC_SRC}/zcnsts.f
)

# CMC-owned Tausworth RNG sources. taus113-ran3.f provides ran3 (overrides
# COSMIC's NR ran3 via linker symbol resolution) plus RandomNormal and
# RandomTruncatedNormal (which the rest of COSMIC's Fortran calls).
# tausworth.f / tausworth.h provide the Tausworth-113 helpers.
set(CMC_TAUSWORTH_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/tausworth.f
    ${CMAKE_CURRENT_SOURCE_DIR}/taus113-ran3.f
)

set(SSE_SOURCES
    ${COSMIC_SRC}/SSE/SSE_deltat.f
    ${COSMIC_SRC}/SSE/SSE_mlwind.f
    ${COSMIC_SRC}/SSE/SSE_hrdiag.f
    ${COSMIC_SRC}/SSE/SSE_star.f
    ${COSMIC_SRC}/SSE/SSE_zcnsts.f
    ${COSMIC_SRC}/SSE/SSE_zfuncs.f
    ${COSMIC_SRC}/SSE/SSE_gntage.f
)

set(METISSE_SOURCES
    ${COSMIC_SRC}/METISSE/src/track_support.f90
    ${COSMIC_SRC}/METISSE/src/c_m_interface.f90
    ${COSMIC_SRC}/METISSE/src/METISSE_gntage.f90
    ${COSMIC_SRC}/METISSE/src/METISSE_deltat.f90
    ${COSMIC_SRC}/METISSE/src/METISSE_mlwind.f90
    ${COSMIC_SRC}/METISSE/src/METISSE_hrdiag.f90
    ${COSMIC_SRC}/METISSE/src/METISSE_star.f90
    ${COSMIC_SRC}/METISSE/src/METISSE_zcnsts.f90
    ${COSMIC_SRC}/METISSE/src/z_support.f90
    ${COSMIC_SRC}/METISSE/src/sse_support.f90
    ${COSMIC_SRC}/METISSE/src/remnant_support.f90
    ${COSMIC_SRC}/METISSE/src/interp_support.f90
    ${COSMIC_SRC}/METISSE/src/comenv_lambda.f90
    ${COSMIC_SRC}/METISSE/src/METISSE_miscellaneous.f90
    ${COSMIC_SRC}/METISSE_utils.f90
)

add_library(bse STATIC
    ${COSMIC_CORE_SOURCES}
    ${CMC_TAUSWORTH_SOURCES}
    ${SSE_SOURCES}
    ${METISSE_SOURCES}
)

# Include paths to headers. CMC's dir first so our tausworth.h overrides
# COSMIC's stale copy; COSMIC's dir still provides int64.h and const_bse.h.
include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${COSMIC_SRC})

install(TARGETS bse DESTINATION lib)

IF(BUILD_COSMIC)
	install(CODE "execute_process(COMMAND ${Python3_EXECUTABLE} -m pip install numpy)")
	install(CODE "execute_process(COMMAND ${Python3_EXECUTABLE} -m pip install . WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/COSMIC/)")
ENDIF(BUILD_COSMIC)
