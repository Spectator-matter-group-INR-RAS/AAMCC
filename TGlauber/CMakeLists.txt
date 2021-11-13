project(TGlauberLib)

set(Source
	TGlauNucleon.cc
	TGlauNucleus.cc
	TGlauberMC.cc)

add_library(TGlauberLib
	${Source})
target_link_libraries(TGlauberLib
	${ROOT_LIBRARIES})
