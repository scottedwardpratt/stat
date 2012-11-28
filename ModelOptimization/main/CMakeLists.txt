include_directories( 
	../include 
	/opt/local/include
	~/local/include
)

set( LIBRARIES 	
  IHU 
  ${CORELIBS}
)

add_executable(        mcmc
                       mcmc.cc )
target_link_libraries( mcmc ${LIBRARIES} )
