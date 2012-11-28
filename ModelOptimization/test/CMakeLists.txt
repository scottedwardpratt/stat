include_directories( 
  ../include 
  /opt/local/include
  ~/local/include
)

set( LIBRARIES 
  IHU
  ${CORELIBS} 
)

add_executable(        Gaussian2DModelTest
                       Gaussian2DModelTest.cxx )
target_link_libraries( Gaussian2DModelTest ${LIBRARIES} )
add_test(              Gaussian2DModelTest
                       Gaussian2DModelTest )

add_executable(        RegularStepGradientDescentOptimizerTest
                       RegularStepGradientDescentOptimizerTest.cxx )
target_link_libraries( RegularStepGradientDescentOptimizerTest ${LIBRARIES} )
add_test(              RegularStepGradientDescentOptimizerTest
                       RegularStepGradientDescentOptimizerTest )

add_executable(        external_model_plus_mhmcmc_test
                       external_model_plus_mhmcmc_test.cxx )
target_link_libraries( external_model_plus_mhmcmc_test ${LIBRARIES} )
add_test(              external_model_plus_mhmcmc_test
                       external_model_plus_mhmcmc_test )
