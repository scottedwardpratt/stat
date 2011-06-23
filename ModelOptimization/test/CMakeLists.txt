include_directories( ../include )

set( LIBRARIES ModelOptimization )

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