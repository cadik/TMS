set(PALLAS_SOURCES
    basinhopping.cc
    brute.cc
    cooling_schedule.cc
    differential_evolution.cc
    history_concept.cc
    simulated_annealing.cc
    step_function.cc
    types.cc
    internal/crossover_strategy.cc
    internal/mutation_strategy.cc
    internal/state.cc
    internal/stringprintf.cc
    internal/wall_time.cc)

add_library(pallas ${PALLAS_SOURCES})
target_link_libraries(pallas ${CERES_LIBRARIES} ${GLOG_LIBRARIES})

if(BUILD_PALLAS_TESTS)
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    # force this option to ON so that Google Test will use /MD instead of /MT
    # /MD is now the default for Visual Studio, so it should be our default, too
    option(gtest_force_shared_crt
           "Use shared (DLL) run-time lib even when Google Test is built as static lib."
           ON)
    endif()

    include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

    MACRO (PALLAS_TEST NAME)
        ADD_EXECUTABLE(${NAME}_test ${NAME}_test.cc)
        TARGET_LINK_LIBRARIES(${NAME}_test pallas gtest)
        ADD_TEST(NAME ${NAME}_test
                 COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${NAME}_test
                 --test_srcdir
                 ${CMAKE_SOURCE_DIR}/data)
    ENDMACRO (PALLAS_TEST)

    PALLAS_TEST(pallas)

endif(BUILD_PALLAS_TESTS)



