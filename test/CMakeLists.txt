include(FetchContent)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG "release-1.10.0"
)
FetchContent_MakeAvailable(googletest)

include(GoogleTest)
foreach (t test-${PROJECT_NAME} test-bfgs)
    add_executable(${t} ${t}.cpp)
    target_link_libraries(${t} PUBLIC gmock_main molpro::${PROJECT_NAME})
    if (NOT CMAKE_CROSSCOMPILING)
        gtest_discover_tests(${t})
        add_test(${t} ${t})
    endif ()
endforeach ()

foreach (molecule benzene ferrocene h2o h2o-nosym allene buckminsterfullerene hexamethylbenzene p4 sulfur-hexafluoride cyclohexane s8 methane ch4 allene45 n2)
    configure_file(${molecule}.xyz . COPYONLY)
endforeach ()
