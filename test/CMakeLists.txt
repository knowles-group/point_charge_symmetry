include(FetchContent)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG "release-1.11.0"
)
FetchContent_MakeAvailable(googletest)

include(GoogleTest)
foreach (t test-${PROJECT_NAME})
    add_executable(${t} ${t}.cpp)
    target_link_libraries(${t} PUBLIC gmock_main molpro::${PROJECT_NAME})
    if (NOT CMAKE_CROSSCOMPILING)
        gtest_discover_tests(${t})
        add_test(${t} ${t})
    endif ()
endforeach ()

foreach (molecule pentagonx c60_from_klaus co2 adamantane n2-displaced Icosahedron c60_clean c60 phloroglucinol ethane benzene ferrocene ferrocene-2 dummy ar h2o h2o-nosym allene buckminsterfullerene hexamethylbenzene p4 sulfur-hexafluoride cyclohexane s8 methane ch4 allene45 n2 waterinf 18-crown ethene co2-bent co)
    configure_file(${molecule}.xyz . COPYONLY)
endforeach ()
