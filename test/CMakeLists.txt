include(FetchContent)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG "release-1.10.0"
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

#file(COPY h2o.xyz DESTINATION .)
#add_custom_target(h2o.xyz DEPENDS h2o.xyz)
#add_dependencies(test-${PROJECT_NAME} h2o.xyz)
configure_file(h2o.xyz h2o.xyz COPYONLY)