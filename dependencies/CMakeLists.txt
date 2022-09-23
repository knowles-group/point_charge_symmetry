include(FetchContent)
FetchContent_Declare(
        dependency_manager
        GIT_REPOSITORY https://github.com/molpro/dependency-manager.git
        GIT_TAG 0.2.0
)
FetchContent_MakeAvailable(dependency_manager)
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" PARENT_SCOPE)

DependencyManager_Declare(profiler https://github.com/molpro/profiler.git)

DependencyManager_Declare(iterative-solver https://github.com/molpro/iterative-solver.git)

DependencyManager_Declare(utilities https://github.com/molpro/utilities.git)
