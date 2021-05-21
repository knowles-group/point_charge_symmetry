include(FetchContent)
FetchContent_Declare(
        dependency_manager
        GIT_REPOSITORY https://gitlab.com/dependencymanager/dependency-manager.git
        GIT_TAG 0.2.0
)
FetchContent_MakeAvailable(dependency_manager)
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" PARENT_SCOPE)

DependencyManager_Declare(Profiler https://gitlab+deploy-token-74792:Bj6FL_cWqz_8y2PNDycx@gitlab.com/molpro/profiler.git)

DependencyManager_Declare(LinearAlgebra https://gitlab+deploy-token-186248:BVoGyrUG_o7zpBfCcsLD@gitlab.com/molpro/linearalgebra.git)

DependencyManager_Declare(utilities https://gitlab+deploy-token-181030:whWoM4qkenzjGjWsFTxj@gitlab.com/molpro/utilities.git)
