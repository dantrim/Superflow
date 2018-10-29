###################################################################
# Package: Superflow
###################################################################

# declare the package name
atlas_subdir( Superflow )

# extra deps based on the build environment
set( extra_deps )
set( extra_libs )

# declare the Superflow package dependencies
atlas_depends_on_subdirs(
    PUBLIC
    SusyNtuple
    ${extra_deps}
)

# external dependencies
find_package( ROOT COMPONENTS Gpad Graf Graf3d Core Tree MathCore Hist RIO )

# libraries
atlas_add_library( SuperflowLib
    Superflow/*.h Root/*.cxx
    PUBLIC_HEADERS Superflow
    PRIVATE_INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
    LINK_LIBRARIES SusyNtupleLib
    PRIVATE_LINK_LIBRARIES ${ROOT_LIBRARIES}
)

# exectuables
set( extra_libs )

# define a function to call the add_executable macro on our executables
function(SuperflowExec filename)
    set(execname)
    get_filename_component(execname ${filename} NAME_WE)
    atlas_add_executable( ${execname} "util/${execname}.cxx"
        INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
        LINK_LIBRARIES ${ROOT_LIBRARIES} SuperflowLib ${extra_libs}
    )
endfunction(SuperflowExec)

# build everything in util/ as an executable
file(GLOB files "util/*.cxx")
foreach(file ${files})
    SuperflowExec(${file})
endforeach()

atlas_install_data( data/* )


