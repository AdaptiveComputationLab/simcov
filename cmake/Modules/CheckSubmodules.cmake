find_package(Git QUIET)
function(CHECK_SUBMODULES)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
    # Update submodule if needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    set(EXPECTED_SUBMODULES ${ARGN})
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update of ${EXPECTED_SUBMODULES}")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT
                        OUTPUT_VARIABLE GIT_SUBMOD_OUTPUT
                        )
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()

        string(REPLACE "\n" ";" SUBMOD_LIST ${GIT_SUBMOD_OUTPUT})
        set(UPDATE_SUBMODULES "")
        foreach(tmp ${EXPECTED_SUBMODULES})
           set(IS_SUB_OK FALSE)
           foreach(tmp2 ${SUBMOD_LIST})
             string(REGEX MATCH "^ .* ${tmp} " IS_OK ${tmp2})
             if(IS_OK)
               set(IS_SUB_OK TRUE)
               message(STATUS "${tmp2} matched '^ .* ${tmp}'")
               break()
             endif()
           endforeach()
           if(NOT IS_SUB_OK)
             message(STATUS "${tmp} submodule needs update")
             set(UPDATE_SUBMODULES ${UPDATE_SUBMODULES} ${tmp})
           endif()
        endforeach()
        foreach(tmp ${UPDATE_SUBMODULES})
          message("Executing git submodule update for ${tmp}")
          execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive ${tmp}
                          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                          RESULT_VARIABLE GIT_SUBMOD_RESULT)
          if(NOT GIT_SUBMOD_RESULT EQUAL "0")
              message(WARNING "'git submodule update --init --recursive ${tmp}' failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
          endif()
        endforeach()
    endif()
endif()
endfunction()

