### Pre-process: .fpp -> .f90 via Fypp

# Create a list of the files to be preprocessed
#set(fppFiles file1.fpp file2.fpp file3.fpp)

# Pre-process
macro(PreProcessFyppFiles)
	message(STATUS "PreProcessFyppFiles: ${FYPP_FLAGS} ${ARGN}")
foreach(infileName ${ARGN})

    # Generate output file name
    string(REGEX REPLACE ".fpp\$"  ".F90" outfileName "${infileName}")
    string(REGEX REPLACE ".fypp\$" ".F90" outfileName "${infileName}")

    message(STATUS ": ${infileName} -> ${outfileName}")
    # Create the full path for the new file
    set(outfile "${CMAKE_CURRENT_BINARY_DIR}/${outfileName}")

    # Generate input file name
    set(infile "${CMAKE_CURRENT_SOURCE_DIR}/${infileName}")

    # Custom command to do the processing
    add_custom_command(
        OUTPUT "${outfile}"
	COMMAND ${CMAKE_SOURCE_DIR}/bin/fypp ${FYPP_FLAGS} "${infile}" "${outfile}"
        MAIN_DEPENDENCY "${infile}"
        VERBATIM)

## Finally add output file to a list
#set(outFiles ${outFiles} "${outfile}")

endforeach(infileName)
endmacro(PreProcessFyppFiles)
