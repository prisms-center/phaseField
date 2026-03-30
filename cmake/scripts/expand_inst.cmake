# Run expand_template_instantiations_exe on the file, piping
# output to a tmp output file in case the command fails. This
# is necessary so that we don't try and compile with incomplete
# instantiations.
execute_process(
  COMMAND
    ${EXE} ${TEMPLATES}
  INPUT_FILE ${INPUT}
  OUTPUT_FILE ${OUTPUT_TMP}
  RESULT_VARIABLE _result
)
if(NOT _result EQUAL 0)
  message(FATAL_ERROR "expand_template_instantiations failed for ${INPUT}")
endif()
file(RENAME ${OUTPUT_TMP} ${OUTPUT})
