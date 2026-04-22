nextflow.enable.dsl=2

include { ECHO } from './workflows/echo.nf'

workflow {
  ECHO()
}

