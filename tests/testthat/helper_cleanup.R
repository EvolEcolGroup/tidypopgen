# Delete, when the calling test file ends, the files that it left in
# tempdir() (mostly the .bk, .rds and .gt files of its gen_tibbles).
#
# Call it at the top level of a test file. A new tempfile() name can match
# the name of a leftover .bk file (tempfile() checks the name without the
# extension, and its random part has 31 bits), and the backing file can then
# not be created ("File '...bk' already exists"): the fewer files left in
# tempdir() over the suite, the lower that chance. Only files are deleted
# (not directories), and only those created while the test file ran, which
# no other test file uses.
#
# @param envir the frame whose exit triggers the clean up (the test file)
local_tempdir_cleanup <- function(envir = parent.frame()) {
  files_before <- list.files(tempdir(), full.names = TRUE)
  withr::defer(
    {
      new_files <- setdiff(
        list.files(tempdir(), full.names = TRUE),
        files_before
      )
      unlink(new_files[!dir.exists(new_files)])
    },
    envir = envir
  )
}
