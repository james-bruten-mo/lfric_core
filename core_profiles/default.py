fortran_compile_flags_common =[
    "-std=f2008",
    "-ffree-line-length-none",
]

compile_profile = {
    "min_compiler_version": 12,
    "fortran_compile_flags_common": fortran_compile_flags_common,
    "openmp": True,
}
