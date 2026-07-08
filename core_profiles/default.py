fortran_compile_flags_common =[
    "-std=f2008",
    "-ffree-line-length-none",
    "-fallow-argument-mismatch",
]

compile_profile = {
    "min_compiler_version": 12,
    "fortran_compile_flags_common": fortran_compile_flags_common,
    "openmp": True,
    "fortran_pp_flags": ["-traditional-cpp", "-P"],
}
