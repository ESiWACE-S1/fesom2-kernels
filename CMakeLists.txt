cmake_minimum_required(VERSION 3.8...3.19)

project(fesom2-benchmarks VERSION 1.0.0 LANGUAGES C Fortran)

option(WITH_ACC "compile with OpenACC" OFF)
# https://cmake.org/cmake/help/latest/manual/cmake-generator-expressions.7.html
# https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_ID.html#variable:CMAKE_%3CLANG%3E_COMPILER_ID
add_library(options_flags INTERFACE)
target_compile_options(options_flags INTERFACE
  # or use $<COMPILE_LANG_AND_ID:language,compiler_ids>
  $<$<Fortran_COMPILER_ID:Intel>:
    -O3
   >

   $<$<Fortran_COMPILER_ID:GNU>:
    -O3 -finit-local-zero  -finline-functions -march=native -fimplicit-none -ffree-line-length-none
    $<$<VERSION_GREATER_EQUAL:${CMAKE_Fortran_COMPILER_VERSION},10>:
      -fallow-argument-mismatch # gfortran v10 is strict about erroneous API calls: "Rank mismatch between actual argument at (1) and actual argument at (2) (scalar and rank-1)"
     >
   >

  $<$<Fortran_COMPILER_ID:PGI>:
    -fast -fastsse -O3 -Minline -Munroll
    $<$<BOOL:${WITH_ACC}>:
      -acc -ta=tesla:cc80
      -Minfo=accel
     >
   >
  )

target_link_options(options_flags INTERFACE
  $<$<Fortran_COMPILER_ID:PGI>:
    $<$<BOOL:${WITH_ACC}>:
      -acc -ta=tesla:cc80
    >
  >
)
target_link_libraries(options_flags INTERFACE wallclock)

add_library(wallclock wallclock.c)

add_executable(oce_adv_tra_fct_loop_a1            wallclock.F90 oce_adv_tra_fct_loop_a1.F90)
add_executable(oce_adv_tra_fct_loop_a1_vlimit1    wallclock.F90 oce_adv_tra_fct_loop_a1_vlimit1.F90)
add_executable(oce_adv_tra_fct_loop_a2            wallclock.F90 oce_adv_tra_fct_loop_a2.F90)
add_executable(oce_adv_tra_fct_loop_b1_horizontal wallclock.F90 oce_adv_tra_fct_loop_b1_horizontal.F90)

target_link_libraries(oce_adv_tra_fct_loop_a1            options_flags)
target_link_libraries(oce_adv_tra_fct_loop_a1_vlimit1    options_flags)
target_link_libraries(oce_adv_tra_fct_loop_a2            options_flags)
target_link_libraries(oce_adv_tra_fct_loop_b1_horizontal options_flags)