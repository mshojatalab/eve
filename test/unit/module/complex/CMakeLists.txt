##==================================================================================================
##  EVE - Expressive Vector Engine
##  Copyright : EVE Contributors & Maintainers
##  SPDX-License-Identifier: MIT
##==================================================================================================

add_custom_target(unit.complex.exe          )
add_dependencies(unit.exe unit.complex.exe  )

make_unit("unit.complex"  constructor.cpp   )
make_unit("unit.complex"  operators.cpp     )
make_unit("unit.complex"  abs.cpp           )
make_unit("unit.complex"  arg.cpp           )
make_unit("unit.complex"  conj.cpp          )
make_unit("unit.complex"  cosh.cpp          )
make_unit("unit.complex"  div.cpp           )
make_unit("unit.complex"  i.cpp             )
make_unit("unit.complex"  if_else.cpp       )
make_unit("unit.complex"  mul.cpp           )
