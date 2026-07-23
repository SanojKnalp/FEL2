# Register LinearAlgebra test sources

list(APPEND UNIT_TEST_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/ShapeFunctions_QP_test.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Jacobian_test.cpp
)

# Export list to parent scope
set(UNIT_TEST_SOURCES ${UNIT_TEST_SOURCES} PARENT_SCOPE)
