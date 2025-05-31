set(myfiles AARotamer.cpp AtomResInfo.cpp ClashCollection.cpp defs.cpp MutInfo.cpp RotamerCluster.cpp Rotamer.cpp RotamerGeneration.cpp RotamerNeighborList.cpp RotConnInfo.cpp Rotlib.cpp RotlibCollection.cpp sc_AABackBone.cpp sc_AAChain.cpp sc_AASideChain.cpp sc_AminoAcid.cpp sc_BackBone.cpp sc_bgf_handler.cpp sc_Hetatm.cpp sc_Ligand.cpp sc_ProteinComponent.cpp sc_Protein.cpp scream_atom.cpp scream_coulomb_EE.cpp scream_ctl_reader.cpp scream_EE.cpp scream_E_functionals_coulomb.cpp scream_E_functionals_hb.cpp scream_E_functionals_vdw.cpp scream_hb_EE.cpp scream_helper_classes.cpp scream_matrix.cpp scream_model.cpp scream_Rotation.cpp scream_rtf.cpp scream_tools.cpp scream_vdw_EE.cpp scream_vdw_hb_exclusion_EE.cpp scream_vector.cpp scream_VolumeOverlap.cpp sc_SideChain.cpp sc_Water.cpp )

foreach(file ${myfiles})
  set(myfileswithpath ${myfileswithpath} ${SCREAM_SOURCE_DIR}/src/${file})
endforeach()
include_directories(${SCREAM_SOURCE_DIR}/include)

FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})

FIND_PACKAGE(PythonLibs)
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(CMAKE_SWIG_FLAGS "-c++")

SET_SOURCE_FILES_PROPERTIES(py_scream_ee.i PROPERTIES CPLUSPLUS ON)
# SET_SOURCE_FILES_PROPERTIES(py_scream_ee.i PROPERTIES SWIG_FLAGS "-includeall")
SWIG_ADD_MODULE(py_scream_ee python py_scream_ee.i ${myfileswithpath})
SWIG_LINK_LIBRARIES(py_scream_ee ${PYTHON_LIBRARIES})

