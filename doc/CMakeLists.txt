set(KMC_INSTALL_DOCS TRUE CACHE BOOL "Indicates whether to install the documentation")

if (KMC_INSTALL_DOCS)
  install(DIRECTORY html DESTINATION doc)
  install(DIRECTORY example-code DESTINATION doc
    USE_SOURCE_PERMISSIONS)
  install(FILES latex/refman.pdf DESTINATION doc RENAME manual.pdf)
endif()