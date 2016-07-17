INCLUDE_DIRECTORIES( BEFORE
  ${AO-LLSM_SOURCE_DIR}/ImageStitching
  ${AO-LLSM_BINARY_DIR}/ImageStitching
  ${AO-LLSM_SOURCE_DIR}/Options
)

ADD_EXECUTABLE(roiExtract ROIExtract.cxx )
TARGET_LINK_LIBRARIES(roiExtract ${ITK_LIBRARIES} )

ADD_EXECUTABLE(maxProj MaxProj.cxx)
TARGET_LINK_LIBRARIES(maxProj ${ITK_LIBRARIES})

ADD_EXECUTABLE(deconvolve RichardsonLucyDeconvolution.cxx)
TARGET_LINK_LIBRARIES(deconvolve ${ITK_LIBRARIES})

ADD_EXECUTABLE( tileReader TileReader.cxx)
TARGET_LINK_LIBRARIES( tileReader ${ITK_LIBRARIES} Option )

ADD_EXECUTABLE( convertFormat ConvertFormat.cxx)
TARGET_LINK_LIBRARIES( convertFormat ${ITK_LIBRARIES} Option )
