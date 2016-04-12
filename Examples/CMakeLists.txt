ADD_EXECUTABLE(roiExtract ROIExtract.cxx )
TARGET_LINK_LIBRARIES(roiExtract ${ITK_LIBRARIES} )

ADD_EXECUTABLE(maxProj MaxProj.cxx)
TARGET_LINK_LIBRARIES(maxProj ${ITK_LIBRARIES})

ADD_EXECUTABLE(deconvolve RichardsonLucyDeconvolution.cxx)
TARGET_LINK_LIBRARIES(deconvolve ${ITK_LIBRARIES})


