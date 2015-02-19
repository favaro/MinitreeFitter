# Source files to compile                                                                                                                                                          
FILES := MiniTreeFitter1D_background
FILES += MiniTreeFitter1D_signal
FILES += MiniTreeFitter1D_simFit
FILES += HiggsCrossSectionReader
FILES += SetupReader
FILES += SpinFitterConfig
FILES += RooPower

DICTFILES :=

PROGRAMS := OptiMvaPhoID
PROGRAMS += MiniTreeFitter
PROGRAMS += VBFFitter

NEEDS_ROOT  := yes
NEEDS_BOOST := yes