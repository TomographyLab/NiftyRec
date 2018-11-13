
#-----------------------------------------------------------------------------
# SET(NAME tclap)
# ADD_LIBRARY(${NAME} ${LIB_TYPE} nifti1.h nifti1_io.c nifti1_io.h znzlib.c znzlib.h)
# INSTALL(TARGETS ${NAME} DESTINATION ${LIB_INSTALL_DIR})
INSTALL(FILES CmdLineInterface.h ArgException.h CmdLine.h XorHandler.h MultiArg.h UnlabeledMultiArg.h ValueArg.h UnlabeledValueArg.h 
Visitor.h Arg.h HelpVisitor.h SwitchArg.h MultiSwitchArg.h VersionVisitor.h IgnoreRestVisitor.h CmdLineOutput.h StdOutput.h DocBookOutput.h 
ZshCompletionOutput.h OptionalUnlabeledTracker.h Constraint.h ValuesConstraint.h ArgTraits.h StandardTraits.h DESTINATION ${INC_INSTALL_DIR}/tclap/)
#-----------------------------------------------------------------------------
