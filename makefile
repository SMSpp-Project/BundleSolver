##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of BundleSolver and ParallelBundleSolver                        #
#                                                                            #
#   The makefile takes in input the -I directives for all the external       #
#   libraries needed by BundleSolver, i.e., core SMS++.                      #
#                                                                            #
#   Note that, conversely, $(SMS++INC) is also assumed to include any        #
#   -I directive corresponding to external libraries needed by SMS++, at     #
#   least to the extent in which they are needed by the parts of SMS++       #
#   used by BundleSolver.                                                    #
#                                                                            #
#   Input:  $(CC)          = compiler command                                #
#           $(SW)          = compiler options                                #
#           $(SMS++INC)    = the -I$( core SMS++ include directory )         #
#           $(SMS++OBJ)    = the core SMS++ library                          #
#           $(libNDOINC)   = the -I$( libNDO include directory )             #
#           $(MILPSINC)    = the -I$( MILPSolver include directory )         #
#           $(MILPSH)      = the .h files to include for MILPSolver          #
#           $(BNDSLVSDR)   = the directory where the source is               #
#                                                                            #
#   Output: $(BNDSLVOBJ)   = the final object(s) / library                   #
#           $(BNDSLVH)     = the .h files to include                         #
#           $(BNDSLVINC)   = the -I$( source directory )                     #
#                                                                            #
#                              Antonio Frangioni                             #
#                               Enrico Gorgone                               #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

# macros to be exported - - - - - - - - - - - - - - - - - - - - - - - - - - -

BNDSLVOBJ = $(BNDSLVSDR)/obj/BundleSolver.o \
	$(BNDSLVSDR)/obj/ParallelBundleSolver.o \
	$(BNDSLVSDR)/obj/MasterProblemBlock.o \
	$(BNDSLVSDR)/obj/BundleSolver.o \
	$(BNDSLVSDR)/obj/MILPMPSolver.o

BNDSLVINC = -I$(BNDSLVSDR)/include -I$(BNDSLVSDR)/MILPMPSolver

BNDSLVH   = $(BNDSLVSDR)/include/BundleSolver.h \
	$(BNDSLVSDR)/include/ParallelBundleSolver.h \
	$(BNDSLVSDR)/include/MasterProblemBlock.h \
	$(BNDSLVSDR)/include/BundleSolver.h \
	$(BNDSLVSDR)/MILPMPSolver/MILPMPSolver.h

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(BNDSLVOBJ) $(BNDSLVSDR)/*~

# dependencies: every .o from its .cpp + every recursively included .h- - - -

$(BNDSLVSDR)/obj/LegacyBundleSolver.o: $(BNDSLVSDR)/src/LegacyBundleSolver.cpp \
	$(BNDSLVSDR)/include/LegacyBundleSolver.h $(SMS++OBJ) $(MILPSH) $(libNDOOBJ)
	$(CC) -c $(BNDSLVSDR)/src/LegacyBundleSolver.cpp -o $@ $(BNDSLVINC) \
	$(SMS++INC) $(MILPSINC) $(libNDOINC) $(SW)

$(BNDSLVSDR)/obj/ParallelBundleSolver.o: $(BNDSLVSDR)/src/ParallelBundleSolver.cpp \
	$(BNDSLVH) $(SMS++OBJ) $(MILPSH) $(libNDOOBJ)
	$(CC) -c $(BNDSLVSDR)/src/ParallelBundleSolver.cpp -o $@ \
	$(BNDSLVINC) $(SMS++INC) $(MILPSINC) $(libNDOINC) $(SW)

$(BNDSLVSDR)/obj/MasterProblemBlock.o: $(BNDSLVSDR)/src/MasterProblemBlock.cpp \
	$(BNDSLVSDR)/include/MasterProblemBlock.h $(SMS++OBJ)
	$(CC) -c $(BNDSLVSDR)/src/MasterProblemBlock.cpp -o $@ \
	$(BNDSLVINC) $(SMS++INC) $(SW)

$(BNDSLVSDR)/obj/BundleSolver.o: $(BNDSLVSDR)/src/BundleSolver.cpp \
	$(BNDSLVSDR)/include/BundleSolver.h \
	$(BNDSLVSDR)/include/MasterProblemBlock.h $(SMS++OBJ) $(MILPSH)
	$(CC) -c $(BNDSLVSDR)/src/BundleSolver.cpp -o $@ \
	$(BNDSLVINC) $(SMS++INC) $(MILPSINC) $(SW)

$(BNDSLVSDR)/obj/MILPMPSolver.o: $(BNDSLVSDR)/MILPMPSolver/MILPMPSolver.cpp \
	$(BNDSLVSDR)/MILPMPSolver/MILPMPSolver.h $(SMS++OBJ) $(MILPSH) $(libNDOOBJ)
	$(CC) -c $(BNDSLVSDR)/MILPMPSolver/MILPMPSolver.cpp -o $@ \
	$(BNDSLVINC) $(SMS++INC) $(MILPSINC) $(libNDOINC) $(SW)

########################## End of makefile ###################################
