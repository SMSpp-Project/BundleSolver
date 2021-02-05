##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of BundleSolver                                                 #
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
#           $(SMS++INC)    = the -I$( core SMS++ directory )                 #
#           $(SMS++OBJ)    = the core SMS++ library                          #
#           $(BNDSLVSDR)   = the directory where the source is               #
#                                                                            #
#   Output: $(BNDSLVOBJ)   = the final object(s) / library                   #
#           $(BNDSLVH)     = the .h files to include                         #
#           $(BNDSLVINC)   = the -I$( source directory )                     #
#                                                                            #
#                                VERSION 1.01                                #
#                               30 - 12 - 2020                               #
#                                                                            #
#                              Antonio Frangioni                             #
#                               Enrico Gorgone                               #
#                          Operations Research Group                         #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

# macros to be exported - - - - - - - - - - - - - - - - - - - - - - - - - - -

BNDSLVOBJ = $(BNDSLVSDR)BundleSolver.o 

BNDSLVINC = -I$(BNDSLVSDR)

BNDSLVH   = $(BNDSLVSDR)BundleSolver.h 

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(BNDSLVOBJ) $(BNDSLVSDR)*~

# dependencies: every .o from its .cpp + every recursively included .h- - - -

$(BNDSLVSDR)BundleSolver.o: $(BNDSLVSDR)BundleSolver.cpp $(BNDSLVH) \
	$(SMS++OBJ) $(libNDOOBJ)
	$(CC) -c $*.cpp -o $@ $(BNDSLVINC) $(SMS++INC) $(libNDOINC) $(SW)

########################## End of makefile ###################################
