# --------------------------------------------------------------------------- #
#    Dockerfile for CI/CD                                                     #
#                                                                             #
#    This file contains the commands to build a Docker image containing       #
#    all the packages needed to build and test the project.                   #
#    Once built and uploaded in the repository's container registry           #
#    (See: https://gitlab.com/smspp/bundlesolver/container_registry),         #
#    the image can be fetched and used by the GitLab Runner.                  #
#                                                                             #
#    Note: Once built, this image will contain parts of the                   #
#          IBM ILOG CPLEX Optimization Studio suite. For this reason,         #
#          the image cannot be publicly distributed and can be used only      #
#          by team members that own a valid license for that software.        #
#    ---------------------------------------------------------------------    #
#    DISCLAIMER: The author of this file is not affiliated, associated,       #
#    authorized, endorsed by, or in any way officially connected with IBM,    #
#    or any of its subsidiaries or its affiliates. The names IBM, ILOG and    #
#    CPLEX as well as related names, marks, emblems and images are            #
#    registered trademarks of their respective owners.                        #
#    ---------------------------------------------------------------------    #
#                                                                             #
#    Build the image with:                                                    #
#                                                                             #
#        $ docker build -t registry.gitlab.com/smspp/bundlesolver .           #
#                                                                             #
#    Upload with:                                                             #
#                                                                             #
#        $ docker push registry.gitlab.com/smspp/bundlesolver                 #
#                                                                             #
#    Run (locally) with:                                                      #
#                                                                             #
#        $ docker run --rm -it registry.gitlab.com/smspp/bundlesolver:latest  #
#                                                                             #
#    Note: you need to rebuild and upload the image only when this file       #
#          changes, not when BundleSolver changes.                            #
#                                                                             #
#                              Niccolo' Iardella                              #
#                          Operations Research Group                          #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #

# An image containing everything needed to build MILPSolver.
# Note: CPLEX is inside this image.
FROM registry.gitlab.com/smspp/milpsolver

## Install required packages
RUN set -ex; \
		apt update; \
		apt install -y --no-install-recommends libbz2-dev; \
		rm -rf /var/lib/apt/lists/*;

# Install Coinutils
ARG COINUTILS_VER="2.11.4"
RUN set -ex; \
		mkdir -p "coinutils"; \
		curl -SsL "https://github.com/coin-or/CoinUtils/archive/releases/$COINUTILS_VER.tar.gz" | \
		tar -xzC "coinutils" --strip-components 1; \
		cd "coinutils"; \
		./configure --prefix=/usr/local --enable-static; \
		make; \
		make install; \
		cd .. && rm -r "coinutils";

# Install Osi
ARG OSI_VER="0.108.6"
RUN set -ex; \
		CPLEX_DIR=`ls -bd1 /opt/ibm/ILOG/CPLEX_Studio* | tail -n1`; \
		CPLEX_LIB_DIR=`ls -bd1 $CPLEX_DIR/cplex/lib/*/static_pic | tail -n1`; \
		mkdir -p "osi"; \
		curl -SsL "https://github.com/coin-or/Osi/archive/releases/$OSI_VER.tar.gz" | \
		tar -xzC "osi" --strip-components 1; \
		cd "osi"; \
		./configure --prefix=/usr/local --enable-static \
			--with-cplex-incdir="$CPLEX_DIR/cplex/include/ilcplex" \
			--with-cplex-lib="-L$CPLEX_LIB_DIR -lcplex -lpthread -lm -ldl"; \
		make; \
		make install; \
		cd .. && rm -r "osi";

# Install Clp
ARG CLP_VER="1.17.6"
RUN set -ex; \
		mkdir -p "clp"; \
		curl -SsL "https://github.com/coin-or/Clp/archive/releases/$CLP_VER.tar.gz" | \
		tar -xzC "clp" --strip-components 1; \
		cd "clp"; \
		./configure --prefix=/usr/local --enable-static; \
		make; \
		make install; \
		cd .. && rm -r "clp";
