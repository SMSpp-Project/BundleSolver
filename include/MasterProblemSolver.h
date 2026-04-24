/*--------------------------------------------------------------------------*/
/*-------------------- File MasterProblemBlock.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the MPBlock class, which implements the master problem
 * interface, using a "Generalized Bundle" algorithm for the solution of 
 * convex nondifferentiable problems.
 *
 * This interface is tought to be always used with the BundleSolver class.
 * TBD
 *
 * This class represents the master problem block within the SMS++ 
 * decomposition framework. The current version provides only the 
 * structural skeleton and does not implement any functionality.
 *
 * The user is assumed to be familiar with the algorithm: refer to
 *
 *  A. Frangioni "Generalized Bundle Methods"
 *  SIAM Journal on Optimization 13(1), 117--156, 2002
 *
 * available at
 *
 * \link
 *  http://www.di.unipi.it/~frangio/abstracts.html#SIOPT02
 * \endlink
 *
 * or
 *
 *  A. Frangioni "Standard Bundle Methods: Untrusted Models and Duality"
 *  in Numerical Nonsmooth Optimization: State of the Art Algorithms,
 *  A.M. Bagirov, M. Gaudioso, N. Karmitsa, M. Mäkelä, S. Taheri (Eds.),
 *  61 - 116, Springer, 2020
 *
 * available at
 *
 * \link
 *  http://www.di.unipi.it/~frangio/abstracts.html#NDOB18
 * \endlink
 *
 * In particular, BundleSolver implements the Incremental version of the
 * (Generalised) Proximal Bundle approach using upper models (for all the
 * components that provide a Lipschitz constant) described in
 *
 *  W. van Ackooij, A. Frangioni "Incremental Bundle Methods Using Upper
 *  Models" SIAM Journal on Optimization 28(1), 379 – 410, 2018
 *
 * available at
 *
 * \link
 *  http://www.di.unipi.it/~frangio/abstracts.html#SIOPT16
 * \endlink
 *
 * BundleSolver is capable of solving any Block such that:
 *
 * - only has "continuous" ColVariable (is_integer() == false);
 *
 * - has no Constraint, except possibly BoxConstraint (upper and/or lower
 *   bounds on the ColVariable);
 *
 * - either its Objective is a FRealObjective containing a C05Function, and
 *   it has no sub-Block;
 *
 * - or its Objective is a is a FRealObjective containing a LinearFunction,
 *   and each of its sub-Block has no Constraint and Variable, and its
 *   Objective is a FRealObjective containing a C05Function.
 *
 * A special treatment is given to the case where some of the C05Function
 * actually are LagBFunction whose inner Block only contains ColVariable,
 * whose Objective is linear (a FRealObjective containing a LinearFunction)
 * and whose Constraint are linear (either FRowConstraint containing a
 * LinearFunction, or BoxConstraint). These can be passed to the Master
 * Problem of the bundle algorithm as "easy components", see
 *
 *   A. Frangioni, E. Gorgone "Generalized Bundle Methods for Sum-Functions
 *   with ``Easy'' Components: Applications to Multicommodity Network Design"
 *   Mathematical Programming 145(1), 133 – 161, 2014
 *
 * available at
 *
 * \link
 *  http://www.di.unipi.it/~frangio/abstracts.html#MP11c
 * \endlink
 *
 * In that case, the LagBFunction is never evaluated, which means that there
 * is no need for a Solver to be attached to the inner Block.
 *
 * If the Block has multiple Objective (that is, it has sub-Block whose
 * Objective is a FRealObjective containing a C05Function), a very strong
 * assumption is required on them:
 *
 *     ALL THE Function IN THE Objective HAVE EXACTLY THE SAME SET OF
 *     "ACTIVE" Variable, ORDERED IN THE SAME WAY, AT ALL TIMES; THIS MEANS
 *     THAT IF THE SET OF "ACTIVE" Variable IS MODIFIED FOR ONE OF THE
 *     Function, IT MUST BE MODIFIED FOR ALL OF THEM AT THE SAME TIME
 *
 * The only exception is that
 *
 *     THE Objective OF THE Block CAN BE EMPTY, I.E., EITHER THERE IS NO
 *     Objective, OR THE FRealObjective HAS NO Function, OR THE
 *     LinearFunction IN THE FRealObjective HAS EXACTLY ZERO "ACTIVE"
 *     Variable; IN THE LATTER CASE, THE SET OF "ACTIVE" Variable IN THE
 *     LinearFunction MUST NEVER CHANGE
 *
 * To ensure that the rule about the list of "ACTIVE" Variable in the
 * (multiple) Objective is respected, an analogous very strong assumption is
 * made on the Modification that change that:
 *
 *     ALL THE Modification THAT CHANGE THE "ACTIVE" Variable MUST BE
 *     BUNCHED TOGETHER IN A SINGLE GroupModification. THIS MUST CONTAIN
 *     EXACTLY AS MANY Modification AS THERE ARE sub-Block (AND, THEREFORE,
 *     DIFFERENT OBJECTIVE), PLUS ONE IF THE (LinearFunction IN THE)
 *     Objective OF THE Block IS NOT EMPTY. ALL Modification MUST BE OF
 *     THE VERY SAME TYPE, I.E., EITHER ALL C05FunctionModVarsAddd, OR ALL
 *     C05FunctionModVarsRngd, OR ALL C05FunctionModVarsSbst, AND THEY MUST
 *     CHANGE THE "ACTIVE" Variable IN PRECISELY THE SAME WAY.
 *
 * \author Enrico Calandrini \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Enrico Calandrini, Enrico Gorgone
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __MasterProblemBlock
 #define __MasterProblemBlock

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

//TBD

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE & USING -----------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/
/*---------------------- CLASS MasterProblemBlock --------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
//TBD

 class MasterProblemBlock : public Block {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
 *
 * "Import" basic types from Function, C05Function.
 *
*  @{ */
/*----------------------------- CONSTANTS ----------------------------------*/


/*--------------------------------------------------------------------------*/
 /// public enum for the int algorithmic parameters TBD
 /** Public enum describing the different algorithmic parameters of int type
  * that MasterProblemBlock has in addition to these of ?. */

 enum int_par_type_MPBlock {

 };  // end( int_par_type_MPBlock )

/*--------------------------------------------------------------------------*/
 /// public enum for the double algorithmic parameters YBD
 /** Public enum describing the different algorithmic parameters of double
  * type that MasterProblemBlock has in addition to these of ?. */

 enum dbl_par_type_MPBlock {
  
 };  // end( dbl_par_type_MPBlock )

/*--------------------------------------------------------------------------*/
 /// public enum for the string algorithmic parameters
 /** Public enum describing the different algorithmic parameters of string
  * type that MasterProblemBlock has in addition to these of ?. */

 enum str_par_type_MPBlock {

  };  // end( str_par_type_MPBlock )

/*--------------------------------------------------------------------------*/
 /// public enum for the vector-of-int parameters
 /** Public enum describing the different algorithmic parameters of
  * vector-of-int type that MasterProblemBlock has in addition to these of ?. */

 enum vint_par_type_MPBlock {

  };  // end( vint_par_type_MPBlock )

/*--------------------------------------------------------------------------*/
 /// public enum for the vector-of-string parameters
 /** Public enum describing the different parameters of vector-of-string type
  * that MasterProblemBlock has in addition to these of ?. */

 enum vstr_par_type_MPBlock {

  };  // end( vstr_par_type_MPBlock )

/** @} ---------------------------------------------------------------------*/
/*----------- CONSTRUCTING AND DESTRUCTING MasterProblemBlock --------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing MasterProblemBlock
 *  @{ */

 /// constructor: ensure every field is initialized

 /** Default constructor */
  MasterProblemBlock( void ) : NoEasyCmps( 0 ) , NoHardCmps( 0 )
  { }

/*--------------------------------------------------------------------------*/
 /// destructor

 virtual ~MasterProblemBlock( ) {  }

/*--------------------------------------------------------------------------*/
 /// destructor  
 void clear();

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations TBD
 *
 *  @{ */

/** Provides the MPBlock with basic information about the "size" of the
 * Master problem:
 *
 * - MxBSz is the maximum number of different items (subgradients and
 *   constraints) that can be managed by the class, i.e. the maximum size of
 *   the bundle: the protected field MaxBSize is offered by  the base class
 *   to store this information;
 *
 * - NVars is the number of variables in in the Primal Master Problem (PMP).
 *   The protected field NumVars is offered by the base class to store this 
 *   information. \note in the current implementation the number of variables 
 *   canno't change. Hence, this number is also equal to the maximum length of 
 *   any item.
 *
 * - NrFi is the number of "components". The function Fi to be 
 *   minimized may be "decomposable", i.e., the sum of k functions; 
 *   in this case, its model Fi_{B,Lambda} is usually decomposable as well 
 *   (for instance, the Cutting Plane model is). For decomposable functions, 
 *   it is possible (although not necessary) to keep k "separate bundles", 
 *   each one describing a model Fi[ k ]_{B,Lambda} of each component Fi[ k ] 
 *   of Fi, which may be beneficial to improve the "exactness" of the model.
 *   Special support is offered for the case when one of the components is a 
 *   linear function, by considering it as the "0-th component" of Fi [see
 *   ReadFiBLambda() and SetItem() and GetItem()].
 * 
 * - NrFiEasy is the number of "easy components". A special treatment is 
 *   given to the case where some of the components actually are Lagrangian 
 *   subproblem where the domain is an "easy" polyhedron", see
 *
 *      A. Frangioni, E. Gorgone "Generalized Bundle Methods for Sum-Functions
 *      with ``Easy'' Components: Applications to Multicommodity Network Design"
 *      Mathematical Programming 145(1), 133 – 161, 2014
 * 
 *   available at
 *
 *     \link
 *     http://www.di.unipi.it/~frangio/abstracts.html#MP11c
 *     \endlink
 *  
 *   for further details.
 * 
 * TBD from here
 * \note "easy" components of Fi() have basically to be dealt with by the
 *   MPSolver, by inserting their description in the Master Problem;
 *	 if this is not possible, the Master Problem has to signal it (e.g.
 *	 by throwing an exception) because the Bundle algorithm relies on
 *	 this and the FiOracle is not going to provide "ordinary" black-box
 *	 information (function values, subgradients, ...) for these
 *	 components.
 *
 * \note *Important*: "easy" components of Fi() in the Master Problem are
 *       *not* translated *by value* (but they are by argument), meaning that
 *	 the function that is included in the MP is
 *
 *	   Fi_{Lambda}[ k ]( d ) = Fi[ k ]( Lambda + d )
 *
 *	 and *not*
 *
 *	   Fi_{Lambda}[ k ]( d ) = Fi[ k ]( Lambda + d ) - Fi[ k ]( Lambda )
 *
 *	 as one would expect by analogy with the ordinary "difficult"
 *	 components. This is done to spare the Bundle algorithm with the
 *	 need to compute Fi[ k ]( Lambda ) for all "easy" components k and
 *	 each current point Lambda, which may be problematic especially for
 *	 the *initial* current point---consider the case where Lambda *is
 *	 not feasible*, so Fi[ k ]( Lambda ) = +INF! However, the value of
 *	 Fi[ k ] is actually computed by the MPSolver at Lambda + d* (d*
 *	 being the optimal solution of the Primal MP), so the value of
 *	 Fi[ k ] is known for the current point at least after every Serious
 *	 Step with step 1 along d*. Yet, this choice has some impact on the
 *	 "output" methods [see ReadFiBLambda() and ReadSigma() below].
 *
 * This method can be called more than once to modify the settings, but expect
 * the implementation to be quite expensive in time and/or memory. There are
 *  essentially three different ways for calling SetDim():
 *
 * - SetDim( 0 , ... ) makes the MPSolver to deallocate all its memory and to
 *   quietly wait for new instructions.
 *
 * - SetDim( n , 0 , ASV ) with n != 0 sets the max bundle size to n and
 *   activate/deactivate the Active Set Mechanism without changing anything
 *   else. The existing items in the bundle (if any) are all kept if n is
 *   larger than the previous setting, but a smaller value will force deletion
 *   of all the items with "name" [see [Get/Set]Item() and RmvItem() below]
 *   greater than or equal to n. Also, if the active set technique is
 *   initialized (ASV == true while it was false previously), the set of
 *   "active" variables is set to be *empty*.
 *
 * - SetDim( n , Orcl , ASV ) with n != 0 and Orcl != 0 discards all the
 *   previous settings and re-allocates everything; all the existing items in
 *   the bundle are lost. Note that such a call also resets every parameter
 *   of the algorithm, such as the starting point (which is set to 0).
 *
 * In general, calling this method with a non-empty bundle could be costly, so
 * if the items are to be discarded anyway, this should be done *before* the
 * call to SetDim() [see RmvItems() below]. */

 void SetDim( int MxBSz , int NVars , int NrFi , int NrFiEasy );

/*--------------------------------------------------------------------------*/
 /// set the Solver of MasterProblemBlock
 /** This method let MasterProblemBlock register its own Solver. This method
  *  expect that solv_cfg_filename corresponds to an appropriate 
  *  BlockSolverConfig *. If empty, MasterProblemBlock will attempt to use 
  *  a default Solver (i.e. GRBMILPSolver).
  * 
  * \note: if Gurobi is not properly installed on the private machine and 
  *        no Solver is provided, this will end up in an error.
  * 
  */
void register_Solver( std::string && solv_cfg_filename );

/*--------------------------------------------------------------------------*/
 /// prepare the default Solver of MasterProblemBlock
 /** This method is used only when a BlockSolverConfiguration is not 
  *  provided for setting the inner Solver of the class. The method outputs
  *  a simple BlockSolverConfig specifying GRBMILPSolver as inner Solver
  *  of the class.
  * 
  * \note: if Gurobi is not properly installed on the private machine and 
  *        no Solver is provided, this will end up in an error.
  * 
  */
BlockSolverConfig * use_default_Solver( void );

/*--------------------------------------------------------------------------*/
 /// set the int parameters of MasterProblemBlock
 /** Set the int parameters specific of MasterProblemBlock, together with the
  * parameters of ? that MasterProblemBlock actually "listens to":
  * 
  */

 void set_par( idx_type par , int value ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// set the double parameters of MasterProblemBlock
 /** Set the double parameters specific of MasterProblemBlock, together with the
  * parameters of ? that MasterProblemBlock actually "listens to":
  * 
  */

 void set_par( idx_type par , double value ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// move the string parameters of MasterProblemBlock
 /** Move in the string parameters specific of MasterProblemBlock, together with
  * the parameters of ? that MasterProblemBlock actually "listens to"
  * 
 */
  
 void set_par( idx_type par , std::string && value ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// move in the vector-of-int parameters of MasterProblemBlock
 /** Move in the given vector-of-int parameters specific of MasterProblemBlock,
  * together with the parameters of ? that MasterProblemBlock actually
  * "listens to"
  * 
  */
  
 void set_par( idx_type par , std::string && value ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// move in the vector-of-int parameters of MasterProblemBlock
 /** Move in the given vector-of-int parameters specific of MasterProblemBlock,
  * together with the parameters of ? that MasterProblemBlock actually
  * "listens to"
  * 
  */

 void set_par( idx_type par , std::vector< int > && value ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// move in the vector-of-string parameters of MasterProblemBlock
 /** Move in the given vector-of-string parameters specific of MasterProblemBlock,
  * together with the parameters of ? that MasterProblemBlock actually
  * "listens to"
  * 
  */

 void set_par( idx_type par , std::vector< std::string > && value ) override;

/*--------------------------------------------------------------------------*/
 /// set the ostream for the MasterProblemBlock log

 void set_log( std::ostream *log_stream = nullptr ) override;

/** @} ---------------------------------------------------------------------*/
/*----------------- METHODS FOR ACCESSING THE DATA OF THE Block ------------*/
/*--------------------------------------------------------------------------*/
/** @name Accessing the data of the Block TBD
 *
 * These methods provide convenient shortcuts for directly asking to the
 * MasterProblemBlock some relevant data about the Block it is solving.
 *  @{ */

/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the parameters of the MasterProblemBlock TBD
 *
 *  @{ */

 [[nodiscard]] idx_type get_num_int_par( void ) const override {
  return(  );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type get_num_dbl_par( void ) const override {
  return(  );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type get_num_str_par( void ) const override {
  return(  );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type get_num_vint_par( void ) const override {
  return(  );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type get_num_vstr_par( void ) const override {
  return( );
  }

/*--------------------------------------------------------------------------*/
 
 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override {

  return(  );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
 [[nodiscard]] double get_dflt_dbl_par( idx_type par ) const override {

  return(  );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const std::string & get_dflt_str_par( idx_type par ) const override {
  
  return(  );
  }

/*--------------------------------------------------------------------------*/
 
 [[nodiscard]] int get_int_par( idx_type par ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
 [[nodiscard]] double get_dbl_par( idx_type par ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
 [[nodiscard]] const std::string & get_str_par( idx_type par )
  const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
 [[nodiscard]] const std::vector< int > & get_vint_par( idx_type par )
  const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
 [[nodiscard]] const std::vector< std::string > & get_vstr_par(
					       idx_type par ) const override;

/*--------------------------------------------------------------------------*/

 [[nodiscard]] idx_type int_par_str2idx( const std::string & name )
  const override {

   return(  );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type dbl_par_str2idx( const std::string & name )
  const override {

  return(  );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type str_par_str2idx( const std::string & name )
  const override {

  return(  );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type vint_par_str2idx( const std::string & name )
  const override {

  return(  );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type vstr_par_str2idx( const std::string & name )
  const override {

  return(  );
  }

/*--------------------------------------------------------------------------*/

 [[nodiscard]] const std::string & int_par_idx2str( idx_type idx )
  const override {

  return(  );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::string & dbl_par_idx2str( idx_type idx )
  const override {

  return(  );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::string & str_par_idx2str( idx_type idx )
  const override {

  return(  );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::string & vint_par_idx2str( idx_type idx )
  const override {

  return(  );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::string & vstr_par_idx2str( idx_type idx )
  const override {

  return(  );
  }

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED TYPES ------------------------------*/
/*--------------------------------------------------------------------------*/

// TBD

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

// TBD

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS  ---------------------------*/
/*--------------------------------------------------------------------------*/

 // algorithmic parameters - - - - - - - - - - - - - - - - - - - - - - - - - -

 int MaxBSize; // maximum bundle size

 int MaxSGLen; // maximum subgradient length

 /* This is the Solver that will be responsable of reading each Sub-Block and 
  * provide the physical representation of the problem (e.g. any *MILPSolver).
 */
 CDASolver InnerSolver; 

 int NoEasyCmps; // Number of easy components

 int NoHardCmps; // Number of easy components

 /* List of sub-blocks identified to be easy components.
 */
 std::vector< Block* > EasyCmps;

 /* List of sub-blocks to be treated as hard components.
 */
 std::vector< Block* > HardCmps; 

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ PRIVATE FIELDS  ---------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class MasterProblemBlock )

} // namespace SMSpp

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* MasterProblemBlock.h included */

/*--------------------------------------------------------------------------*/
/*--------------------- End File MasterProblemBlock.h ----------------------*/
/*--------------------------------------------------------------------------*/
