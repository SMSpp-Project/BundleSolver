/*--------------------------------------------------------------------------*/
/*-------------------- File MasterProblemBlock.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the MPBlock class, implementing the master problem
 * interface for the BundleSolver class.
 * 
 * This class implements and solves the Master Problem for a generic Bundle
 * Method, with the possibility of initializing the problem either in its
 * primal or dual formulation.
 * 
 * It is particularly suitable for problems where the objective function is
 * expressed as the sum of independent components coupled through
 * (possibly) complicating constraints, although it can also handle the
 * case where only a single component is present. In the former setting,
 * each component is treated separately through a dedicated sub-block
 * reproducing the corresponding sub-problem, while the current class is
 * responsible for reproducing the coupling constraints. For a more
 * detailed overview of the general algorithm, please refer to
 * GeneralizedBundleSolver.h.
 * 
 * In its most general form, we assume that the objective function f can
 * be expressed as the sum of N sub-functions f^k, each of which is only
 * accessible through an oracle providing the function value and a
 * subgradient at a given point. Assuming that a bundle
 * B^k = {(g^k_i, α^k_i)}_i stores information related to subgradients and
 * linearization errors, each sub-problem can be written as:
 * 
 * \min v^k                                               (1)
 * s.t  v^k \ge g^k_i * d + α^k_i       \forall i         (2)
 *      α^k_h \ge g^k_h * d             \forall h         (3)
 *      v^k \ge LB^k                                      (4)
 * 
 * where constraints (2) and (3) are commonly referred to as subgradient
 * cuts and feasibility cuts, respectively, while constraint (4)
 * represents a possible known lower bound on the value of f^k. In the
 * formulation above, d denotes the step to be performed at the current
 * iteration and will be a responsability of MPB of initializing such
 * variable and providing it to each sub-block.
 * 
 * This formulation is naturally represented by the "primal
 * version" of a generic PolyhedralFunctionBlock. For this reason, one
 * such block is created and registered as a sub-block of the current
 * class for each sub-function of the original function f.
 * 
 * The lower approximation model of the complete function f is then
 * obtained by introducing an appropriate variable v and the constraint
 * 
 *  v \ge d * b + \sum_k v_k                               (5)
 * 
 * where b represents the constant part of the objective function f, when
 * present.
 * 
 * Finally, an appropriate stabilization mechanism must be considered in
 * the master problem. Currently, MPB supports:
 *  - Proximal Stabilization;
 *  - Level Stabilization;
 *  - Doubly Stabilized Bundle Methods.
 * 
 * The Master Problem obtained through the minimization of v together with
 * the selected stabilization term corresponds to the formulation
 * implemented in this class under the name "Primal Form".
 * 
 * In addition, MasterProblemBlock is also able to reproduce the dual
 * formulation of the model described above. This is particularly useful
 * when some of the sub-functions are so-called "easy components", i.e.,
 * components that can be evaluated exactly instead of relying on an
 * approximated model (see BundleSolver.h for further details).
 * 
 * In this case, each sub-block introduces dual multipliers θ^k_i for
 * each cut in the bundle B^k and a multiplier \gamma^k associated with
 * the lower bound LB^k. The corresponding contribution
 * 
 * - \sum_i \theta^k_i \alpha^k_i - \sum_h \theta^k_h \alpha^k_h
 * 
 * is then added to the objective function, together with the
 * "normalization constraints"
 * 
 * \sum_i \theta^k_i + \gamma^k = \lambda
 * 
 * where \lambda is the dual multiplier associated with constraint (5).
 * 
 * MasterProblemBlock is instead responsible for coupling all terms
 * \theta^k_i g^k_i coming from each sub-block into a single vector R,
 * which is then used in the objective function of the maximization
 * problem.
 * 
 * As in the primal case, specific approaches are adopted depending on
 * the stabilization mechanism being used.
 * 
 * Also in this setting, the dual formulation described above is naturally
 * represented by the "dual version" of a generic
 * PolyhedralFunctionBlock. Hence, one such block is created for each
 * component, while specialized approaches are adopted for the treatment
 * of the "easy components".
 * 
 * NOTE: for a more extensive description of the specific master problems
 * constructed in both their primal and dual formulations, please refer to
 * the accompanying notes.
 * 
 * TODO: possible link to the related notes/documentation.
 * 
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

#include "Block.h"

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
/// A generic Block containing the Master Problem generated by a BundleSolver
/** The MasterProblemBlock class implements the interface within
 *  the SMS++ framework for representing and solving a generic master problem
 *  for the BundleSolver class.
 * 
 * */

 class MasterProblemBlock : public Block {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
 * @{ */

 /// Definition of the possible type of stabilization
 /** The enum stabilization_type defines all possible "types" of stabilization
  *  that can be inserted into the Master Problem. Currently, 
  *  MasterProblemBlock supports 3 types of stabilization:
  *
  *  - Proximal;
  *
  *  - Level;
  *
  *  - Doubly-Stabilized. */

 typedef enum {
  kProximal         =  0 ,  ///< proximal
  kLevel            =  1 ,  ///< level
  kDoublyStabilized =  2    ///< doubly-stabilized
  } stabilization_type;

/*----------------------------- CONSTANTS ----------------------------------*/


/** @} ---------------------------------------------------------------------*/
/*----------- CONSTRUCTING AND DESTRUCTING MasterProblemBlock --------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing MasterProblemBlock
 *  @{ */

 /// constructor: ensure every field is initialized

 /** Default constructor */
  MasterProblemBlock( void ) : NoEasyCmps( 0 ) , NoHardCmps( 0 ) , 
    StblType( kDoublyStabilized ) , 
  { }

/*--------------------------------------------------------------------------*/
 /// destructor: cleanly detaches the MasterProblemBlock from the Block

 virtual ~MasterProblemBlock( ) { }

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
 *   MPBlock, by inserting their description in the Master Problem;
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
 *	 Fi[ k ] is actually computed by the MPBlock at Lambda + d* (d*
 *	 being the optimal solution of the Primal MP), so the value of
 *	 Fi[ k ] is known for the current point at least after every Serious
 *	 Step with step 1 along d*. Yet, this choice has some impact on the
 *	 "output" methods [see ReadFiBLambda() and ReadSigma() below].
 *
 * This method can be called more than once to modify the settings, but expect
 * the implementation to be quite expensive in time and/or memory. There are
 *  essentially three different ways for calling SetDim():
 *
 * - SetDim( 0 , ... ) makes the MPBlock to deallocate all its memory and to
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
 /// Initialize the primal version of the Master Problem
 /** This method TBD
  * 
  */
void CreatePrimalMP( stabilization_type Stbl ,  );

/*--------------------------------------------------------------------------*/
 /// Initialize the dual version of the Master Problem
 /** This method TBD
  * 
  */
void CreateDualMP( stabilization_type Stbl ,  );

/*--------------------------------------------------------------------------*/
 /// Load the probelm into the Block
 /** This method tells the registered solver to load the master problem. 
  *  If the problem type is 0, then a simple Solver->load_problem() is enough
  *  to upload information from the unique sub-block. Otherwise, we have to
  *  distinguish on where information is coming from (i.e. an hard or easy
  *  component). TBD
  * 
  */
void load_problem( void );

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

 bool IsPrimal; // wether we are solving the primal or dual version of the MP

 stabilization_type StblType; // type of stabilization

 int MaxBSize; // maximum bundle size

 int MaxSGLen; // maximum subgradient length

 int NoTotCmps; // Number of total components

 int NoEasyCmps; // Number of easy components

 std::vector< Bool > IsEasyCmp; // which component is easy

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
