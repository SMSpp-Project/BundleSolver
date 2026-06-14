/*--------------------------------------------------------------------------*/
/*--------------------- File BundleSolverML.h ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the BundleSolverML class, which derives from BundleSolver
 * and replaces the classical rule-based step-size heuristics with a small
 * feedforward neural network (implemented with the Torch C++ API) that
 * predicts the proximal regularization parameter t out of a set of features
 * describing the current state of the algorithm.
 *
 * The network can be trained online: after each compute() the caller can
 * invoke Backward(), which computes a loss out of the recorded iterations,
 * back-propagates it through the network and applies an Adam update to its
 * weights. A single network can be shared among multiple BundleSolverML
 * objects (one per instance of the problem) so that the gradient updates of
 * all of them accumulate into the same weights, which is the standard
 * pattern for multi-instance training. The weights can be saved to / loaded
 * from TorchScript archives that are cross-compatible with Python.
 *
 * \author Francesca Demelas \n
 *         Laboratoire d'Informatique de Paris Nord \n
 *         Universite' Sorbonne Paris Nord \n
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Francesca Demelas, Antonio Frangioni, Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __BundleSolverML
 #define __BundleSolverML
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <iostream>
#include <memory>
#include <string>

#include <torch/torch.h>
#include <torch/optim.h>

#include "BundleSolver.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- MACROS -----------------------------------*/
/*--------------------------------------------------------------------------*/
/* BML_LOG (so named to avoid clashing with Torch's c10 LOG( n ) macro)
 * prints debug output of BundleSolverML on std::cout; it is governed by the
 * VERBOSE macro, set to 1 to enable the output and to 0 to disable it. */

#ifndef VERBOSE
 #define VERBOSE 0
#endif

#if VERBOSE
 #define BML_LOG( x ) std::cout << x
#else
 #define BML_LOG( x )
#endif

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE & USING -----------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------------- CLASS Net ----------------------------------*/
/*--------------------------------------------------------------------------*/
/// two-layer feedforward neural network predicting the step-size t
/** Net is a simple two-layer feedforward neural network used by
 * BundleSolverML to predict the step-size t:
 *
 * - input layer: 20 features describing the current state of the bundle
 *   algorithm;
 *
 * - hidden layer: 16 units with Softplus activation;
 *
 * - output layer: 1 unit with Softplus activation, producing t > 0.
 *
 * The Softplus output ensures that the predicted step-size is always
 * strictly positive, as required for the convergence of the bundle method.
 *
 * A single Net object can be shared among multiple BundleSolverML objects
 * (one per instance of the problem) via BundleSolverML::set_shared_net() /
 * BundleSolverML::get_shared_net(). This is the recommended pattern for
 * training: all the solvers read from and write to the same parameter
 * tensors, so that every Backward() call accumulates gradients into the
 * shared weights.
 *
 * The weights are saved / loaded via BundleSolverML::SaveModel() and
 * BundleSolverML::LoadModel() in the TorchScript .pt format, which is
 * cross-compatible with Python's torch.load(). */

struct Net : torch::nn::Module
{
 /// first fully-connected layer: 20 inputs --> 16 hidden units
 torch::nn::Linear fc1{ nullptr };

 /// second fully-connected layer: 16 hidden units --> 1 output (step-size)
 torch::nn::Linear fc2{ nullptr };

 /// constructor: registers the two layers with the Torch module system
 /** Constructor. The output layer is initialized with small weights and a
  * bias such that softplus( bias ) is about 1, so that the first
  * predictions of the untrained network are close to the classical
  * default initial step-size t = 1 (instead of some arbitrary value
  * dictated by the random initialization, which can be pathological
  * enough to stall the bundle algorithm) and the training only has to
  * refine them. */
 Net( void ) {
  fc1 = register_module( "fc1" , torch::nn::Linear( 20 , 16 ) );
  fc2 = register_module( "fc2" , torch::nn::Linear( 16 , 1 ) );

  torch::NoGradGuard no_grad;
  fc2->weight *= 0.1;
  fc2->bias.fill_( 0.5413 );  // softplus( 0.5413 ) ~= 1
  }

 /// forward pass: maps a feature vector into a positive step-size
 /** Forward pass of the network: maps the input feature tensor \p x, of
  * shape { 20 } (or { batch , 20 }), into a scalar tensor representing the
  * predicted step-size t > 0. */
 torch::Tensor forward( torch::Tensor x ) {
  x = torch::softplus( fc1->forward( x ) );
  x = fc2->forward( x );
  return( torch::softplus( x ) + 1.0e-8 );  // ensure strictly positive output
  }

 };  // end( struct( Net ) )

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS BundleSolverML ----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// BundleSolver with a machine-learning heuristic for the step-size t
/** BundleSolverML derives from BundleSolver and replaces the classical
 * rule-based step-size heuristics (Heuristic1() ... Heuristic4()) with a
 * neural network (Net) that predicts the proximal regularization parameter
 * t out of a 20-dimensional feature vector built from the current state of
 * the algorithm.
 *
 * By default each BundleSolverML owns a private Net, allocated at
 * construction time. For training across multiple instances the caller
 * should rather create a single shared Net and inject it into every solver:
 *
 *     auto shared_net = std::make_shared< Net >();
 *
 *     // ... for each instance:
 *     auto bml = dynamic_cast< BundleSolverML * >( get_solver_for( i ) );
 *     bml->set_shared_net( shared_net );  // inject the shared weights
 *     bml->compute( false );
 *     bml->Backward();                    // update the shared weights
 *
 * After set_shared_net() all forward and backward passes operate on the
 * same parameter tensors, so the gradient updates of every instance
 * accumulate into the same network, which is exactly the online SGD
 * pattern needed for multi-instance training. To revert to the private
 * network call clear_shared_net(). The weights of whichever network is
 * currently active are saved / loaded with SaveModel() / LoadModel().
 *
 * Note that copy and move operations are deleted to protect the Torch
 * nn::Module internals, which manage reference-counted parameter tensors;
 * all the tensors needed for differentiation are stored in per-iteration
 * vectors (phi_vecs, w_vecs, Gs, ...) and consumed during Backward(). */

class BundleSolverML : public BundleSolver
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/

 using Index = Function::Index;

 typedef std::vector< double > VecMem;
 ///< a vector of floating-point values (e.g., a single iterate)

 using MemMultiVector = std::vector< VecMem >;
 ///< a collection of VecMem vectors (e.g., previous iterates)

/*--------------------------------------------------------------------------*/
/*------------------- CONSTRUCTING AND DESTRUCTING --------------------------*/
/*--------------------------------------------------------------------------*/

 /// constructor: allocates the private Net
 /** Constructor: it allocates the private network and makes nn point to
  * it; the feature vector is sized to the input dimension of Net. */

 BundleSolverML( void ) : BundleSolver() {
  f_owned_net = std::make_shared< Net >();
  nn = f_owned_net.get();         // point to the private network by default
  size_features = 20;
  features.resize( size_features );

  /* The BundleSolver constructor caches the parameter defaults before the
   * BundleSolverML overrides of get_dflt_*_par() are active, since virtual
   * dispatch does not reach the derived class during base construction;
   * hence, re-apply here the defaults that differ [see
   * get_dflt_int_par() / get_dflt_dbl_par()]. */
  set_par( intMnSSC , get_dflt_int_par( intMnSSC ) );
  set_par( intMnNSC , get_dflt_int_par( intMnNSC ) );
  set_par( inttSPar1 , get_dflt_int_par( inttSPar1 ) );
  set_par( dblmnIncr , get_dflt_dbl_par( dblmnIncr ) );
  set_par( dblmnDecr , get_dflt_dbl_par( dblmnDecr ) );
  }

/*--------------------------------------------------------------------------*/
 // copy and move constructor and assignment are deleted: copying or moving
 // would invalidate the parameter registry of the Torch nn::Module

 BundleSolverML( const BundleSolverML & ) = delete;

 BundleSolverML & operator=( const BundleSolverML & ) = delete;

 BundleSolverML( BundleSolverML && ) = delete;

 BundleSolverML & operator=( BundleSolverML && ) = delete;

/*--------------------------------------------------------------------------*/
 /// destructor: cleanly detaches the BundleSolverML from the Block
 /** Destructor. The shared_ptr to the shared network (if any) is released
  * here; the Net is destroyed only if no other BundleSolverML still holds
  * a reference to it. */

 virtual ~BundleSolverML() { set_Block( nullptr ); }

/*--------------------------------------------------------------------------*/
/*-------------------- SHARED NETWORK MANAGEMENT ----------------------------*/
/*--------------------------------------------------------------------------*/
 /// injects an externally-owned shared network into this solver
 /** Injects the externally-owned shared network \p net into this solver:
  * after this call all the forward passes (in Heuristic()) and backward
  * passes (in Backward()) operate on the parameters of \p net instead of
  * those of the private network allocated at construction time.
  *
  * Since multiple BundleSolverML objects can point to the same Net, the
  * gradient updates of all of them accumulate into the same tensors, which
  * enables true multi-instance online training without any weight-file I/O
  * between instances. The shared_ptr keeps the Net alive as long as at
  * least one solver holds a reference, regardless of the order in which
  * the solvers are destroyed.
  *
  * @param net  shared pointer to the network to be used, must not be
  *             nullptr (exception thrown otherwise) */

 void set_shared_net( std::shared_ptr< Net > net ) {
  if( ! net )
   throw( std::invalid_argument(
	     "BundleSolverML::set_shared_net: net must not be nullptr" ) );
  f_shared_net = std::move( net );
  nn = f_shared_net.get();    // redirect the raw pointer used internally
  }

/*--------------------------------------------------------------------------*/
 /// returns the currently active network as a shared_ptr
 /** Returns the currently active network: the shared one if
  * set_shared_net() has been called, the private one otherwise (still
  * usable as a shared_ptr, e.g., to pass it to another solver). The
  * returned pointer is never nullptr. */

 std::shared_ptr< Net > get_shared_net( void ) const {
  return( f_shared_net ? f_shared_net : f_owned_net );
  }

/*--------------------------------------------------------------------------*/
 /// reverts to the private network, releasing the shared reference
 /** Reverts to the private network: after this call nn points back to the
  * network allocated at construction time and the shared_ptr to the
  * external network is released (the Net is destroyed only if no other
  * solver holds a reference to it). */

 void clear_shared_net( void ) {
  f_shared_net.reset();
  nn = f_owned_net.get();
  }

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR HANDLING THE MODEL -----------------------*/
/*--------------------------------------------------------------------------*/
 /// saves the active network weights to a TorchScript archive (.pt)
 /** Saves the weights of the currently active network (shared or private)
  * to the TorchScript archive \p filepath, whose parent directory must
  * exist (the file is overwritten if present). The format is
  * cross-compatible with Python's torch.load(). An std::runtime_error is
  * thrown on any I/O or Torch error. */

 void SaveModel( const std::string & filepath );

/*--------------------------------------------------------------------------*/
 /// loads weights from a TorchScript archive into the active network
 /** Loads the weights found in the TorchScript archive \p filepath into
  * the currently active network (shared or private); the architecture
  * must match exactly (same layer names and tensor shapes). The network
  * is left in eval() mode; nn->train() is called automatically at the
  * next Heuristic() or Backward(). An std::runtime_error is thrown if
  * the file is missing or malformed, or the architecture mismatches. */

 void LoadModel( const std::string & filepath );

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE Block -----------------------*/
/*--------------------------------------------------------------------------*/
 /// solve the problem, optionally training the network online
 /** Runs BundleSolver::compute() and, if intMLTrainOnline is nonzero,
  * trains the network online on the just-completed solve by calling
  * Backward() and then ClearBuffers(). This makes online training fully
  * transparent to the driver: configuring a BundleSolverML with
  * intMLTrainOnline == 1 is enough, no ML-specific calls are needed.
  *
  * @param changedvars  forwarded to BundleSolver::compute()
  *
  * @return the status returned by BundleSolver::compute() */

 int compute( bool changedvars = true ) override;

/*--------------------------------------------------------------------------*/
/*----------------------- ML-SPECIFIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
 /// step-size heuristic overriding the base class rule
 /** Method overriding the rule-based short-term t-strategies of the base
  * BundleSolver. It is called at each iteration where the corresponding
  * bits of inttSPar1 are set [see get_dflt_int_par()]: it builds the
  * feature vector out of the current state of the algorithm, runs a
  * forward pass through the active network (shared or private) and
  * returns the predicted step-size t > 0. It also records all the tensors
  * (G, Q, alpha, w, ...) needed by the subsequent Backward() call.
  *
  * @param whch  unused, kept for interface compatibility with BundleSolver
  *
  * @return the predicted step-size t > 0 */

 HpNum Heuristic( Index whch ) override;

/*--------------------------------------------------------------------------*/
 /// computes the gradients and updates the active network parameters
 /** Iterates over all the recorded iterations, accumulates a scalar loss
  * (weighted by the step type and discounted over time), back-propagates
  * it through the active network and applies an Adam update to its
  * parameters. If a shared network is active the update affects all the
  * solvers sharing the same Net. This is meant to be called once per
  * solve, after compute() returns. */

 void Backward( void );

/*--------------------------------------------------------------------------*/
 /// clears all the per-solve tensor buffers
 /** Empties all the per-iteration tensor buffers recorded by Heuristic().
  * Call this between successive solves on the same BundleSolverML object
  * to free memory and avoid stale data; it is not needed when a fresh
  * BundleSolverML is created for each instance. */

 void ClearBuffers( void ) {
  phi_vecs.clear();
  w_vecs.clear();
  coeff_vecs.clear();
  Gs.clear();
  Gs_aggreg.clear();
  Qs.clear();
  alphaS.clear();
  thetaS.clear();
  phiS.clear();
  lastIndexS.clear();
  tS.clear();
  FiS.clear();
  }

/*--------------------------------------------------------------------------*/
 /// constructs the search direction tensor for iteration f with step t
 /** Wraps the custom autograd function computing the search direction so
  * that the computation is differentiable w.r.t. the output of the
  * network (used internally by Backward()).
  *
  * @param f  iteration index into the recorded tensor vectors
  *
  * @param t  step-size value
  *
  * @return a scalar tensor representing sum( w ), connected to the
  *         autograd graph */

 torch::Tensor w( size_t f , double t );

/*--------------------------------------------------------------------------*/
/*---------------------- PARAMETERS OF THE SOLVER ----------------------------*/
/*--------------------------------------------------------------------------*/
 /// public enum of the int parameters added by BundleSolverML
 /** Public enum of the additional int algorithmic parameters of
  * BundleSolverML, on top of those of BundleSolver:
  *
  * - intMLTrainOnline: if nonzero, at the end of each compute() the network
  *   is trained online, i.e. Backward() and then ClearBuffers() are called
  *   automatically. With the default 0 the solver only predicts t and the
  *   training, if any, is left to an external driver (as the standalone
  *   BundleSolverML tester does);
  *
  * - intMLSeed: if >= 0, torch::manual_seed() is called with this value and
  *   the (privately-owned) network is re-initialized, so that the initial
  *   weights, and hence the whole training trajectory, are reproducible;
  *   the default -1 leaves the global Torch RNG untouched;
  *
  * - intNTrainRounds: number of times the same instance should be (re-)solved
  *   for online training. It is not used by the solver itself (which cannot
  *   reload the Block), but it is a hint that the training driver can read to
  *   set its re-solve loop; the default is 1, i.e. a single solve. */

 enum int_par_type_BndSlvML {
  intMLTrainOnline = intLastBndSlvPar ,  ///< auto-train at end of compute()

  intMLSeed ,        ///< Torch manual seed (< 0 = leave RNG untouched)

  intNTrainRounds ,  ///< number of training re-solves (read by the driver)

  intLastBndSlvMLPar ///< first allowed new int parameter for derived classes
  };

/*--------------------------------------------------------------------------*/
 /// get the number of int parameters

 [[nodiscard]] idx_type get_num_int_par( void ) const override {
  return( idx_type( intLastBndSlvMLPar ) );
  }

/*--------------------------------------------------------------------------*/
 /// get the int parameter par

 [[nodiscard]] int get_int_par( idx_type par ) const override {
  switch( par ) {
   case( intMLTrainOnline ): return( f_train_online );
   case( intMLSeed ):        return( f_ML_seed );
   case( intNTrainRounds ):  return( f_n_train_rounds );
   default:                  return( BundleSolver::get_int_par( par ) );
   }
  }

/*--------------------------------------------------------------------------*/
 /// set the int parameter par
 /** Sets the BundleSolverML-specific int parameters and delegates all the
  * others to BundleSolver. Setting intMLSeed to a nonnegative value seeds
  * the global Torch RNG and re-initializes the privately-owned network. */

 using BundleSolver::set_par;  // keep the other set_par() overloads visible

 void set_par( idx_type par , int value ) override {
  switch( par ) {
   case( intMLTrainOnline ):
    f_train_online = value;
    break;
   case( intMLSeed ):
    f_ML_seed = value;
    if( value >= 0 ) {
     torch::manual_seed( value );
     reset_net();  // re-init the owned net so the weights are reproducible
     }
    break;
   case( intNTrainRounds ):
    if( value < 1 )
     throw( std::invalid_argument(
		   "BundleSolverML::set_par: intNTrainRounds must be >= 1" ) );
    f_n_train_rounds = value;
    break;
   default:
    BundleSolver::set_par( par , value );
   }
  }

/*--------------------------------------------------------------------------*/
 /// string-to-index map for the BundleSolverML int parameters

 [[nodiscard]] idx_type int_par_str2idx( const std::string & name )
  const override {
  if( name == "intMLTrainOnline" ) return( intMLTrainOnline );
  if( name == "intMLSeed" )        return( intMLSeed );
  if( name == "intNTrainRounds" )  return( intNTrainRounds );

  return( BundleSolver::int_par_str2idx( name ) );
  }

/*--------------------------------------------------------------------------*/
 /// index-to-string map for the BundleSolverML int parameters

 [[nodiscard]] const std::string & int_par_idx2str( idx_type idx )
  const override {
  static const std::array< std::string , 3 > int_pars_str_ML = {
   "intMLTrainOnline" , "intMLSeed" , "intNTrainRounds" };

  if( ( idx >= intLastBndSlvPar ) && ( idx < intLastBndSlvMLPar ) )
   return( int_pars_str_ML[ idx - intLastBndSlvPar ] );

  return( BundleSolver::int_par_idx2str( idx ) );
  }

/*--------------------------------------------------------------------------*/
 /// get the default value of an int parameter
 /** Returns the default value of the int parameter \p par. The defaults
  * differ from the BundleSolver ones where needed to ensure that the
  * network actually drives the step-size:
  *
  * - inttSPar1 == 3: Heuristic() is called after every SS (bit 0) and
  *   after every NS (bit 1), with no long-term t-strategy interfering;
  *
  * - intMnSSC == intMnNSC == 0: t is allowed to change at every
  *   iteration, for otherwise the [ tm , tp ] window collapses onto
  *   { t } and the prediction of the network is computed but discarded. */

 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override {
  static const std::array< int , 22 > dflt_int_par = {
   10 ,   // intBPar1
   100 ,  // intBPar2
   1 ,    // intBPar3
   1 ,    // intBPar4
   0 ,    // intBPar6
   3 ,    // intBPar7
   0 ,    // intMnSSC
   0 ,    // intMnNSC
   3 ,    // inttSPar1
   2 ,    // intMaxNrEvls
   1 ,    // intDoEasy
   2 ,    // intWZNorm
   0 ,    // intFrcLstSS
   0 ,    // intTrgtMng
   0 ,    // intMPName
   0 ,    // intMPlvl
   0 ,    // intQPmp1
   0 ,    // intQPmp2
   4 ,    // intOSImp1
   0 ,    // intOSImp2
   1 ,    // intOSImp3
   2      // intRstAlg
   };

  static const std::array< int , 3 > dflt_int_par_ML = {
   0 ,    // intMLTrainOnline
   -1 ,   // intMLSeed
   1      // intNTrainRounds
   };

  if( ( par >= intLastBndSlvPar ) && ( par < intLastBndSlvMLPar ) )
   return( dflt_int_par_ML[ par - intLastBndSlvPar ] );

  if( ( par >= intLastParCDAS ) && ( par < intLastBndSlvPar ) )
   return( dflt_int_par[ par - intLastParCDAS ] );

  return( CDASolver::get_dflt_int_par( par ) );
  }

/*--------------------------------------------------------------------------*/
 /// get the default value of a double parameter
 /** Returns the default value of the double parameter \p par. Same as
  * BundleSolver except dblmnIncr and dblmnDecr, both (essentially) == 1:
  * the minimum "significant" change of t is none, so that after a SS the
  * network can pick any t in [ t , t * mxIncr ] and after a NS any t in
  * [ t * mxDecr , t ] (including keeping it essentially unchanged)
  * instead of being forced to move it. Note that mnIncr is required by
  * BundleSolver::set_par() to be strictly > 1, hence the tiny offset. */

 [[nodiscard]] double get_dflt_dbl_par( idx_type par ) const override {
  static const std::array< double , 17 > dflt_dbl_par = {
   0 ,         // dblNZEps
   1e+2 ,      // dbltStar
   0 ,         // dblMinNrEvls
   30 ,        // dblBPar5
   0.01 ,      // dblm1
   0.99 ,      // dblm2
   0.99 ,      // dblm3
   10 ,        // dblmxIncr
   1.000001 ,  // dblmnIncr
   0.1 ,       // dblmxDecr
   1 ,         // dblmnDecr
   1e+6 ,      // dbltMaior
   1e-6 ,      // dbltMinor
   1 ,         // dbltInit
   1e-3 ,      // dbltSPar2
   0 ,         // dbltSPar3
   1e-1        // dblCtOff
   };

  if( ( par >= dblLastParCDAS ) && ( par < dblLastBndSlvPar ) )
   return( dflt_dbl_par[ par - dblLastParCDAS ] );

  return( CDASolver::get_dflt_dbl_par( par ) );
  }

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC FIELDS --------------------------------*/
/*--------------------------------------------------------------------------*/
 // the fields are public to allow direct access by training drivers

 Net * nn = nullptr;
 ///< raw pointer to the active network (shared or private), never nullptr
 /**< Raw pointer to the active network; use set_shared_net() /
  * clear_shared_net() to redirect it, do not reassign it directly. */

 VecMem features;     ///< feature vector fed to the network at each iteration

 int size_features;   ///< number of features (the input dimension of Net)

 Vec_VarValue Lambda0;  ///< starting point (initial stability center)

 MemMultiVector pw;     ///< history of previous iterates (currently unused)

 // per-iteration tensors recorded by Heuristic() for the backward pass

 std::vector< torch::Tensor > phi_vecs;    ///< feature tensors (net inputs)
 std::vector< torch::Tensor > w_vecs;      ///< search direction tensors w
 std::vector< torch::Tensor > coeff_vecs;  ///< step-type signs: +1 SS, -1 NS
 std::vector< torch::Tensor > Gs;          ///< subgradient matrices G
 std::vector< torch::Tensor > Gs_aggreg;   ///< aggregated subgradient vectors
 std::vector< torch::Tensor > Qs;          ///< Gram matrices Q = G G^T
 std::vector< torch::Tensor > alphaS;      ///< linearization error vectors
 std::vector< torch::Tensor > thetaS;      ///< currently unused
 std::vector< torch::Tensor > phiS;        ///< currently unused
 std::vector< torch::Tensor > lastIndexS;  ///< currently unused
 std::vector< double > tS;                 ///< predicted step-sizes
 std::vector< double > FiS;                ///< best upper bounds on f*

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE FIELDS -------------------------------*/
/*--------------------------------------------------------------------------*/

 std::shared_ptr< Net > f_owned_net;
 ///< private network allocated at construction time

 std::shared_ptr< Net > f_shared_net;
 ///< shared network injected by set_shared_net(), nullptr if none

 std::unique_ptr< torch::optim::Adam > f_optimizer;
 ///< Adam optimizer, lazily initialized at the first Backward() call;
 ///< kept as a field so that the moment estimates persist across calls

 int f_train_online = 0;     ///< value of intMLTrainOnline [0]

 int f_ML_seed = -1;         ///< value of intMLSeed [-1]

 int f_n_train_rounds = 1;   ///< value of intNTrainRounds [1]

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

 /// re-initializes the privately-owned network
 /** Re-allocates the private network (and drops the stale Adam optimizer,
  * which is bound to the parameters of the old one); used by set_par() when
  * intMLSeed is set, so that the new initial weights are drawn right after
  * the seeding. If a shared network is active it is left untouched, since
  * it is externally managed and possibly shared with other solvers. */

 void reset_net( void ) {
  if( f_shared_net )  // shared net is externally managed: do not touch it
   return;
  f_owned_net = std::make_shared< Net >();
  nn = f_owned_net.get();
  f_optimizer.reset();
  }

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class( BundleSolverML ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* BundleSolverML.h included */

/*--------------------------------------------------------------------------*/
/*--------------------- End File BundleSolverML.h --------------------------*/
/*--------------------------------------------------------------------------*/
