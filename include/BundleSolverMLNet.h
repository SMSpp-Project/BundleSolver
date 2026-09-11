/*--------------------------------------------------------------------------*/
/*----------------------- File BundleSolverMLNet.h -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * PROPOSAL / PROTOTYPE (not yet compiled against the real repo).
 *
 * A *runtime-parametric* replacement for the hard-coded `Net` struct
 * currently found in BundleSolverML.h.
 *
 * NOTE: the struct below is deliberately named `Net`, in the same namespace,
 * because it is meant to take the place of the existing one. Nothing includes
 * this file yet, so there is no clash today; but this header and
 * BundleSolverML.h cannot both be included in the same translation unit until
 * the substitution is actually made. Wiring it in is the next step, and is
 * left out of this commit on purpose.
 *
 * The aim is to be able to change the underlying ML model *without* editing
 * the code and recompiling. Ideally everything is driven by algorithmic
 * parameters (i.e. the BSPar-ML.txt file); at worst by a few macros.
 *
 * Findings on how flexible Torch actually is, summarised:
 *
 *  (1) FULLY RUNTIME-CONFIGURABLE
 *      - layer sizes, number of layers, hidden size:
 *        every module takes an *Options* object built at run time,
 *        e.g. torch::nn::LinearOptions( in , out ),
 *             torch::nn::LSTMOptions( in , hid ).num_layers( L ).
 *      - depth of a feed-forward stack: torch::nn::Sequential is a *container*
 *        that can be push_back()-ed in a loop, so an arbitrary number of
 *        layers can be assembled from an integer parameter.
 *      - choice of activation: torch::nn::AnyModule type-erases the module,
 *        so ReLU / Softplus / Tanh can be selected by an int at run time
 *        and pushed into the same Sequential.
 *      - choice of recurrent cell (RNN / GRU / LSTM): all three exist as
 *        separate module types; we hold all three as (null) members and
 *        instantiate exactly one.
 *
 *  (2) THE ONE AWKWARD POINT
 *      torch::nn::Sequential can only chain modules whose forward() takes a
 *      single Tensor and returns a single Tensor. The recurrent modules
 *      return std::tuple< Tensor , std::tuple< Tensor , Tensor > > (output +
 *      hidden state), so an LSTM cannot be pushed into a Sequential *as is*.
 *      The route taken below is therefore: keep the recurrent core as a
 *      separate member with a hand-written forward(), and use a Sequential
 *      only for the MLP *head*.
 *
 *      UPDATE: this is not the only route. A thin wrapper module whose
 *      forward() is Tensor -> Tensor, holding ( h , c ) internally, does let
 *      an LSTM sit inside a Sequential; see LSTMBlock.h, where it is
 *      implemented and tested. Both routes give the same run-time
 *      configurability, and neither is forced on us -- the wrapper keeps the
 *      whole network as one uniform container, the explicit member below
 *      keeps the state handling visible at the call site. The choice can be
 *      made deliberately rather than by default.
 *
 *  (3) CONSEQUENCE FOR State / netCDF (relevant to the get_State/set_State
 *      idea): if the architecture is a run-time parameter, then the *shape*
 *      of the parameter vector is no longer known at compile time. Hence a
 *      BundleSolverMLState must serialise the architecture descriptor
 *      (model type, n_layers, hidden size, head sizes) *alongside* the flat
 *      tensor data, otherwise set_State() cannot rebuild the modules before
 *      loading the weights. named_parameters() gives a stable, ordered
 *      (name -> Tensor) map that flattens cleanly into netCDF doubles.
 *
 *  (4) NOTE: the optimiser (Adam) must be constructed *after* the modules,
 *      since it captures parameters(). Rebuilding the net invalidates it.
 */

#ifndef __BundleSolverMLNet
#define __BundleSolverMLNet

#include <torch/torch.h>
#include <vector>
#include <stdexcept>

namespace SMSpp_di_unipi_it {

/*--------------------------------------------------------------------------*/
/*------------------------------ NetOptions --------------------------------*/
/*--------------------------------------------------------------------------*/
/// everything that used to be hard-coded, now a plain data object
/** Filled from the algorithmic parameters of BundleSolverML (see the
 * set_par() sketch at the bottom of this file), so that the whole
 * architecture is decided at run time. */

struct NetOptions {

 /// which core the net uses
 enum ModelType {
  eMLP  = 0 ,   ///< no recurrence (the current behaviour)
  eRNN  = 1 ,   ///< vanilla Elman RNN
  eGRU  = 2 ,   ///< gated recurrent unit
  eLSTM = 3     ///< long short-term memory
  };

 /// which nonlinearity the head uses
 enum ActType {
  eReLU     = 0 ,
  eSoftplus = 1 ,
  eTanh     = 2
  };

 int  model_type  = eMLP;   ///< ModelType
 int  input_size  = 20;     ///< size of the feature vector phi_t
 int  hidden_size = 16;     ///< hidden size of the recurrent core
 int  num_layers  = 1;      ///< number of stacked recurrent layers
 int  activation  = eSoftplus;  ///< ActType used in the head

 /// sizes of the hidden layers of the MLP head; empty == single Linear
 std::vector< int > head_sizes = { 16 };

 /// if true, the core output is turned into a Gaussian (mu,sigma) and
 /// sampled with the reparametrization trick (cf. thesis, Fig. 5.1)
 bool stochastic = false;

 /// the predicted t is clamped into [ t_min , t_max ] (thesis: 1e-5..1e4)
 double t_min = 1e-5;
 double t_max = 1e+4;

 };  // end( struct NetOptions )

/*--------------------------------------------------------------------------*/
/*--------------------------------- Net ------------------------------------*/
/*--------------------------------------------------------------------------*/
/// a Net whose architecture is decided at construction time from NetOptions
/** Replaces the fixed { Linear(20,16) -> Softplus -> Linear(16,1) } network.
 *
 * Structure (mirrors Figure 5.1 of Demelas' thesis):
 *
 *     phi_t --> [ recurrent core ] --> [ (mu,sigma) sampler ] --> [ MLP head ]
 *                     |                     (optional)                |
 *                 h_t , c_t                                           v
 *              (carried across                                  softplus + clamp
 *               bundle iterations)                                    |
 *                                                                     v
 *                                                                     t > 0
 *
 * When model_type == eMLP the recurrent core degenerates into the identity
 * and the class behaves exactly like the current implementation, so the old
 * results remain reproducible with a single parameter value. */

struct Net : torch::nn::Module {

 NetOptions opt;   ///< the architecture descriptor (needed by set_State())

 // --- the recurrent core: exactly one of these is instantiated -----------
 torch::nn::RNN  rnn { nullptr };
 torch::nn::GRU  gru { nullptr };
 torch::nn::LSTM lstm{ nullptr };

 // --- the (optional) stochastic bottleneck -------------------------------
 torch::nn::Linear fc_mu    { nullptr };
 torch::nn::Linear fc_logvar{ nullptr };

 // --- the head: a dynamically assembled stack of Linear + activation -----
 torch::nn::Sequential head{ nullptr };

 // --- the hidden state, carried across bundle iterations -----------------
 torch::Tensor h , c;

/*--------------------------------------------------------------------------*/

 explicit Net( const NetOptions & o ) : opt( o ) { build(); }

/*--------------------------------------------------------------------------*/
 /// assemble every module from opt; no size is hard-coded

 void build( void ) {

  const int in = opt.input_size;
  const int hs = opt.hidden_size;

  // ---- 1) recurrent core ------------------------------------------------
  // note: an Options object is an ordinary C++ value, so all of this is
  // decided at run time. Nothing here needs a macro.

  switch( opt.model_type ) {

   case( NetOptions::eRNN ):
    rnn = register_module( "rnn" , torch::nn::RNN(
     torch::nn::RNNOptions( in , hs ).num_layers( opt.num_layers )
                                     .batch_first( true ) ) );
    break;

   case( NetOptions::eGRU ):
    gru = register_module( "gru" , torch::nn::GRU(
     torch::nn::GRUOptions( in , hs ).num_layers( opt.num_layers )
                                     .batch_first( true ) ) );
    break;

   case( NetOptions::eLSTM ):
    lstm = register_module( "lstm" , torch::nn::LSTM(
     torch::nn::LSTMOptions( in , hs ).num_layers( opt.num_layers )
                                      .batch_first( true ) ) );
    break;

   case( NetOptions::eMLP ):
    break;   // no core: the head reads phi_t directly

   default:
    throw( std::invalid_argument( "Net: unknown model_type" ) );
   }

  // width of whatever reaches the head
  const int core_out = ( opt.model_type == NetOptions::eMLP ) ? in : hs;

  // ---- 2) optional Gaussian bottleneck ----------------------------------

  int head_in = core_out;
  if( opt.stochastic ) {
   fc_mu     = register_module( "fc_mu" ,
                                torch::nn::Linear( core_out , hs ) );
   fc_logvar = register_module( "fc_logvar" ,
                                torch::nn::Linear( core_out , hs ) );
   head_in = hs;
   }

  // ---- 3) the head: an arbitrary-depth stack, built in a loop -----------
  // this is where Sequential shines: push_back() in a loop over a vector
  // read from the parameter file. AnyModule type-erases the activation so
  // that its *type* can also be a run-time choice.

  head = torch::nn::Sequential();

  int prev = head_in;
  for( int width : opt.head_sizes ) {
   head->push_back( torch::nn::Linear( prev , width ) );
   head->push_back( make_activation() );
   prev = width;
   }

  head->push_back( torch::nn::Linear( prev , 1 ) );   // scalar output: t
  register_module( "head" , head );

  reset_state();
  }

/*--------------------------------------------------------------------------*/
 /// type-erased activation, chosen by an integer parameter

 torch::nn::AnyModule make_activation( void ) const {
  switch( opt.activation ) {
   case( NetOptions::eReLU ):
    return( torch::nn::AnyModule( torch::nn::ReLU() ) );
   case( NetOptions::eSoftplus ):
    return( torch::nn::AnyModule( torch::nn::Softplus() ) );
   case( NetOptions::eTanh ):
    return( torch::nn::AnyModule( torch::nn::Tanh() ) );
   default:
    throw( std::invalid_argument( "Net: unknown activation" ) );
   }
  }

/*--------------------------------------------------------------------------*/
 /// zero the hidden state; MUST be called at the start of each new instance
 /** The RNN is applied *dynamically to the bundle iterations*: one call of
  * forward() == one time step. Hence the state has to survive between calls
  * but be cleared between different problems. */

 void reset_state( void ) {
  if( opt.model_type == NetOptions::eMLP )
   return;
  h = torch::zeros( { opt.num_layers , 1 , opt.hidden_size } );
  if( opt.model_type == NetOptions::eLSTM )
   c = torch::zeros( { opt.num_layers , 1 , opt.hidden_size } );
  }

/*--------------------------------------------------------------------------*/
 /// one bundle iteration == one time step; returns the predicted t > 0
 /** x is the feature tensor phi_t of shape { input_size }. The hidden state
  * is updated in place, so the time dependency is handled implicitly and no
  * external sequence buffer is needed. */

 torch::Tensor forward( torch::Tensor x ) {

  torch::Tensor y;

  // ---- 1) recurrent core -----------------------------------------------
  // NOTE: this is the part that cannot live inside a Sequential, because
  // forward() here returns a tuple rather than a single Tensor.

  switch( opt.model_type ) {

   case( NetOptions::eMLP ):
    y = x;                                   // identity core
    break;

   case( NetOptions::eRNN ): {
    auto inp = x.view( { 1 , 1 , -1 } );     // { batch , time , feature }
    auto out = rnn->forward( inp , h );
    y = std::get< 0 >( out ).view( { -1 } ); // last (only) time step
    h = std::get< 1 >( out );                // carry the state forward
    break;
    }

   case( NetOptions::eGRU ): {
    auto inp = x.view( { 1 , 1 , -1 } );
    auto out = gru->forward( inp , h );
    y = std::get< 0 >( out ).view( { -1 } );
    h = std::get< 1 >( out );
    break;
    }

   case( NetOptions::eLSTM ): {
    auto inp = x.view( { 1 , 1 , -1 } );
    auto out = lstm->forward( inp , std::make_tuple( h , c ) );
    y = std::get< 0 >( out ).view( { -1 } );
    auto st = std::get< 1 >( out );
    h = std::get< 0 >( st );
    c = std::get< 1 >( st );
    break;
    }
   }

  // ---- 2) optional reparametrization trick ------------------------------

  if( opt.stochastic ) {
   auto mu     = fc_mu->forward( y );
   auto logvar = fc_logvar->forward( y );
   auto sigma  = torch::exp( 0.5 * logvar );
   // sample = mu + sigma * eps, eps ~ N(0,I): differentiable wrt mu,sigma
   y = is_training() ? mu + sigma * torch::randn_like( sigma ) : mu;
   }

  // ---- 3) head + positivity + clamping ----------------------------------

  auto t = head->forward( y );
  t = torch::softplus( t ) + 1.0e-8;         // strictly positive
  t = torch::clamp( t , opt.t_min , opt.t_max );

  return( t );
  }

 };  // end( struct Net )

/*--------------------------------------------------------------------------*/
/*------------------- SKETCH: wiring into the parameters --------------------*/
/*--------------------------------------------------------------------------*/
/* Inside BundleSolverML, following the usual SMS++ convention, we would add
 * the new algorithmic parameters right after the last BundleSolver one (the
 * exact names of the "last" enumerators must be taken from BundleSolver.h):
 *
 *   enum int_par_type_BSML {
 *    intMLModelType = intLastBSlvPar ,   ///< NetOptions::ModelType
 *    intMLHiddenSize ,                   ///< hidden size of the core
 *    intMLNumLayers ,                    ///< stacked recurrent layers
 *    intMLActivation ,                   ///< NetOptions::ActType
 *    intMLHeadDepth ,                    ///< how many hidden layers in head
 *    intMLHeadWidth ,                    ///< width of each of them
 *    intMLStochastic ,                   ///< 0/1: reparametrization on/off
 *    intLastBSMLPar
 *    };
 *
 *   enum dbl_par_type_BSML {
 *    dblMLtMin = dblLastBSlvPar ,
 *    dblMLtMax ,
 *    dblMLLearnRate ,
 *    dblLastBSMLPar
 *    };
 *
 * and then, in set_par(), simply record the value and mark the net as dirty;
 * the Net is (re)built lazily on the first compute() after any change:
 *
 *   void BundleSolverML::set_par( idx_type par , int value ) {
 *    switch( par ) {
 *     case( intMLModelType ):  f_opt.model_type  = value; f_dirty = true; break;
 *     case( intMLHiddenSize ): f_opt.hidden_size = value; f_dirty = true; break;
 *     ...
 *     default: BundleSolver::set_par( par , value );
 *     }
 *    }
 *
 *   void BundleSolverML::rebuild_net_if_needed( void ) {
 *    if( ! f_dirty ) return;
 *    f_opt.head_sizes.assign( f_head_depth , f_head_width );
 *    f_owned_net = std::make_shared< Net >( f_opt );
 *    nn = f_owned_net.get();
 *    // the optimiser captures parameters(), so it must be rebuilt too
 *    f_optimizer = std::make_unique< torch::optim::Adam >(
 *                   nn->parameters() , torch::optim::AdamOptions( f_lr ) );
 *    f_dirty = false;
 *    }
 *
 * With this, a BSPar-ML.txt entry such as
 *
 *   intMLModelType   3     # 0=MLP 1=RNN 2=GRU 3=LSTM
 *   intMLHiddenSize  64
 *   intMLNumLayers   2
 *   intMLHeadDepth   1
 *   intMLHeadWidth   32
 *   intMLStochastic  1
 *
 * switches the whole architecture with no recompilation. Only two things
 * would still need compile-time work: adding a genuinely *new* module type
 * (a new case in build()), and any change to the shape of the feature
 * vector phi_t itself.
 */

}  // end( namespace SMSpp_di_unipi_it )

#endif  /* BundleSolverMLNet.h included */

/*--------------------------------------------------------------------------*/
/*--------------------- End File BundleSolverMLNet.h -----------------------*/
/*--------------------------------------------------------------------------*/
