/*--------------------------------------------------------------------------*/
/*----------------------- File BundleSolverMLNet.h -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Definition of Net, the network with which BundleSolverML predicts the
 * step-size t at each iteration, and of NetOptions, the plain data object
 * that describes its architecture.
 *
 * The architecture is decided at run time: a recurrent core (none, i.e. a
 * feed-forward network, or an RNN, a GRU or an LSTM, of given hidden size
 * and number of layers), an optional Gaussian bottleneck, and a head made of
 * a stack of Linear layers with a chosen activation, ending in a single
 * strictly positive output clamped to [ t_min , t_max ]. BundleSolverML
 * builds it from intMLModel and intMLHidden; the remaining fields of
 * NetOptions keep their defaults.
 *
 * The recurrent core is applied to the bundle iterations, one call of
 * forward() being one time step, and hence it carries its hidden state from
 * a call to the next; reset_state() clears it at the start of each instance
 * and detach_state() cuts it off the autograd graph. The core is a separate
 * member with its own forward() rather than an element of a
 * torch::nn::Sequential, since the recurrent modules of Torch return their
 * state alongside the output (LSTMBlock.h gives the alternative of a wrapper
 * that does fit in a Sequential).
 *
 * Since the shape of the parameters depends on NetOptions, whatever saves
 * the weights has to save NetOptions too (see BundleSolverMLCheckpoint.h),
 * and the Adam optimizer, which is bound to parameters(), has to be rebuilt
 * whenever the network is.
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
/// the whole architecture of the network as a plain data object
/** model_type and hidden_size are set from the intMLModel and intMLHidden
 * parameters of BundleSolverML, the other fields keep their defaults. */

struct NetOptions {

 /// which core the net uses
 enum ModelType {
  eMLP  = 0 ,   ///< no recurrence: the head reads phi_t directly
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
/** Structure (mirrors Figure 5.1 of Demelas' thesis):
 *
 *     phi_t --> [ recurrent core ] --> [ (mu,sigma) sampler ] --> [ MLP head ]
 *                     |                     (optional)                |
 *                 h_t , c_t                                           v
 *              (carried across                                  softplus + clamp
 *               bundle iterations)                                    |
 *                                                                     v
 *                                                                     t > 0
 *
 * When model_type == eMLP the recurrent core is the identity, and with the
 * default NetOptions the network is { Linear(20,16) -> Softplus ->
 * Linear(16,1) }. */

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
 /// cut the hidden state off the autograd graph that produced it
 /** The state keeps its value, but the gradient no longer flows back into
  * the time steps that computed it; this is what truncates the
  * back-propagation through time at the start of each training window. */

 void detach_state( void ) {
  if( h.defined() )
   h = h.detach();
  if( c.defined() )
   c = c.detach();
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

}  // end( namespace SMSpp_di_unipi_it )

#endif  /* BundleSolverMLNet.h included */

/*--------------------------------------------------------------------------*/
/*--------------------- End File BundleSolverMLNet.h -----------------------*/
/*--------------------------------------------------------------------------*/
