/*--------------------------------------------------------------------------*/
/*--------------------------- File LSTMBlock.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * A thin wrapper that makes a torch::nn::LSTM usable inside a
 * torch::nn::Sequential.
 *
 * THE PROBLEM
 *   torch::nn::Sequential can only chain modules whose forward() takes a
 *   single Tensor and returns a single Tensor. torch::nn::LSTM returns
 *   std::tuple< Tensor , std::tuple< Tensor , Tensor > > (output plus hidden
 *   and cell state), so it cannot be pushed into a Sequential directly.
 *
 * THE SOLUTION
 *   Wrap the LSTM in a small nn::Module whose forward() is Tensor -> Tensor,
 *   keeping ( h , c ) as internal members. From the outside the wrapper looks
 *   like any ordinary one-input / one-output layer, so Sequential accepts it.
 *
 *   The recurrence is preserved: successive forward() calls share the same
 *   ( h , c ), so the network genuinely accumulates memory across bundle
 *   iterations. reset_state() clears it when moving to a new instance.
 *
 * RELATION TO BundleSolverMLNet.h
 *   BundleSolverMLNet takes the other available route: it holds the recurrent
 *   core as a separate member with a hand-written forward(), and uses a
 *   Sequential only for the MLP head. Both routes give a run-time selectable
 *   architecture; this file exists so that the Sequential route is available
 *   as well, and so the choice between the two can be made deliberately.
 *
 * NOTE ON GRADIENTS
 *   ( h , c ) are detached on entry, i.e. truncated BPTT with a window of one
 *   step: gradients do not flow across bundle iterations through the
 *   recurrent state. Removing the two detach() calls below gives full BPTT
 *   over the unrolled horizon, at the cost of a graph that grows with the
 *   number of unrolled steps. This is a design choice, not a constraint.
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __LSTMBlock
 #define __LSTMBlock

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <torch/torch.h>

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it {

/*--------------------------------------------------------------------------*/
/*------------------------------ LSTMBlock ---------------------------------*/
/*--------------------------------------------------------------------------*/
/// an LSTM presented as a Tensor -> Tensor module
/** Holds a torch::nn::LSTM together with its hidden and cell state, and
 * exposes a forward() with the signature that torch::nn::Sequential
 * requires. */

class LSTMBlockImpl : public torch::nn::Module {

 public:

/*--------------------------------------------------------------------------*/
 /// construct the block
 /** @param input_size  size of the feature vector fed at each step
  *  @param hidden_size size of the hidden and cell state
  *  @param num_layers  number of stacked LSTM layers */

 LSTMBlockImpl( int64_t input_size ,
                int64_t hidden_size ,
                int64_t num_layers = 1 )
  : lstm( torch::nn::LSTMOptions( input_size , hidden_size )
           .num_layers( num_layers )
           .batch_first( true ) ) ,
    hidden_size_( hidden_size ) ,
    num_layers_( num_layers )
 {
  register_module( "lstm" , lstm );
 }

/*--------------------------------------------------------------------------*/
 /// Tensor in -> Tensor out; hidden and cell state live inside the module

 torch::Tensor forward( torch::Tensor x ) {

  // add batch and time dimensions if the caller passes a plain feature
  // vector, so the same wrapper works from either shape
  if( x.dim() == 1 )
   x = x.unsqueeze( 0 ).unsqueeze( 0 );        // [ F ] -> [ 1 , 1 , F ]
  else
   if( x.dim() == 2 )
    x = x.unsqueeze( 1 );                      // [ B , F ] -> [ B , 1 , F ]

  // lazily create the state on the first call, then keep it across calls
  if( ! state_initialised_ ) {
   auto opts = torch::TensorOptions().device( x.device() ).dtype( x.dtype() );
   h_ = torch::zeros( { num_layers_ , x.size( 0 ) , hidden_size_ } , opts );
   c_ = torch::zeros( { num_layers_ , x.size( 0 ) , hidden_size_ } , opts );
   state_initialised_ = true;
   }

  // detach the state on entry: truncated BPTT with a window of one step,
  // see the note at the top of this file
  auto h_in = h_.detach();
  auto c_in = c_.detach();

  auto out_tuple = lstm->forward( x , std::make_tuple( h_in , c_in ) );
  auto y  = std::get< 0 >( out_tuple );        // [ B , T , H ]
  auto hc = std::get< 1 >( out_tuple );
  h_ = std::get< 0 >( hc );
  c_ = std::get< 1 >( hc );

  // take the last time step and drop the time dimension -> [ B , H ]
  return( y.select( 1 , y.size( 1 ) - 1 ) );
  }

/*--------------------------------------------------------------------------*/
 /// clear the hidden and cell state
 /** To be called at the start of a new problem instance, so that memory
  * accumulated on the previous one does not leak into it. */

 void reset_state( void ) {
  state_initialised_ = false;
  h_ = torch::Tensor();
  c_ = torch::Tensor();
  }

/*--------------------------------------------------------------------------*/

 private:

 torch::nn::LSTM lstm{ nullptr };
 torch::Tensor   h_ , c_;
 int64_t         hidden_size_;
 int64_t         num_layers_;
 bool            state_initialised_ = false;

 };  // end( class LSTMBlockImpl )

/*--------------------------------------------------------------------------*/

TORCH_MODULE( LSTMBlock );   // usable as: LSTMBlock block( 20 , 32 );

/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*-------------------------------- USAGE -----------------------------------*/
/*--------------------------------------------------------------------------*/
/*
 * Building a Sequential from run-time parameters:
 *
 *   torch::nn::Sequential net;
 *   net->push_back( LSTMBlock( input_size , hidden_size , num_layers ) );
 *   net->push_back( torch::nn::Linear( hidden_size , head_hidden ) );
 *   net->push_back( torch::nn::ReLU() );
 *   net->push_back( torch::nn::Linear( head_hidden , 1 ) );
 *
 * Unrolling over bundle iterations:
 *
 *   for( int step = 0 ; step < unrolled_steps ; ++step ) {
 *    torch::Tensor prediction = net->forward( features_at_step[ step ] );
 *    ...
 *    }
 *
 * Moving to a new instance:
 *
 *   dynamic_cast< LSTMBlockImpl * >( net[ 0 ].get() )->reset_state();
 *
 * The dynamic_cast above is the one rough edge of this approach: the caller
 * has to know which element of the Sequential is stateful. A small
 * StatefulModule interface with a virtual reset_state(), which Sequential
 * could walk over, would remove it.
 */

#endif  /* LSTMBlock.h included */

/*--------------------------------------------------------------------------*/
/*--------------------------- End File LSTMBlock.h -------------------------*/
/*--------------------------------------------------------------------------*/
