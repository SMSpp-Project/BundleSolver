/*--------------------------------------------------------------------------*/
/*------------------------ File LSTMBlock_test.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Checks that LSTMBlock behaves the way BundleSolverML needs it to:
 *
 *  (1) torch::nn::Sequential accepts it as an ordinary layer, which a bare
 *      torch::nn::LSTM cannot be, since its forward() returns a tuple;
 *
 *  (2) the recurrent state survives between calls, so that feeding the same
 *      input twice does *not* give the same answer -- this is what makes the
 *      network able to accumulate information over bundle iterations;
 *
 *  (3) reset_state() really clears that state, so a new instance starts from
 *      a clean slate rather than inheriting the previous one's memory.
 *
 * Returns 0 if all three hold, 1 otherwise.
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "LSTMBlock.h"

#include <torch/torch.h>

#include <cmath>
#include <iostream>
#include <vector>

/*--------------------------------------------------------------------------*/
/*------------------------------ NAMESPACE ---------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------------- MAIN ------------------------------------*/
/*--------------------------------------------------------------------------*/

int main( void )
{
 torch::manual_seed( 0 );

 constexpr int64_t input_size  = 20;
 constexpr int64_t hidden_size = 32;
 constexpr int64_t head_hidden = 16;
 constexpr double  tolerance   = 1e-9;

 bool ok = true;

 /*------------------------- (1) build the network -------------------------*/
 // exactly the shape a bare LSTM could not take, because Sequential requires
 // Tensor -> Tensor at every stage

 torch::nn::Sequential net;
 net->push_back( LSTMBlock( input_size , hidden_size ) );
 net->push_back( torch::nn::Linear( hidden_size , head_hidden ) );
 net->push_back( torch::nn::ReLU() );
 net->push_back( torch::nn::Linear( head_hidden , 1 ) );

 std::cout << "Sequential built with " << net->size()
           << " layers ( LSTMBlock + Linear + ReLU + Linear )" << std::endl
           << std::endl;

 if( net->size() != 4 ) {
  std::cout << "FAILED: the Sequential did not accept all four layers"
            << std::endl;
  return( 1 );
  }

 /*--------------------- (2) state carried across calls --------------------*/
 // same input every time; if the state were not being kept, every output
 // would be identical

 std::cout << "--- successive forward calls on the same input ---"
           << std::endl;

 torch::Tensor x = torch::ones( { input_size } );
 std::vector< double > first_run;

 for( int step = 1 ; step <= 5 ; ++step ) {
  double y = net->forward( x ).item< double >();
  first_run.push_back( y );
  std::cout << " iteration " << step << "   output = " << y << std::endl;
  }

 bool all_equal = true;
 for( std::size_t i = 1 ; i < first_run.size() ; ++i )
  if( std::abs( first_run[ i ] - first_run[ 0 ] ) > tolerance ) {
   all_equal = false;
   break;
   }

 if( all_equal ) {
  std::cout << std::endl
            << "FAILED: every output is the same, so no state is being kept"
            << std::endl;
  ok = false;
  }

 /*------------------------- (3) reset_state() -----------------------------*/
 // after clearing the state the network must reproduce the very first output

 std::cout << std::endl << "--- reset_state() , then run again ---"
           << std::endl;

 auto wrapper =
  std::dynamic_pointer_cast< LSTMBlockImpl >( net->children()[ 0 ] );

 if( ! wrapper ) {
  std::cout << "FAILED: could not recover the LSTMBlock from the Sequential"
            << std::endl;
  return( 1 );
  }

 wrapper->reset_state();

 for( int step = 1 ; step <= 3 ; ++step ) {
  double y = net->forward( x ).item< double >();
  std::cout << " iteration " << step << "   output = " << y << std::endl;

  if( std::abs( y - first_run[ step - 1 ] ) > tolerance ) {
   std::cout << std::endl << "FAILED: output " << step
             << " after reset differs from the same step before reset ( "
             << y << " against " << first_run[ step - 1 ] << " )"
             << std::endl;
   ok = false;
   }
  }

 /*------------------------------ verdict ----------------------------------*/

 std::cout << std::endl;

 if( ok )
  std::cout << "All tests passed!" << std::endl
            << "An LSTM lives inside a Sequential, its state is carried "
               "across calls, and reset_state() clears it." << std::endl;
 else
  std::cout << "Some tests failed." << std::endl;

 return( ok ? 0 : 1 );

 }  // end( main )

/*--------------------------------------------------------------------------*/
/*--------------------- End File LSTMBlock_test.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
