/*--------------------------------------------------------------------------*/
/*------------- File BundleSolverMLCheckpoint_test.cpp ---------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Checks that a checkpoint written by save_checkpoint() comes back as the
 * same network:
 *
 *  (1) the architecture descriptor survives the round trip, so that a
 *      checkpoint is self-describing and does not rely on the caller
 *      remembering which configuration produced it;
 *
 *  (2) an LSTM net reloaded with restore_state = true takes up where the
 *      original left off, rather than starting the sequence again -- this
 *      is the case that failed before ( h , c ) were written;
 *
 *  (3) the same net reloaded with restore_state = false starts the sequence
 *      from the beginning, as a new instance should;
 *
 *  (4) the MLP path, which has no recurrent state, still round-trips
 *      exactly, so the fix has not disturbed the existing behaviour.
 *
 * Returns 0 if all four hold, 1 otherwise.
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "BundleSolverMLCheckpoint.h"

#include <torch/torch.h>

#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*------------------------------- CONSTANTS --------------------------------*/
/*--------------------------------------------------------------------------*/

// exact equality is the right bar here: the same weights fed the same input
// on the same machine should give bit-identical results, and anything looser
// would have accepted the bug this test exists to catch
static constexpr double Tolerance = 1e-12;

static const std::string CkptPath = "checkpoint_test.ckpt";

/*--------------------------------------------------------------------------*/
/*------------------------------- HELPERS ----------------------------------*/
/*--------------------------------------------------------------------------*/

static NetOptions lstm_options( void )
{
 NetOptions o;
 o.model_type  = NetOptions::eLSTM;
 o.input_size  = 20;
 o.hidden_size = 16;
 o.num_layers  = 1;
 o.activation  = NetOptions::eSoftplus;
 o.head_sizes  = { 16 , 8 };     // two hidden layers: exercises the vector
 return( o );
 }

/*--------------------------------------------------------------------------*/

static bool same_options( const NetOptions & a , const NetOptions & b )
{
 return( a.model_type  == b.model_type  &&
         a.input_size  == b.input_size  &&
         a.hidden_size == b.hidden_size &&
         a.num_layers  == b.num_layers  &&
         a.activation  == b.activation  &&
         a.stochastic  == b.stochastic  &&
         a.head_sizes  == b.head_sizes  &&
         std::abs( a.t_min - b.t_min ) < Tolerance &&
         std::abs( a.t_max - b.t_max ) < Tolerance );
 }

/*--------------------------------------------------------------------------*/
 /// run the net forward n times on the same input, collecting the outputs
 /** Feeding a constant input is deliberate: with a working recurrent state
  * the outputs still differ from one call to the next, so the sequence is a
  * fingerprint of the state as much as of the weights. */

static std::vector< double > run( std::shared_ptr< Net > net ,
                                  const torch::Tensor & x , int n )
{
 std::vector< double > out;
 torch::NoGradGuard no_grad;
 for( int i = 0 ; i < n ; ++i )
  out.push_back( net->forward( x ).item< double >() );
 return( out );
 }

/*--------------------------------------------------------------------------*/

static bool matches( const std::vector< double > & a ,
                     const std::vector< double > & b )
{
 if( a.size() != b.size() )
  return( false );
 for( std::size_t i = 0 ; i < a.size() ; ++i )
  if( std::abs( a[ i ] - b[ i ] ) > Tolerance )
   return( false );
 return( true );
 }

/*--------------------------------------------------------------------------*/

static void report( const std::vector< double > & a ,
                    const std::vector< double > & b )
{
 for( std::size_t i = 0 ; i < a.size() && i < b.size() ; ++i )
  std::cout << "   step " << ( i + 1 ) << " : " << a[ i ]
            << "   against " << b[ i ]
            << "   difference " << std::abs( a[ i ] - b[ i ] ) << std::endl;
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------- MAIN -----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( void )
{
 torch::manual_seed( 0 );

 bool ok = true;
 torch::Tensor x = torch::ones( { 20 } );

 /*---------------------- (1) architecture round trip ---------------------*/

 std::cout << "--- architecture descriptor ---" << std::endl;

 auto original = std::make_shared< Net >( lstm_options() );
 original->reset_state();

 // steps 1..5 from a cleared state: what a *new* instance should see
 auto from_clean = run( original , x , 5 );

 // save mid-instance, with ( h , c ) carrying the memory of those 5 steps
 save_checkpoint( original , CkptPath );

 // steps 6..10 on the original: what *continuing* should look like
 auto continued = run( original , x , 5 );

 const NetOptions recovered = peek_checkpoint_options( CkptPath );

 if( same_options( original->opt , recovered ) )
  std::cout << " descriptor survives the round trip" << std::endl;
 else {
  std::cout << " FAILED: descriptor differs after reload" << std::endl;
  ok = false;
  }

 /*------------------- (2) reload keeping ( h , c ) -----------------------*/
 // the reloaded net must take up at step 6, not start over at step 1;
 // this is the case that failed before ( h , c ) were written

 std::cout << std::endl
           << "--- reload with restore_state = true ---" << std::endl;

 auto resumed     = load_checkpoint( CkptPath , true );
 auto resumed_out = run( resumed , x , 5 );

 if( matches( continued , resumed_out ) )
  std::cout << " predictions continue exactly where they left off"
            << std::endl;
 else {
  std::cout << " FAILED: resumed net does not match the original's"
            << " steps 6..10" << std::endl;
  report( continued , resumed_out );
  ok = false;
  }

 if( matches( from_clean , resumed_out ) ) {
  std::cout << " FAILED: resumed net reproduced steps 1..5, so ( h , c )"
            << " were not restored" << std::endl;
  ok = false;
  }

 /*------------------- (3) reload discarding ( h , c ) --------------------*/
 // the same weights from a cleared state: what a new instance should see,
 // i.e. steps 1..5 again

 std::cout << std::endl
           << "--- reload with restore_state = false ---" << std::endl;

 auto fresh     = load_checkpoint( CkptPath , false );
 auto fresh_out = run( fresh , x , 5 );

 if( matches( from_clean , fresh_out ) )
  std::cout << " a fresh load starts from a cleared state" << std::endl;
 else {
  std::cout << " FAILED: fresh load does not reproduce steps 1..5"
            << std::endl;
  report( from_clean , fresh_out );
  ok = false;
  }

 /*------------------------- (4) the MLP path -----------------------------*/

 std::cout << std::endl << "--- MLP, which has no recurrent state ---"
           << std::endl;

 NetOptions mlp;                       // defaults: eMLP, 20 -> 16 -> 1
 auto mlp_net = std::make_shared< Net >( mlp );
 auto mlp_before = run( mlp_net , x , 3 );

 save_checkpoint( mlp_net , CkptPath );
 auto mlp_reloaded = load_checkpoint( CkptPath );
 auto mlp_after    = run( mlp_reloaded , x , 3 );

 if( matches( mlp_before , mlp_after ) )
  std::cout << " round-trips exactly" << std::endl;
 else {
  std::cout << " FAILED: MLP predictions differ after reload" << std::endl;
  report( mlp_before , mlp_after );
  ok = false;
  }

 /*------------------------------ verdict ---------------------------------*/

 std::remove( CkptPath.c_str() );

 std::cout << std::endl;

 if( ok )
  std::cout << "All tests passed!" << std::endl
            << "Architecture, weights and recurrent state all survive a "
               "checkpoint, and a fresh load starts clean." << std::endl;
 else
  std::cout << "Some tests failed." << std::endl;

 return( ok ? 0 : 1 );

 }  // end( main )

/*--------------------------------------------------------------------------*/
/*----------- End File BundleSolverMLCheckpoint_test.cpp -------------------*/
/*--------------------------------------------------------------------------*/
