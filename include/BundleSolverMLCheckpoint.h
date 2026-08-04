/*--------------------------------------------------------------------------*/
/*-------------------- File BundleSolverMLCheckpoint.h ---------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Saving and restoring the network of BundleSolverML across *different*
 * Blocks.
 *
 * WHY THIS IS NOT A Solver::State
 *   A State is meant to be saved and restored while solving one given Block,
 *   or an identical one; it is not meant to be carried onto a different
 *   Block. A trained network is the opposite case: it is trained on many
 *   instances and then applied to instances it has never seen. The two
 *   concepts serve different purposes, so they are kept separate.
 *
 * WHAT GOES INTO A CHECKPOINT
 *   (1) the architecture descriptor, i.e. NetOptions;
 *   (2) the learned parameters of every module;
 *   (3) the recurrent state ( h , c ), when the model has one.
 *
 *   (1) is there because the network is parametric at run time: torch::save
 *   serialises the parameters, not the shape of the module that held them,
 *   so LSTM weights cannot be loaded into a net that was built as an MLP.
 *   Storing the descriptor next to the weights means a checkpoint is
 *   self-describing -- load_checkpoint() rebuilds the right architecture and
 *   only then fills it in, with no need to remember which configuration
 *   produced which file.
 *
 *   (3) is there because ( h , c ) are ordinary members rather than
 *   registered buffers, so Module::save() does not see them. Whether they
 *   *should* be restored depends on what is being resumed: continuing on the
 *   same instance needs them, starting a new instance does not and should
 *   call reset_state() instead. Both are available -- see the flag on
 *   load_checkpoint() -- but the information is always written, since a
 *   checkpoint that has dropped it can never be made whole again.
 *
 * WHAT DOES *NOT* GO IN
 *   The optimizer. Its state is only meaningful while training continues,
 *   and it is roughly the size of the model again; inference does not need
 *   it. torch::save( optimizer , ... ) alongside, as the training loop
 *   already does, keeps the common case small.
 *
 * USAGE
 *   save_checkpoint( net , "experiments/run_001/best.ckpt" );
 *
 *   auto net = load_checkpoint( "experiments/run_001/best.ckpt" );   // fresh
 *   auto net = load_checkpoint( path , true );   // resume: keep ( h , c )
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __BundleSolverMLCheckpoint
 #define __BundleSolverMLCheckpoint

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "BundleSolverMLNet.h"

#include <torch/torch.h>

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

/*--------------------------------------------------------------------------*/
/*------------------------------ NAMESPACE ---------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it {

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONSTANTS -------------------------------*/
/*--------------------------------------------------------------------------*/

/// bumped whenever the layout below changes in a way older files cannot
/// satisfy; load_checkpoint() refuses anything it does not recognise rather
/// than reading garbage into a net

static constexpr int64_t CheckpointVersion = 1;

/*--------------------------------------------------------------------------*/
/*--------------------------- HELPER FUNCTIONS -----------------------------*/
/*--------------------------------------------------------------------------*/

namespace {   // implementation details

/*--------------------------------------------------------------------------*/
 /// write one scalar as a 0-dim tensor, the archive's only currency

 inline void put_scalar( torch::serialize::OutputArchive & ar ,
                         const std::string & key , int64_t value )
 {
  ar.write( key , torch::tensor( value , torch::kInt64 ) );
  }

/*--------------------------------------------------------------------------*/

 inline void put_scalar( torch::serialize::OutputArchive & ar ,
                         const std::string & key , double value )
 {
  ar.write( key , torch::tensor( value , torch::kFloat64 ) );
  }

/*--------------------------------------------------------------------------*/

 inline int64_t get_int( torch::serialize::InputArchive & ar ,
                         const std::string & key )
 {
  torch::Tensor t;
  ar.read( key , t );
  return( t.item< int64_t >() );
  }

/*--------------------------------------------------------------------------*/

 inline double get_double( torch::serialize::InputArchive & ar ,
                           const std::string & key )
 {
  torch::Tensor t;
  ar.read( key , t );
  return( t.item< double >() );
  }

/*--------------------------------------------------------------------------*/
 /// the architecture descriptor, field by field
 /** Written out one field at a time rather than as an opaque blob, so that a
  * checkpoint stays readable by anything that knows the key names, and so
  * that adding a field later does not invalidate the files already on
  * disk. */

 inline void write_options( torch::serialize::OutputArchive & ar ,
                            const NetOptions & o )
 {
  put_scalar( ar , "opt/model_type"  , int64_t( o.model_type  ) );
  put_scalar( ar , "opt/input_size"  , int64_t( o.input_size  ) );
  put_scalar( ar , "opt/hidden_size" , int64_t( o.hidden_size ) );
  put_scalar( ar , "opt/num_layers"  , int64_t( o.num_layers  ) );
  put_scalar( ar , "opt/activation"  , int64_t( o.activation  ) );
  put_scalar( ar , "opt/stochastic"  , int64_t( o.stochastic ? 1 : 0 ) );
  put_scalar( ar , "opt/t_min"       , o.t_min );
  put_scalar( ar , "opt/t_max"       , o.t_max );

  // the head is a variable-length list, so it goes in as one tensor plus
  // its length; an empty head is a legitimate configuration (a single
  // Linear), which is why the count is written even when it is zero
  const int64_t n = int64_t( o.head_sizes.size() );
  put_scalar( ar , "opt/head_count" , n );

  if( n > 0 ) {
   std::vector< int64_t > widths( o.head_sizes.begin() , o.head_sizes.end() );
   ar.write( "opt/head_sizes" ,
             torch::tensor( widths , torch::kInt64 ).clone() );
   }
  }

/*--------------------------------------------------------------------------*/

 inline NetOptions read_options( torch::serialize::InputArchive & ar )
 {
  NetOptions o;

  o.model_type  = int( get_int( ar , "opt/model_type"  ) );
  o.input_size  = int( get_int( ar , "opt/input_size"  ) );
  o.hidden_size = int( get_int( ar , "opt/hidden_size" ) );
  o.num_layers  = int( get_int( ar , "opt/num_layers"  ) );
  o.activation  = int( get_int( ar , "opt/activation"  ) );
  o.stochastic  = ( get_int( ar , "opt/stochastic" ) != 0 );
  o.t_min       = get_double( ar , "opt/t_min" );
  o.t_max       = get_double( ar , "opt/t_max" );

  const int64_t n = get_int( ar , "opt/head_count" );

  o.head_sizes.clear();
  if( n > 0 ) {
   torch::Tensor widths;
   ar.read( "opt/head_sizes" , widths );
   auto flat = widths.contiguous().view( { -1 } );
   for( int64_t i = 0 ; i < n ; ++i )
    o.head_sizes.push_back( int( flat[ i ].item< int64_t >() ) );
   }

  return( o );
  }

/*--------------------------------------------------------------------------*/
 /// one recurrent tensor, guarded by a presence flag
 /** h is undefined for an MLP and before the first reset_state(), so the
  * flag distinguishes "not written" from "written as zeros" -- which are
  * different situations on the way back in. */

 inline void write_state_tensor( torch::serialize::OutputArchive & ar ,
                                 const std::string & key ,
                                 const torch::Tensor & t )
 {
  const bool present = t.defined() && t.numel() > 0;
  put_scalar( ar , key + "_present" , int64_t( present ? 1 : 0 ) );

  if( present )
   ar.write( key , t.detach().to( torch::kCPU ).clone() );
  }

/*--------------------------------------------------------------------------*/

 inline void read_state_tensor( torch::serialize::InputArchive & ar ,
                                const std::string & key ,
                                torch::Tensor & t )
 {
  if( get_int( ar , key + "_present" ) == 0 ) {
   t = torch::Tensor();
   return;
   }

  torch::Tensor tmp;
  ar.read( key , tmp );
  t = tmp;
  }

/*--------------------------------------------------------------------------*/

 }  // end( anonymous namespace )

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
 /// write net to path: architecture, weights and recurrent state
 /** The three parts go into one file so that a checkpoint cannot be
  * separated from the description of what it is. */

inline void save_checkpoint( const std::shared_ptr< Net > & net ,
                             const std::string & path )
{
 if( ! net )
  throw( std::invalid_argument( "save_checkpoint: null net" ) );

 torch::serialize::OutputArchive ar;

 put_scalar( ar , "version" , CheckpointVersion );

 // (1) architecture
 write_options( ar , net->opt );

 // (2) learned parameters, in a nested archive of their own so that the
 //     module's own naming cannot collide with the keys used here
 torch::serialize::OutputArchive weights;
 net->save( weights );
 ar.write( "weights" , weights );

 // (3) recurrent state
 write_state_tensor( ar , "state/h" , net->h );
 write_state_tensor( ar , "state/c" , net->c );

 ar.save_to( path );
 }

/*--------------------------------------------------------------------------*/
 /// rebuild a net from path
 /** @param path            the file written by save_checkpoint()
  *  @param restore_state   if true, ( h , c ) are restored as they were;
  *                         if false, the net starts from reset_state().
  *
  * The default is false: a checkpoint is usually being loaded to run on a
  * new instance, and carrying over the memory of whichever instance
  * happened to be in progress when it was written would be wrong. Pass true
  * to continue on the same instance. */

inline std::shared_ptr< Net > load_checkpoint( const std::string & path ,
                                               bool restore_state = false )
{
 torch::serialize::InputArchive ar;
 ar.load_from( path );

 const int64_t version = get_int( ar , "version" );
 if( version != CheckpointVersion )
  throw( std::runtime_error(
   "load_checkpoint: file version " + std::to_string( version ) +
   ", expected " + std::to_string( CheckpointVersion ) ) );

 // (1) architecture first: the net has to have the right shape before
 //     anything can be loaded into it
 auto net = std::make_shared< Net >( read_options( ar ) );

 // (2) learned parameters
 torch::serialize::InputArchive weights;
 ar.read( "weights" , weights );
 net->load( weights );

 // (3) recurrent state
 if( restore_state ) {
  read_state_tensor( ar , "state/h" , net->h );
  read_state_tensor( ar , "state/c" , net->c );
  }
 else
  net->reset_state();

 return( net );
 }

/*--------------------------------------------------------------------------*/
 /// the NetOptions stored in a checkpoint, without building the net
 /** Useful for reporting on a directory of experiments, or for checking that
  * a checkpoint matches the instance family it is about to be run on. */

inline NetOptions peek_checkpoint_options( const std::string & path )
{
 torch::serialize::InputArchive ar;
 ar.load_from( path );

 const int64_t version = get_int( ar , "version" );
 if( version != CheckpointVersion )
  throw( std::runtime_error( "peek_checkpoint_options: unknown version " +
                             std::to_string( version ) ) );

 return( read_options( ar ) );
 }

/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* BundleSolverMLCheckpoint.h included */

/*--------------------------------------------------------------------------*/
/*------------- End File BundleSolverMLCheckpoint.h ------------------------*/
/*--------------------------------------------------------------------------*/
