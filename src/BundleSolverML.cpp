/*--------------------------------------------------------------------------*/
/*--------------------- File BundleSolverML.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the BundleSolverML class.
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
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include "BundleSolverML.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

using Index = Function::Index;

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// maximum L2 norm for gradient clipping applied after loss.backward();
/// if the norm of the gradients exceeds this value, all the gradients are
/// rescaled by GRAD_CLIP_NORM / norm (typical values are 1 ... 10)
static constexpr float GRAD_CLIP_NORM = 1.0f;

/*--------------------------------------------------------------------------*/
/*------------------------- FACTORY MANAGEMENT -----------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it
{
 SMSpp_insert_in_factory_cpp_0( BundleSolverML );
 }

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE Block -----------------------*/
/*--------------------------------------------------------------------------*/

int BundleSolverML::compute( bool changedvars )
{
 // the threads of Torch, for the predictions and for the training that
 // follows [see intMaxThread]
 torch::set_num_threads( std::max( 1 , f_max_thread ) );

 // a new instance: the recurrent core of the network starts from a null
 // state (nothing happens with a feed-forward network)
 nn->reset_state();

 const int ret = BundleSolver::compute( changedvars );

 if( f_train_online ) {  // online training: learn from the just-ended solve
  Backward();
  ClearBuffers();
  }

 return( ret );
 }

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR HANDLING THE MODEL -----------------------*/
/*--------------------------------------------------------------------------*/

void BundleSolverML::SaveModel( const std::string & filepath )
{
 BML_LOG( "SaveModel: saving to \"" << filepath << "\"" << std::endl );
 try {
  nn->eval();
  torch::serialize::OutputArchive archive;
  nn->save( archive );
  archive.save_to( filepath );
  nn->train();
  }
 catch( const c10::Error & e ) {
  nn->train();
  throw( std::runtime_error(
	  std::string( "BundleSolverML::SaveModel: Torch error: " ) +
	  e.what() ) );
  }
 catch( const std::exception & e ) {
  nn->train();
  throw( std::runtime_error(
	  std::string( "BundleSolverML::SaveModel: " ) + e.what() ) );
  }
 }

/*--------------------------------------------------------------------------*/

void BundleSolverML::LoadModel( const std::string & filepath )
{
 BML_LOG( "LoadModel: loading from \"" << filepath << "\"" << std::endl );
 try {
  torch::serialize::InputArchive archive;
  archive.load_from( filepath );
  nn->load( archive );
  nn->eval();
  }
 catch( const c10::Error & e ) {
  throw( std::runtime_error(
	  std::string( "BundleSolverML::LoadModel: Torch error: " ) +
	  e.what() ) );
  }
 catch( const std::exception & e ) {
  throw( std::runtime_error(
	  std::string( "BundleSolverML::LoadModel: " ) + e.what() ) );
  }
 }

/*--------------------------------------------------------------------------*/
/*----------------------- ML-SPECIFIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

double BundleSolverML::Heuristic( Index whch )
{
 if( ! nn->is_training() )
  nn->train();

 /* Feature normalization: the raw values span many orders of magnitude,
  * hence each feature is divided by a per-feature scale that keeps the
  * inputs O( 1 ) at initialization, preventing the Softplus output from
  * being huge with random weights. The normalized tensor is what is stored
  * in phi_vecs, so that Backward() re-runs the forward pass with the same
  * scaled inputs that produced the original prediction. */
 features[ 0 ] = float( tHasChgd );             // binary { 0 , 1 } 
if( G1Norm == INFshift ) {    // Heuristic4() computes this lazily, but
 double n2 = 0;               // BundleSolverML replaces Heuristic4()
 for( auto gi : G1 )
  n2 += double( gi ) * double( gi );
G1Norm = std::sqrt( n2 );
  }
 features[ 1 ] = float( G1Norm ) / 1e3f;        // subgradient norm
 features[ 2 ] = float( ScPr1 ) / 1e4f;         // scalar product
 features[ 3 ] = float( Alfa1 ) / 1e4f;         // linearization error
 features[ 4 ] = 0.0f;
 features[ 5 ] = float( t ) / 1e3f;             // current step-size
 features[ 6 ] = float( Sigma ) / 1e4f;         // predicted decrease
 features[ 7 ] = float( DSTS ) / 1e3f;          // distance to center
 features[ 8 ] = 0.0f;
 features[ 9 ] = float( EpsU ) * 1e3f;          // EpsU ~ 1e-3 ==> * 1e3
 features[ 10 ] = float( CSSCntr ) / 100.0f;    // SS counter
 features[ 11 ] = float( CNSCntr ) / 100.0f;    // NS counter
 features[ 12 ] = float( UpTrgt ) / 1e5f;       // upper target
 features[ 13 ] = float( LwTrgt ) / 1e5f;       // lower target
 features[ 14 ] = float( Fi0Lmb ) / 1e5f;       // f at tentative point
 features[ 15 ] = float( Fi0Lmb1 ) / 1e5f;      // f at center
 features[ 16 ] = float( DST ) / 1e3f;          // distance to center (alt)
 features[ 17 ] = float( NrmD ) / 1e3f;         // norm of the direction
 features[ 18 ] = float( NrmZ ) / 1e3f;         // norm of agg. subgradient
 features[ 19 ] = float( NrmZFctr ) / 1e3f;     // normalization factor

 /* Saturation guard: clamp the normalized features to [ -10 , 10 ], as a
  * single wildly out-of-scale feature (e.g., EpsU at the very first
  * iterations) would otherwise saturate the Softplus units, flattening
  * the prediction onto its minimum and killing the gradients, with no way
  * for the training to ever recover. */
 for( auto & fi : features )
  fi = std::min( std::max( fi , -10.0 ) , 10.0 );

 /* The prediction needs no gradient, since Backward() evaluates the network
  * again on the recorded features; this also keeps the state of a recurrent
  * core from chaining the autograd graphs of all the iterations. */
 auto input = torch::tensor( features );
 torch::Tensor output;
 {
  torch::NoGradGuard no_grad;
  output = nn->forward( input );
  }

 double t_pred = output.item< double >();
 if( ( ! std::isfinite( t_pred ) ) || ( t_pred <= 0 ) ) {
  BML_LOG( "Heuristic: invalid t = " << t_pred << ", clamping to 1"
	   << std::endl );
  t_pred = 1;
  }

 BML_LOG( "Heuristic: predicted t = " << t_pred << std::endl );

 /* Nothing below is used unless the trajectory is replayed in Backward(),
  * i.e., unless the training is on and this iteration is in the window
  * [ intMLIterFirst , intMLIterLast ]. */
 if( ( ! f_train_online ) || ( ++f_ML_iter < f_ML_iter_first ) ||
     ( f_ML_iter > f_ML_iter_last ) )
  return( t_pred );

 /* What is recorded is the prediction, the step-size the current trial
  * point was computed with, and the outcome of that step: the cosine of
  * the angle between the new subgradient g (the sum of those of the
  * components evaluated at the trial point) and the step d = \bar{x}_+ -
  * \bar{x} from the center to the trial point. By convexity
  * f( \bar{x}_+ ) - f( \bar{x} ) <= < g , d >, and along the path d( t )
  * of the proximal Master problem d f( \bar{x} + d( t ) ) / d t ~
  * < g , d > / t, so a positive < g , d > says that a smaller t would have
  * done better and a negative one that a larger t would. This holds
  * whatever the sign of f, since g is the subgradient of the convex
  * function the algorithm actually minimizes, i.e., of - f for a concave
  * one [see Backward()]. */
 phi_vecs.push_back( input );
 tS.push_back( t_pred );
 tU.push_back( t );

 double gd = 0;
 double nd = 0;
 for( Index i = 0 ; i < NumVar ; ++i ) {
  const double di = Lambda1[ i ] - Lambda[ i ];
  gd += G1[ i ] * di;
  nd += di * di;
  }
 nd = std::sqrt( nd );

 double cs = 0;
 if( ( G1Norm > 0 ) && ( nd > 0 ) && std::isfinite( gd ) )
  cs = std::min( std::max( gd / ( G1Norm * nd ) , -1.0 ) , 1.0 );
 csS.push_back( cs );

 return( t_pred );

 }  // end( BundleSolverML::Heuristic )

/*--------------------------------------------------------------------------*/
/* The prediction t_f made at the recorded iteration f is the step-size of
 * the trial point of the next recorded iteration f + 1, provided nothing
 * else has changed t in between (a noise reduction step, the bounds
 * [ tMinor , tMaior ], the counters intMnSSC / intMnNSC, ...), which is
 * checked by comparing t_f with the step-size actually used there. The
 * outcome of that step is the cosine c_{f+1} recorded at f + 1 [see
 * Heuristic()], and the loss is
 *
 *   L = ( 1 / n ) \sum_f c_{f+1} \log t_f
 *
 * over the n iterations f whose prediction was used. Since c_{f+1} is a
 * constant, the gradient of L w.r.t. the parameters of the network is the
 * mean of c_{f+1} d \log t_f, hence an Adam step decreases t_f where the
 * step went too far ( c_{f+1} > 0 ) and increases it where it did not go
 * far enough ( c_{f+1} < 0 ); the logarithm makes the update independent
 * of the scale of t, and the cosine of the scale of f and of its
 * subgradients, so that one network can be trained on instances of
 * different size. The stationary point is where, on average, the trial
 * points are minima of f along the path d( t ) (i.e., < g , d > = 0), which
 * is the step-size an exact line search along that path would give.
 *
 * The recorded iterations are split in windows of intMLWindow iterations,
 * each one giving an Adam step. A recurrent core is replayed on all the
 * recorded iterations, in their order and starting from a null state, so
 * that at each of them it has the state it had when the prediction was
 * made; the state crosses the windows, but the gradient does not
 * (truncated back-propagation through time). After loss.backward() the
 * gradients are clipped to GRAD_CLIP_NORM. The optimizer is kept as a field
 * so that its moment estimates persist across successive Backward() calls
 * on the same (possibly shared) network. */

void BundleSolverML::Backward( void )
{
 BML_LOG( "Backward: starting (" << phi_vecs.size()
	  << " iterations recorded)" << std::endl );

 try {
  /* The Adam optimizer is lazily initialized at the first Backward() call,
   * which means it is always bound to the *active* network (nn), whether
   * it is the private one or a shared network injected via
   * set_shared_net(). */
  if( ! f_optimizer ) {
   auto opts = torch::optim::AdamOptions( 1e-4 ).betas( { 0.9 , 0.999 } )
                          .eps( 1e-8 ).weight_decay( 1e-5 );
   f_optimizer = std::make_unique< torch::optim::Adam >( nn->parameters() ,
							 opts );
   }

  const size_t n_rec = phi_vecs.size();
  if( n_rec < 2 ) {
   BML_LOG( "Backward: no prediction with a known outcome, skipping"
	    << std::endl );
   return;
   }

  if( ! nn->is_training() )
   nn->train();

  const size_t k_win = ( f_ML_window >= int( n_rec ) )
                       ? n_rec : size_t( f_ML_window );

  nn->reset_state();

  for( size_t wstart = 0 ; wstart < n_rec ; wstart += k_win ) {
   const size_t wend = std::min( wstart + k_win , n_rec );
   nn->detach_state();
   f_optimizer->zero_grad();

   torch::Tensor loss = torch::zeros( {} , torch::kFloat32 );
   size_t n_used = 0;

   for( size_t f = wstart ; f < wend ; ++f ) {
    // the network is run on every iteration, as each one moves the state
    // of a recurrent core
    auto t_f = nn->forward( phi_vecs[ f ] ).sum();

    if( f + 1 >= n_rec )  // the last prediction has no outcome
     continue;

    if( std::abs( tU[ f + 1 ] - tS[ f ] ) >
	1e-9 * std::max( 1.0 , std::abs( tS[ f ] ) ) )
     continue;            // t_f has not been used as such

    if( csS[ f + 1 ] == 0 )
     continue;

    loss = loss + float( csS[ f + 1 ] ) * torch::log( t_f + 1e-12f );
    ++n_used;
    }

   if( ! n_used ) {
    BML_LOG( "Backward: window [ " << wstart << " , " << wend
	     << " ) has no used prediction, skipping" << std::endl );
    continue;
    }

   loss = loss / float( n_used );

   BML_LOG( "Backward: window [ " << wstart << " , " << wend
	    << " ), loss = " << loss.item< float >() << " over " << n_used
	    << " predictions" << std::endl );

   loss.backward();

   // clip the gradients so that their global L2 norm is <= GRAD_CLIP_NORM
   auto pre_clip_norm = torch::nn::utils::clip_grad_norm_(
					     nn->parameters() , GRAD_CLIP_NORM );
   BML_LOG( "Backward: gradient norm (pre-clip) = " << pre_clip_norm
	    << std::endl );

   f_optimizer->step();
   }  // end( for( wstart ) )
  }
 catch( const std::exception & e ) {
  std::cerr << "Backward: exception: " << e.what() << std::endl;
  }
 catch( ... ) {
  std::cerr << "Backward: unknown exception" << std::endl;
  }

 BML_LOG( "Backward: done" << std::endl );

 }  // end( BundleSolverML::Backward )

/*--------------------------------------------------------------------------*/
/*------------------- End File BundleSolverML.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
