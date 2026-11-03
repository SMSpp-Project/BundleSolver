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

static constexpr auto InINF = SMSpp_di_unipi_it::Inf< Index >();

/*--------------------------------------------------------------------------*/
/*----------------------------- FUNCTIONS ----------------------------------*/
/*--------------------------------------------------------------------------*/

static void chgsign( double * v , Index n )
{
 for( const auto ev = v + n ; v < ev ; ++v )
  *v = - *v;
 }

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS BundleSolverML_Fi -------------------------*/
/*--------------------------------------------------------------------------*/
/* Custom autograd function used to attach the discounted search-direction
 * accumulator w_cum to the loss computation graph.
 *
 * The signature is apply( input , gs_scalar , fi_val ), where
 *
 * - input is the tensor whose gradient is sought (w_cum, of shape { N });
 *
 * - gs_scalar is a scalar tensor (gs.sum()), used as a constant multiplier
 *   in the backward pass, i.e., it is not differentiated;
 *
 * - fi_val is a double constant, not differentiated.
 *
 * The forward pass returns fi_val as a scalar (shape { }) with the same
 * dtype / device as input.
 *
 * The backward pass returns a null gradient for all the inputs: w_cum does
 * not require_grad() itself, as the gradient path towards the network
 * parameters runs through nn_out (which multiplies the output of this
 * function in the loss term of Backward()), not through w_cum directly.
 * Returning a null gradient is correct and avoids any shape mismatch
 * between w_cum (whose shape is set at the first recorded iteration) and
 * the bundle size at later iterations. */

struct BundleSolverML_Fi
 : public torch::autograd::Function< BundleSolverML_Fi >
{
 static torch::Tensor forward( torch::autograd::AutogradContext * ctx ,
			       torch::Tensor input ,
			       torch::Tensor gs_scalar ,
			       double fi_val ) {
  // save gs_scalar (shape { }) for backward: the shape is independent of
  // the bundle size
  ctx->save_for_backward( { gs_scalar } );
  // return fi_val as a scalar tensor matching input's dtype / device
  return( torch::tensor( fi_val , input.options() ).squeeze() );
  }

 static torch::autograd::tensor_list backward(
				torch::autograd::AutogradContext * ctx ,
				torch::autograd::tensor_list grad_outputs ) {
  return { torch::Tensor() ,   // grad w.r.t. input (w_cum), not needed
	   torch::Tensor() ,   // grad w.r.t. gs_scalar, a constant
	   torch::Tensor() };  // grad w.r.t. fi_val, a constant
  }
 };

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS BundleSolverML_W --------------------------*/
/*--------------------------------------------------------------------------*/
/* Custom autograd function for the projected search direction w.
 *
 * The forward pass returns sum( w ) so as to attach the computation graph;
 * the backward pass computes the gradient w.r.t. the step-size t via the
 * KKT projection
 *
 *   \f[ \lambda^* = ( e^T Q^{-1} \alpha ) / ( e^T Q^{-1} e ) \f]
 *
 * out of which
 *
 *   \f[ dw / dt = ( 1 / t^2 ) G^T Q^{-1} ( \lambda^* e - \alpha ) \f] */

struct BundleSolverML_W
 : public torch::autograd::Function< BundleSolverML_W >
{
 static torch::Tensor forward( torch::autograd::AutogradContext * ctx ,
			       torch::Tensor w , torch::Tensor G ,
			       torch::Tensor Q , torch::Tensor alpha ,
			       torch::Tensor E , double ts ) {
  torch::Tensor ts_tensor = torch::tensor( ts , Q.options() );
  ctx->save_for_backward( { G , Q , alpha , E , ts_tensor } );
  return( w.sum() );
  }

 static torch::autograd::tensor_list backward(
				torch::autograd::AutogradContext * ctx ,
				torch::autograd::tensor_list grad_outputs ) {
  auto saved = ctx->get_saved_variables();
  torch::Tensor G = saved[ 0 ];
  torch::Tensor Q = saved[ 1 ];
  torch::Tensor alpha = saved[ 2 ];
  torch::Tensor E = saved[ 3 ];
  torch::Tensor ts = saved[ 4 ];
  torch::Tensor grad_out = grad_outputs[ 0 ];

  // regularize Q and compute its inverse via Cholesky factorization
  auto I = torch::eye( Q.size( 0 ) , Q.options() );
  auto Qreg = Q + 1e+1 * I;
  auto L = torch::linalg_cholesky( Qreg );
  auto Qinv = torch::cholesky_inverse( L );

  // KKT projection: lambda* = ( e^T Q^{-1} alpha ) / ( e^T Q^{-1} e )
  auto Qinvalpha = torch::matmul( Qinv , alpha );
  auto QinvEt = torch::matmul( Qinv , E.t() );
  auto M = torch::matmul( E , QinvEt );
  auto Ik = torch::eye( M.size( 0 ) , M.options() );
  auto Mreg = M + 1e-8 * Ik;
  auto ek = torch::ones( { M.size( 0 ) } , Q.options() );
  auto rhs = ek + torch::matmul( E , Qinvalpha );
  auto lambda = torch::linalg_solve( Mreg , rhs );
  auto proj = torch::matmul( E.t() , lambda ) - alpha;

  // dw / dt = ( 1 / t^2 ) G^T Q^{-1} proj
  double ts_val = ts.item< double >();
  auto result = torch::matmul( G.transpose( 0 , 1 ) ,
			       torch::matmul( Qinv , proj ) );
  result = ( 1.0 / ( ts_val * ts_val ) ) * result * grad_out;

  return { result , torch::Tensor() , torch::Tensor() , torch::Tensor() ,
	   torch::Tensor() , torch::Tensor() };
  }
 };

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------- FACTORY MANAGEMENT -----------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it
{
 SMSpp_insert_in_factory_cpp_0( BundleSolverML );
 }

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE Block ----------------------*/
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
/*-------------------- METHODS FOR HANDLING THE MODEL ----------------------*/
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
/*----------------------- ML-SPECIFIC METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

double BundleSolverML::Heuristic( Index whch )
{
 if( ! nn->is_training() )
  nn->train();

 if( G1Norm == INFshift ) {  // Heuristic4() computes it lazily, but
  double n2 = 0;             // BundleSolverML replaces Heuristic4()
  for( auto gi : G1 )
   n2 += double( gi ) * double( gi );
  G1Norm = std::sqrt( n2 );
  }

 // < g , d > of the new subgradient and the step from the center to the
 // trial point, and the norm of the step [see the recording below]
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

 if( f_ML_features )
  dimensionless_features( cs , nd );
 else
  scaled_features();

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

 /* With intMLFeatures == 0 the output of the network is t itself. With
  * intMLFeatures == 1 it is tau, dimensionless as its input, and t is the
  * step-size with which the decrease the model predicts, t || z* ||^2, is
  * tau times the current estimate of the gap G = Sigma* + D_t*( z* ):
  * t = tau G / || z* ||^2. Since G and || z* || go to 0 together as the
  * algorithm converges, t does not drift as a product of factors would.
  * tau is the output divided by log( 2 ), the output of the final Softplus
  * at 0, so that tau = 1 for a network with null weights; t is left as it
  * is when G or || z* || is null, i.e., at the first iterations. Either
  * way t = s * output, with the scale s recorded for Backward() (0 when t
  * does not depend on the output). */
 double t_pred = t;
 double scale = 0;
 if( f_ML_features ) {
  const double G = Sigma + DSTS;
  const double out = output.item< double >();
  if( std::isfinite( out ) && ( out > 0 ) && ( G > 0 ) && ( NrmZ > 0 ) &&
      ( G < INFshift ) ) {
   scale = G / ( std::log( 2.0 ) * double( NrmZ ) * double( NrmZ ) );
   t_pred = scale * out;
   }
  else
   BML_LOG( "Heuristic: output = " << out << ", G = " << G
	    << ", || z* || = " << NrmZ << ", t unchanged" << std::endl );
  }
 else {
  t_pred = output.item< double >();
  scale = 1;
  if( ( ! std::isfinite( t_pred ) ) || ( t_pred <= 0 ) ) {
   BML_LOG( "Heuristic: invalid t = " << t_pred << ", clamping to 1"
	    << std::endl );
   t_pred = 1;
   }
  }

 BML_LOG( "Heuristic: predicted t = " << t_pred << std::endl );

 /* Nothing below is used unless the trajectory is replayed in Backward(),
  * i.e., unless the training is on and this iteration is in the window
  * [ intMLIterFirst , intMLIterLast ]. */
 if( ( ! f_train_online ) || ( ++f_ML_iter < f_ML_iter_first ) ||
     ( f_ML_iter > f_ML_iter_last ) )
  return( t_pred );

 if( f_ML_loss ) {
  /* The step loss needs the prediction, the step-size the current trial
   * point was computed with, and the outcome of that step: the cosine cs
   * of the angle between the new subgradient g (the sum of those of the
   * components evaluated at the trial point) and the step d = \bar{x}_+ -
   * \bar{x} from the center to the trial point. By convexity
   * f( \bar{x}_+ ) - f( \bar{x} ) <= < g , d >, and along the path d( t )
   * of the proximal Master problem d f( \bar{x} + d( t ) ) / d t ~
   * < g , d > / t, so a positive < g , d > says that a smaller t would
   * have done better and a negative one that a larger t would. This holds
   * whatever the sign of f, since g is the subgradient of the convex
   * function the algorithm actually minimizes, i.e., of - f for a concave
   * one [see Backward()]. */
  phi_vecs.push_back( input );
  tS.push_back( t_pred );
  sS.push_back( scale );
  tU.push_back( t );
  csS.push_back( cs );
  return( t_pred );
  }

 /* The objective loss needs the tensors out of which the search direction
  * is differentiated w.r.t. t [see BundleSolverML_W]. */
 phi_vecs.push_back( input );

 // search direction w
 auto tZ = MasterPB->get_z_vector();
 w_vecs.push_back( torch::tensor( tZ ).requires_grad_( true ) );

 // how many items the bundle holds now, which is how many rows G has
 int col_num = 0;
 for( Index i = 0 ; i < get_max_name() ; ++i )
  if( is_bundle_item( i ) )
   col_num++;

 // subgradient matrix G, one row for each item in the bundle; the global
 // name of each of them is kept to read its linearization error below
 std::vector< Index > G_names;
 std::vector< std::vector< VarValue > > G_mat;
 std::vector< Index > comp_of_row;
 comp_of_row.reserve( col_num );

 G_mat.reserve( col_num );
 for( Index i = 0 ; i < get_max_name() ; ++i )
  if( is_bundle_item( i ) ) {
   std::vector< VarValue > G( NumVar );
   lin_coeffs( ItemVcblr[ i ].first ,
	       G.data() , Range( 0 , NumVar ) , ItemVcblr[ i ].second );
   if( ! f_convex )
    chgsign( G.data() , NumVar );
   G_mat.push_back( std::move( G ) );
   G_names.push_back( i );
   comp_of_row.push_back( ItemVcblr[ i ].first );
   }

 if( G_mat.empty() ) {
  // no linearization in the bundle yet: nothing to differentiate through,
  // drop the entries recorded so far for this iteration
  phi_vecs.pop_back();
  w_vecs.pop_back();
  return( t_pred );
  }

 int rows = G_mat.size();
 int cols = G_mat[ 0 ].size();
 std::vector< double > flat;
 flat.reserve( rows * cols );
 for( const auto & row : G_mat )
  flat.insert( flat.end() , row.begin() , row.end() );

 torch::Tensor G_tensor = torch::from_blob( flat.data() , { rows , cols } ,
					    torch::kDouble ).clone();
 /* Indicator matrix E of the convexity constraints: one row per component
  * that the bundle has an item of, compressed over the components actually
  * present, which is what the projection in the backward pass solves with
  * [see BundleSolverML_W::backward()]. */
 std::map< Index , long > comp2row;
 for( auto c : comp_of_row )
  if( ! comp2row.count( c ) ) {
   long next = long( comp2row.size() );
   comp2row[ c ] = next;
   }

 auto E_tensor = torch::zeros( { long( comp2row.size() ) , long( rows ) } ,
			       torch::kDouble );
 for( size_t r = 0 ; r < comp_of_row.size() ; ++r )
  E_tensor[ comp2row[ comp_of_row[ r ] ] ][ long( r ) ] = 1.0;

 Es.push_back( E_tensor );

 Gs.push_back( G_tensor );

 // aggregated subgradient (whisG1 indices)
 std::vector< int64_t > valid_indices;
 for( auto k : whisG1 ) {
  if( k == InINF )
   break;
  if( Index( k ) < Index( cols ) )
   valid_indices.push_back( int64_t( k ) );
  }

 torch::Tensor grad_sum;
 if( valid_indices.empty() )
  grad_sum = torch::zeros( { rows } , torch::kDouble );
 else {
  auto idx = torch::tensor( valid_indices , torch::kLong );
  grad_sum = G_tensor.index_select( 1 , idx ).sum( 1 );
  }
 Gs_aggreg.push_back( grad_sum );

 /* Step-type coefficient: +1 SS, -1 NS. Note that SSDone already reflects
  * the type of the current step when Heuristic() is called, while the
  * SS / NS counters do not, as they are updated later. */
 coeff_vecs.push_back( torch::tensor( SSDone ? 1.0f : -1.0f ) );

 // Gram matrix Q and linearization errors alpha
 torch::Tensor Q = torch::matmul( G_tensor , G_tensor.transpose( 0 , 1 ) );
 Qs.push_back( Q );
 std::vector< double > alpha( G_names.size() , 0 );
 for( Index i = 0 ; i < G_names.size() ; ++i )
  alpha[ i ] = read_alpha_global( G_names[ i ] );
 alphaS.push_back( torch::tensor( alpha , torch::kDouble ) );

 tS.push_back( t_pred );
 sS.push_back( scale );
 FiS.push_back( UpFiBest );

 return( t_pred );

 }  // end( BundleSolverML::Heuristic )

/*--------------------------------------------------------------------------*/
/* The features of intMLFeatures == 0: each one is divided by a fixed scale
 * that keeps it O( 1 ) on the instances the network was first tried on,
 * which prevents the Softplus output from being huge with random weights.
 */

void BundleSolverML::scaled_features( void )
{
 features[ 0 ] = float( tHasChgd );             // binary { 0 , 1 }
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
 }

/*--------------------------------------------------------------------------*/
/* The features of intMLFeatures == 1 are dimensionless, so that the network
 * sees the same values on instances of different size and scale: the
 * values of the function are relative to Fs = max( | f( \bar{x} ) | , 1 )
 * (the same normalization of EpsU), the step-size enters through
 * D_t( z* ) / Fs, i.e., the decrease the model predicts relative to f, and
 * the norms only through their ratios; the logarithms are taken of the
 * quantities that span orders of magnitude. */

void BundleSolverML::dimensionless_features( double cs , double nd )
{
 const double fc = UpRifFi.back();
 const double Fs = ( std::abs( fc ) < INFshift ) ?
                   std::max( std::abs( fc ) , 1.0 ) : 1.0;
 const auto rel = [ Fs ]( double v ) {
  return( std::abs( v ) < INFshift ? v / Fs : 0.0 );
  };
 const auto lg = []( double v ) {
  return( ( v > 0 ) && ( v < INFshift ) ? std::log10( v ) : -10.0 );
  };

 features[ 0 ] = float( tHasChgd );                    // t has changed
 features[ 1 ] = NrmZ > 0 ? lg( G1Norm / NrmZ ) : 0;   // || g || / || z* ||
 features[ 2 ] = cs;                                    // cos( g , d )
 features[ 3 ] = rel( Alfa1 );                          // linearization error
 features[ 4 ] = 0;
 features[ 5 ] = lg( DST / Fs );                        // D_t( z* ) / Fs
 features[ 6 ] = rel( Sigma );                          // Sigma* / Fs
 features[ 7 ] = rel( DSTS );                           // D_t*( z* ) / Fs
 features[ 8 ] = 0;
 features[ 9 ] = lg( EpsU );                            // relative gap
 features[ 10 ] = CSSCntr / 100.0;                      // SS counter
 features[ 11 ] = CNSCntr / 100.0;                      // NS counter
 features[ 12 ] = rel( UpTrgt - fc );                   // upper target
 features[ 13 ] = rel( LwTrgt - fc );                   // lower target
 features[ 14 ] = rel( UpFiLmb1.back() - fc );          // change of f
 features[ 15 ] = rel( ScPr1 );                         // < g , d* > / Fs
 features[ 16 ] = ( t > 0 ) && ( NrmZ > 0 ) ?          // || d || / t || z* ||
                  nd / ( t * NrmZ ) : 0;
 features[ 17 ] = NrmD > 0 ? NrmDInf / NrmD : 0;        // shape of d
 features[ 18 ] = 0;
 features[ 19 ] = NrmZFctr;                             // normalization factor
 }

/*--------------------------------------------------------------------------*/

torch::Tensor BundleSolverML::w( size_t f , double t )
{
 return( BundleSolverML_W::apply( w_vecs[ f ] , Gs[ f ] , Qs[ f ] ,
				  alphaS[ f ] , Es[ f ] , t ) );
 }

/*--------------------------------------------------------------------------*/
/* The recorded iterations are split in windows of intMLWindow iterations,
 * each one giving an Adam step on the loss of intMLLoss [see
 * objective_loss() and step_loss()]. A recurrent core is replayed on all
 * the recorded iterations, in their order and starting from a null state,
 * so that at each of them it has the state it had when the prediction was
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
   * set_shared_net(). lr = 1e-4 is a conservative choice for this loss
   * scale. */
  if( ! f_optimizer ) {
   auto opts = torch::optim::AdamOptions( 1e-4 ).betas( { 0.9 , 0.999 } )
                          .eps( 1e-8 ).weight_decay( 1e-5 );
   f_optimizer = std::make_unique< torch::optim::Adam >( nn->parameters() ,
							 opts );
   }

  const size_t n_rec = phi_vecs.size();
  if( n_rec < ( f_ML_loss ? 2 : 1 ) ) {
   BML_LOG( "Backward: nothing to learn from, skipping" << std::endl );
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

   size_t n_used = 0;
   torch::Tensor loss = f_ML_loss ? step_loss( wstart , wend , n_used )
                                  : objective_loss( wstart , wend , n_used );
   if( ! n_used ) {
    BML_LOG( "Backward: window [ " << wstart << " , " << wend
	     << " ) has nothing in the loss, skipping" << std::endl );
    continue;
    }

   // mean reduction: the loss magnitude must not grow with the iterations
   loss = loss / float( n_used );

   BML_LOG( "Backward: window [ " << wstart << " , " << wend
	    << " ), mean loss = " << loss.item< float >() << " over "
	    << n_used << " iterations" << std::endl );

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
/* The loss of intMLLoss == 0 is constructed as follows: for each recorded
 * iteration f that either is a SS or is the last one of the window
 *
 *   t_f    = s_f nn->forward( phi_vecs[ f ] ).sum()    [ re-evaluated ]
 *   w_curr = BundleSolverML_W::apply( ... )            [ scalar ]
 *   w_cum += 0.9^( total - ss_count ) * t_f * w_curr
 *   loss  += coeff * t_f * BundleSolverML_Fi( w_cum , gs_scalar , fi )
 *
 * with s_f the scale of the prediction [see Heuristic()], i.e., 1 with
 * intMLFeatures == 0, coeff = +1 for a SS and -1 for a NS, and fi the best
 * value of the function found so far. Since BundleSolverML_Fi returns the
 * constant fi and a null gradient, and BundleSolverML_W receives t_f as a
 * double, the gradient of the term w.r.t. the parameters of the network is
 * coeff * fi * d t_f, i.e., the loss is linear in t with the sign of fi.
 * The mean (rather than the sum) is taken, and fi is not multiplied in, so
 * that the magnitude of the loss does not grow with the iterations and with
 * the scale of f. */

torch::Tensor BundleSolverML::objective_loss( size_t wstart , size_t wend ,
					      size_t & n_used )
{
 const bool rec = ( nn->opt.model_type != NetOptions::eMLP );
 torch::Tensor loss = torch::zeros( {} , torch::kFloat32 );
 const size_t last_idx = wend - 1;

 // scalar accumulator for the discounted search directions
 torch::Tensor w_cum = torch::zeros( {} , torch::kFloat32 );

 for( size_t f = wstart ; f < wend ; ++f ) {
  const bool is_ss = ( coeff_vecs[ f ].item< float >() > 0 );
  if( ( ( ! is_ss ) && ( f != last_idx ) ) || ( sS[ f ] == 0 ) ) {
   if( rec )  // not in the loss, but it moves the state
    nn->forward( phi_vecs[ f ] );
   continue;
   }
  ++n_used;

  if( ( ! w_vecs[ f ].defined() ) || ( ! Gs[ f ].defined() ) ||
      ( ! Qs[ f ].defined() ) || ( f >= alphaS.size() ) ||
      ( ! alphaS[ f ].defined() ) ) {
   std::cerr << "Backward: missing tensor at iteration " << f << std::endl;
   if( rec )
    nn->forward( phi_vecs[ f ] );
   continue;
   }

  try {
   std::unique_ptr< torch::NoGradGuard > no_grad;
   if( ! f_train_online )
    no_grad = std::make_unique< torch::NoGradGuard >();

   auto nn_out = nn->forward( phi_vecs[ f ] ).sum();  // scalar, float32
   if( ! nn_out.requires_grad() ) {
    std::cerr << "Backward: no grad_fn on nn_out at iteration " << f
	      << std::endl;
    continue;
    }
   if( sS[ f ] != 1 )
    nn_out = float( sS[ f ] ) * nn_out;

   auto w_curr = BundleSolverML_W::apply( w_vecs[ f ] , Gs[ f ] , Qs[ f ] ,
					  alphaS[ f ] , Es[ f ] ,
					  nn_out.item< double >() );

   const double discount = std::pow( 0.9 ,
				     double( ( wend - wstart ) - n_used ) );
   w_cum = w_cum + float( discount ) * nn_out * w_curr;

   auto gs_scalar = Gs_aggreg[ f ].to( torch::kFloat32 ).sum();

   auto term = coeff_vecs[ f ] * nn_out *
               BundleSolverML_Fi::apply( w_cum , gs_scalar , FiS[ f ] );
   loss += term;

   BML_LOG( "Backward: iteration " << f << ( is_ss ? " [SS]" : " [last]" )
	    << ", discount = " << discount << ", t = "
	    << nn_out.item< float >() << ", contribution = "
	    << term.item< float >() << std::endl );
   }
  catch( const std::exception & e ) {
   std::cerr << "Backward: error at iteration " << f << ": " << e.what()
	     << std::endl;
   }
  }  // end( for( f ) )

 return( loss );
 }

/*--------------------------------------------------------------------------*/
/* The loss of intMLLoss == 1 pairs the prediction t_f made at the recorded
 * iteration f, which gives the step-size of the trial point of the next
 * recorded iteration f + 1, with the outcome of that step, i.e., the cosine
 * c_{f+1} recorded there [see Heuristic()]:
 *
 *   L = ( 1 / n ) \sum_f c_{f+1} \log t_f
 *
 * over the n iterations f that take part in it. The step-size actually
 * used may differ from t_f, since the algorithm keeps t within bounds
 * (tMinor, tMaior and the interval allowed by the counters intMnSSC /
 * intMnNSC) and changes it in a noise reduction step: then c_{f+1} says
 * how the t that was used should have changed, and the prediction takes
 * part in L only if c_{f+1} moves it towards that t, as moving it the
 * other way would change nothing in the algorithm. Since c_{f+1} is a
 * constant, the gradient of L w.r.t. the parameters of the network is the
 * mean of c_{f+1} d \log t_f, hence an Adam step decreases t_f where the
 * step went too far ( c_{f+1} > 0 ) and increases it where it did not go
 * far enough ( c_{f+1} < 0 ); the logarithm makes the update independent
 * of the scale of t (and of the scale s_f of the prediction), and the
 * cosine of the scale of f and of its subgradients. The stationary point
 * is where, on average, the trial points are minima of f along the path
 * d( t ) (i.e., < g , d > = 0), which is the step-size an exact line search
 * along that path would give. */

torch::Tensor BundleSolverML::step_loss( size_t wstart , size_t wend ,
					 size_t & n_used )
{
 const size_t n_rec = phi_vecs.size();
 torch::Tensor loss = torch::zeros( {} , torch::kFloat32 );

 for( size_t f = wstart ; f < wend ; ++f ) {
  // the network is run on every iteration, as each one moves the state of
  // a recurrent core; d log t_f = d log of the output, s_f being fixed
  auto t_f = nn->forward( phi_vecs[ f ] ).sum();

  if( f + 1 >= n_rec )  // the last prediction has no outcome
   continue;

  const double c = csS[ f + 1 ];
  if( ( c == 0 ) || ( sS[ f ] == 0 ) )
   continue;

  // a t_f that has not been used as such (the bounds on t, a noise
  // reduction, ...) is only moved towards the t that was used
  const double gap = tU[ f + 1 ] - tS[ f ];
  if( std::abs( gap ) > 1e-9 * std::max( 1.0 , std::abs( tS[ f ] ) ) )
   if( ( gap < 0 ) != ( c > 0 ) )
    continue;

  loss = loss + float( c ) * torch::log( t_f + 1e-12f );
  ++n_used;
  }

 return( loss );
 }

/*--------------------------------------------------------------------------*/
/*------------------- End File BundleSolverML.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
