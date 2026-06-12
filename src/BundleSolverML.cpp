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

#include <cmath>
#include <iostream>
#include <stdexcept>

#include "BundleSolverML.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

using HpNum = NDO_di_unipi_it::HpNum;
using cHpRow = NDO_di_unipi_it::cHpRow;
using Index = Function::Index;

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

static constexpr auto InINF = SMSpp_di_unipi_it::Inf< Index >();

/// maximum L2 norm for gradient clipping applied after loss.backward();
/// if the norm of the gradients exceeds this value, all the gradients are
/// rescaled by GRAD_CLIP_NORM / norm (typical values are 1 ... 10)
static constexpr float GRAD_CLIP_NORM = 1.0f;

/*--------------------------------------------------------------------------*/
/*----------------------------- FUNCTIONS ----------------------------------*/
/*--------------------------------------------------------------------------*/

static void chgsign( double * v , Index n )
{
 for( const auto ev = v + n ; v < ev ; ++v )
  *v = - *v;
 }

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS BundleSolverML_Fi --------------------------*/
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
/*------------------------ CLASS BundleSolverML_W ---------------------------*/
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
			       double ts ) {
  torch::Tensor ts_tensor = torch::tensor( ts , Q.options() );
  ctx->save_for_backward( { G , Q , alpha , ts_tensor } );
  return( w.sum() );
  }

 static torch::autograd::tensor_list backward(
				torch::autograd::AutogradContext * ctx ,
				torch::autograd::tensor_list grad_outputs ) {
  auto saved = ctx->get_saved_variables();
  torch::Tensor G = saved[ 0 ];
  torch::Tensor Q = saved[ 1 ];
  torch::Tensor alpha = saved[ 2 ];
  torch::Tensor ts = saved[ 3 ];
  torch::Tensor grad_out = grad_outputs[ 0 ];

  // regularize Q and compute its inverse via Cholesky factorization
  auto I = torch::eye( Q.size( 0 ) , Q.options() );
  auto Qreg = Q + 1e+1 * I;
  auto L = torch::linalg_cholesky( Qreg );
  auto Qinv = torch::cholesky_inverse( L );

  // KKT projection: lambda* = ( e^T Q^{-1} alpha ) / ( e^T Q^{-1} e )
  auto e = torch::ones( { Q.size( 0 ) } , Q.options() );
  auto Qinvalpha = torch::matmul( Qinv , alpha );
  auto Qinv_e = torch::matmul( Qinv , e );
  auto num = torch::dot( e , Qinvalpha );
  auto den = torch::dot( e , Qinv_e );
  auto proj = ( num / den ) * e - alpha;

  // dw / dt = ( 1 / t^2 ) G^T Q^{-1} proj
  double ts_val = ts.item< double >();
  auto result = torch::matmul( G.transpose( 0 , 1 ) ,
			       torch::matmul( Qinv , proj ) );
  result = ( 1.0 / ( ts_val * ts_val ) ) * result * grad_out;

  return { result , torch::Tensor() , torch::Tensor() , torch::Tensor() ,
	   torch::Tensor() };
  }
 };

/*--------------------------------------------------------------------------*/
/*------------------------- FACTORY MANAGEMENT -----------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it
{
 SMSpp_insert_in_factory_cpp_0( BundleSolverML );
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

HpNum BundleSolverML::Heuristic( Index whch )
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

 auto input = torch::tensor( features );
 auto output = nn->forward( input );

 double t_pred = output.item< double >();
 if( ( ! std::isfinite( t_pred ) ) || ( t_pred <= 0 ) ) {
  BML_LOG( "Heuristic: invalid t = " << t_pred << ", clamping to 1"
	   << std::endl );
  t_pred = 1;
  }

 phi_vecs.push_back( input );
 BML_LOG( "Heuristic: predicted t = " << t_pred << std::endl );

 // search direction w
 Index dim;
 const Index * nms;
 std::vector< double > tZ( NumVar );
 Master->ReadZ( tZ.data() , nms , dim );
 tZ.resize( dim );
 w_vecs.push_back( torch::tensor( tZ ).requires_grad_( true ) );

 // subgradient matrix G
 int col_num = 0;
 for( Index i = 0 ; i < Master->MaxName() ; ++i )
  if( ItemVcblr[ i ].second < vBPar2[ ItemVcblr[ i ].first ] )
   col_num++;

 std::vector< std::vector< VarValue > > G_mat;
 G_mat.reserve( col_num );
 for( Index i = 0 ; i < Master->MaxName() ; ++i )
  if( ItemVcblr[ i ].second < vBPar2[ ItemVcblr[ i ].first ] ) {
   std::vector< VarValue > G( NumVar );
   v_c05f[ ItemVcblr[ i ].first ]->get_linearization_coefficients(
		    G.data() , Range( 0 , NumVar ) , ItemVcblr[ i ].second );
   if( ! f_convex )
    chgsign( G.data() , NumVar );
   G_mat.push_back( std::move( G ) );
   }

 if( G_mat.empty() ) {
  // no linearization in the bundle yet: nothing to differentiate through,
  // drop the entries recorded so far for this iteration
  phi_vecs.pop_back();
  w_vecs.pop_back();
  return( HpNum( t_pred ) );
  }

 int rows = G_mat.size();
 int cols = G_mat[ 0 ].size();
 std::vector< double > flat;
 flat.reserve( rows * cols );
 for( const auto & row : G_mat )
  flat.insert( flat.end() , row.begin() , row.end() );

 torch::Tensor G_tensor = torch::from_blob( flat.data() , { rows , cols } ,
					    torch::kDouble ).clone();
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
 cHpRow tA = Master->ReadLinErr();
 coeff_vecs.push_back( torch::tensor( SSDone ? 1.0f : -1.0f ) );

 // Gram matrix Q and linearization errors alpha
 torch::Tensor Q = torch::matmul( G_tensor , G_tensor.transpose( 0 , 1 ) );
 Qs.push_back( Q );
 std::vector< double > alpha( Q.sizes()[ 0 ] , 0 );
 for( Index i = 0 ; i < Q.sizes()[ 0 ] ; ++i )
  alpha[ i ] = tA[ i ];
 alphaS.push_back( torch::tensor( alpha , torch::kDouble ) );

 tS.push_back( t_pred );
 FiS.push_back( UpFiBest );

 return( HpNum( t_pred ) );

 }  // end( BundleSolverML::Heuristic )

/*--------------------------------------------------------------------------*/

torch::Tensor BundleSolverML::w( size_t f , double t )
{
 return( BundleSolverML_W::apply( w_vecs[ f ] , Gs[ f ] , Qs[ f ] ,
				  alphaS[ f ] , t ) );
 }

/*--------------------------------------------------------------------------*/
/* The loss is constructed as follows: for each recorded iteration f that
 * either is a SS or is the last one
 *
 *   nn_out = nn->forward( phi_vecs[ f ] ).sum()        [ re-evaluated ]
 *   w_curr = BundleSolverML_W::apply( ... )            [ scalar ]
 *   w_cum += 0.9^( total - ss_count ) * nn_out * w_curr
 *   loss  += coeff * nn_out * BundleSolverML_Fi( w_cum , gs_scalar , fi )
 *
 * and finally loss /= ss_count (mean reduction), so that the magnitude of
 * the loss does not scale with the number of iterations. After
 * loss.backward() the gradients are clipped to GRAD_CLIP_NORM and the Adam
 * update is applied.
 *
 * Three measures guard against exploding gradients:
 *
 * 1. feature normalization [see Heuristic()]: phi_vecs stores the already
 *    normalized tensors, so the forward pass re-run here uses the same
 *    scaled inputs that produced the original prediction;
 *
 * 2. mean (rather than sum) loss reduction, and no fi_val multiplier in
 *    the loss term, since it is a constant w.r.t. the network parameters
 *    whose ~ 1e5 magnitude only inflates the gradients;
 *
 * 3. gradient clipping + Adam: the optimizer is kept as a field so that
 *    its moment estimates persist across successive Backward() calls on
 *    the same (possibly shared) network. */

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

  f_optimizer->zero_grad();

  if( phi_vecs.empty() ) {
   BML_LOG( "Backward: no iterations, skipping" << std::endl );
   return;
   }

  if( ! nn->is_training() )
   nn->train();

  torch::Tensor loss = torch::zeros( {} , torch::kFloat32 );
  size_t last_idx = phi_vecs.size() - 1;
  int ss_count = 0;

  // scalar accumulator for the discounted search directions
  torch::Tensor w_cum = torch::zeros( {} , torch::kFloat32 );

  for( size_t f = 0 ; f < phi_vecs.size() ; ++f ) {
   float coeff_val = coeff_vecs[ f ].item< float >();
   bool is_ss = ( coeff_val > 0 );
   if( ( ! is_ss ) && ( f != last_idx ) )
    continue;
   ss_count++;

   if( ( ! w_vecs[ f ].defined() ) || ( ! Gs[ f ].defined() ) ||
       ( ! Qs[ f ].defined() ) || ( f >= alphaS.size() ) ||
       ( ! alphaS[ f ].defined() ) ) {
    std::cerr << "Backward: missing tensor at iteration " << f << std::endl;
    continue;
    }

   try {
    auto nn_out = nn->forward( phi_vecs[ f ] ).sum();  // scalar, float32

    if( ! nn_out.requires_grad() ) {
     std::cerr << "Backward: no grad_fn on nn_out at iteration " << f
	       << std::endl;
     continue;
     }

    auto w_curr = BundleSolverML_W::apply( w_vecs[ f ] , Gs[ f ] , Qs[ f ] ,
					   alphaS[ f ] ,
					   nn_out.item< double >() );

    double discount = std::pow( 0.9 ,
				double( phi_vecs.size() - ss_count ) );
    w_cum = w_cum + float( discount ) * nn_out * w_curr;

    auto gs_scalar = Gs_aggreg[ f ].to( torch::kFloat32 ).sum();

    auto term = coeff_vecs[ f ] * nn_out *
                BundleSolverML_Fi::apply( w_cum , gs_scalar , FiS[ f ] );
    loss += term;

    BML_LOG( "Backward: iteration " << f << ( is_ss ? " [SS]" : " [last]" )
	     << ", discount = " << discount << ", nn_out = "
	     << nn_out.item< float >() << ", contribution = "
	     << term.item< float >() << std::endl );
    }
   catch( const std::exception & e ) {
    std::cerr << "Backward: error at iteration " << f << ": " << e.what()
	      << std::endl;
    }
   }  // end( for( f ) )

  if( ! ss_count ) {
   BML_LOG( "Backward: no SS / last iteration, skipping" << std::endl );
   return;
   }

  // mean reduction: the loss magnitude must not grow with the iterations
  loss = loss / float( ss_count );

  BML_LOG( "Backward: mean loss = " << loss.item< float >()
	   << ", ss_count = " << ss_count << " / " << phi_vecs.size()
	   << std::endl );

  loss.backward();

  // clip the gradients so that their global L2 norm is <= GRAD_CLIP_NORM
  auto params = nn->parameters();
  auto pre_clip_norm = torch::nn::utils::clip_grad_norm_( params ,
							  GRAD_CLIP_NORM );
  BML_LOG( "Backward: gradient norm (pre-clip) = " << pre_clip_norm
	   << std::endl );

  f_optimizer->step();
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
