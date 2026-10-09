/*--------------------------------------------------------------------------*/
/*------------------------ File test_osbdo_mcf.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The multicommodity flow instances of
 *
 *   T. Parshakova, F. Zhang, S. Boyd "Implementation of an Oracle-Structured
 *   Bundle Method for Distributed Optimization", Optimization and
 *   Engineering, 2023
 *
 * solved in the same resource-directive form: the master variables are the
 * capacities x_i assigned to each commodity, each agent is the value function
 *
 *   f_i( x_i ) = min { - b_i t : A z + t F_i = 0 , 0 <= z <= x_i , t >= 0 }
 *
 * (a BendersBFunction whose rows z_e <= x_ie take the right-hand side from
 * the master), and the coupling sum_i x_i <= cap is the indicator
 *
 *   g( x ) = max { mu' ( sum_i x_i - cap ) : mu >= 0 }
 *
 * written as a LagBFunction that BundleSolver handles as an "easy"
 * component, i.e., exactly in the master problem. The same instance is also
 * solved as a single LP, which gives the reference value; the test fails if
 * the two differ by more than tol relative.
 *
 * With scale = 1 the master variables are x_ie = cap_e * y_ie with
 * y_ie in [ 0 , 1 ], i.e., the diagonal preconditioning D = diag( u - l ) of
 * the paper done in the model: the mapping of the BendersBFunction and the
 * dual pairs of the LagBFunction get the coefficient cap_e instead of 1.
 *
 * With nseq > 0 the instance is then re-solved nseq more times, each time
 * with the capacities cap_e of the coupling (only them: the agents do not
 * change) drawn uniformly in [ ( 1 - spread ) cap0_e , cap0_e ], cap0 those
 * of the file, with the given seed; the box of the master variables stays
 * the one of cap0, which contains the feasible region of every problem of
 * the sequence. The changes reach BundleSolver as Modification of the
 * LagBFunction, and each problem is checked against its own LP. With the
 * intRstAlg of BSPar one chooses to keep the bundle across the problems
 * (re-optimization) or to empty it at every call (from scratch, starting
 * from the last stability center); with cold = 1 the master variables are
 * also set to 0 before every call, so that with intRstAlg 6 every problem is
 * solved exactly as the first one.
 *
 * Usage: osbdo_mcf_test instance BSPar LPPar LPPar_inner [ scale [ tol
 *                       [ nseq [ spread [ seed [ cold ] ] ] ] ] ]
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/

#include "AbstractBlock.h"
#include "BendersBFunction.h"
#include "BlockSolverConfig.h"
#include "FRealObjective.h"
#include "FRowConstraint.h"
#include "LagBFunction.h"
#include "LinearFunction.h"
#include "OneVarConstraint.h"
#include "Solver.h"

#include <boost/multi_array.hpp>

#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/

struct Instance {
 int V , E , M;
 std::vector< int > tail , head , src , dst;
 std::vector< double > cap , b;
 };

static Instance read_instance( const std::string & fn )
{
 std::ifstream f( fn );
 if( ! f.is_open() )
  throw( std::invalid_argument( "read_instance: cannot open " + fn ) );
 Instance I;
 f >> I.V >> I.E >> I.M;
 I.tail.resize( I.E ); I.head.resize( I.E ); I.cap.resize( I.E );
 for( int e = 0 ; e < I.E ; ++e )
  f >> I.tail[ e ] >> I.head[ e ] >> I.cap[ e ];
 I.src.resize( I.M ); I.dst.resize( I.M ); I.b.resize( I.M );
 for( int i = 0 ; i < I.M ; ++i )
  f >> I.src[ i ] >> I.dst[ i ] >> I.b[ i ];
 if( ! f )
  throw( std::invalid_argument( "read_instance: malformed " + fn ) );
 return( I );
 }

/*--------------------------------------------------------------------------*/
// the rows A z + t F_i = 0 of one commodity, A[ tail ][ e ] = +1 and
// A[ head ][ e ] = -1, F_i[ src ] = +1 and F_i[ dst ] = -1

static std::vector< FRowConstraint > * flow_rows( const Instance & I , int i ,
                                                  ColVariable * z ,
                                                  ColVariable * t )
{
 std::vector< LinearFunction::v_coeff_pair > rows( I.V );
 for( int e = 0 ; e < I.E ; ++e ) {
  rows[ I.tail[ e ] ].emplace_back( & z[ e ] , 1.0 );
  rows[ I.head[ e ] ].emplace_back( & z[ e ] , -1.0 );
  }
 rows[ I.src[ i ] ].emplace_back( t , 1.0 );
 rows[ I.dst[ i ] ].emplace_back( t , -1.0 );

 auto fc = new std::vector< FRowConstraint >( I.V );
 for( int v = 0 ; v < I.V ; ++v ) {
  ( *fc )[ v ].set_function( new LinearFunction( std::move( rows[ v ] ) ) );
  ( *fc )[ v ].set_both( 0 );
  }
 return( fc );
 }

/*--------------------------------------------------------------------------*/

static Solver * attach( Block * block , const std::string & fn )
{
 auto c = Configuration::deserialize( fn );
 auto bsc = dynamic_cast< BlockSolverConfig * >( c );
 if( ! bsc )
  throw( std::invalid_argument( "attach: not a BlockSolverConfig: " + fn ) );
 bsc->apply( block );
 delete( bsc );
 if( block->get_registered_solvers().empty() )
  throw( std::logic_error( "attach: no Solver from " + fn ) );
 return( block->get_registered_solvers().front() );
 }

/*--------------------------------------------------------------------------*/
// the whole problem as one LP: max sum_i b_i t_i, i.e., min - sum_i b_i t_i

static AbstractBlock * build_lp( const Instance & I ,
                                 std::vector< FRowConstraint > * & cc )
{
 auto blk = new AbstractBlock();

 using array2 = boost::multi_array< ColVariable , 2 >;
 auto z = new array2( boost::extents[ I.M ][ I.E ] );
 auto t = new std::vector< ColVariable >( I.M );
 for( auto p = z->data() ; p != z->data() + z->num_elements() ; ++p )
  p->is_positive( true );
 for( auto & ti : *t )
  ti.is_positive( true );
 blk->add_static_variable( *z );
 blk->add_static_variable( *t );

 for( int i = 0 ; i < I.M ; ++i )
  blk->add_static_constraint( *flow_rows( I , i , & ( *z )[ i ][ 0 ] ,
                                          & ( *t )[ i ] ) );

 cc = new std::vector< FRowConstraint >( I.E );
 for( int e = 0 ; e < I.E ; ++e ) {
  LinearFunction::v_coeff_pair cf;
  for( int i = 0 ; i < I.M ; ++i )
   cf.emplace_back( & ( *z )[ i ][ e ] , 1.0 );
  ( *cc )[ e ].set_function( new LinearFunction( std::move( cf ) ) );
  ( *cc )[ e ].set_lhs( -Inf< double >() );
  ( *cc )[ e ].set_rhs( I.cap[ e ] );
  }
 blk->add_static_constraint( *cc );

 LinearFunction::v_coeff_pair of;
 for( int i = 0 ; i < I.M ; ++i )
  of.emplace_back( & ( *t )[ i ] , - I.b[ i ] );
 auto obj = new FRealObjective( blk , new LinearFunction( std::move( of ) ) );
 obj->set_sense( Objective::eMin );
 blk->set_objective( obj );
 return( blk );
 }

/*--------------------------------------------------------------------------*/
// the resource-directive form, see the header

static AbstractBlock * build_decomposed( const Instance & I ,
                                         const std::string & inner_cfg ,
                                         bool scale ,
                                         boost::multi_array< ColVariable , 2 > * & x ,
                                         LinearFunction * & cplobj )
{
 // the coefficient of master variable x_ie wherever it appears
 auto sc = [ & ]( int e ) { return( scale ? I.cap[ e ] : 1.0 ); };

 auto root = new AbstractBlock();

 using array2 = boost::multi_array< ColVariable , 2 >;
 x = new array2( boost::extents[ I.M ][ I.E ] );
 for( auto p = x->data() ; p != x->data() + x->num_elements() ; ++p )
  p->is_positive( true );
 root->add_static_variable( *x );

 auto box = new std::vector< BoxConstraint >( I.M * I.E );
 for( int i = 0 ; i < I.M ; ++i )
  for( int e = 0 ; e < I.E ; ++e ) {
   auto & bc = ( *box )[ i * I.E + e ];
   bc.set_variable( & ( *x )[ i ][ e ] );
   bc.set_lhs( 0 );
   bc.set_rhs( scale ? 1.0 : I.cap[ e ] );
   }
 root->add_static_constraint( *box );

 auto robj = new FRealObjective( root , new LinearFunction() );
 robj->set_sense( Objective::eMin );
 root->set_objective( robj );

 auto & nested = root->access_nested_Blocks();

 // the agents
 for( int i = 0 ; i < I.M ; ++i ) {
  auto inner = new AbstractBlock();
  auto z = new std::vector< ColVariable >( I.E );
  auto t = new std::vector< ColVariable >( 1 );
  for( auto & ze : *z )
   ze.is_positive( true );
  ( *t )[ 0 ].is_positive( true );
  inner->add_static_variable( *z );
  inner->add_static_variable( *t );
  inner->add_static_constraint( *flow_rows( I , i , z->data() ,
                                            t->data() ) );

  auto cr = new std::vector< FRowConstraint >( I.E );
  for( int e = 0 ; e < I.E ; ++e ) {
   ( *cr )[ e ].set_function( new LinearFunction(
                         LinearFunction::v_coeff_pair{ { & ( *z )[ e ] , 1.0 } } ) );
   ( *cr )[ e ].set_lhs( -Inf< double >() );
   ( *cr )[ e ].set_rhs( 0 );
   }
  inner->add_static_constraint( *cr );

  auto iobj = new FRealObjective( inner , new LinearFunction(
                   LinearFunction::v_coeff_pair{ { t->data() , - I.b[ i ] } } ) );
  iobj->set_sense( Objective::eMin );
  inner->set_objective( iobj );
  attach( inner , inner_cfg );

  auto bf = new BendersBFunction();
  bf->set_inner_block( inner );
  BendersBFunction::VarVector xi( I.E );
  for( int e = 0 ; e < I.E ; ++e )
   xi[ e ] = & ( *x )[ i ][ e ];
  bf->set_variables( std::move( xi ) );
  for( int e = 0 ; e < I.E ; ++e ) {
   BendersBFunction::RealVector Ae( I.E , 0 );
   Ae[ e ] = sc( e );
   bf->add_row( std::move( Ae ) , 0 , & ( *cr )[ e ] ,
                BendersBFunction::eRHS );
   }

  auto sub = new AbstractBlock( root );
  auto sobj = new FRealObjective( sub , bf );
  sobj->set_sense( Objective::eMin );
  sub->set_objective( sobj );
  nested.push_back( sub );
  }

 // the coupling, an easy LagBFunction: max { - cap' mu + sum_ie x_ie mu_e }
 {
  auto inner = new AbstractBlock();
  auto mu = new std::vector< ColVariable >( I.E );
  for( auto & m : *mu )
   m.is_positive( true );
  inner->add_static_variable( *mu );
  LinearFunction::v_coeff_pair of;
  for( int e = 0 ; e < I.E ; ++e )
   of.emplace_back( & ( *mu )[ e ] , - I.cap[ e ] );
  cplobj = new LinearFunction( std::move( of ) );
  auto iobj = new FRealObjective( inner , cplobj );
  iobj->set_sense( Objective::eMax );   // convex LagBFunction
  inner->set_objective( iobj );

  auto lbf = new LagBFunction( inner );
  LagBFunction::v_dual_pair dp;
  dp.reserve( I.M * I.E );
  for( int i = 0 ; i < I.M ; ++i )
   for( int e = 0 ; e < I.E ; ++e )
    dp.emplace_back( & ( *x )[ i ][ e ] , new LinearFunction(
                        LinearFunction::v_coeff_pair{ { & ( *mu )[ e ] , sc( e ) } } ) );
  lbf->set_dual_pairs( std::move( dp ) );

  auto sub = new AbstractBlock( root );
  auto sobj = new FRealObjective( sub , lbf );
  sobj->set_sense( Objective::eMin );
  sub->set_objective( sobj );
  nested.push_back( sub );
  }

 return( root );
 }

/*--------------------------------------------------------------------------*/

int main( int argc , char ** argv )
{
 if( argc < 5 ) {
  std::cerr << "Usage: " << argv[ 0 ]
            << " instance BSPar LPPar LPPar_inner [ scale [ tol [ nseq"
            << " [ spread [ seed [ cold ] ] ] ] ] ]" << std::endl;
  return( 1 );
  }

 const auto I = read_instance( argv[ 1 ] );
 const bool scale = ( argc > 5 ) && ( std::atoi( argv[ 5 ] ) != 0 );
 const double tol = ( argc > 6 ) ? std::atof( argv[ 6 ] ) : 1e-6;
 const int nseq = ( argc > 7 ) ? std::atoi( argv[ 7 ] ) : 0;
 const double spread = ( argc > 8 ) ? std::atof( argv[ 8 ] ) : 0.2;
 const int seed = ( argc > 9 ) ? std::atoi( argv[ 9 ] ) : 0;
 const bool cold = ( argc > 10 ) && ( std::atoi( argv[ 10 ] ) != 0 );

 std::cout << std::setprecision( 12 );
 std::cout << "instance " << argv[ 1 ] << ": V = " << I.V << ", E = " << I.E
           << ", M = " << I.M << std::endl;

 using clk = std::chrono::steady_clock;
 auto since = []( clk::time_point t0 ) {
  return( std::chrono::duration< double >( clk::now() - t0 ).count() );
  };

 std::vector< FRowConstraint > * cc;
 auto lp = build_lp( I , cc );
 auto slp = attach( lp , argv[ 3 ] );

 boost::multi_array< ColVariable , 2 > * x;
 LinearFunction * cplobj;
 auto t0 = clk::now();
 auto blk = build_decomposed( I , argv[ 4 ] , scale , x , cplobj );
 double tb = since( t0 );
 auto s = attach( blk , argv[ 2 ] );
 s->set_log( & std::cout );

 std::mt19937 rng( seed );
 std::uniform_real_distribution< double > U( 1 - spread , 1 );

 bool ok = true;
 double sum_t = 0;
 long sum_it = 0;
 for( int k = 0 ; k <= nseq ; ++k ) {
  if( k > 0 ) {  // new capacities of the coupling, in the LP and the bundle
   Function::Vec_FunctionValue nc( I.E );
   for( int e = 0 ; e < I.E ; ++e ) {
    const double c = I.cap[ e ] * U( rng );
    ( *cc )[ e ].set_rhs( c );
    nc[ e ] = - c;
    }
   cplobj->modify_coefficients( std::move( nc ) , Block::Range( 0 , I.E ) );
   if( cold )
    for( auto p = x->data() ; p != x->data() + x->num_elements() ; ++p )
     p->set_value( 0 );
   }

  t0 = clk::now();
  auto stlp = slp->compute();
  double ref = slp->get_var_value();
  double tlp = since( t0 );

  t0 = clk::now();
  auto st = s->compute();
  double tt = since( t0 );
  double v = s->get_var_value();
  long it = s->get_elapsed_iterations();
  const double gap = std::abs( v - ref ) / std::max( 1.0 , std::abs( ref ) );
  const bool good = ( stlp == Solver::kOK ) && ( st == Solver::kOK ) &&
                    ( gap <= tol );
  ok = ok && good;
  if( k > 0 ) {
   sum_t += tt;
   sum_it += it;
   }
  std::cout << "problem " << k << ": LP status " << stlp << ", value " << ref
            << ", time " << tlp << " s | bundle"
            << ( scale ? " (scaled)" : "" ) << ": status " << st
            << ", value " << v << ", iterations " << it << ", time " << tt
            << " s" << ( k == 0 ? " (build " + std::to_string( tb ) + " s)"
                                : std::string() )
            << ", rel. gap " << gap << ( good ? "" : " KO" ) << std::endl;
  }

 if( nseq > 0 )
  std::cout << "sequence  : " << nseq << " problems after the first, "
            << sum_it << " iterations, " << sum_t << " s" << std::endl;

 lp->unregister_Solvers();
 delete( lp );

 std::cout << ( ok ? "OK" : "KO" ) << std::endl;
 return( ok ? 0 : 1 );
 }

/*--------------------------------------------------------------------------*/
/*---------------------- End File test_osbdo_mcf.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
