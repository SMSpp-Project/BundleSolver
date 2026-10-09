/*--------------------------------------------------------------------------*/
/*------------------------- File test_osbdo_sc.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The supply chain instances of
 *
 *   T. Parshakova, F. Zhang, S. Boyd "Implementation of an Oracle-Structured
 *   Bundle Method for Distributed Optimization", Optimization and
 *   Engineering, 2023
 *
 * solved in the same form: a chain of N agents, agent i with m_i input and
 * n_i output nodes and public variable x_i = ( inputs , outputs ) in
 * [ 0 , u_i ], whose function is the value of its transport QP
 *
 *   f_i( x_i ) = min lin_i . X + 0.5 quad_i . X^2 + 0.5 mu || r ||_1
 *                s.t. ( colsum X , rowsum X ) = w , sum w_in = sum w_out ,
 *                     0 <= X <= cap_i , 0 <= w <= u_i , w - r = x_i ,
 *
 * a BendersBFunction whose rows w - r+ + r- = x_i take both sides from the
 * master (r = r+ - r-, so the penalty is linear). The coupling is the cost
 * sale' x_0,in - retail' x_N-1,out, the linear part of the objective, plus
 * the equalities x_i,out = x_i+1,in and sum x_i,in = sum x_i,out, i.e., the
 * indicator
 *
 *   g( x ) = max { mu' E x : mu free }
 *
 * written as a LagBFunction that BundleSolver handles as an "easy"
 * component, i.e., exactly in the master problem, as the paper does with its
 * structured function. The same instance is also solved as one QP, which
 * gives the reference value; the test fails if the two differ by more than
 * tol, relative to max( 1 , |ref| ).
 *
 * Usage: osbdo_sc_test instance BSPar QPPar QPPar_inner [ tol ]
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
#include "DQuadFunction.h"
#include "FRealObjective.h"
#include "FRowConstraint.h"
#include "LagBFunction.h"
#include "LinearFunction.h"
#include "OneVarConstraint.h"
#include "Solver.h"

#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/

struct Agent {
 int m , n;
 double mu;
 std::vector< double > ub;                           // m + n
 std::vector< std::vector< double > > lin , quad , cap;  // n x m
 };

struct Instance {
 std::vector< double > sale , retail;
 std::vector< Agent > ag;
 };

static Instance read_instance( const std::string & fn )
{
 std::ifstream f( fn );
 if( ! f.is_open() )
  throw( std::invalid_argument( "read_instance: cannot open " + fn ) );
 Instance I;
 int N;
 f >> N;
 I.ag.resize( N );
 // the sizes of sale and retail are those of the first and last agent,
 // which come later: read the two lines as they are
 std::string line;
 std::getline( f , line );
 auto read_line = [ & ]( std::vector< double > & v ) {
  std::getline( f , line );
  std::istringstream is( line );
  double a;
  while( is >> a )
   v.push_back( a );
  };
 read_line( I.sale );
 read_line( I.retail );
 for( auto & a : I.ag ) {
  f >> a.m >> a.n >> a.mu;
  a.ub.resize( a.m + a.n );
  for( auto & u : a.ub )
   f >> u;
  for( auto p : { & a.lin , & a.quad , & a.cap } ) {
   p->assign( a.n , std::vector< double >( a.m ) );
   for( auto & row : *p )
    for( auto & e : row )
     f >> e;
   }
  }
 if( ( ! f ) || ( int( I.sale.size() ) != I.ag.front().m ) ||
     ( int( I.retail.size() ) != I.ag.back().n ) )
  throw( std::invalid_argument( "read_instance: malformed " + fn ) );
 for( int i = 0 ; i + 1 < N ; ++i )
  if( I.ag[ i ].n != I.ag[ i + 1 ].m )
   throw( std::invalid_argument( "read_instance: the chain does not match" ) );
 return( I );
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

static std::vector< BoxConstraint > * boxes( std::vector< ColVariable > & v ,
                                            const std::vector< double > & u )
{
 auto b = new std::vector< BoxConstraint >( v.size() );
 for( std::size_t j = 0 ; j < v.size() ; ++j ) {
  ( *b )[ j ].set_variable( & v[ j ] );
  ( *b )[ j ].set_lhs( 0 );
  ( *b )[ j ].set_rhs( u[ j ] );
  }
 return( b );
 }

/*--------------------------------------------------------------------------*/
// the transport problem of agent a in blk: returns the rows
// w_j - r+_j + r-_j - x_j = 0 (x == nullptr: = 0, x comes from the master)
// and adds the terms of the objective to of

static std::vector< FRowConstraint > * agent( AbstractBlock * blk ,
                                              const Agent & a ,
                                              ColVariable * x ,
                                              DQuadFunction::v_coeff_triple & of )
{
 const int d = a.m + a.n;
 auto X = new std::vector< ColVariable >( a.n * a.m );
 auto w = new std::vector< ColVariable >( d );
 auto rp = new std::vector< ColVariable >( d );
 auto rm = new std::vector< ColVariable >( d );
 for( auto v : { X , w , rp , rm } ) {
  for( auto & c : *v )
   c.is_positive( true );
  blk->add_static_variable( *v );
  }

 std::vector< double > capv( a.n * a.m );
 for( int k = 0 ; k < a.n ; ++k )
  for( int j = 0 ; j < a.m ; ++j ) {
   capv[ k * a.m + j ] = a.cap[ k ][ j ];
   of.emplace_back( & ( *X )[ k * a.m + j ] , a.lin[ k ][ j ] ,
                    0.5 * a.quad[ k ][ j ] );
   }
 blk->add_static_constraint( *boxes( *X , capv ) );
 blk->add_static_constraint( *boxes( *w , a.ub ) );
 for( int j = 0 ; j < d ; ++j ) {
  of.emplace_back( & ( *rp )[ j ] , 0.5 * a.mu , 0.0 );
  of.emplace_back( & ( *rm )[ j ] , 0.5 * a.mu , 0.0 );
  }

 // w = ( colsum X , rowsum X ) and sum w_in = sum w_out
 auto bal = new std::vector< FRowConstraint >( d + 1 );
 for( int j = 0 ; j < d ; ++j ) {
  LinearFunction::v_coeff_pair cf{ { & ( *w )[ j ] , 1.0 } };
  if( j < a.m )
   for( int k = 0 ; k < a.n ; ++k )
    cf.emplace_back( & ( *X )[ k * a.m + j ] , -1.0 );
  else
   for( int h = 0 ; h < a.m ; ++h )
    cf.emplace_back( & ( *X )[ ( j - a.m ) * a.m + h ] , -1.0 );
  ( *bal )[ j ].set_function( new LinearFunction( std::move( cf ) ) );
  ( *bal )[ j ].set_both( 0 );
  }
 {
  LinearFunction::v_coeff_pair cf;
  for( int j = 0 ; j < d ; ++j )
   cf.emplace_back( & ( *w )[ j ] , j < a.m ? 1.0 : -1.0 );
  ( *bal )[ d ].set_function( new LinearFunction( std::move( cf ) ) );
  ( *bal )[ d ].set_both( 0 );
  }
 blk->add_static_constraint( *bal );

 auto lr = new std::vector< FRowConstraint >( d );
 for( int j = 0 ; j < d ; ++j ) {
  LinearFunction::v_coeff_pair cf{ { & ( *w )[ j ] , 1.0 } ,
                                   { & ( *rp )[ j ] , -1.0 } ,
                                   { & ( *rm )[ j ] , 1.0 } };
  if( x )
   cf.emplace_back( & x[ j ] , -1.0 );
  ( *lr )[ j ].set_function( new LinearFunction( std::move( cf ) ) );
  ( *lr )[ j ].set_both( 0 );
  }
 blk->add_static_constraint( *lr );
 return( lr );
 }

/*--------------------------------------------------------------------------*/
// the rows of the coupling: for each of them the pairs ( agent , position
// in x_agent , coefficient )

using CRow = std::vector< std::tuple< int , int , double > >;

static std::vector< CRow > coupling_rows( const Instance & I )
{
 std::vector< CRow > rows;
 const int N = I.ag.size();
 for( int i = 0 ; i < N ; ++i ) {
  const auto & a = I.ag[ i ];
  if( i + 1 < N )  // x_i,out = x_i+1,in
   for( int k = 0 ; k < a.n ; ++k )
    rows.push_back( { { i , a.m + k , 1.0 } , { i + 1 , k , -1.0 } } );
  CRow cons;       // sum x_i,in = sum x_i,out
  for( int j = 0 ; j < a.m + a.n ; ++j )
   cons.emplace_back( i , j , j < a.m ? 1.0 : -1.0 );
  rows.push_back( std::move( cons ) );
  }
 return( rows );
 }

/*--------------------------------------------------------------------------*/
// the linear cost: sale' x_0,in - retail' x_N-1,out

static LinearFunction::v_coeff_pair cost( const Instance & I ,
                                          std::vector< std::vector< ColVariable > * > & x )
{
 LinearFunction::v_coeff_pair of;
 for( std::size_t j = 0 ; j < I.sale.size() ; ++j )
  of.emplace_back( & ( *x.front() )[ j ] , I.sale[ j ] );
 const auto & l = I.ag.back();
 for( std::size_t k = 0 ; k < I.retail.size() ; ++k )
  of.emplace_back( & ( *x.back() )[ l.m + k ] , - I.retail[ k ] );
 return( of );
 }

/*--------------------------------------------------------------------------*/
// the whole problem as one QP

static AbstractBlock * build_qp( const Instance & I )
{
 auto blk = new AbstractBlock();
 std::vector< std::vector< ColVariable > * > x;
 for( auto & a : I.ag ) {
  x.push_back( new std::vector< ColVariable >( a.m + a.n ) );
  blk->add_static_variable( *x.back() );
  blk->add_static_constraint( *boxes( *x.back() , a.ub ) );
  }

 DQuadFunction::v_coeff_triple of;
 for( std::size_t i = 0 ; i < I.ag.size() ; ++i )
  agent( blk , I.ag[ i ] , x[ i ]->data() , of );
 for( auto & [ v , c ] : cost( I , x ) )
  of.emplace_back( v , c , 0.0 );

 const auto cr = coupling_rows( I );
 auto cc = new std::vector< FRowConstraint >( cr.size() );
 for( std::size_t r = 0 ; r < cr.size() ; ++r ) {
  LinearFunction::v_coeff_pair cf;
  for( auto & [ i , j , c ] : cr[ r ] )
   cf.emplace_back( & ( *x[ i ] )[ j ] , c );
  ( *cc )[ r ].set_function( new LinearFunction( std::move( cf ) ) );
  ( *cc )[ r ].set_both( 0 );
  }
 blk->add_static_constraint( *cc );

 auto obj = new FRealObjective( blk , new DQuadFunction( std::move( of ) ) );
 obj->set_sense( Objective::eMin );
 blk->set_objective( obj );
 return( blk );
 }

/*--------------------------------------------------------------------------*/
// the decomposed form, see the header

static AbstractBlock * build_decomposed( const Instance & I ,
                                         const std::string & inner_cfg )
{
 auto root = new AbstractBlock();
 std::vector< std::vector< ColVariable > * > x;
 for( auto & a : I.ag ) {
  x.push_back( new std::vector< ColVariable >( a.m + a.n ) );
  root->add_static_variable( *x.back() );
  root->add_static_constraint( *boxes( *x.back() , a.ub ) );
  }
 // the linear part has every master variable, most of them with 0 cost,
 // since BundleSolver wants it to cover all those of the components
 LinearFunction::v_coeff_pair rc;
 {
  std::map< ColVariable * , double > c;
  for( auto & [ v , cv ] : cost( I , x ) )
   c[ v ] = cv;
  for( auto xi : x )
   for( auto & v : *xi )
    rc.emplace_back( & v , c.count( & v ) ? c[ & v ] : 0.0 );
  }
 auto robj = new FRealObjective( root , new LinearFunction( std::move( rc ) ) );
 robj->set_sense( Objective::eMin );
 root->set_objective( robj );

 auto & nested = root->access_nested_Blocks();

 // the agents
 for( std::size_t i = 0 ; i < I.ag.size() ; ++i ) {
  const auto & a = I.ag[ i ];
  const int d = a.m + a.n;
  auto inner = new AbstractBlock();
  DQuadFunction::v_coeff_triple of;
  auto lr = agent( inner , a , nullptr , of );
  auto iobj = new FRealObjective( inner , new DQuadFunction( std::move( of ) ) );
  iobj->set_sense( Objective::eMin );
  inner->set_objective( iobj );
  attach( inner , inner_cfg );

  auto bf = new BendersBFunction();
  bf->set_inner_block( inner );
  BendersBFunction::VarVector xi( d );
  for( int j = 0 ; j < d ; ++j )
   xi[ j ] = & ( *x[ i ] )[ j ];
  bf->set_variables( std::move( xi ) );
  for( int j = 0 ; j < d ; ++j ) {
   BendersBFunction::RealVector Aj( d , 0 );
   Aj[ j ] = 1;
   bf->add_row( std::move( Aj ) , 0 , & ( *lr )[ j ] ,
                BendersBFunction::eBoth );
   }

  auto sub = new AbstractBlock( root );
  auto sobj = new FRealObjective( sub , bf );
  sobj->set_sense( Objective::eMin );
  sub->set_objective( sobj );
  nested.push_back( sub );
  }

 // the coupling, an easy LagBFunction: max { sum_r mu_r E_r x : mu free }
 {
  const auto cr = coupling_rows( I );
  auto inner = new AbstractBlock();
  auto mu = new std::vector< ColVariable >( cr.size() );
  inner->add_static_variable( *mu );
  auto iobj = new FRealObjective( inner , new LinearFunction() );
  iobj->set_sense( Objective::eMax );   // convex LagBFunction
  inner->set_objective( iobj );

  // the coefficients of each x_ij in the rows
  std::vector< std::vector< LinearFunction::v_coeff_pair > > col( I.ag.size() );
  for( std::size_t i = 0 ; i < I.ag.size() ; ++i )
   col[ i ].resize( I.ag[ i ].m + I.ag[ i ].n );
  for( std::size_t r = 0 ; r < cr.size() ; ++r )
   for( auto & [ i , j , c ] : cr[ r ] )
    col[ i ][ j ].emplace_back( & ( *mu )[ r ] , c );

  auto lbf = new LagBFunction( inner );
  LagBFunction::v_dual_pair dp;
  for( std::size_t i = 0 ; i < I.ag.size() ; ++i )
   for( std::size_t j = 0 ; j < col[ i ].size() ; ++j )
    dp.emplace_back( & ( *x[ i ] )[ j ] ,
                     new LinearFunction( std::move( col[ i ][ j ] ) ) );
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
            << " instance BSPar QPPar QPPar_inner [ tol ]" << std::endl;
  return( 1 );
  }

 const auto I = read_instance( argv[ 1 ] );
 const double tol = ( argc > 5 ) ? std::atof( argv[ 5 ] ) : 1e-6;
 std::cout << std::setprecision( 12 );
 std::cout << "instance " << argv[ 1 ] << ": N = " << I.ag.size()
           << std::endl;

 using clk = std::chrono::steady_clock;
 auto since = []( clk::time_point t0 ) {
  return( std::chrono::duration< double >( clk::now() - t0 ).count() );
  };

 double ref = 0;
 {
  auto t0 = clk::now();
  auto qp = build_qp( I );
  auto s = attach( qp , argv[ 3 ] );
  auto st = s->compute();
  ref = s->get_var_value();
  std::cout << "QP        : status " << st << ", value " << ref
            << ", time " << since( t0 ) << " s" << std::endl;
  qp->unregister_Solvers();
  delete( qp );
  }

 auto t0 = clk::now();
 auto blk = build_decomposed( I , argv[ 4 ] );
 double tb = since( t0 );
 auto s = attach( blk , argv[ 2 ] );
 s->set_log( & std::cout );
 t0 = clk::now();
 auto st = s->compute();
 double tt = since( t0 );
 double v = s->get_var_value();
 const double gap = std::abs( v - ref ) / std::max( 1.0 , std::abs( ref ) );
 std::cout << "bundle    : status " << st << ", value " << v
           << ", iterations " << s->get_elapsed_iterations() << ", time "
           << tt << " s (build " << tb << " s)" << std::endl;
 std::cout << "gap       : " << gap << std::endl;
 const bool ok = ( st == Solver::kOK ) && ( gap <= tol );
 std::cout << ( ok ? "OK" : "KO" ) << std::endl;
 return( ok ? 0 : 1 );
 }

/*--------------------------------------------------------------------------*/
/*----------------------- End File test_osbdo_sc.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
