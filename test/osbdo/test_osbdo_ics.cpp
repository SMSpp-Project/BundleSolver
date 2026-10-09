/*--------------------------------------------------------------------------*/
/*------------------------ File test_osbdo_ics.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The intersection-of-convex-sets instances of
 *
 *   T. Parshakova, F. Zhang, S. Boyd "Implementation of an Oracle-Structured
 *   Bundle Method for Distributed Optimization", Optimization and
 *   Engineering, 2023
 *
 * i.e., find a point in the intersection of M polyhedra { y : A_i y <= b_i }
 * by minimizing the sum of the squared distances from them,
 *
 *   min sum_i f_i( x ) ,  f_i( x ) = min { || x - y ||^2 : A_i y <= b_i } ,
 *
 * whose optimal value is 0 when the intersection is not empty, as it is for
 * the generated instances. The paper gives each agent its own copy x_i with
 * the coupling x_i = x_0; here all the agents share the same master
 * variables x, which is the same problem without the coupling. Each f_i is a
 * BendersBFunction whose inner problem is the QP
 *
 *   min { || d ||^2 : d + y = x , A_i y <= b_i }
 *
 * with the rows d + y = x taking both sides from the master. The same
 * instance is also solved as one QP, which gives the reference value; the
 * test fails if the two differ by more than tol, relative to max( 1 , |ref| ).
 *
 * Usage: osbdo_ics_test instance BSPar QPPar QPPar_inner [ tol ]
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
#include "LinearFunction.h"
#include "Solver.h"

#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/

struct Instance {
 int R , C , M;
 std::vector< std::vector< double > > A;
 std::vector< double > b;
 };

static Instance read_instance( const std::string & fn )
{
 std::ifstream f( fn );
 if( ! f.is_open() )
  throw( std::invalid_argument( "read_instance: cannot open " + fn ) );
 Instance I;
 f >> I.R >> I.C >> I.M;
 I.A.assign( I.R , std::vector< double >( I.C ) );
 I.b.resize( I.R );
 for( int r = 0 ; r < I.R ; ++r ) {
  for( auto & a : I.A[ r ] )
   f >> a;
  f >> I.b[ r ];
  }
 if( ( ! f ) || ( I.R % I.M ) )
  throw( std::invalid_argument( "read_instance: malformed " + fn ) );
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
// the rows A_i y <= b_i of agent i

static std::vector< FRowConstraint > * set_rows( const Instance & I , int i ,
                                                 ColVariable * y )
{
 const int nr = I.R / I.M;
 auto rows = new std::vector< FRowConstraint >( nr );
 for( int h = 0 ; h < nr ; ++h ) {
  const int r = i * nr + h;
  LinearFunction::v_coeff_pair cf;
  cf.reserve( I.C );
  for( int j = 0 ; j < I.C ; ++j )
   cf.emplace_back( & y[ j ] , I.A[ r ][ j ] );
  ( *rows )[ h ].set_function( new LinearFunction( std::move( cf ) ) );
  ( *rows )[ h ].set_lhs( -Inf< double >() );
  ( *rows )[ h ].set_rhs( I.b[ r ] );
  }
 return( rows );
 }

/*--------------------------------------------------------------------------*/
// the rows d + y - x = 0 (x == nullptr: d + y = 0, x comes from the master)

static std::vector< FRowConstraint > * link_rows( const Instance & I ,
                                                  ColVariable * d ,
                                                  ColVariable * y ,
                                                  ColVariable * x )
{
 auto rows = new std::vector< FRowConstraint >( I.C );
 for( int j = 0 ; j < I.C ; ++j ) {
  LinearFunction::v_coeff_pair cf{ { & d[ j ] , 1.0 } , { & y[ j ] , 1.0 } };
  if( x )
   cf.emplace_back( & x[ j ] , -1.0 );
  ( *rows )[ j ].set_function( new LinearFunction( std::move( cf ) ) );
  ( *rows )[ j ].set_both( 0 );
  }
 return( rows );
 }

/*--------------------------------------------------------------------------*/
// the whole problem as one QP

static AbstractBlock * build_qp( const Instance & I )
{
 auto blk = new AbstractBlock();
 auto x = new std::vector< ColVariable >( I.C );
 blk->add_static_variable( *x );
 DQuadFunction::v_coeff_triple of;
 for( int i = 0 ; i < I.M ; ++i ) {
  auto y = new std::vector< ColVariable >( I.C );
  auto d = new std::vector< ColVariable >( I.C );
  blk->add_static_variable( *y );
  blk->add_static_variable( *d );
  blk->add_static_constraint( *set_rows( I , i , y->data() ) );
  blk->add_static_constraint( *link_rows( I , d->data() , y->data() ,
                                          x->data() ) );
  for( auto & dj : *d )
   of.emplace_back( & dj , 0.0 , 1.0 );
  }
 auto obj = new FRealObjective( blk , new DQuadFunction( std::move( of ) ) );
 obj->set_sense( Objective::eMin );
 blk->set_objective( obj );
 return( blk );
 }

/*--------------------------------------------------------------------------*/
// one BendersBFunction per agent, all on the same master variables

static AbstractBlock * build_decomposed( const Instance & I ,
                                         const std::string & inner_cfg )
{
 auto root = new AbstractBlock();
 auto x = new std::vector< ColVariable >( I.C );
 root->add_static_variable( *x );
 auto robj = new FRealObjective( root , new LinearFunction() );
 robj->set_sense( Objective::eMin );
 root->set_objective( robj );

 auto & nested = root->access_nested_Blocks();
 for( int i = 0 ; i < I.M ; ++i ) {
  auto inner = new AbstractBlock();
  auto y = new std::vector< ColVariable >( I.C );
  auto d = new std::vector< ColVariable >( I.C );
  inner->add_static_variable( *y );
  inner->add_static_variable( *d );
  inner->add_static_constraint( *set_rows( I , i , y->data() ) );
  auto lr = link_rows( I , d->data() , y->data() , nullptr );
  inner->add_static_constraint( *lr );
  DQuadFunction::v_coeff_triple of;
  for( auto & dj : *d )
   of.emplace_back( & dj , 0.0 , 1.0 );
  auto iobj = new FRealObjective( inner , new DQuadFunction( std::move( of ) ) );
  iobj->set_sense( Objective::eMin );
  inner->set_objective( iobj );
  attach( inner , inner_cfg );

  auto bf = new BendersBFunction();
  bf->set_inner_block( inner );
  BendersBFunction::VarVector xi( I.C );
  for( int j = 0 ; j < I.C ; ++j )
   xi[ j ] = & ( *x )[ j ];
  bf->set_variables( std::move( xi ) );
  for( int j = 0 ; j < I.C ; ++j ) {
   BendersBFunction::RealVector Aj( I.C , 0 );
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
 std::cout << "instance " << argv[ 1 ] << ": R = " << I.R << ", C = " << I.C
           << ", M = " << I.M << std::endl;

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
/*---------------------- End File test_osbdo_ics.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
