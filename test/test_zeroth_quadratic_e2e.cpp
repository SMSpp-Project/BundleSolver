/*--------------------------------------------------------------------------*/
/*----------------- File test_zeroth_quadratic_e2e.cpp ---------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * End-to-end test of the quadratic "0-th" component of the bundle
 * [see MasterProblemBlock::set_zeroth_quadratic()].
 *
 * The problem is
 * \f[
 *   \min_x \; \frac{1}{2} \sum_j \rho_j x_j^2 + b^\top x +
 *             \max_{i} \{ a_i^\top x + c_i \} \; ,
 * \f]
 * i.e. a regularised risk, which is what the Benders master of a support
 * vector machine and any other strongly convex nonsmooth problem look like.
 * It is solved twice: by the bundle, which sees the quadratic term as the
 * 0-th component of the sum-function and hands it to its master, and as one
 * monolithic quadratic program in ( x , v ) with the epigraph constraints
 * written out, which any general-purpose solver handles. The two optimal
 * values have to agree.
 *
 * With rho = 0 the same instance is solved again, which is the case that has
 * to keep behaving exactly as it did before the quadratic component existed.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cmath>
#include <iostream>
#include <random>

#include "AbstractBlock.h"

#include "BlockSolverConfig.h"

#include "ColVariable.h"

#include "DQuadFunction.h"

#include "FRealObjective.h"

#include "FRowConstraint.h"

#include "LinearFunction.h"

#include "PolyhedralFunction.h"

#include "BundleSolver.h"

/*--------------------------------------------------------------------------*/

using namespace std;
using namespace SMSpp_di_unipi_it;

using Index = Block::Index;
using MultiVector = PolyhedralFunction::MultiVector;
using RealVector = PolyhedralFunction::RealVector;

#define RED( x )   "\x1B[31m" #x "\033[0m"
#define GREEN( x ) "\x1B[32m" #x "\033[0m"

const double BND = 1e+6;   ///< global lower bound of the polyhedral component

static int failed = 0;

/*--------------------------------------------------------------------------*/

static void check_close( double got , double expected , double tol ,
                         const std::string & what )
{
 const bool ok = std::abs( got - expected ) <=
                 tol * std::max( 1.0 , std::abs( expected ) );
 cout << ( ok ? GREEN( ok ) : RED( KO ) ) << "   " << what
      << ": " << got << " against " << expected << endl;
 if( ! ok )
  ++failed;
 }

/*--------------------------------------------------------------------------*/

static BlockSolverConfig * apply_bsc( Block * block , const std::string & fn )
{
 auto * bsc = dynamic_cast< BlockSolverConfig * >(
                                     Configuration::deserialize( fn ) );
 if( ! bsc ) {
  cerr << "Error: " << fn << " does not hold a BlockSolverConfig" << endl;
  exit( 1 );
  }
 bsc->apply( block );
 return( bsc );
 }

/*--------------------------------------------------------------------------*/

/// the instance: rows a_i and constants c_i of the polyhedral component,
/// plus the linear part b of the 0-th one

static MultiVector A;
static RealVector c;
static std::vector< double > b_lin;

static void make_instance( Index n , Index m , unsigned seed )
{
 std::mt19937 rg( seed );
 std::uniform_real_distribution<> dis( -1.0 , 1.0 );

 A.resize( m );
 for( auto & Ai : A ) {
  Ai.resize( n );
  for( auto & aij : Ai )
   aij = 10 * dis( rg );
  }

 c.resize( m );
 for( auto & ci : c )
  ci = 10 * dis( rg );

 b_lin.resize( n );
 for( auto & bi : b_lin )
  bi = dis( rg );
 }

/*--------------------------------------------------------------------------*/

/// the problem as the bundle sees it: the quadratic 0-th component in the
/// Objective of the root, the polyhedral one in a sub-Block

static double solve_with_bundle( Index n , const std::vector< double > & rho ,
                                 const char * bsc_fn )
{
 auto root = new AbstractBlock();

 auto x = new std::vector< ColVariable >( n );
 for( auto & xi : *x )
  xi.is_unitary( false , eNoMod );
 root->add_static_variable( *x , "x" );

 // the 0-th component: ( 1 / 2 ) sum_j rho_j x_j^2 + b . x, i.e. a
 // DQuadFunction whose quadratic coefficient is rho_j / 2
 DQuadFunction::v_coeff_triple triples( n );
 for( Index i = 0 ; i < n ; ++i )
  triples[ i ] = std::make_tuple( &(*x)[ i ] , b_lin[ i ] ,
                                  rho[ i ] / 2.0 );

 auto obj0 = new FRealObjective( root , new DQuadFunction(
                                              std::move( triples ) ) );
 obj0->set_sense( Objective::eMin , eNoMod );
 root->set_objective( obj0 , eNoMod );

 // the hard component, in a sub-Block of its own
 auto sub = new AbstractBlock( root );
 PolyhedralFunction::VarVector vars( n );
 for( Index i = 0 ; i < n ; ++i )
  vars[ i ] = &(*x)[ i ];

 auto PF = new PolyhedralFunction( std::move( vars ) , MultiVector( A ) ,
                                   RealVector( c ) , - BND );
 auto objs = new FRealObjective( sub , PF );
 objs->set_sense( Objective::eMin , eNoMod );
 sub->set_objective( objs , eNoMod );
 root->add_nested_Block( sub );

 auto bsc = apply_bsc( root , bsc_fn );
 auto slvr = root->get_registered_solvers().front();
 slvr->compute( false );
 const double value = slvr->get_var_value();

 bsc->clear();
 bsc->apply( root );
 delete bsc;
 delete root;
 return( value );
 }

/*--------------------------------------------------------------------------*/

/// the same problem written out as one quadratic program in ( x , v )

static double solve_monolithic( Index n , const std::vector< double > & rho ,
                                const char * bsc_fn )
{
 const Index m = A.size();

 auto blk = new AbstractBlock();

 auto x = new std::vector< ColVariable >( n );
 auto v = new ColVariable();
 v->is_unitary( false , eNoMod );
 v->is_positive( false , eNoMod );
 blk->add_static_variable( *x , "x" );
 blk->add_static_variable( *v , "v" );

 // v >= a_i . x + c_i for every i
 auto cons = new std::vector< FRowConstraint >( m );
 for( Index i = 0 ; i < m ; ++i ) {
  LinearFunction::v_coeff_pair cf;
  cf.reserve( n + 1 );
  for( Index j = 0 ; j < n ; ++j )
   cf.emplace_back( &(*x)[ j ] , A[ i ][ j ] );
  cf.emplace_back( v , -1.0 );
  (*cons)[ i ].set_function( new LinearFunction( std::move( cf ) ) , eNoMod );
  (*cons)[ i ].set_rhs( - c[ i ] , eNoMod );
  (*cons)[ i ].set_lhs( - Inf< double >() , eNoMod );
  }
 blk->add_static_constraint( *cons , "epi" );

 // ( 1 / 2 ) sum_j rho_j x_j^2 + b . x + v
 DQuadFunction::v_coeff_triple triples;
 triples.reserve( n + 1 );
 for( Index j = 0 ; j < n ; ++j )
  triples.emplace_back( &(*x)[ j ] , b_lin[ j ] , rho[ j ] / 2.0 );
 triples.emplace_back( v , 1.0 , 0.0 );

 auto obj = new FRealObjective( blk , new DQuadFunction(
                                            std::move( triples ) ) );
 obj->set_sense( Objective::eMin , eNoMod );
 blk->set_objective( obj , eNoMod );

 auto bsc = apply_bsc( blk , bsc_fn );
 auto slvr = blk->get_registered_solvers().front();
 slvr->compute( false );
 const double value = slvr->get_var_value();

 bsc->clear();
 bsc->apply( blk );
 delete bsc;
 delete blk;
 return( value );
 }

/*--------------------------------------------------------------------------*/

int main( int argc , char ** argv )
{
 const Index n = ( argc > 1 ) ? std::stoi( argv[ 1 ] ) : 5;
 const Index m = ( argc > 2 ) ? std::stoi( argv[ 2 ] ) : 8;
 const unsigned seed = ( argc > 3 ) ? std::stoi( argv[ 3 ] ) : 1;

 make_instance( n , m , seed );

 cout << n << " variables, " << m << " pieces" << endl;

 /* rho = 0 is the instance without the quadratic component, which here is
  * unbounded below: what is checked there is that both say so, the value
  * comparison being meaningless. */

 for( double rho : { 1.0 , 3.0 , 10.0 } ) {
  const std::vector< double > rhos( n , rho );
  const double bundle = solve_with_bundle( n , rhos , "BSPar-rho.txt" );
  const double mono = solve_monolithic( n , rhos , "MPBCfg-rho.txt" );
  check_close( bundle , mono , 1e-6 ,
               "rho = " + std::to_string( rho ) );
  }

 // the coefficients need not be all equal: a Benders master where only some
 // of the variables are regularised is the reason why

 { std::vector< double > rhos( n );
  for( Index j = 0 ; j < n ; ++j )
   rhos[ j ] = 1.0 + j;
  const double bundle = solve_with_bundle( n , rhos , "BSPar-rho.txt" );
  const double mono = solve_monolithic( n , rhos , "MPBCfg-rho.txt" );
  check_close( bundle , mono , 1e-6 , "a rho per coordinate" );
  }

 cout << ( failed ? RED( Shit happened!! ) : GREEN( All tests passed!! ) )
      << endl;

 return( failed ? 1 : 0 );

 }  // end( main )

/*--------------------------------------------------------------------------*/
/*--------------- End File test_zeroth_quadratic_e2e.cpp -------------------*/
/*--------------------------------------------------------------------------*/
