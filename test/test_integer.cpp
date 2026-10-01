/*--------------------------------------------------------------------------*/
/*------------------------- File test_integer.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * End-to-end test of the bundle on integer Variable [see intIntVars].
 *
 * The problem is
 * \f[
 *   \min_x \; b^\top x + \max_{i} \{ a_i^\top x + c_i \} \; , \qquad
 *   x \in [ -R , R ]^n \, , \; x_j \in \mathbb{Z} \; \forall j \in J \; ,
 * \f]
 * with J either all the coordinates or half of them. It is solved twice: by
 * the bundle, whose master is then a mixed-integer linear problem
 * [see BSPar-int.txt], and as one mixed-integer program in ( x , v ) with
 * the epigraph constraints written out. The two optimal values have to
 * agree.
 *
 * The same instance is then solved by the bundle with intIntVars == 0
 * [see BSPar-relax.txt], which has to ignore the integrality and give the
 * value of the continuous relaxation, i.e., what the bundle gave before the
 * integer Variable were dealt with.
 *
 * Both the master and the mixed-integer program are solved by HiGHS, so
 * that the test needs no licence.
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

#include "FRealObjective.h"

#include "FRowConstraint.h"

#include "LinearFunction.h"

#include "OneVarConstraint.h"

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

const double R = 3;        ///< the box is [ -R , R ]^n

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

/// the Variable x, integer where integer[ j ] is, in the box [ -R , R ]^n

static std::vector< ColVariable > * add_x( AbstractBlock * blk ,
                                           const std::vector< bool > & integer )
{
 const Index n = integer.size();

 auto x = new std::vector< ColVariable >( n );
 for( Index j = 0 ; j < n ; ++j ) {
  (*x)[ j ].is_positive( false , eNoMod );
  (*x)[ j ].is_unitary( false , eNoMod );
  (*x)[ j ].is_integer( integer[ j ] , eNoMod );
  }
 blk->add_static_variable( *x , "x" );

 auto box = new std::vector< BoxConstraint >( n );
 for( Index j = 0 ; j < n ; ++j ) {
  (*box)[ j ].set_variable( &(*x)[ j ] , eNoMod );
  (*box)[ j ].set_lhs( - R , eNoMod );
  (*box)[ j ].set_rhs( R , eNoMod );
  }
 blk->add_static_constraint( *box , "box" );

 return( x );
 }

/*--------------------------------------------------------------------------*/

/// the problem as the bundle sees it: the linear 0-th component in the
/// Objective of the root, the polyhedral one in a sub-Block

static double solve_with_bundle( const std::vector< bool > & integer ,
                                 const char * bsc_fn )
{
 const Index n = integer.size();

 auto root = new AbstractBlock();
 auto x = add_x( root , integer );

 LinearFunction::v_coeff_pair cf;
 cf.reserve( n );
 for( Index j = 0 ; j < n ; ++j )
  cf.emplace_back( &(*x)[ j ] , b_lin[ j ] );

 auto obj0 = new FRealObjective( root , new LinearFunction( std::move( cf ) ) );
 obj0->set_sense( Objective::eMin , eNoMod );
 root->set_objective( obj0 , eNoMod );

 // the hard component, in a sub-Block of its own
 auto sub = new AbstractBlock( root );
 PolyhedralFunction::VarVector vars( n );
 for( Index j = 0 ; j < n ; ++j )
  vars[ j ] = &(*x)[ j ];

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

/// the same problem written out as one mixed-integer program in ( x , v )

static double solve_monolithic( const std::vector< bool > & integer ,
                                const char * bsc_fn )
{
 const Index n = integer.size();
 const Index m = A.size();

 auto blk = new AbstractBlock();
 auto x = add_x( blk , integer );

 auto v = new ColVariable();
 v->is_unitary( false , eNoMod );
 v->is_positive( false , eNoMod );
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

 // b . x + v
 LinearFunction::v_coeff_pair cf;
 cf.reserve( n + 1 );
 for( Index j = 0 ; j < n ; ++j )
  cf.emplace_back( &(*x)[ j ] , b_lin[ j ] );
 cf.emplace_back( v , 1.0 );

 auto obj = new FRealObjective( blk , new LinearFunction( std::move( cf ) ) );
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
 const Index n = ( argc > 1 ) ? std::stoi( argv[ 1 ] ) : 6;
 const Index m = ( argc > 2 ) ? std::stoi( argv[ 2 ] ) : 12;
 const unsigned seed = ( argc > 3 ) ? std::stoi( argv[ 3 ] ) : 1;

 cout << n << " variables, " << m << " pieces" << endl;

 const std::vector< bool > all( n , true );
 std::vector< bool > half( n , false );
 for( Index j = 0 ; j < n ; j += 2 )
  half[ j ] = true;
 const std::vector< bool > none( n , false );

 for( unsigned s = seed ; s < seed + 3 ; ++s ) {
  make_instance( n , m , s );
  const auto sd = " (seed " + std::to_string( s ) + ")";

  check_close( solve_with_bundle( all , "BSPar-int.txt" ) ,
               solve_monolithic( all , "MPBCfg-int.txt" ) , 1e-6 ,
               "all integer" + sd );

  check_close( solve_with_bundle( half , "BSPar-int.txt" ) ,
               solve_monolithic( half , "MPBCfg-int.txt" ) , 1e-6 ,
               "half integer" + sd );

  // with intIntVars == 0 the integrality is ignored
  check_close( solve_with_bundle( all , "BSPar-relax.txt" ) ,
               solve_monolithic( none , "MPBCfg-int.txt" ) , 1e-6 ,
               "relaxation" + sd );
  }

 cout << ( failed ? RED( Shit happened!! ) : GREEN( All tests passed!! ) )
      << endl;

 return( failed ? 1 : 0 );

 }  // end( main )

/*--------------------------------------------------------------------------*/
/*-------------------------- End File test_integer.cpp ---------------------*/
/*--------------------------------------------------------------------------*/
