/*--------------------------------------------------------------------------*/
/*------------------------- File test_integer.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * End-to-end test of the bundle on integer Variable [see intIntVars].
 *
 * The problem is
 * \f[
 *   \min_x \; b^\top x + \sum_{k=1}^K \max_{i \in I_k}
 *   \{ a_i^\top x + c_i \} \; , \qquad
 *   x \in [ -R , R ]^n \, , \; x_j \in \mathbb{Z} \; \forall j \in J \; ,
 * \f]
 * the pieces being split among the K polyhedral components (I_k being the i
 * with i % K == k), each in a sub-Block of its own, with J either all the
 * coordinates or half of them, and with the box of some integer coordinates
 * [ 0 , 1 ] in place of [ -R , R ], so that they are binary and the trust
 * region leaves them to the local branching [see
 * MasterProblemBlock::set_local_branching()]: all of them, half of them, or
 * a third of them with a third general integer and a third continuous. It is
 * solved by every Solver of BSPar-int.txt, i.e., by the BundleSolver, whose
 * master is then a mixed-integer linear problem, and by the
 * ParallelBundleSolver, which evaluates the K components on its threads, and
 * as one mixed-integer program in ( x , v ) with the epigraph constraints
 * written out. The optimal values have to agree, and the
 * ParallelBundleSolver has to take as many iterations as the BundleSolver,
 * since it inserts the linearizations in the bundle in the same order.
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

static void check_all( const std::vector< std::pair< double , long > > & got ,
                       double expected , double tol ,
                       const std::string & what )
{
 for( std::size_t h = 0 ; h < got.size() ; ++h )
  check_close( got[ h ].first , expected , tol ,
               what + ", Solver " + std::to_string( h ) );

 // the Solver after the first evaluate the same components and insert the
 // linearizations in the same order, hence they follow the same path
 for( std::size_t h = 1 ; h < got.size() ; ++h ) {
  const bool ok = ( got[ h ].second == got[ 0 ].second );
  cout << ( ok ? GREEN( ok ) : RED( KO ) ) << "   " << what << ", Solver "
       << h << ": " << got[ h ].second << " iterations against "
       << got[ 0 ].second << endl;
  if( ! ok )
   ++failed;
  }
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

static Index K = 3;        ///< the number of polyhedral components

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
/// save where binary[ j ] is, where the box is [ 0 , 1 ]

static std::vector< ColVariable > * add_x(
                                     AbstractBlock * blk ,
                                     const std::vector< bool > & integer ,
                                     const std::vector< bool > & binary )
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
  const bool bin = ( j < binary.size() ) && binary[ j ];
  (*box)[ j ].set_lhs( bin ? 0 : - R , eNoMod );
  (*box)[ j ].set_rhs( bin ? 1 : R , eNoMod );
  }
 blk->add_static_constraint( *box , "box" );

 return( x );
 }

/*--------------------------------------------------------------------------*/

/// the problem as the bundle sees it: the linear 0-th component in the
/// Objective of the root, each polyhedral one in a sub-Block; returns the
/// value and the iterations of each Solver of the BlockSolverConfig

static std::vector< std::pair< double , long > > solve_with_bundle(
                                      const std::vector< bool > & integer ,
                                      const char * bsc_fn ,
                                      const std::vector< bool > & binary = {} )
{
 const Index n = integer.size();

 auto root = new AbstractBlock();
 auto x = add_x( root , integer , binary );

 LinearFunction::v_coeff_pair cf;
 cf.reserve( n );
 for( Index j = 0 ; j < n ; ++j )
  cf.emplace_back( &(*x)[ j ] , b_lin[ j ] );

 auto obj0 = new FRealObjective( root , new LinearFunction( std::move( cf ) ) );
 obj0->set_sense( Objective::eMin , eNoMod );
 root->set_objective( obj0 , eNoMod );

 // the hard components, each in a sub-Block of its own
 for( Index k = 0 ; k < K ; ++k ) {
  auto sub = new AbstractBlock( root );
  PolyhedralFunction::VarVector vars( n );
  for( Index j = 0 ; j < n ; ++j )
   vars[ j ] = &(*x)[ j ];

  MultiVector Ak;
  RealVector ck;
  for( Index i = k ; i < A.size() ; i += K ) {
   Ak.push_back( A[ i ] );
   ck.push_back( c[ i ] );
   }

  auto PF = new PolyhedralFunction( std::move( vars ) , std::move( Ak ) ,
                                    std::move( ck ) , - BND );
  auto objs = new FRealObjective( sub , PF );
  objs->set_sense( Objective::eMin , eNoMod );
  sub->set_objective( objs , eNoMod );
  root->add_nested_Block( sub );
  }

 auto bsc = apply_bsc( root , bsc_fn );
 std::vector< std::pair< double , long > > res;
 for( auto slvr : root->get_registered_solvers() ) {
  slvr->compute( false );
  res.emplace_back( slvr->get_var_value() ,
                    static_cast< CDASolver * >( slvr )->
                                                get_elapsed_iterations() );
  }

 bsc->clear();
 bsc->apply( root );
 delete bsc;
 delete root;
 return( res );
 }

/*--------------------------------------------------------------------------*/

/// the same problem written out as one mixed-integer program in ( x , v )

static double solve_monolithic( const std::vector< bool > & integer ,
                                const char * bsc_fn ,
                                const std::vector< bool > & binary = {} )
{
 const Index n = integer.size();
 const Index m = A.size();

 auto blk = new AbstractBlock();
 auto x = add_x( blk , integer , binary );

 auto v = new std::vector< ColVariable >( K );
 for( auto & vk : *v ) {
  vk.is_unitary( false , eNoMod );
  vk.is_positive( false , eNoMod );
  }
 blk->add_static_variable( *v , "v" );

 // v_k >= a_i . x + c_i for every i of component k
 auto cons = new std::vector< FRowConstraint >( m );
 for( Index i = 0 ; i < m ; ++i ) {
  LinearFunction::v_coeff_pair cf;
  cf.reserve( n + 1 );
  for( Index j = 0 ; j < n ; ++j )
   cf.emplace_back( &(*x)[ j ] , A[ i ][ j ] );
  cf.emplace_back( &(*v)[ i % K ] , -1.0 );
  (*cons)[ i ].set_function( new LinearFunction( std::move( cf ) ) , eNoMod );
  (*cons)[ i ].set_rhs( - c[ i ] , eNoMod );
  (*cons)[ i ].set_lhs( - Inf< double >() , eNoMod );
  }
 blk->add_static_constraint( *cons , "epi" );

 // b . x + sum_k v_k
 LinearFunction::v_coeff_pair cf;
 cf.reserve( n + K );
 for( Index j = 0 ; j < n ; ++j )
  cf.emplace_back( &(*x)[ j ] , b_lin[ j ] );
 for( auto & vk : *v )
  cf.emplace_back( &vk , 1.0 );

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
 if( argc > 4 )
  K = std::stoi( argv[ 4 ] );

 cout << n << " variables, " << m << " pieces in " << K << " components"
      << endl;

 const std::vector< bool > all( n , true );
 std::vector< bool > half( n , false );
 for( Index j = 0 ; j < n ; j += 2 )
  half[ j ] = true;
 const std::vector< bool > none( n , false );

 // binary coordinates: all, the even ones, and a third of them with a
 // third general integer and a third continuous
 std::vector< bool > third_int( n , true );
 std::vector< bool > third_bin( n , false );
 for( Index j = 0 ; j < n ; ++j ) {
  third_int[ j ] = ( j % 3 != 2 );
  third_bin[ j ] = ( j % 3 == 0 );
  }

 for( unsigned s = seed ; s < seed + 3 ; ++s ) {
  make_instance( n , m , s );
  const auto sd = " (seed " + std::to_string( s ) + ")";

  check_all( solve_with_bundle( all , "BSPar-int.txt" ) ,
               solve_monolithic( all , "MPBCfg-int.txt" ) , 1e-6 ,
               "all integer" + sd );

  check_all( solve_with_bundle( half , "BSPar-int.txt" ) ,
               solve_monolithic( half , "MPBCfg-int.txt" ) , 1e-6 ,
               "half integer" + sd );

  check_all( solve_with_bundle( all , "BSPar-int.txt" , all ) ,
               solve_monolithic( all , "MPBCfg-int.txt" , all ) , 1e-6 ,
               "all binary" + sd );

  check_all( solve_with_bundle( all , "BSPar-int.txt" , half ) ,
               solve_monolithic( all , "MPBCfg-int.txt" , half ) , 1e-6 ,
               "half binary" + sd );

  check_all( solve_with_bundle( third_int , "BSPar-int.txt" , third_bin ) ,
               solve_monolithic( third_int , "MPBCfg-int.txt" , third_bin ) ,
               1e-6 , "a third binary" + sd );

  // with intIntVars == 0 the integrality is ignored
  check_all( solve_with_bundle( all , "BSPar-relax.txt" ) ,
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
