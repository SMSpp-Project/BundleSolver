/*--------------------------------------------------------------------------*/
/*----------------------- File test_reattach.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Test of a BundleSolver detached from a Block and attached to another one
 * (or to the same one again), as a LagrangianDualSolver does with its inner
 * Solver whenever it is moved to another Block.
 *
 * The problem is
 * \f[
 *   \min_x \; \frac{1}{2} \sum_j \rho_j x_j^2 + b^\top x +
 *             \max_{i} \{ a_i^\top x + c_i \} \; ,
 * \f]
 * whose 0-th component is either linear (rho = 0, with -1 <= x <= 1 so that
 * the problem is bounded) or quadratic (rho > 0, x free). One BundleSolver
 * is carried across a sequence of such Block, and at each of them its value
 * has to be that of a BundleSolver attached to that Block only: the master
 * problem of the moved Solver must not keep anything of the previous Block,
 * in particular the linear part b of its 0-th component, which a stale
 * value turns into a wrong direction and a wrong bound.
 *
 * The sequence covers the same Block attached again, another b, a zero b
 * after a nonzero one and the other way round, a Block with another number
 * of variables, and the move from a linear 0-th component to a quadratic
 * one and back.
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
#include <string>
#include <vector>

#include "AbstractBlock.h"

#include "BlockSolverConfig.h"

#include "ColVariable.h"

#include "DQuadFunction.h"

#include "FRealObjective.h"

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

/// an instance: n variables, m pieces, the linear part b scaled by bscale
/// (0 gives b = 0) and rho the quadratic coefficient (0 gives the linear
/// 0-th component and the box -1 <= x <= 1)

static AbstractBlock * make_block( Index n , Index m , unsigned seed ,
                                   double bscale , double rho )
{
 std::mt19937 rg( seed );
 std::uniform_real_distribution<> dis( -1.0 , 1.0 );

 MultiVector A( m , RealVector( n ) );
 for( auto & Ai : A )
  for( auto & aij : Ai )
   aij = 10 * dis( rg );
 RealVector c( m );
 for( auto & ci : c )
  ci = 10 * dis( rg );
 std::vector< double > b( n );
 for( auto & bi : b )
  bi = bscale * dis( rg );

 auto root = new AbstractBlock();

 auto x = new std::vector< ColVariable >( n );
 root->add_static_variable( *x , "x" );

 // the box, which the bundle reads from the OneVarConstraint of a Variable
 if( rho == 0 ) {
  auto box = new std::vector< BoxConstraint >( n );
  for( Index i = 0 ; i < n ; ++i ) {
   ( *box )[ i ].set_variable( &(*x)[ i ] );
   ( *box )[ i ].set_lhs( -1 );
   ( *box )[ i ].set_rhs( 1 );
   }
  root->add_static_constraint( *box , "box" );
  }

 // the 0-th component, linear or quadratic
 Function * f0;
 if( rho == 0 ) {
  LinearFunction::v_coeff_pair pairs( n );
  for( Index i = 0 ; i < n ; ++i )
   pairs[ i ] = std::make_pair( &(*x)[ i ] , b[ i ] );
  f0 = new LinearFunction( std::move( pairs ) );
  }
 else {
  DQuadFunction::v_coeff_triple triples( n );
  for( Index i = 0 ; i < n ; ++i )
   triples[ i ] = std::make_tuple( &(*x)[ i ] , b[ i ] , rho / 2.0 );
  f0 = new DQuadFunction( std::move( triples ) );
  }

 auto obj0 = new FRealObjective( root , f0 );
 obj0->set_sense( Objective::eMin , eNoMod );
 root->set_objective( obj0 , eNoMod );

 // the hard component, in a sub-Block of its own
 auto sub = new AbstractBlock( root );
 PolyhedralFunction::VarVector vars( n );
 for( Index i = 0 ; i < n ; ++i )
  vars[ i ] = &(*x)[ i ];

 auto PF = new PolyhedralFunction( std::move( vars ) , std::move( A ) ,
                                   std::move( c ) , - BND );
 auto objs = new FRealObjective( sub , PF );
 objs->set_sense( Objective::eMin , eNoMod );
 sub->set_objective( objs , eNoMod );
 root->add_nested_Block( sub );

 return( root );
 }

/*--------------------------------------------------------------------------*/

/// the value of a BundleSolver attached to a Block of its own

static double solve_fresh( Index n , Index m , unsigned seed , double bscale ,
                           double rho )
{
 auto root = make_block( n , m , seed , bscale , rho );
 auto bsc = apply_bsc( root , "BSPar-rho.txt" );
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

int main( void )
{
 // the Block of the sequence: n , m , seed , scale of b , rho , whether it
 // is the previous Block attached again, and what the step tests
 struct Step {
  Index n , m;
  unsigned seed;
  double bscale , rho;
  bool same;
  std::string what;
  };

 const std::vector< Step > steps = {
  { 5 , 8 , 1 , 1 , 0 , false , "a linear 0-th component" } ,
  { 5 , 8 , 1 , 1 , 0 , true , "the same Block attached again" } ,
  { 5 , 8 , 2 , 1 , 0 , false , "another b" } ,
  { 5 , 8 , 3 , 0 , 0 , false , "b = 0 after a nonzero b" } ,
  { 5 , 8 , 4 , 1 , 0 , false , "a nonzero b after b = 0" } ,
  { 7 , 6 , 5 , 1 , 0 , false , "another number of variables" } ,
  { 5 , 8 , 6 , 1 , 1 , false ,
    "a quadratic 0-th component after a linear one" } ,
  { 5 , 8 , 7 , 1 , 3 , false , "another quadratic 0-th component" } ,
  { 5 , 8 , 8 , 1 , 0 , false ,
    "a linear 0-th component after a quadratic one" }
  };

 // the Solver carried along: that of the first Block, then moved
 AbstractBlock * held = nullptr;
 Solver * slvr = nullptr;
 BlockSolverConfig * bsc = nullptr;

 for( const auto & s : steps ) {
  AbstractBlock * root = ( held && s.same ) ? held
                         : make_block( s.n , s.m , s.seed , s.bscale , s.rho );
  if( ! slvr ) {
   bsc = apply_bsc( root , "BSPar-rho.txt" );
   slvr = root->get_registered_solvers().front();
   }
  else {
   held->unregister_Solver( slvr );
   if( held != root )
    delete held;
   root->register_Solver( slvr );
   }
  held = root;

  slvr->compute( false );
  check_close( slvr->get_var_value() ,
               solve_fresh( s.n , s.m , s.seed , s.bscale , s.rho ) , 1e-6 ,
               s.what );
  }

 held->unregister_Solver( slvr , true );
 delete held;
 delete bsc;

 cout << ( failed ? RED( Shit happened!! ) : GREEN( All tests passed!! ) )
      << endl;

 return( failed ? 1 : 0 );

 }  // end( main )

/*--------------------------------------------------------------------------*/
/*-------------------- End File test_reattach.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
