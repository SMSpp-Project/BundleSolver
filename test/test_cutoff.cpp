/*--------------------------------------------------------------------------*/
/*------------------------ File test_cutoff.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Test of the cutoffs of BundleSolver [see dblUpCutOff and dblLwCutOff in
 * BundleSolver.h], in particular of the one that a point reaches being the
 * global lower bound of the dual master, which is only conditionally valid.
 *
 * The problem is
 * \f[
 *   \min_x \{ \, b^\top x + \max_{i} \{ a_i^\top x + c_i \} \, : \,
 *            -1 \leq x \leq 1 \, \} \; ,
 * \f]
 * and its concave mirror, the maximum of the opposite function, whose value
 * is the opposite one. Each instance is solved by each BundleSolver of
 * BSPar-cutoff.txt, one for each stabilization of the dual master, first with
 * no cutoff, which gives the optimal value v*, and then with the cutoff that
 * a point reaches (dblLwCutOff for the minimization, dblUpCutOff for the
 * maximization) put
 *
 * - beyond v* by a gap, so that it is never reached: the result has to be
 *   the same as with no cutoff, kOK and v*;
 *
 * - on the other side of v* by a gap, so that a point reaches it: the
 *   result has to be kCutOff, with a point at least as good as asked;
 *
 * - on the other side of v* by less than the accuracy, so that only a
 *   point within the accuracy of v* reaches it, which the master with the
 *   cutoff as its global lower bound finds optimal at the bound: the
 *   result has to be kCutOff, with a point as good as asked up to the
 *   accuracy.
 *
 * In all cases the bound the solver reports (get_lb() for the minimization,
 * get_ub() for the maximization) has to be valid, i.e., never beyond v*,
 * whatever the cutoff is. The number of iterations of each run is printed.
 *
 *       ./cutoff_test [seed n m]
 *       seed: pseudo-random generator seed of the first instance [0]
 *       n:    number of variables [10]
 *       m:    number of pieces of the max [40]
 *
 * Three instances are solved, from seed on.
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
#include <iomanip>
#include <iostream>
#include <iterator>
#include <random>
#include <string>
#include <vector>

#include "AbstractBlock.h"

#include "BlockSolverConfig.h"

#include "ColVariable.h"

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

const double BND = 1e+6;     ///< bound of the polyhedral component
const double ACC = 1e-6;     ///< dblRelAcc of BSCfg-cut.txt
const double GAP = 1e-2;     ///< relative gap of the cutoffs off v*

static const char * const STBL[] = { "proximal" , "level" , "doubly" };

static int failed = 0;

/*--------------------------------------------------------------------------*/

static void check( bool ok , const std::string & what )
{
 cout << ( ok ? GREEN( ok ) : RED( KO ) ) << "   " << what << endl;
 if( ! ok )
  ++failed;
 }

/*--------------------------------------------------------------------------*/

/// an instance: n variables, m pieces, the minimization of the convex
/// function or the maximization of the opposite, concave, one

static AbstractBlock * make_block( Index n , Index m , unsigned seed ,
                                   bool concave )
{
 std::mt19937 rg( seed );
 std::uniform_real_distribution<> dis( -1.0 , 1.0 );

 const double sgn = concave ? -1 : 1;
 MultiVector A( m , RealVector( n ) );
 for( auto & Ai : A )
  for( auto & aij : Ai )
   aij = sgn * 10 * dis( rg );
 RealVector c( m );
 for( auto & ci : c )
  ci = sgn * 10 * dis( rg );
 std::vector< double > b( n );
 for( auto & bi : b )
  bi = sgn * dis( rg );

 auto root = new AbstractBlock();

 auto x = new std::vector< ColVariable >( n );
 root->add_static_variable( *x , "x" );

 auto box = new std::vector< BoxConstraint >( n );
 for( Index i = 0 ; i < n ; ++i ) {
  ( *box )[ i ].set_variable( &(*x)[ i ] );
  ( *box )[ i ].set_lhs( -1 );
  ( *box )[ i ].set_rhs( 1 );
  }
 root->add_static_constraint( *box , "box" );

 LinearFunction::v_coeff_pair pairs( n );
 for( Index i = 0 ; i < n ; ++i )
  pairs[ i ] = std::make_pair( &(*x)[ i ] , b[ i ] );
 auto obj0 = new FRealObjective( root ,
                                 new LinearFunction( std::move( pairs ) ) );
 obj0->set_sense( concave ? Objective::eMax : Objective::eMin , eNoMod );
 root->set_objective( obj0 , eNoMod );

 auto sub = new AbstractBlock( root );
 PolyhedralFunction::VarVector vars( n );
 for( Index i = 0 ; i < n ; ++i )
  vars[ i ] = &(*x)[ i ];

 auto PF = new PolyhedralFunction( std::move( vars ) , std::move( A ) ,
                                   std::move( c ) , concave ? BND : - BND ,
                                   ! concave );
 auto objs = new FRealObjective( sub , PF );
 objs->set_sense( concave ? Objective::eMax : Objective::eMin , eNoMod );
 sub->set_objective( objs , eNoMod );
 root->add_nested_Block( sub );

 return( root );
 }

/*--------------------------------------------------------------------------*/

/// what a run gives

struct Run {
 int status;     ///< the return value of compute()
 double value;   ///< the value of the solution
 double bound;   ///< the bound: get_lb() if min, get_ub() if max
 long iters;     ///< the number of iterations
 };

/*--------------------------------------------------------------------------*/

/// a run of the j-th BundleSolver of BSPar-cutoff.txt on a fresh instance,
/// with the cutoff that a point reaches set to cut

static Run solve( Index n , Index m , unsigned seed , bool concave ,
                  Index j , double cut )
{
 auto root = make_block( n , m , seed , concave );
 auto bsc = dynamic_cast< BlockSolverConfig * >(
                          Configuration::deserialize( "BSPar-cutoff.txt" ) );
 if( ! bsc ) {
  cerr << "Error: BSPar-cutoff.txt does not hold a BlockSolverConfig"
       << endl;
  exit( 1 );
  }
 bsc->apply( root );

 auto slvr = dynamic_cast< BundleSolver * >(
             *std::next( root->get_registered_solvers().begin() , j ) );
 slvr->set_par( concave ? Solver::dblUpCutOff : Solver::dblLwCutOff , cut );

 Run r;
 r.status = slvr->compute( false );
 r.value = slvr->get_var_value();
 r.bound = concave ? slvr->get_ub() : slvr->get_lb();
 r.iters = slvr->get_elapsed_iterations();

 bsc->clear();
 bsc->apply( root );
 delete bsc;
 delete root;
 return( r );
 }

/*--------------------------------------------------------------------------*/

static std::string status_name( int s )
{
 switch( s ) {
  case( Solver::kOK ):           return( "kOK" );
  case( Solver::kCutOff ):       return( "kCutOff" );
  case( Solver::kLowPrecision ): return( "kLowPrecision" );
  case( Solver::kStopIter ):     return( "kStopIter" );
  default:                       return( std::to_string( s ) );
  }
 }

/*--------------------------------------------------------------------------*/

int main( int argc , char ** argv )
{
 unsigned seed = 0;
 Index n = 10;
 Index m = 40;
 if( argc > 1 )
  seed = std::stoul( argv[ 1 ] );
 if( argc > 2 )
  n = std::stoul( argv[ 2 ] );
 if( argc > 3 )
  m = std::stoul( argv[ 3 ] );

 cout << std::setprecision( 10 );

 for( unsigned s = seed ; s < seed + 3 ; ++s )
  for( bool concave : { false , true } )
   for( Index j = 0 ; j < 3 ; ++j ) {
    // in the minimization the solver does, the cutoff is a lower bound:
    // "beyond" v* is below it, the "other side" above it
    const double sgn = concave ? -1 : 1;
    const std::string inst = "seed " + std::to_string( s ) +
                             ( concave ? " max " : " min " ) + STBL[ j ];

    const auto r0 = solve( n , m , s , concave , j ,
                           concave ? Inf< double >() : - Inf< double >() );
    const double v = r0.value;
    const double tol = 10 * ACC * std::max( 1.0 , std::abs( v ) );
    const double gap = GAP * std::max( 1.0 , std::abs( v ) );
    cout << inst << ": v* = " << v << ", " << status_name( r0.status )
         << ", " << r0.iters << " iterations" << endl;
    check( r0.status == Solver::kOK , inst + ": kOK with no cutoff" );

    // a bound is valid if it is not beyond v* (by more than the accuracy)
    auto valid = [ & ]( const Run & r ) {
     return( sgn * ( r.bound - v ) <= tol );
     };

    // never reached: the same as with no cutoff
    const auto r1 = solve( n , m , s , concave , j , v - sgn * gap );
    cout << "  cutoff never reached: " << status_name( r1.status )
         << ", value " << r1.value << ", bound " << r1.bound << ", "
         << r1.iters << " iterations" << endl;
    check( r1.status == Solver::kOK , inst + ": kOK, cutoff never reached" );
    check( std::abs( r1.value - v ) <= tol ,
           inst + ": v*, cutoff never reached" );
    check( valid( r1 ) , inst + ": valid bound, cutoff never reached" );

    // reached by a point at a gap from v*
    const double cut = v + sgn * gap;
    const auto r2 = solve( n , m , s , concave , j , cut );
    cout << "  cutoff reached: " << status_name( r2.status ) << ", value "
         << r2.value << ", bound " << r2.bound << ", " << r2.iters
         << " iterations" << endl;
    check( r2.status == Solver::kCutOff , inst + ": kCutOff, reached" );
    check( sgn * ( r2.value - cut ) <= tol ,
           inst + ": as good as asked, reached" );
    check( valid( r2 ) , inst + ": valid bound, reached" );

    // reached only within the accuracy of v*
    const double cuta = v + sgn * ACC * std::max( 1.0 , std::abs( v ) ) / 10;
    const auto r3 = solve( n , m , s , concave , j , cuta );
    cout << "  cutoff at the optimum: " << status_name( r3.status )
         << ", value " << r3.value << ", bound " << r3.bound << ", "
         << r3.iters << " iterations" << endl;
    check( ( r3.status == Solver::kCutOff ) || ( r3.status == Solver::kOK ) ,
           inst + ": kCutOff or kOK, at the optimum" );
    check( sgn * ( r3.value - cuta ) <= tol ,
           inst + ": as good as asked up to the accuracy, at the optimum" );
    check( valid( r3 ) , inst + ": valid bound, at the optimum" );
    }

 if( failed )
  cout << RED( KO ) << "   " << failed << " checks failed" << endl;
 else
  cout << GREEN( ok ) << "   all checks passed" << endl;

 return( failed ? 1 : 0 );
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- End File test_cutoff.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
