/*--------------------------------------------------------------------------*/
/*------------------- File test_linear_constraints.cpp ---------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Test of the linear Constraint of the Block of BundleSolver, i.e., the
 * FRowConstraint with a LinearFunction of its ColVariable, which the master
 * problem has as rows (the primal one) or as multipliers (the dual one), and
 * of their Modification between two calls of compute().
 *
 * The problem is
 * \f[
 *   \min_x \{ \, b^\top x + \max_{i} \{ a_i^\top x + c_i \} \, : \,
 *            l \leq R x \leq u \, , \, -1 \leq x \leq 1 \, \} \; ,
 * \f]
 * the rows of R being a dynamic group of FRowConstraint of the Block. Each
 * BundleSolver of BSPar-lc.txt is attached to the same instance, and all of
 * them solve it, and solve it again after each of the following changes:
 *
 * - none: the rows are a minimum on the sum of x, which the starting point
 *   x = 0 violates, and a ramp | x_0 - x_1 | <= 0.2;
 *
 * - the minimum on the sum of x is raised (a change of the left-hand side);
 *
 * - a coefficient of the ramp changes (a change of its LinearFunction);
 *
 * - a row x_2 + x_3 <= -0.5 is added;
 *
 * - the minimum on the sum of x is relaxed;
 *
 * - the ramp is removed.
 *
 * Each time the value has to be the one of the same problem written as a
 * Linear Program, with the epigraph variable of the max, which a
 * :MILPSolver of MPBCfg.txt solves on a Block of its own. The number of
 * iterations of each run is printed.
 *
 *       ./linear_constraints_test [seed n m]
 *       seed: pseudo-random generator seed of the first instance [0]
 *       n:    number of variables [10]
 *       m:    number of pieces of the max [40]
 *
 * Two instances are solved, from seed on.
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
#include <list>
#include <random>
#include <string>
#include <vector>

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

const double BND = 1e+6;     ///< bound of the polyhedral component
const double ACC = 1e-6;     ///< dblRelAcc of BSCfg-cut.txt

static const char * const NAME[] = {
 "proximal, dual" , "proximal, primal" , "trust region, dual" ,
 "trust region, primal, iterate" , "level, dual, iterate" };

static const Index NSLV = 5;  ///< the BundleSolver of BSPar-lc.txt

static int failed = 0;

/*--------------------------------------------------------------------------*/

static void check( bool ok , const std::string & what )
{
 cout << ( ok ? GREEN( ok ) : RED( KO ) ) << "   " << what << endl;
 if( ! ok )
  ++failed;
 }

/*--------------------------------------------------------------------------*/

/// the data of an instance

struct Data {
 MultiVector A;               ///< the pieces of the max
 RealVector c;                ///< their constants
 std::vector< double > b;     ///< the linear part
 };

static Data make_data( Index n , Index m , unsigned seed )
{
 std::mt19937 rg( seed );
 std::uniform_real_distribution<> dis( -1.0 , 1.0 );

 Data d;
 d.A.assign( m , RealVector( n ) );
 for( auto & Ai : d.A )
  for( auto & aij : Ai )
   aij = 10 * dis( rg );
 d.c.resize( m );
 for( auto & ci : d.c )
  ci = 10 * dis( rg );
 d.b.resize( n );
 for( auto & bi : d.b )
  bi = dis( rg );
 return( d );
 }

/*--------------------------------------------------------------------------*/

/// the instance that the BundleSolver solve: x, its box and the rows in the
/// root, whose objective is b x, and the max in its only sub-Block

struct Instance {
 AbstractBlock * root;
 std::vector< ColVariable > * x;
 std::list< FRowConstraint > * rows;
 };

/// a new row l <= sum_j t_j x_j <= u at the end of rows (a Constraint
/// cannot be copied, hence it is built where it stays)

static void add_row( std::list< FRowConstraint > & rows ,
                     std::vector< ColVariable > & x ,
                     const std::vector< std::pair< Index , double > > & t ,
                     double l , double u )
{
 LinearFunction::v_coeff_pair pairs;
 for( const auto & [ j , a ] : t )
  pairs.emplace_back( & x[ j ] , a );
 auto & row = rows.emplace_back();
 row.set_function( new LinearFunction( std::move( pairs ) ) , eNoMod );
 row.set_lhs( l , eNoMod );
 row.set_rhs( u , eNoMod );
 }

static Instance make_instance( const Data & d , Index n )
{
 Instance I;
 I.root = new AbstractBlock();

 I.x = new std::vector< ColVariable >( n );
 I.root->add_static_variable( *I.x , "x" );

 auto box = new std::vector< BoxConstraint >( n );
 for( Index i = 0 ; i < n ; ++i ) {
  ( *box )[ i ].set_variable( &(*I.x)[ i ] );
  ( *box )[ i ].set_lhs( -1 );
  ( *box )[ i ].set_rhs( 1 );
  }
 I.root->add_static_constraint( *box , "box" );

 // the minimum on the sum, which x = 0 violates, and the ramp
 std::vector< std::pair< Index , double > > all( n );
 for( Index i = 0 ; i < n ; ++i )
  all[ i ] = { i , 1.0 };
 I.rows = new std::list< FRowConstraint >;
 add_row( *I.rows , *I.x , all , 0.2 * n , Inf< double >() );
 add_row( *I.rows , *I.x , { { 0 , 1.0 } , { 1 , -1.0 } } , -0.2 , 0.2 );
 I.root->add_dynamic_constraint( *I.rows , "rows" );

 LinearFunction::v_coeff_pair pairs( n );
 for( Index i = 0 ; i < n ; ++i )
  pairs[ i ] = std::make_pair( &(*I.x)[ i ] , d.b[ i ] );
 auto obj0 = new FRealObjective( I.root ,
                                 new LinearFunction( std::move( pairs ) ) );
 obj0->set_sense( Objective::eMin , eNoMod );
 I.root->set_objective( obj0 , eNoMod );

 auto sub = new AbstractBlock( I.root );
 PolyhedralFunction::VarVector vars( n );
 for( Index i = 0 ; i < n ; ++i )
  vars[ i ] = &(*I.x)[ i ];
 auto A = d.A;
 auto c = d.c;
 auto PF = new PolyhedralFunction( std::move( vars ) , std::move( A ) ,
                                   std::move( c ) , - BND , true );
 auto objs = new FRealObjective( sub , PF );
 objs->set_sense( Objective::eMin , eNoMod );
 sub->set_objective( objs , eNoMod );
 I.root->add_nested_Block( sub );

 return( I );
 }

/*--------------------------------------------------------------------------*/

/// the optimal value of the same problem as a Linear Program, with the rows
/// that are in I now: min b x + v : v >= a_i x + c_i, the rows, the box

static double reference( const Data & d , Index n , const Instance & I )
{
 auto lp = new AbstractBlock();
 auto x = new std::vector< ColVariable >( n + 1 );  // x and v
 lp->add_static_variable( *x , "xv" );

 auto box = new std::vector< BoxConstraint >( n );
 for( Index i = 0 ; i < n ; ++i ) {
  ( *box )[ i ].set_variable( &(*x)[ i ] );
  ( *box )[ i ].set_lhs( -1 );
  ( *box )[ i ].set_rhs( 1 );
  }
 lp->add_static_constraint( *box , "box" );

 auto cuts = new std::vector< FRowConstraint >( d.A.size() );
 for( Index i = 0 ; i < d.A.size() ; ++i ) {
  LinearFunction::v_coeff_pair pairs;
  pairs.emplace_back( &(*x)[ n ] , 1.0 );
  for( Index j = 0 ; j < n ; ++j )
   pairs.emplace_back( &(*x)[ j ] , - d.A[ i ][ j ] );
  ( *cuts )[ i ].set_function( new LinearFunction( std::move( pairs ) ) ,
                               eNoMod );
  ( *cuts )[ i ].set_lhs( d.c[ i ] , eNoMod );
  ( *cuts )[ i ].set_rhs( Inf< double >() , eNoMod );
  }
 lp->add_static_constraint( *cuts , "cuts" );

 // the rows of I that are not relaxed, on the same coordinates
 auto rows = new std::list< FRowConstraint >;
 for( auto & r : *I.rows ) {
  if( r.is_relaxed() )
   continue;
  auto lf = static_cast< LinearFunction * >( r.get_function() );
  std::vector< std::pair< Index , double > > t;
  for( Index h = 0 ; h < lf->get_num_active_var() ; ++h ) {
   const auto v = static_cast< ColVariable * >( lf->get_active_var( h ) );
   t.emplace_back( Index( v - I.x->data() ) , lf->get_coefficient( h ) );
   }
  add_row( *rows , *x , t , r.get_lhs() , r.get_rhs() );
  }
 lp->add_dynamic_constraint( *rows , "rows" );

 LinearFunction::v_coeff_pair pairs( n + 1 );
 for( Index i = 0 ; i < n ; ++i )
  pairs[ i ] = std::make_pair( &(*x)[ i ] , d.b[ i ] );
 pairs[ n ] = std::make_pair( &(*x)[ n ] , 1.0 );
 auto obj = new FRealObjective( lp ,
                                new LinearFunction( std::move( pairs ) ) );
 obj->set_sense( Objective::eMin , eNoMod );
 lp->set_objective( obj , eNoMod );

 auto bsc = dynamic_cast< BlockSolverConfig * >(
                                Configuration::deserialize( "MPBCfg.txt" ) );
 if( ! bsc ) {
  cerr << "Error: MPBCfg.txt does not hold a BlockSolverConfig" << endl;
  exit( 1 );
  }
 bsc->apply( lp );
 auto slvr = lp->get_registered_solvers().front();
 const int st = slvr->compute( false );
 const double v = ( st == Solver::kOK ) ? slvr->get_var_value()
                                        : Inf< double >();
 bsc->clear();
 bsc->apply( lp );
 delete bsc;
 delete lp;
 return( v );
 }

/*--------------------------------------------------------------------------*/

/// all the BundleSolver of I solve it, and each has to find the value of
/// the Linear Program

static void solve_all( const Data & d , Index n , const Instance & I ,
                       const std::string & step )
{
 const double ref = reference( d , n , I );
 const double tol = 10 * ACC * std::max( 1.0 , std::abs( ref ) );
 cout << step << ": LP = " << ref << endl;

 auto it = I.root->get_registered_solvers().begin();
 for( Index j = 0 ; j < NSLV ; ++j , ++it ) {
  auto slvr = dynamic_cast< BundleSolver * >( *it );
  const int st = slvr->compute( false );
  const double v = slvr->get_var_value();
  cout << "  " << NAME[ j ] << ": " << v << ", status " << st << ", "
       << slvr->get_elapsed_iterations() << " iterations" << endl;
  check( ( st == Solver::kOK ) && ( std::abs( v - ref ) <= tol ) ,
         step + ", " + NAME[ j ] );
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

 for( unsigned s = seed ; s < seed + 2 ; ++s ) {
  const auto d = make_data( n , m , s );
  auto I = make_instance( d , n );
  const std::string inst = "seed " + std::to_string( s );

  auto bsc = dynamic_cast< BlockSolverConfig * >(
                              Configuration::deserialize( "BSPar-lc.txt" ) );
  if( ! bsc ) {
   cerr << "Error: BSPar-lc.txt does not hold a BlockSolverConfig" << endl;
   exit( 1 );
   }
  bsc->apply( I.root );

  solve_all( d , n , I , inst + ", the first rows" );

  auto sum = I.rows->begin();
  auto ramp = std::next( sum );

  sum->set_lhs( 0.4 * n );
  solve_all( d , n , I , inst + ", the minimum raised" );

  auto ramp_lf = static_cast< LinearFunction * >( ramp->get_function() );
  ramp_lf->modify_coefficient( 1 , -2 );
  solve_all( d , n , I , inst + ", a coefficient of the ramp changed" );

  std::list< FRowConstraint > added;
  add_row( added , *I.x , { { 2 , 1.0 } , { 3 , 1.0 } } , - Inf< double >() ,
           -0.5 );
  I.root->add_dynamic_constraints( *I.rows , added );
  solve_all( d , n , I , inst + ", a row added" );

  sum->relax( true );
  solve_all( d , n , I , inst + ", the minimum relaxed" );

  I.root->remove_dynamic_constraint( *I.rows , ramp );
  solve_all( d , n , I , inst + ", the ramp removed" );

  bsc->clear();
  bsc->apply( I.root );
  delete bsc;
  delete I.root;
  }

 if( failed )
  cout << RED( KO ) << "   " << failed << " checks failed" << endl;
 else
  cout << GREEN( ok ) << "   all checks passed" << endl;

 return( failed ? 1 : 0 );
 }

/*--------------------------------------------------------------------------*/
/*-------------------- End File test_linear_constraints.cpp ----------------*/
/*--------------------------------------------------------------------------*/
