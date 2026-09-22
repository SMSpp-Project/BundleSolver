/*--------------------------------------------------------------------------*/
/*------------------------- File test_osbdo_fl.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The federated learning instances of
 *
 *   T. Parshakova, F. Zhang, S. Boyd "Implementation of an Oracle-Structured
 *   Bundle Method for Distributed Optimization", Optimization and
 *   Engineering, 2023
 *
 * i.e., L1-regularized logistic regression with the samples split among M
 * agents,
 *
 *   min sum_i f_i( x ) + lambda || x ||_1 ,
 *   f_i( x ) = sum_{s in S_i} log( 1 + exp( - y_s a_s x ) ) .
 *
 * The paper gives each agent its own copy x_i with the coupling x_i = x_0;
 * here all the agents share the same master variables x, which is the same
 * problem without the coupling. Each f_i is a LogisticFunction, a
 * C05Function whose oracle is the closed form of the value and the gradient,
 * with a global pool; the regularization is
 *
 *   lambda || x ||_1 = max { lambda w' x : -1 <= w <= 1 } ,
 *
 * a LagBFunction that BundleSolver handles as an "easy" component, i.e.,
 * exactly in the master problem, as the paper does with its structured
 * function. There is no LP to compare with: the value is printed, and if a
 * reference value is given the test fails if the two differ by more than tol,
 * relative to max( 1 , |ref| ).
 *
 * Usage: osbdo_fl_test instance BSPar [ ref [ tol ] ]
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/

#include "AbstractBlock.h"
#include "BlockSolverConfig.h"
#include "C05Function.h"
#include "FRealObjective.h"
#include "LagBFunction.h"
#include "LinearFunction.h"
#include "OneVarConstraint.h"
#include "Solver.h"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <unordered_map>

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS LogisticFunction ------------------------*/
/*--------------------------------------------------------------------------*/
/// the logistic loss of a set of samples, with a global pool
/** f( x ) = sum_s log( 1 + exp( - y_s a_s x ) ), with a_s the rows of a
 * dense matrix A and y_s in { -1 , 1 }; it is finite and continuously
 * differentiable everywhere, so its only linearization at x is the gradient
 *
 *   g = sum_s - y_s sigma( - y_s a_s x ) a_s ,  sigma( t ) = 1 / ( 1 + e^-t )
 *
 * and alpha = f( x ) - g x. The global pool keeps the ( g , alpha ) that
 * are stored and their non-negative combinations. */

class LogisticFunction : public C05Function {

/*--------------------------------------------------------------------------*/

 using v_col_var = std::vector< ColVariable * >;

 class v_iterator : public ThinVarDepInterface::v_iterator {
  public:
   explicit v_iterator( v_col_var::iterator itr ) : itr_( itr ) {}
   v_iterator * clone( void ) override { return( new v_iterator( itr_ ) ); }
   void operator++( void ) override { ++itr_; }
   reference operator*( void ) const override { return( **itr_ ); }
   pointer operator->( void ) const override { return( *itr_ ); }
   bool operator==( const ThinVarDepInterface::v_iterator & rhs )
    const override {
    auto t = dynamic_cast< const v_iterator * >( & rhs );
    return( t ? itr_ == t->itr_ : false );
    }
   bool operator!=( const ThinVarDepInterface::v_iterator & rhs )
    const override { return( ! ( *this == rhs ) ); }
  private:
   v_col_var::iterator itr_;
  };

 class v_const_iterator : public ThinVarDepInterface::v_const_iterator {
  public:
   explicit v_const_iterator( v_col_var::const_iterator itr ) : itr_( itr ) {}
   v_const_iterator * clone( void ) override {
    return( new v_const_iterator( itr_ ) );
    }
   void operator++( void ) override { ++itr_; }
   reference operator*( void ) const override { return( **itr_ ); }
   pointer operator->( void ) const override { return( *itr_ ); }
   bool operator==( const ThinVarDepInterface::v_const_iterator & rhs )
    const override {
    auto t = dynamic_cast< const v_const_iterator * >( & rhs );
    return( t ? itr_ == t->itr_ : false );
    }
   bool operator!=( const ThinVarDepInterface::v_const_iterator & rhs )
    const override { return( ! ( *this == rhs ) ); }
  private:
   v_col_var::const_iterator itr_;
  };

/*--------------------------------------------------------------------------*/

 public:

 /// A is n x p by rows, y has n entries in { -1 , 1 }, vars has p entries
 LogisticFunction( v_col_var && vars , std::vector< double > && A ,
                   std::vector< double > && y )
  : C05Function() , v_vars( std::move( vars ) ) , f_A( std::move( A ) ) ,
    f_y( std::move( y ) ) , f_val( Inf< FunctionValue >() ) ,
    f_g( v_vars.size() , 0 ) {
  for( Index i = 0 ; i < v_vars.size() ; ++i )
   f_pos.emplace( v_vars[ i ] , i );
  }

 ~LogisticFunction() override = default;

 void clear( void ) override {}

/*--------------------------------------------------------------------------*/

 void set_par( idx_type par , int value ) override {
  if( par == intGPMaxSz ) {
   if( value < 0 )
    throw( std::invalid_argument( "LogisticFunction::set_par: negative "
                                  "intGPMaxSz" ) );
   if( Index( value ) < v_pool.size() )
    v_pool.resize( value );
   f_gp_max = value;
   return;
   }
  if( par == intLPMaxSz )
   return;  // there is only one linearization anyway
  C05Function::set_par( par , value );
  }

 using C05Function::set_par;

 int get_int_par( idx_type par ) const override {
  if( par == intGPMaxSz )
   return( f_gp_max );
  return( C05Function::get_int_par( par ) );
  }

/*--------------------------------------------------------------------------*/

 void remove_variable( Index i , ModParam issueMod = eModBlck ) override {
  throw( std::logic_error( "LogisticFunction::remove_variable: not "
                           "allowed" ) );
  }

/*--------------------------------------------------------------------------*/

 int compute( bool changedvars = true ) override {
  const Index p = v_vars.size();
  const Index n = f_y.size();
  f_x.resize( p );
  for( Index j = 0 ; j < p ; ++j )
   f_x[ j ] = v_vars[ j ]->get_value();

  f_val = 0;
  std::fill( f_g.begin() , f_g.end() , 0 );
  for( Index s = 0 ; s < n ; ++s ) {
   const double * a = f_A.data() + s * p;
   double z = 0;
   for( Index j = 0 ; j < p ; ++j )
    z += a[ j ] * f_x[ j ];
   const double t = - f_y[ s ] * z;  // loss = log( 1 + e^t )
   f_val += ( t > 0 ) ? t + std::log1p( std::exp( - t ) )
                      : std::log1p( std::exp( t ) );
   const double sg = ( t >= 0 ) ? 1 / ( 1 + std::exp( - t ) )
                                : std::exp( t ) / ( 1 + std::exp( t ) );
   const double c = - f_y[ s ] * sg;
   for( Index j = 0 ; j < p ; ++j )
    f_g[ j ] += c * a[ j ];
   }
  f_alpha = f_val;
  for( Index j = 0 ; j < p ; ++j )
   f_alpha -= f_g[ j ] * f_x[ j ];
  return( kOK );
  }

 FunctionValue get_value( void ) override { return( f_val ); }

 bool is_convex( void ) override { return( true ); }

 bool is_continuously_differentiable( void ) const override {
  return( true );
  }

/*--------------------------------------------------------------------------*/

 void store_linearization( Index name , ModParam issueMod = eModBlck )
  override {
  if( name >= f_gp_max )
   throw( std::invalid_argument( "LogisticFunction::store_linearization: "
                                 "invalid name" ) );
  if( name >= v_pool.size() )
   v_pool.resize( name + 1 );
  v_pool[ name ].g = f_g;
  v_pool[ name ].alpha = f_alpha;
  v_pool[ name ].there = true;
  }

 bool is_linearization_there( Index name ) const override {
  return( ( name < v_pool.size() ) && v_pool[ name ].there );
  }

 void store_combination_of_linearizations(
                    c_LinearCombination & coefficients , Index name ,
                    ModParam issueMod = eModBlck ) override {
  if( name >= f_gp_max )
   throw( std::invalid_argument( "LogisticFunction::store_combination_of_"
                                 "linearizations: invalid name" ) );
  std::vector< double > g( v_vars.size() , 0 );
  double alpha = 0;
  for( auto & [ i , c ] : coefficients ) {
   const auto & e = entry( i );
   for( Index j = 0 ; j < g.size() ; ++j )
    g[ j ] += c * e.g[ j ];
   alpha += c * e.alpha;
   }
  if( name >= v_pool.size() )
   v_pool.resize( name + 1 );
  v_pool[ name ].g = std::move( g );
  v_pool[ name ].alpha = alpha;
  v_pool[ name ].there = true;
  }

 void delete_linearization( Index name , ModParam issueMod = eModBlck )
  override {
  if( name < v_pool.size() ) {
   v_pool[ name ].there = false;
   v_pool[ name ].g.clear();
   }
  }

/*--------------------------------------------------------------------------*/

 void get_linearization_coefficients( FunctionValue * g ,
                                      Range range = INFRange ,
                                      Index name = Inf< Index >() )
  override {
  const auto & src = ( name == Inf< Index >() ) ? f_g : entry( name ).g;
  range.second = std::min( range.second , Index( src.size() ) );
  for( Index j = range.first ; j < range.second ; ++j )
   *(g++) = src[ j ];
  }

 void get_linearization_coefficients( FunctionValue * g ,
                                      c_Subset & subset ,
                                      bool ordered = false ,
                                      Index name = Inf< Index >() )
  override {
  const auto & src = ( name == Inf< Index >() ) ? f_g : entry( name ).g;
  for( auto j : subset )
   *(g++) = src[ j ];
  }

 FunctionValue get_linearization_constant( Index name = Inf< Index >() )
  override {
  return( ( name == Inf< Index >() ) ? f_alpha : entry( name ).alpha );
  }

/*--------------------------------------------------------------------------*/

 Index get_num_active_var( void ) const override { return( v_vars.size() ); }

 Index is_active( const Variable * var ) const override {
  auto it = f_pos.find( var );
  return( it == f_pos.end() ? Inf< Index >() : it->second );
  }

 void map_active( c_Vec_p_Var & vars , Subset & map , bool ordered = false )
  const override {
  map.resize( vars.size() );
  for( Index i = 0 ; i < vars.size() ; ++i ) {
   map[ i ] = is_active( vars[ i ] );
   if( map[ i ] == Inf< Index >() )
    throw( std::invalid_argument( "LogisticFunction::map_active: not an "
                                  "active Variable" ) );
   }
  }

 Variable * get_active_var( Index i ) const override { return( v_vars[ i ] ); }

 ThinVarDepInterface::v_iterator * v_begin( void ) override {
  return( new v_iterator( v_vars.begin() ) );
  }
 ThinVarDepInterface::v_const_iterator * v_begin( void ) const override {
  return( new v_const_iterator( v_vars.begin() ) );
  }
 ThinVarDepInterface::v_iterator * v_end( void ) override {
  return( new v_iterator( v_vars.end() ) );
  }
 ThinVarDepInterface::v_const_iterator * v_end( void ) const override {
  return( new v_const_iterator( v_vars.end() ) );
  }

/*--------------------------------------------------------------------------*/

 protected:

 void print( std::ostream & output ) override {
  output << "LogisticFunction [" << this << "] with " << v_vars.size()
         << " Variable and " << f_y.size() << " samples";
  }

/*--------------------------------------------------------------------------*/

 private:

 struct Entry {
  std::vector< double > g;
  double alpha = 0;
  bool there = false;
  };

 const Entry & entry( Index name ) const {
  if( ! is_linearization_there( name ) )
   throw( std::invalid_argument( "LogisticFunction: no linearization " +
                                 std::to_string( name ) + " in the pool" ) );
  return( v_pool[ name ] );
  }

 v_col_var v_vars;
 std::unordered_map< const Variable * , Index > f_pos;
 std::vector< double > f_A , f_y , f_x;
 FunctionValue f_val , f_alpha = 0;
 std::vector< double > f_g;
 Index f_gp_max = 0;
 std::vector< Entry > v_pool;

 };  // end( class LogisticFunction )

/*--------------------------------------------------------------------------*/

struct Instance {
 int n , p , M;
 double lambda;
 std::vector< double > A , y;  // A is n x p by rows
 };

static Instance read_instance( const std::string & fn )
{
 std::ifstream f( fn );
 if( ! f.is_open() )
  throw( std::invalid_argument( "read_instance: cannot open " + fn ) );
 Instance I;
 f >> I.n >> I.p >> I.M >> I.lambda;
 I.A.resize( std::size_t( I.n ) * I.p );
 I.y.resize( I.n );
 for( int s = 0 ; s < I.n ; ++s ) {
  for( int j = 0 ; j < I.p ; ++j )
   f >> I.A[ std::size_t( s ) * I.p + j ];
  f >> I.y[ s ];
  }
 if( ( ! f ) || ( I.n % I.M ) )
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

static AbstractBlock * build( const Instance & I )
{
 auto root = new AbstractBlock();
 auto x = new std::vector< ColVariable >( I.p );
 root->add_static_variable( *x );
 LinearFunction::v_coeff_pair zero;
 for( auto & xj : *x )
  zero.emplace_back( & xj , 0.0 );
 auto robj = new FRealObjective( root , new LinearFunction( std::move( zero ) ) );
 robj->set_sense( Objective::eMin );
 root->set_objective( robj );

 auto & nested = root->access_nested_Blocks();
 const int ns = I.n / I.M;
 for( int i = 0 ; i < I.M ; ++i ) {
  std::vector< ColVariable * > vars( I.p );
  for( int j = 0 ; j < I.p ; ++j )
   vars[ j ] = & ( *x )[ j ];
  std::vector< double > A( I.A.begin() + std::size_t( i ) * ns * I.p ,
                           I.A.begin() + std::size_t( i + 1 ) * ns * I.p );
  std::vector< double > y( I.y.begin() + i * ns ,
                           I.y.begin() + ( i + 1 ) * ns );
  auto lf = new LogisticFunction( std::move( vars ) , std::move( A ) ,
                                  std::move( y ) );
  auto sub = new AbstractBlock( root );
  auto sobj = new FRealObjective( sub , lf );
  sobj->set_sense( Objective::eMin );
  sub->set_objective( sobj );
  nested.push_back( sub );
  }

 // lambda || x ||_1 as an easy LagBFunction: max { lambda w' x : |w| <= 1 }
 {
  auto inner = new AbstractBlock();
  auto w = new std::vector< ColVariable >( I.p );
  inner->add_static_variable( *w );
  auto box = new std::vector< BoxConstraint >( I.p );
  for( int j = 0 ; j < I.p ; ++j ) {
   ( *box )[ j ].set_variable( & ( *w )[ j ] );
   ( *box )[ j ].set_lhs( -1 );
   ( *box )[ j ].set_rhs( 1 );
   }
  inner->add_static_constraint( *box );
  auto iobj = new FRealObjective( inner , new LinearFunction() );
  iobj->set_sense( Objective::eMax );   // convex LagBFunction
  inner->set_objective( iobj );

  auto lbf = new LagBFunction( inner );
  LagBFunction::v_dual_pair dp;
  for( int j = 0 ; j < I.p ; ++j )
   dp.emplace_back( & ( *x )[ j ] , new LinearFunction(
                LinearFunction::v_coeff_pair{ { & ( *w )[ j ] , I.lambda } } ) );
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
 if( argc < 3 ) {
  std::cerr << "Usage: " << argv[ 0 ] << " instance BSPar [ ref [ tol ] ]"
            << std::endl;
  return( 1 );
  }

 using clk = std::chrono::steady_clock;
 auto since = []( clk::time_point t0 ) {
  return( std::chrono::duration< double >( clk::now() - t0 ).count() );
  };

 auto t0 = clk::now();
 const auto I = read_instance( argv[ 1 ] );
 std::cout << std::setprecision( 12 );
 std::cout << "instance " << argv[ 1 ] << ": n = " << I.n << ", p = " << I.p
           << ", M = " << I.M << ", lambda = " << I.lambda << ", read in "
           << since( t0 ) << " s" << std::endl;

 t0 = clk::now();
 auto blk = build( I );
 double tb = since( t0 );
 auto s = attach( blk , argv[ 2 ] );
 s->set_log( & std::cout );
 t0 = clk::now();
 auto st = s->compute();
 double tt = since( t0 );
 double v = s->get_var_value();
 std::cout << "bundle    : status " << st << ", value " << v
           << ", iterations " << s->get_elapsed_iterations() << ", time "
           << tt << " s (build " << tb << " s)" << std::endl;

 bool ok = ( st == Solver::kOK );
 if( argc > 3 ) {
  const double ref = std::atof( argv[ 3 ] );
  const double tol = ( argc > 4 ) ? std::atof( argv[ 4 ] ) : 1e-6;
  const double gap = std::abs( v - ref ) / std::max( 1.0 , std::abs( ref ) );
  std::cout << "reference : " << ref << ", gap " << gap << std::endl;
  ok = ok && ( gap <= tol );
  }
 std::cout << ( ok ? "OK" : "KO" ) << std::endl;
 return( ok ? 0 : 1 );
 }

/*--------------------------------------------------------------------------*/
/*----------------------- End File test_osbdo_fl.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
