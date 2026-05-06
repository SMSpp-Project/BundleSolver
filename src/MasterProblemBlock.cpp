/*--------------------------------------------------------------------------*/
/*-------------------- File MasterProblemBlock.cpp ------------------------*/
/*--------------------------------------------------------------------------*/

#include "MasterProblemBlock.h"

#include <algorithm>
#include <stdexcept>
#include <utility>
#include <ostream>
#include <string>
#include <vector>
#include <sstream>

namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*-------------------------- utility constants -----------------------------*/
/*--------------------------------------------------------------------------*/

namespace
{
 static const std::string empty_string;
 static const std::vector< int > empty_vint;
 static const std::vector< std::string > empty_vstr;
}

/*--------------------------------------------------------------------------*/
/*------------------------- clear / reinitialize ---------------------------*/
/*--------------------------------------------------------------------------*/

void MasterProblemBlock::clear()
{
 /* Reset the structural information owned by the master problem block. */

 MaxBSize = 0;

 NoEasyCmps = 0;
 NoHardCmps = 0;

 EasyCmps.clear();
 HardCmps.clear();
}

void MasterProblemBlock::SetDim( int MxBSz , int NVars , 
    int NrFi , int NrFiEasy ){
 
 //TBD
}

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::register_Solver( std::string solv_cfg_filename )
{
 // Check if a configuration file has been provided for the internal solver
 // of MPBlock
 BlockSolverConfig * MPBSC = nullptr;
 if( ! solv_cfg_filename.empty() ) {
  auto cfg = Configuration::deserialize( solv_cfg_filename );
  if( ! ( MPBSC = dynamic_cast< BlockSolverConfig * >( cfg ) ) ){
   delete cfg;
   
   // throw an error because the Configuration file is not a BlockSolverConfig
   throw( std::invalid_argument( "The provided Configuration file for the "
            "solver of the MaterProblemBlock is not a BlockSolverConfig" ) );
  }
 }
 else{
  BLOG( 1 , std::endl << " Warning: a Configuration file for the inner Solver "
    "of the Master Problem Block has not been provided. Using the default "
    "GRBMILPSolver." );

  MPBSC = use_default_Solver();
 }

 // Set the Solver in the Master Problem Block
 MPBSC->apply( this ); 

 MPBSC->clear();

} // end( MasterProblemBlock::register_Solver )

/*--------------------------------------------------------------------------*/

BlockSolverConfig * MasterProblemBlock::use_default_Solver( void )
{
 // prepare the istream for the BlockSolverConfig
 std::istringstream input(
    "1 1 GRBMILPSolver 1 ComputeConfig 1 0 0 0 0 0 0 *" );

 // Generate the simple BlockSolverConfig
 return new BlockSolverConfig(input);

} // end( MasterProblemBlock::use_default_Solver )

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::CreateEmptyMP( stabilization_type Stbl , int NoCmps,
  int DoEasyCmp, int NoEasy , std::vector< Bool > IsEasy , )
{
 // Clear all the data structures
 clear();

 // Set the class specific fields
 StblType = Stbl;
 NoTotCmps = NoCmps;
 DoEasy = DoEasyCmp;
 NoEasyCmps = NoEasy;
 IsEasyCmp = IsEasy;

 // Immediately check if we treat differentl. In this case the MP
 // will be initialized in its primal version. Otherwise, the dual one will 
 // be used
 if( NoEasyCmps == 0 )
  CreatePrimalMP();
 else
  CreateDualMP();

} // end( MasterProblemBlock::CreateEmptyMP )

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::CreateDualMP( )
{
 // Generate the dual

 // Set number of variables

 // Generate the static variables - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // Firstly we add the non-negative "dual" multiplier corresponding to the 
 // model, global lower bound and level constraints.
 auto Var_lambda = new ColVariable( this );
 auto Var_r = new ColVariable( this );
 auto Var_omega = new ColVariable( this );

 Var_lambda.is_positive( true );
 Var_r.is_positive( true );

 // Here we have to check if the level stabilization is currently being
 // used. If not, the dual multiplier is simply fixed to be 0.
 if( Stbl == kLevel || Stbl == kDoublyStabilized )
  Var_omega.is_positive( true );
 else{
  Var_omega.is_fixed( true );
  Var_omega.set_value( 0 );
 }

 // Add the static variables to the block
 add_static_variable( Var_lambda, "lambda" );
 add_static_variable( Var_r, "r" );
 add_static_variable( Var_omega, "omega" );

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // Now we introduce the static variables γ^k associated with the
 // lower bounds on each hard component. Note that if LB^k = -∞,
 // then the dual objective function will force γ^k = 0.
 // We add them as a group of static variables having dimension 
 // NoHardCmps = NoTotCmps - NoEasyCmps.
 for( k = 0 ; k < NoTotCmps - NoEasyCmps ; k ++ ){
  auto Var_lambda = new ColVariable( this );
 }

 // Generate the static constraints - - - - - - - - - - - - - - - - - 
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

 // The first constraint regards the bounds on the model, involving 
 // the level constraint and the global lower model, i.e. 
 // \lambda + r - \omega = 1

 // Define pairs of (var, coefficient)
 v_coeff_pair lin_terms = {
  { Var_lambda , 1},
  { Var_r , 1},
  { Var_omega , 1}
 };

 // Create the constraint as an equality one with fixed RHS = LHS = 1
 // and add it to the block
 auto first_cons = new FRowConstraint( this , 1 , 1 , 
                        new LinearFunction( lin_terms ) );
 add_static_constraint( first_cons , "global_bounds" );

 // add the component-specific rows - - - - - - - - - - - - - - - - - 
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

 // NOTE: The easy component constraints 
 //   u^k E^k + λ e^k = 0, ∀k ∈ E ,
 // are not added at this level, as they should be already contained 
 // in the inner block of each LagbFunction, which is registered as
 // a sub-block of MPB. Hence, each solver will specifically construct
 // them when loading the sub-blocks.

 // add the hard-component rows - - - - - - - - - - - - - - - - - - - 
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 for( int k = 0 ; k < NoTotCmps ; k++ ){
  if( ! IsEasyCmp[ k ] ){
    // Hard component

    // For the k-th hard component we just have to add the simplex 
    // constraint:   λ = \sum_{i ∈ β^k_SG} θ_i^k + γ^k

  }
 }

 // Generate the dynamic constraints (?)

 // Start generating the objective function to be maximized
 Objective * MPOF = Objective( this );

 // Set the problem as a maximization one
 MPOF->set_sense( eMax ); // eMax = 1
}

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::load_problem( void )
{
 // Check if any solver has been registered to the Block
 auto solvers_list = this->get_registered_solvers();
 if( solvers_list.empty() )
   throw( std::logic_error(
	      "No solver has been registered to the Master Problem Block" ) );
 else{
   for( auto solver : *solvers_list ){
    // Load in the Solver the representation of the Master Problem using all 
    // the registered sub-blocks
    solver->load_problem();

    // Now if there are easy components, we have to separately manage them.
    if( NoEasyCmps > 0 ){
      // TBD
    }
   }
 }
}

/*--------------------------------------------------------------------------*/
/*--------------------------- set parameters -------------------------------*/
/*--------------------------------------------------------------------------*/

void MasterProblemBlock::set_par( idx_type par , int value )
{
 /* No MPBlock-specific integer parameters are currently defined in the
  * header enum, so forward to the base Block interface.
  */
 Block::set_par( par , value );
}

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::set_par( idx_type par , double value )
{
 /* No MPBlock-specific double parameters are currently defined. */
 Block::set_par( par , value );
}

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::set_par( idx_type par , std::string && value )
{
 /* No MPBlock-specific string parameters are currently defined. */
 Block::set_par( par , std::move( value ) );
}

/*--------------------------------------------------------------------------*/
/* NOTE:
 * The header contains a duplicated signature:
 *
 *   void set_par( idx_type par , std::string && value ) override;
 *
 * twice in a row, while the comment says "vector-of-int".
 * The implementation below assumes the intended declaration is:
 *
 *   void set_par( idx_type par , std::vector< int > && value ) override;
 *
 * Therefore there is no second std::string&& overload here.
 */
/*--------------------------------------------------------------------------*/

void MasterProblemBlock::set_par( idx_type par , std::vector< int > && value )
{
 /* No MPBlock-specific vector<int> parameters are currently defined. */
 Block::set_par( par , std::move( value ) );
}

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::set_par(
 idx_type par , std::vector< std::string > && value )
{
 /* No MPBlock-specific vector<string> parameters are currently defined. */
 Block::set_par( par , std::move( value ) );
}

/*--------------------------------------------------------------------------*/

void MasterProblemBlock::set_log( std::ostream * log_stream )
{
 Block::set_log( log_stream );
}

/*--------------------------------------------------------------------------*/
/*---------------------------- get parameters ------------------------------*/
/*--------------------------------------------------------------------------*/

int MasterProblemBlock::get_int_par( idx_type par ) const
{
 /* No MPBlock-specific integer parameters currently implemented. */
 return Block::get_int_par( par );
}

/*--------------------------------------------------------------------------*/

double MasterProblemBlock::get_dbl_par( idx_type par ) const
{
 /* No MPBlock-specific double parameters currently implemented. */
 return Block::get_dbl_par( par );
}

/*--------------------------------------------------------------------------*/

const std::string & MasterProblemBlock::get_str_par( idx_type par ) const
{
 /* No MPBlock-specific string parameters currently implemented. */
 return Block::get_str_par( par );
}

/*--------------------------------------------------------------------------*/

const std::vector< int > & MasterProblemBlock::get_vint_par(
 idx_type par ) const
{
 /* No MPBlock-specific vector<int> parameters currently implemented. */
 return Block::get_vint_par( par );
}

/*--------------------------------------------------------------------------*/

const std::vector< std::string > & MasterProblemBlock::get_vstr_par(
 idx_type par ) const
{
 /* No MPBlock-specific vector<string> parameters currently implemented. */
 return Block::get_vstr_par( par );
}

}  // end( namespace SMSpp_di_unipi_it )