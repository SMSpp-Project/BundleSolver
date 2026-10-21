/*--------------------------------------------------------------------------*/
/*------------------------- File BundleSolver.cppm -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Interface of the C++20 module smspp.bundlesolver, which exports every name
 * that the headers of BundleSolver declare at namespace scope, together with
 * the modules it depends on; a translation unit can then write
 *
 *     import smspp.bundlesolver;
 *
 * in place of the #include of the headers. The macros (e.g., those of the
 * factory in SMSFactory.h) are not exported by a module, so a file that uses
 * them still includes the header that defines them.
 *
 * This file is written by SMS++/cmake/smspp_module_exports.py, which reads
 * the names from the AST of the headers: it is not to be edited by hand.
 *
 * The headers BundleSolverML.h, BundleSolverMLCheckpoint.h,
 * BundleSolverMLNet.h, LSTMBlock.h are not in this module: they need a
 * library of their own, or only a source includes them.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Enrico Calandrini \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Enrico Calandrini, Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*------------------------- GLOBAL MODULE FRAGMENT -------------------------*/
/*--------------------------------------------------------------------------*/

module;

#include "BundleSolver.h"
#include "ParallelBundleSolver.h"

/*--------------------------------------------------------------------------*/
/*----------------------------- MODULE PURVIEW -----------------------------*/
/*--------------------------------------------------------------------------*/

export module smspp.bundlesolver;

export import smspp.milpsolver;

export namespace SMSpp_di_unipi_it {

 // BundleSolver.h
 using SMSpp_di_unipi_it::BundleSolver;
 using SMSpp_di_unipi_it::BundleSolverState;

 // ParallelBundleSolver.h
 using SMSpp_di_unipi_it::ParallelBundleSolver;
}

/*--------------------------------------------------------------------------*/
/*----------------------- End File BundleSolver.cppm -----------------------*/
/*--------------------------------------------------------------------------*/
