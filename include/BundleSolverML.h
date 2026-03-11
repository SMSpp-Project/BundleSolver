#pragma once

#include <ctime>
#include <memory>
#include <queue>
#include <string>
#include <torch/torch.h>
#include "BundleSolver.h"
#include <torch/optim.h>

// ============================================================================
// VERBOSE macro: set to 1 to enable debug output, 0 to disable.
// Must be consistent with the definition in BundleSolverML.cpp.
// ============================================================================
#ifndef VERBOSE
  #define VERBOSE 1
#endif

#if VERBOSE
#define BML_LOG(x) std::cout << x
#else
#define BML_LOG(x) do {} while(0)
#endif

/*--------------------------------------------------------------------------*/
/*-------------------------- A SIMPLE NEURAL NETWORK -----------------------*/
/*--------------------------------------------------------------------------*/
/**
 * @struct Net
 * @brief Two-layer feedforward neural network used to predict the step-size t.
 *
 * Architecture:
 *   - Input layer  : 20 features describing the current bundle solver state.
 *   - Hidden layer : 16 units with Softplus activation.
 *   - Output layer : 1 unit with Softplus activation, producing t > 0.
 *
 * The Softplus output ensures the predicted step-size is always strictly
 * positive, which is required for the bundle method's convergence.
 *
 * Sharing across instances
 * ------------------------
 * A single Net instance can be shared across multiple BundleSolverML objects
 * (one per MMCF instance) via BundleSolverML::set_shared_net() /
 * BundleSolverML::get_shared_net(). This is the recommended pattern for
 * training: all solvers read from and write to the same parameter tensors,
 * so every Backward() call accumulates gradients into the shared weights.
 *
 * Persistence:
 *   Weights are saved/loaded via BundleSolverML::SaveModel() and
 *   BundleSolverML::LoadModel() (TorchScript .pt format, cross-compatible
 *   with Python's torch.load()).
 */
struct Net : torch::nn::Module {
  /// First fully-connected layer: 20 inputs → 16 hidden units
  torch::nn::Linear fc1{nullptr};
  /// Second fully-connected layer: 16 hidden units → 1 output (step-size)
  torch::nn::Linear fc2{nullptr};

  /// Minimum admissible step-size (used by the commented sigmoid clamping)
  float min_val = 0.000001f;
  /// Maximum admissible step-size (used by the commented sigmoid clamping)
  float max_val = 10000.0f;

  /// Constructor: registers the two linear layers with the LibTorch module system.
  Net() {
    fc1 = register_module("fc1", torch::nn::Linear(20, 16));
    fc2 = register_module("fc2", torch::nn::Linear(16, 1));
  }

  /**
   * @brief Forward pass: maps a feature vector to a positive step-size.
   * @param x  Input feature tensor of shape {20} (or {batch, 20}).
   * @return   Scalar tensor representing the predicted step-size t > 0.
   */
  torch::Tensor forward(torch::Tensor x) {
    x = torch::softplus(fc1->forward(x)); // hidden layer with Softplus activation
    x = fc2->forward(x);                  // linear output layer
    // Alternative clamped output (disabled):
    // x = min_val + (max_val - min_val) * torch::sigmoid(x);
    return torch::softplus(x)+1.0e-8;            // ensure strictly positive output
  }
};

/*--------------------------------------------------------------------------*/
/*--------------------- SMS++ NAMESPACE DECLARATION ------------------------*/
/*--------------------------------------------------------------------------*/
/// Namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it {

/**
 * @class BundleSolverML
 * @brief Bundle solver with a machine-learning heuristic for step-size selection.
 *
 * Extends BundleSolver by replacing the classical rule-based step-size
 * heuristic with a neural network (Net) that predicts the proximal
 * regularization parameter t from a 20-dimensional feature vector built from
 * the current solver state.
 *
 * Shared network model
 * --------------------
 * By default each BundleSolverML owns a private Net instance (nn_owned_),
 * allocated at construction time. For training across multiple instances the
 * caller should create a single shared Net and inject it into every solver:
 *
 * @code
 *   // Create one shared network (lives as long as training continues)
 *   auto shared_net = std::make_shared<Net>();
 *
 *   // ... for each instance:
 *   auto* bml = dynamic_cast<BundleSolverML*>(Slvr2->get_inner_Solver());
 *   bml->set_shared_net(shared_net);  // inject shared weights
 *   Slvr2->compute(false);
 *   bml->Backward();                  // updates shared_net's parameters
 * @endcode
 *
 * After set_shared_net(), the solver's nn reference points to the shared
 * object. All forward and backward passes operate on the same parameter
 * tensors, so gradient updates from every instance accumulate into the same
 * network — exactly the online SGD pattern needed for multi-instance training.
 *
 * To revert to the private network, call clear_shared_net().
 *
 * Model persistence
 * -----------------
 * @code
 *   bml->SaveModel("checkpoints/model.pt");
 *   bml->LoadModel("checkpoints/model.pt");
 * @endcode
 * Both methods operate on whichever network is currently active (shared or
 * private). The .pt format is readable from Python via torch.load().
 *
 * Key design choices:
 *   - Copy and move operations are deleted to protect the LibTorch nn::Module
 *     internals, which manage reference-counted parameter tensors.
 *   - All tensors needed for differentiation are stored in per-iteration
 *     vectors (phi_vecs, w_vecs, Gs, …) and consumed during Backward().
 */
class BundleSolverML : public BundleSolver {
public:
  using Index = Function::Index;

  /// A vector of floating-point values (e.g., a single iterate)
  typedef std::vector<double> VecMem;
  /// A collection of VecMem vectors (e.g., previous iterates)
  using MemMultiVector = std::vector<VecMem>;

  // --------------------------------------------------------------------------
  // Constructor / destructor
  // --------------------------------------------------------------------------

  /**
   * @brief Default constructor.
   *
   * Allocates the private Net (nn_owned_) and sets nn to point to it.
   * The feature vector is resized to size_features = 20.
   */
  BundleSolverML(void) : BundleSolver() {
    nn_owned_ = std::make_shared<Net>();
    nn = nn_owned_.get();             // point to the private network by default
    size_features = 20;
    features.resize(size_features);
  }

  /// Copy constructor deleted: copying would invalidate nn's parameter registry
  BundleSolverML(const BundleSolverML&) = delete;
  /// Copy assignment deleted: same reason as copy constructor
  BundleSolverML& operator=(const BundleSolverML&) = delete;

  /// Move constructor deleted: moving might invalidate internal nn pointers
  BundleSolverML(BundleSolverML&&) = delete;
  /// Move assignment deleted: same reason as move constructor
  BundleSolverML& operator=(BundleSolverML&&) = delete;

  /**
   * @brief Destructor.
   *
   * Detaches the Block before destruction (inherited pattern from BundleSolver).
   * The shared_net_ shared_ptr is released here; the Net is destroyed only if
   * no other BundleSolverML still holds a reference to it.
   */
  virtual ~BundleSolverML() { set_Block(nullptr); }

  // --------------------------------------------------------------------------
  // Shared network management
  // --------------------------------------------------------------------------

  /**
   * @brief Injects an externally-owned shared network into this solver.
   *
   * After this call, all forward passes (in Heuristic()) and backward passes
   * (in Backward()) operate on @p net's parameters instead of the private
   * nn_owned_ instance.  Because multiple BundleSolverML objects can point
   * to the same Net, gradient updates from all of them accumulate into the
   * same tensors — enabling true multi-instance online training without any
   * weight-file I/O between instances.
   *
   * The shared_ptr keeps the Net alive as long as at least one solver holds
   * a reference, regardless of the order in which solvers are destroyed.
   *
   * @param net  Shared pointer to the network to use. Must not be nullptr.
   *
   * Example (training loop):
   * @code
   *   auto shared_net = std::make_shared<Net>();
   *   for (auto& path : train_paths) {
   *     auto* bml = get_bml_for(path);      // creates BundleSolverML
   *     bml->set_shared_net(shared_net);    // inject
   *     solve_and_backward(bml);            // updates shared_net in-place
   *   }
   *   bml->SaveModel("model.pt");           // save final weights
   * @endcode
   */
  void set_shared_net(std::shared_ptr<Net> net) {
    if (!net)
      throw std::invalid_argument("set_shared_net: net must not be nullptr");
    shared_net_ = std::move(net);
    nn = shared_net_.get();           // redirect the raw pointer used internally
  }

  /**
   * @brief Returns the currently active network as a shared_ptr.
   *
   * - If set_shared_net() was called, returns the shared network.
   * - Otherwise returns the private nn_owned_ instance (still usable as a
   *   shared_ptr, e.g. to pass to another solver).
   *
   * @return Shared pointer to the active Net (never nullptr).
   */
  std::shared_ptr<Net> get_shared_net() const {
    return shared_net_ ? shared_net_ : nn_owned_;
  }

  /**
   * @brief Reverts to the private network, releasing the shared reference.
   *
   * After this call, nn points back to nn_owned_ and the solver's
   * shared_ptr to the external network is released (the Net is destroyed
   * only if no other solver holds a reference).
   */
  void clear_shared_net() {
    shared_net_.reset();
    nn = nn_owned_.get();
  }

  // --------------------------------------------------------------------------
  // Core ML-specific methods
  // --------------------------------------------------------------------------

  /**
   * @brief Step-size heuristic overriding the base class rule.
   *
   * Called at each bundle iteration. Builds the feature vector from the
   * current solver state, runs a forward pass through the active network
   * (shared or private), and returns the predicted step-size t. Also
   * collects and stores all tensors (G, Q, alpha, w, gs_aggreg, …) needed
   * by the subsequent Backward() call.
   *
   * @param whch  Unused; kept for interface compatibility with BundleSolver.
   * @return      The predicted step-size t > 0.
   */
  virtual HpNum Heuristic(Index whch) override;

  /**
   * @brief Computes gradients and updates the active network parameters.
   *
   * Iterates over all stored iterations, accumulates a scalar loss
   * (weighted by step type and function value), then calls loss.backward()
   * to propagate gradients into the active network's parameters.
   * If a shared network is active, the gradient update affects all solvers
   * that share the same Net.
   *
   * Should be called once per solve, after compute() returns.
   */
  void Backward();

  /**
   * @brief Clears all per-solve tensor buffers.
   *
   * Empties phi_vecs, w_vecs, Gs, Qs, alphaS, Gs_aggreg, coeff_vecs,
   * thetaS, tS, and FiS. Call this between successive solves on the same
   * BundleSolverML instance to free memory and avoid stale data.
   *
   * Not needed when a fresh BundleSolverML is created per instance.
   */
  void ClearBuffers() {
    phi_vecs.clear();
    w_vecs.clear();
    coeff_vecs.clear();
    Gs.clear();
    Gs_aggreg.clear();
    Qs.clear();
    alphaS.clear();
    thetaS.clear();
    phiS.clear();
    lastIndexS.clear();
    tS.clear();
    FiS.clear();
  }

  /**
   * @brief Constructs the search direction tensor for iteration f with step t.
   *
   * Wraps BundleSolverML_W::apply so the computation is differentiable
   * w.r.t. the network's output (used internally by Backward()).
   *
   * @param f  Iteration index into the stored tensor vectors.
   * @param t  Step-size value (scalar double).
   * @return   Scalar tensor representing sum(w), connected to the autograd graph.
   */
  torch::Tensor w(size_t f, double t);

  // --------------------------------------------------------------------------
  // Model persistence
  // --------------------------------------------------------------------------

  /**
   * @brief Saves the active network weights to a TorchScript archive (.pt).
   *
   * Works on whichever network is currently active (shared or private).
   * The file is cross-compatible with Python:
   * @code
   *   import torch; sd = torch.load("model.pt")
   * @endcode
   *
   * @param filepath  Output path (directory must exist; file is overwritten).
   * @throws std::runtime_error on I/O or LibTorch error.
   */
  void SaveModel(const std::string& filepath);

  /**
   * @brief Loads weights from a TorchScript archive into the active network.
   *
   * Works on whichever network is currently active (shared or private).
   * Architecture must match exactly (same layer names and tensor shapes).
   * Leaves the network in eval() mode; call nn->train() before Backward().
   *
   * @param filepath  Path to the .pt file (must exist).
   * @throws std::runtime_error if the file is missing, malformed, or the
   *         architecture does not match.
   */
  void LoadModel(const std::string& filepath);

  // --------------------------------------------------------------------------
  // Default integer parameter overrides
  // --------------------------------------------------------------------------

  /**
   * @brief Returns the default value of integer parameter @p par.
   *
   * See full parameter table in the implementation for details.
   */
  [[nodiscard]] int get_dflt_int_par(idx_type par) const override {
    static const std::array<int, 22> dflt_int_par = {
        0,    // intBPar1
        100,  // intBPar2
        1,    // intBPar3
        1,    // intBPar4
        0,    // intBPar6
        3,    // intBPar7
        1000, // intMnSSC
        1000, // intMnNSC
        12,   // inttSPar1
        2,    // intMaxNrEvls
        1,    // intDoEasy
        2,    // intWZNorm
        0,    // intFrcLstSS
        0,    // intTrgtMng
        0,    // intMPName
        0,    // intMPlvl
        0,    // intQPmp1
        0,    // intQPmp2
        4,    // intOSImp1
        0,    // intOSImp2
        1,    // intOSImp3
        2     // intRstAlg
    };
    if ((par >= intLastParCDAS) && (par < intLastBndSlvPar))
      return dflt_int_par[par - intLastParCDAS];
    return CDASolver::get_dflt_int_par(par);
  }

  /// Macro that registers BundleSolverML in the SMS++ solver factory
  SMSpp_insert_in_factory_h;

  // --------------------------------------------------------------------------
  // Data members (public for direct access by train_and_test.cpp)
  // --------------------------------------------------------------------------

  /// Raw pointer to the active network (shared or private). Never nullptr.
  /// Use set_shared_net() / clear_shared_net() to redirect it; do not
  /// reassign this pointer directly.
  Net* nn = nullptr;

  /// Feature vector built at each iteration and fed to the network (size 20)
  VecMem features;

  /// Number of features (always 20, matching the network's input dimension)
  int size_features;

  /// Starting point Lambda0 (stability center at the beginning of the solve)
  Vec_VarValue Lambda0;

  /// History of previous iterates (reserved for future use)
  MemMultiVector pw;

  // Per-iteration tensors stored for the backward pass
  std::vector<torch::Tensor> phi_vecs;    ///< Feature tensors (network inputs)
  std::vector<torch::Tensor> w_vecs;      ///< Search direction tensors w
  std::vector<torch::Tensor> coeff_vecs;  ///< Step-type signs (+1 SS / -1 NS)
  std::vector<torch::Tensor> Gs;          ///< Subgradient matrices G
  std::vector<torch::Tensor> Gs_aggreg;  ///< Aggregated subgradient vectors
  std::vector<torch::Tensor> Qs;          ///< Gram matrices Q = G @ G^T
  std::vector<torch::Tensor> alphaS;      ///< Linearization error vectors
  std::vector<torch::Tensor> thetaS;      ///< Dual multiplier vectors
  std::vector<torch::Tensor> phiS;        ///< Reserved for future use
  std::vector<torch::Tensor> lastIndexS;  ///< Reserved for future use
  std::vector<double>        tS;          ///< Predicted step-sizes
  std::vector<double>        FiS;         ///< Best upper bounds on f*

private:
  // --------------------------------------------------------------------------
  // Network ownership
  // --------------------------------------------------------------------------

  /// Private network allocated at construction (used when no shared net is set)
  std::shared_ptr<Net> nn_owned_;

  /// Shared network injected by set_shared_net() (nullptr when using nn_owned_)
  std::shared_ptr<Net> shared_net_;
  
  /// Adam optimizer, lazily initialized on first Backward() call.
  /// Held as a member so moment estimates persist across successive Backward() calls.
  std::unique_ptr<torch::optim::Adam> optimizer_;


}; // class BundleSolverML

} // namespace SMSpp_di_unipi_it
