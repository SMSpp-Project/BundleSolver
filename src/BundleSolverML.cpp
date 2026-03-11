#include <iostream>
#include <stdexcept>
#include <cmath>

#include "BundleSolverML.h"

// ============================================================================
// BML_LOG: renamed from LOG to avoid clash with LibTorch's c10 LOG(n) macro.
// ============================================================================
#define VERBOSE 1
#if VERBOSE
  #define BML_LOG(x) std::cout << x
#else
  #define BML_LOG(x) do {} while(0)
#endif

// ============================================================================
// BundleSolverML_Fi
//
// Custom autograd Function used to attach the discounted search-direction
// accumulator w_cum to the loss computation graph.
//
// Signature: apply(input, gs_aggreg_scalar, fi_val)
//   input           : the tensor whose gradient we want (w_cum, shape {N})
//   gs_aggreg_scalar: a SCALAR tensor (gs.sum()), used as a constant multiplier
//                     in backward — it is NOT differentiated.
//   fi_val          : a double constant, not differentiated.
//
// Forward:  returns fi_val as a scalar (shape {}) with the same dtype/device
//           as input.
//
// Backward: the upstream gradient dL/d(output) is a scalar {}.
//           The gradient w.r.t. input must have EXACTLY the same shape as
//           input, i.e. {N}.
//           Chain rule: dL/d(input) = dL/d(output) * gs_aggreg_scalar * ones(N)
//           Because output = fi_val (a constant), and the loss uses this output
//           as a multiplicative factor, the gradient flows back as a uniform
//           scaling of all elements of input by gs_aggreg_scalar.
//
// PREVIOUS BUG: backward returned grad_out * gs_aggreg where gs_aggreg had
// shape {bundle_size_at_iteration_f}.  w_cum accumulates across iterations
// so its shape (set at iteration 0) != bundle_size at later iterations.
// Fix: gs_aggreg is saved as a SCALAR (gs.sum()) so backward always returns
// a tensor of shape {N} = shape of input, matching LibTorch's requirement.
// ============================================================================
struct BundleSolverML_Fi : public torch::autograd::Function<BundleSolverML_Fi> {
  static torch::Tensor forward(torch::autograd::AutogradContext* ctx,
                               torch::Tensor input,        // w_cum, shape {N}
                               torch::Tensor gs_scalar,    // gs.sum(), shape {}
                               double        fi_val) {
    // Save gs_scalar (shape {}) for backward — shape is independent of bundle size
    ctx->save_for_backward({gs_scalar});
    // Return fi_val as a scalar tensor matching input's dtype/device
    return torch::tensor(fi_val, input.options()).squeeze();
  }

  static torch::autograd::tensor_list backward(
      torch::autograd::AutogradContext* ctx,
      torch::autograd::tensor_list grad_outputs) {
    auto saved      = ctx->get_saved_variables();
    auto gs_scalar  = saved[0];       // shape {}
    auto grad_out   = grad_outputs[0]; // shape {} (upstream scalar gradient)

    // gradient w.r.t. input (w_cum): the autograd engine expects a tensor
    // with EXACTLY the same shape as the original input tensor passed to
    // forward().  We return torch::Tensor() (undefined / null) here and
    // instead rely on the outer multiplication in Backward() to propagate
    // gradients through nn_out, which IS connected to the network parameters.
    //
    // Rationale: w_cum itself does not require_grad() — it is built from
    // w_curr (output of BundleSolverML_W::apply) and nn_out.  The gradient
    // path to the network parameters runs through nn_out, not through w_cum
    // directly.  Returning a null gradient for input is correct and avoids
    // the shape mismatch error.
    return {
      torch::Tensor(),   // grad w.r.t. input (w_cum) — not needed
      torch::Tensor(),   // grad w.r.t. gs_scalar — constant, not differentiated
      torch::Tensor()    // grad w.r.t. fi_val (double) — not differentiated
    };
  }
};

// ============================================================================
// BundleSolverML_W
// Custom autograd Function for the projected search direction w.
// Forward:  returns sum(w) to attach the computation graph.
// Backward: gradient w.r.t. t via KKT projection.
// ============================================================================
struct BundleSolverML_W : public torch::autograd::Function<BundleSolverML_W> {
  static torch::Tensor forward(torch::autograd::AutogradContext* ctx,
                               torch::Tensor w, torch::Tensor G,
                               torch::Tensor Q, torch::Tensor alpha,
                               double ts) {
    torch::Tensor ts_tensor = torch::tensor(ts, Q.options());
    ctx->save_for_backward({G, Q, alpha, ts_tensor});
    return w.sum();
  }

  static torch::autograd::tensor_list backward(
      torch::autograd::AutogradContext* ctx,
      torch::autograd::tensor_list grad_outputs) {
    auto saved          = ctx->get_saved_variables();
    torch::Tensor G     = saved[0];
    torch::Tensor Q     = saved[1];
    torch::Tensor alpha = saved[2];
    torch::Tensor ts    = saved[3];
    torch::Tensor grad_out = grad_outputs[0];

    // Regularize Q and compute inverse via Cholesky
    auto I    = torch::eye(Q.size(0), Q.options());
    auto Qreg = Q + 1e+1 * I;
    auto L    = torch::linalg_cholesky(Qreg);
    auto Qinv = torch::cholesky_inverse(L);

    // KKT projection: lambda* = (e^T Q^{-1} alpha) / (e^T Q^{-1} e)
    auto e         = torch::ones({Q.size(0)}, Q.options());
    auto Qinvalpha = torch::matmul(Qinv, alpha);
    auto Qinv_e    = torch::matmul(Qinv, e);
    auto num       = torch::dot(e, Qinvalpha);
    auto den       = torch::dot(e, Qinv_e);
    auto proj      = (num / den) * e - alpha;

    // dw/dt = (1/t^2) G^T Q^{-1} proj
    double ts_val = ts.item<double>();
    auto result   = torch::matmul(G.transpose(0, 1), torch::matmul(Qinv, proj));
    result        = (1.0 / (ts_val * ts_val)) * result * grad_out;

    return {result, torch::Tensor(), torch::Tensor(), torch::Tensor(),
            torch::Tensor()};
  }
};

namespace SMSpp_di_unipi_it {
using HpNum  = NDO_di_unipi_it::HpNum;
using cHpRow = NDO_di_unipi_it::cHpRow;
using Index  = Function::Index;

static constexpr auto InINF = SMSpp_di_unipi_it::Inf<Index>();

static void chgsign(double* v, Index n) {
  for (const auto ev = v + n; v < ev; ++v) *v = -*v;
}

SMSpp_insert_in_factory_cpp_0(BundleSolverML);

// ============================================================================
// SaveModel / LoadModel
// ============================================================================
void BundleSolverML::SaveModel(const std::string& filepath) {
  BML_LOG("SaveModel: saving to \"" << filepath << "\"\n");
  try {
    nn->eval();
    torch::serialize::OutputArchive archive;
    nn->save(archive);
    archive.save_to(filepath);
    nn->train();
    BML_LOG("SaveModel: success\n");
  } catch (const c10::Error& e) {
    nn->train();
    throw std::runtime_error(std::string("SaveModel: LibTorch error: ") + e.what());
  } catch (const std::exception& e) {
    nn->train();
    throw std::runtime_error(std::string("SaveModel: ") + e.what());
  }
}

void BundleSolverML::LoadModel(const std::string& filepath) {
  BML_LOG("LoadModel: loading from \"" << filepath << "\"\n");
  try {
    torch::serialize::InputArchive archive;
    archive.load_from(filepath);
    nn->load(archive);
    nn->eval();
    BML_LOG("LoadModel: success\n");
  } catch (const c10::Error& e) {
    throw std::runtime_error(std::string("LoadModel: LibTorch error: ") + e.what());
  } catch (const std::exception& e) {
    throw std::runtime_error(std::string("LoadModel: ") + e.what());
  }
}

// ============================================================================
// BundleSolverML::Heuristic
// ============================================================================
HpNum BundleSolverML::Heuristic(Index whch) {

  if (!nn->is_training()) {
    BML_LOG("WARNING: Heuristic in eval() mode — switching to train()\n");
    nn->train();
  }

  // Feature normalization: raw values span many orders of magnitude.
  // Dividing by per-feature scales keeps inputs O(1) at initialization,
  // preventing Softplus from outputting t ~1e10 with random weights.
  // Adjust scales if your instances have very different magnitudes.
  features[0]  = static_cast<float>(tHasChgd);           // binary {0,1}
  features[1]  = static_cast<float>(G1Norm)   / 1e3f;    // subgradient norm
  features[2]  = static_cast<float>(ScPr1)    / 1e4f;    // scalar product
  features[3]  = static_cast<float>(Alfa1)    / 1e4f;    // linearization error
  features[4]  = 0.0f;
  features[5]  = static_cast<float>(t)        / 1e3f;    // current step-size
  features[6]  = static_cast<float>(Sigma)    / 1e4f;    // predicted decrease
  features[7]  = static_cast<float>(DSTS)     / 1e3f;    // dist to center
  features[8]  = 0.0f;
  features[9]  = static_cast<float>(EpsU)     * 1e3f;    // EpsU ~1e-3 → *1e3
  features[10] = static_cast<float>(CSSCntr)  / 100.0f;  // SS counter
  features[11] = static_cast<float>(CNSCntr)  / 100.0f;  // NS counter
  features[12] = static_cast<float>(UpTrgt)   / 1e5f;    // upper target
  features[13] = static_cast<float>(LwTrgt)   / 1e5f;    // lower target
  features[14] = static_cast<float>(Fi0Lmb)   / 1e5f;    // f at tentative pt
  features[15] = static_cast<float>(Fi0Lmb1)  / 1e5f;    // f at center
  features[16] = static_cast<float>(DST)      / 1e3f;    // dist to center (alt)
  features[17] = static_cast<float>(NrmD)     / 1e3f;    // norm of direction
  features[18] = static_cast<float>(NrmZ)     / 1e3f;    // norm of agg. subgrad
  features[19] = static_cast<float>(NrmZFctr) / 1e3f;    // normalisation factor

  // phi_vecs stores the normalized tensor so Backward() re-runs forward
  // with the same scaled inputs that produced the original prediction.
  auto input  = torch::tensor(features);
  auto output = nn->forward(input);

  double t_pred = output.item<double>();
  if (!std::isfinite(t_pred) || t_pred <= 0.0) {
    BML_LOG("WARNING: invalid t=" << t_pred << " — clamping to 1.0\n");
    t_pred = 1.0;
  }

  phi_vecs.push_back(input);
  BML_LOG("Heuristic: predicted t=" << t_pred << "\n");

  // Search direction w
  Index dim;
  const Index* nms;
  std::vector<double> tZ(NumVar);
  Master->ReadZ(tZ.data(), nms, dim);
  tZ.resize(dim);
  w_vecs.push_back(torch::tensor(tZ).requires_grad_(true));

  // Subgradient matrix G
  int col_num = 0;
  for (Index i = 0; i < Master->MaxName(); ++i)
    if (ItemVcblr[i].second < vBPar2[ItemVcblr[i].first]) col_num++;

  std::vector<std::vector<VarValue>> G_mat;
  G_mat.reserve(col_num);
  for (Index i = 0; i < Master->MaxName(); ++i) {
    std::vector<VarValue> G(NumVar);
    if (ItemVcblr[i].second < vBPar2[ItemVcblr[i].first]) {
      v_c05f[ItemVcblr[i].first]->get_linearization_coefficients(
          G.data(), Range(0, NumVar), ItemVcblr[i].second);
      if (!f_convex) chgsign(G.data(), NumVar);
      G_mat.push_back(G);
    }
  }

  int rows = G_mat.size();
  int cols = G_mat[0].size();
  std::vector<double> flat;
  flat.reserve(rows * cols);
  for (const auto& row : G_mat) flat.insert(flat.end(), row.begin(), row.end());

  torch::Tensor G_tensor =
      torch::from_blob(flat.data(), {rows, cols}, torch::kDouble).clone();
  Gs.push_back(G_tensor);

  // Aggregated subgradient (whisG1 indices)
  std::vector<int64_t> valid_indices;
  for (auto k : whisG1) {
    if (k == InINF) break;
    if (k >= 0 && (int)k < cols) valid_indices.push_back(static_cast<int64_t>(k));
  }
  torch::Tensor grad_sum;
  if (valid_indices.empty()) {
    grad_sum = torch::zeros({rows}, torch::kDouble);
  } else {
    auto idx = torch::tensor(valid_indices, torch::kLong);
    grad_sum = G_tensor.index_select(1, idx).sum(1);
  }
  Gs_aggreg.push_back(grad_sum);

  // Step-type coefficient: +1 SS, -1 NS
  int n     = this->FakeFi.GetNumVar();
  cHpRow tA = Master->ReadLinErr();
  coeff_vecs.push_back(torch::tensor((CNSCntr >= 1) ? -1.0f : 1.0f));

  // Variable bounds
  std::vector<LMNum> ubs(n);
  c_Vec_VarValue Lambda1 = this->get_tentative_point();
  for (int i = 0; i < n; i++) ubs[i] = this->FakeFi.GetUB(i);

  // Master QP multipliers
  Index MDBm;
  cIndex_Set MBse;
  cHpRow Mlt = Master->ReadMult(MBse, MDBm);
  std::vector<double> mult(n, 0.0);
  if (MBse) {
    Index h = 0;
    for (Index i = 0; i < n; ++i) {
      if (MBse[h] == InINF) break;
      if (MBse[h] == i) mult[i] = Mlt[h++];
    }
  } else {
    for (Index i = 0; i < n; ++i) mult[i] = Mlt[i];
  }
  thetaS.push_back(torch::tensor(mult, torch::kDouble));

  // Gram matrix Q and alpha
  torch::Tensor Q = torch::matmul(G_tensor, G_tensor.transpose(0, 1));
  Qs.push_back(Q);
  std::vector<double> alpha(Q.sizes()[0], 0.0);
  for (Index i = 0; i < Q.sizes()[0]; ++i) alpha[i] = tA[i];
  alphaS.push_back(torch::tensor(alpha, torch::kDouble));

  tS.push_back(t_pred);
  FiS.push_back(UpFiBest);

  return static_cast<HpNum>(t_pred);
}

// ============================================================================
// BundleSolverML::w
// ============================================================================
torch::Tensor BundleSolverML::w(size_t f, double t) {
  return BundleSolverML_W::apply(w_vecs[f], Gs[f], Qs[f], alphaS[f], t);
}

// ============================================================================
// BundleSolverML::Backward
//
// Loss design:
//   For each included iteration f (SS or last):
//     nn_out_norm = nn->forward(phi_norm[f]).sum()          [re-evaluated]
//     w_curr      = BundleSolverML_W::apply(...)            [scalar]
//     w_cum      += 0.9^(total - ss_count) * nn_out * w_curr
//     gs_scalar   = Gs_aggreg[f].sum()                      [scalar constant]
//     term        = coeff * nn_out * BundleSolverML_Fi(w_cum, gs_scalar, fi)
//     loss       += term
//   loss /= ss_count                                         [mean reduction]
//
// After loss.backward():
//   - Gradient clipping (max_norm = GRAD_CLIP_NORM) prevents inf/nan grads.
//   - optimizer_.step() applies the Adam update to the network weights.
//
// Three sources of exploding gradients fixed here:
//
//  1. FEATURE NORMALIZATION (in Heuristic, applied before storing phi_vecs):
//     Raw features span many orders of magnitude (fi_val ~1e5, norms ~1e3).
//     Each feature is divided by a hand-tuned scale so all inputs are O(1).
//     This keeps the network's linear outputs small at initialization and
//     prevents Softplus from outputting t ~1e10.
//     The same normalization is applied when re-running forward in Backward()
//     by storing already-normalized tensors in phi_vecs.
//
//  2. LOSS NORMALIZATION:
//     Dividing by ss_count (number of contributing iterations) converts the
//     sum loss into a mean loss, preventing the magnitude from scaling with
//     the number of iterations.  The fi_val multiplier is also removed from
//     the loss term — it was contributing an extra ~1e5 factor that served
//     no gradient-signal purpose (it is a constant w.r.t. nn parameters).
//
//  3. GRADIENT CLIPPING + ADAM OPTIMIZER:
//     torch::nn::utils::clip_grad_norm_(params, GRAD_CLIP_NORM) rescales all
//     gradients so their global L2 norm does not exceed GRAD_CLIP_NORM.
//     The Adam optimizer (lr=1e-4) is held as a member (optimizer_) so its
//     moment estimates persist across successive Backward() calls on the same
//     shared network, providing adaptive per-parameter learning rates.
// ============================================================================

/// Maximum L2 norm for gradient clipping applied after loss.backward().
/// If grad_norm > GRAD_CLIP_NORM all gradients are rescaled by
/// GRAD_CLIP_NORM / grad_norm.  Typical values: 1.0 – 10.0.
static constexpr float GRAD_CLIP_NORM = 1.0f;

void BundleSolverML::Backward() {
  BML_LOG("Backward: starting (" << phi_vecs.size() << " iterations recorded)\n");

  try {
    // Lazy initialisation of the Adam optimizer on first Backward() call.
    // Doing this lazily means the optimizer is always bound to the ACTIVE
    // network (nn), whether it is the private nn_owned_ or a shared network
    // injected via set_shared_net().  lr=1e-4 is a conservative starting
    // point for this loss scale; adjust via set_lr() on the optimizer options.
    if (!optimizer_) {
      auto opts = torch::optim::AdamOptions(/*lr=*/1e-4)
                    .betas({0.9, 0.999})
                    .eps(1e-8)
                    .weight_decay(1e-5);
      optimizer_ = std::make_unique<torch::optim::Adam>(nn->parameters(), opts);
    }

    // Zero gradients before accumulation
    optimizer_->zero_grad();

    if (phi_vecs.empty()) {
      BML_LOG("Backward: no iterations — skipping\n");
      return;
    }

    if (!nn->is_training()) {
      BML_LOG("Backward: WARNING in eval() — switching to train()\n");
      nn->train();
    }

    torch::Tensor loss = torch::zeros({}, torch::kFloat32);
    size_t last_idx    = phi_vecs.size() - 1;
    int    ss_count    = 0;

    // Scalar accumulator for discounted search directions.
    torch::Tensor w_cum = torch::zeros({}, torch::kFloat32);

    for (size_t f = 0; f < phi_vecs.size(); ++f) {
      float coeff_val = coeff_vecs[f].item<float>();
      bool  is_ss     = (coeff_val > 0.0f);
      bool  is_last   = (f == last_idx);
      if (!is_ss && !is_last) continue;
      ss_count++;

      if (!w_vecs[f].defined() || !Gs[f].defined() || !Qs[f].defined()
          || f >= alphaS.size() || !alphaS[f].defined()) {
        std::cerr << "Backward: missing tensor at iter " << f << "\n";
        continue;
      }

      try {
        // phi_vecs[f] already contains the NORMALIZED feature tensor stored
        // by Heuristic(), so re-running forward here uses the same scaled
        // inputs that produced the original prediction.
        auto nn_out = nn->forward(phi_vecs[f]).sum(); // scalar, float32

        if (!nn_out.requires_grad()) {
          std::cerr << "Backward: iter " << f << " — no grad_fn on nn_out\n";
          continue;
        }

        auto w_curr = BundleSolverML_W::apply(
            w_vecs[f], Gs[f], Qs[f], alphaS[f], nn_out.item<double>());

        double discount = std::pow(0.9, static_cast<double>(phi_vecs.size() - ss_count));
        w_cum = w_cum + static_cast<float>(discount) * nn_out * w_curr;

        auto gs_scalar = Gs_aggreg[f].to(torch::kFloat32).sum(); // scalar {}

        // Loss term: coeff * nn_out * Fi(w_cum, gs, fi)
        // fi_val is intentionally NOT included as a multiplier: it is a
        // constant w.r.t. the network parameters and its ~1e5 magnitude
        // was the primary driver of exploding loss and gradients.
        auto term = coeff_vecs[f] * nn_out
                    * BundleSolverML_Fi::apply(w_cum, gs_scalar, FiS[f]);
        loss += term;

        BML_LOG("Backward: iter " << f
            << (is_ss ? " [SS]" : " [last-NS]")
            << "  discount=" << discount
            << "  nn_out="   << nn_out.item<float>()
            << "  contrib="  << term.item<float>() << "\n");

      } catch (const std::exception& e) {
        std::cerr << "Backward: ERROR at iter " << f << ": " << e.what() << "\n";
      }
    }

    if (ss_count == 0) {
      BML_LOG("Backward: no SS/last-iter — skipping\n");
      return;
    }

    // Mean reduction: prevents loss magnitude from growing with iteration count
    loss = loss / static_cast<float>(ss_count);

    BML_LOG("Backward: mean loss=" << loss.item<float>()
            << "  ss_count=" << ss_count << "/" << phi_vecs.size() << "\n");

    loss.backward();

    // Gradient clipping: rescale all gradients so global L2 norm <= GRAD_CLIP_NORM.
    // This is the primary defence against inf gradients when loss magnitudes
    // are large.  clip_grad_norm_ returns the pre-clipping norm for logging.
    auto params = nn->parameters();
    float pre_clip_norm = torch::nn::utils::clip_grad_norm_(params, GRAD_CLIP_NORM);
    BML_LOG("Backward: grad_norm (pre-clip)=" << pre_clip_norm
            << "  clip_threshold=" << GRAD_CLIP_NORM << "\n");

    // Log post-clipping gradient norms
    int pidx = 0;
    for (auto& p : params) {
      BML_LOG("Backward: param[" << pidx++ << "] grad_norm="
              << (p.grad().defined() ? p.grad().norm().item<float>() : 0.0f) << "\n");
    }

    // Apply the Adam update.  optimizer_ holds the network's parameters and
    // persists across calls, so moment estimates accumulate correctly.
    optimizer_->step();
    BML_LOG("Backward: optimizer step done\n");

  } catch (const std::exception& e) {
    std::cerr << "Backward: exception: " << e.what() << "\n";
  } catch (...) {
    std::cerr << "Backward: unknown exception\n";
  }

  BML_LOG("Backward: done\n");
}

} // namespace SMSpp_di_unipi_it
