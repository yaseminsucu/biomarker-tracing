import numpy as np
import matplotlib.pyplot as plt
from sklearn.base import BaseEstimator, TransformerMixin, clone
from sklearn.utils import check_random_state, check_X_y
from joblib import Parallel, delayed

class StabilitySelection(BaseEstimator, TransformerMixin):
    def __init__(self, base_estimator, lambda_name, lambda_grid,
                 n_bootstrap=100, sample_fraction=0.5, 
                 threshold=0.6, n_jobs=1, verbose=0, random_state=None):
        self.base_estimator = base_estimator
        self.lambda_name = lambda_name
        self.lambda_grid = lambda_grid
        self.n_bootstrap = n_bootstrap
        self.sample_fraction = sample_fraction
        self.threshold = threshold
        self.n_jobs = n_jobs
        self.verbose = verbose
        self.random_state = random_state

    def fit(self, X, y):
        X, y = check_X_y(X, y, accept_sparse=['csr', 'csc', 'coo'])
        n_samples, n_features = X.shape
        n_lambdas = len(self.lambda_grid)
        rng = check_random_state(self.random_state)

        def _fit_bootstrap(X, y, lambda_name, lambda_value):
            # 1. Create a clone of the model and set the alpha (lambda)
            model = clone(self.base_estimator)
            model.set_params(**{lambda_name: lambda_value})
            
            # 2. Subsample the data (bootstrap)
            indices = rng.choice(n_samples, int(self.sample_fraction * n_samples), replace=False)
            model.fit(X[indices], y[indices])
            
            # 3. Extract coefficients (handle both raw models and Pipelines)
            if hasattr(model, 'coef_'):
                coef = model.coef_
            elif hasattr(model, 'steps'):
                coef = model.steps[-1][1].coef_
            else:
                coef = model.feature_importances_
                
            # 4. Return binary mask (was the feature selected?)
            return (np.abs(coef) > 1e-10).astype(int).flatten()

        stability_scores = np.zeros((n_features, n_lambdas))

        # Iterate through the grid of penalty values
        for idx, lv in enumerate(self.lambda_grid):
            print(f"Fitting lambda {idx+1}/{n_lambdas}: {lv:.4f}")
            
            # Run bootstraps in parallel
            out = Parallel(n_jobs=self.n_jobs)(
                delayed(_fit_bootstrap)(X, y, self.lambda_name, lv)
                for _ in range(self.n_bootstrap)
            )
            
            # Average the binary masks to get the stability score for this lambda
            stability_scores[:, idx] = np.mean(np.vstack(out), axis=0)

        self.stability_scores_ = stability_scores
        return self

    def get_support(self, indices=False, threshold=None):
        if threshold is None: threshold = self.threshold
        # A feature is selected if its stability score exceeds the threshold at ANY lambda
        mask = np.any(self.stability_scores_ >= threshold, axis=1)
        return np.where(mask)[0] if indices else mask

def plot_stability_path(selector, threshold_highlight=0.6):
    fig, ax = plt.subplots(figsize=(10, 6))
    lambdas = selector.lambda_grid
    scores = selector.stability_scores_
    
    # Plot a line for every cell type
    for i in range(scores.shape[0]):
        ax.plot(lambdas, scores[i, :], color='grey', alpha=0.3)
    
    if threshold_highlight:
        ax.axhline(threshold_highlight, color='red', linestyle='--', label=f'Threshold {threshold_highlight}')
    
    ax.set_xscale('log')
    ax.set_xlabel('Regularization Parameter (Alpha/Lambda)')
    ax.set_ylabel('Stability Score (Frequency of Selection)')
    ax.set_title('Lasso Stability Selection Path')
    return fig, ax