from collections.abc import Mapping
from typing import Any, Optional, Union

import numpy as np
import pandas as pd
from formulaic import Formula
from formulaic.utils.context import capture_context as _capture_context

from pyfixest.utils.dev_utils import _create_rng


def ssc(
    adj: bool = True,
    fixef_k: str = "nested",
    cluster_adj: bool = True,
    cluster_df: str = "min",
) -> dict[str, Union[str, bool]]:
    """
    Set the small sample correction factor applied in `get_ssc()`.

    Parameters
    ----------
        adj : bool, default True
            If True, applies a small sample correction of (N-1) / (N-k) where N
            is the number of observations and k is the number of estimated
            coefficients excluding any fixed effects projected out by either
            `feols()` or `fepois()`.
        fixef_k : str, default "none"
            Equal to 'none': the fixed effects parameters are discarded when
            calculating k in (N-1) / (N-k).
        cluster_adj : bool, default True
            If True, a cluster correction G/(G-1) is performed, with G the number
            of clusters.
        cluster_df : str, default "conventional"
            Controls how "G" is computed for multiway clustering if cluster_adj = True.
            Note that the covariance matrix in the multiway clustering case is of
            the form V = V_1 + V_2 - V_12. If "conventional", then each summand G_i
            is multiplied with a small sample adjustment G_i / (G_i - 1). If "min",
            all summands are multiplied with the same value, min(G) / (min(G) - 1)

    Details
    -------
    The small sample correction choices mimic fixest's behavior. For details, see
    https://cran.r-project.org/web/packages/fixest/vignettes/standard_errors.html.

    In general, if adj = True, we multiply the variance covariance matrix V with a
    small sample correction factor of (N-1) / (N-k), where N is the number of
    observations and k is the number of estimated coefficients.

    If fixef_k = "none", the fixed effects parameters are discarded when
    calculating k. This is the default behavior and currently the only
    option. Note that it is not r-fixest's default behavior.

    Hence if adj = True, the covariance matrix is computed as
    V = V x (N-1) / (N-k) for iid and heteroskedastic errors.

    If adj = False, no small sample correction is applied of the type
    above is applied.

    If cluster_adj = True, a cluster correction of G/(G-1) is performed,
    with G the number of clusters.

    If adj = True and cluster_adj = True, V = V x (N - 1) / N - k) x G/(G-1)
    for cluster robust errors where G is the number of clusters.

    If adj = False and cluster_adj = True, V = V x G/(G-1) for cluster robust
    errors, i.e. we drop the (N-1) / (N-k) factor. And if cluster_adj = False,
    no cluster correction is applied.

    Things are slightly more complicated for multiway clustering. In this
    case, we compute the variance covariance matrix as V = V1 + V2 - V_12.

    If cluster_adj = True and cluster_df = "conventional", then
    V += [V x G_i / (G_i - 1) for i in [1, 2, 12]], i.e. each separate
    covariance matrix G_i is multiplied with a small sample adjustment
    G_i / (G_i - 1) corresponding to the number of clusters in the
    respective covariance matrix. This is the default behavior
    for clustered errors.

    If cluster_df = "min", then
    V += [V x min(G) / (min(G) - 1) for i in [1, 2, 12]].

    Returns
    -------
    dict
        A dictionary with encoded info on how to form small sample corrections
    """
    if adj not in [True, False]:
        raise ValueError("adj must be True or False.")
    if fixef_k not in ["none", "full", "nested"]:
        raise ValueError(
            f"fixef_k must be 'none', 'full', or 'nested' but it is {fixef_k}."
        )
    if cluster_adj not in [True, False]:
        raise ValueError("cluster_adj must be True or False.")
    if cluster_df not in ["conventional", "min"]:
        raise ValueError("cluster_df must be 'conventional' or 'min'.")

    return {
        "adj": adj,
        "fixef_k": fixef_k,
        "cluster_adj": cluster_adj,
        "cluster_df": cluster_df,
    }


def get_ssc(
    *,
    ssc_dict: dict[str, Union[str, str | bool]],
    N: int,
    k: int,
    k_fe: int,
    k_fe_nested: int,
    n_fe: int,
    n_fe_fully_nested: int,
    G: int,
    vcov_sign: int,
    vcov_type: "str",
) -> tuple[np.ndarray, int, int]:
    """
    Compute small sample adjustment factors.

    Parameters
    ----------
    ssc_dict : dict
        A dictionary created via the ssc() function.
    N : int
        The number of observations.
    k : int
        The number of estimated parameters (as in the first part of the model formula)
    k_fe : int
        The number of estimated fixed effects (as specified in the second part of the model formula).
    k_fe_nested : int
        The number of estimated fixed effects nested within clusters.
    n_fe : int
        The number of fixed effects in the model. I.e. 'Y ~ X1  | f1 + f2' has 2 fixed effects.
    n_fe_fully_nested : int
        The number of fixed effects that are fully nested within clusters.
    G : int
        The number of clusters.
    vcov_sign : array-like
        A vector that helps create the covariance matrix.
    vcov_type : str
        The type of covariance matrix. Must be one of "iid", "hetero", or "CRV".

    Returns
    -------
    tuple of np.ndarray and int
        A small sample adjustment factor and the effective number of coefficients k used in the adjustment.

    Raises
    ------
    ValueError
        If vcov_type is not "iid", "hetero", or "CRV", or if cluster_df is neither
        "conventional" nor "min".
    """
    adj = ssc_dict["adj"]
    fixef_k = ssc_dict["fixef_k"]
    cluster_adj = ssc_dict["cluster_adj"]
    cluster_df = ssc_dict["cluster_df"]

    cluster_adj_value = 1.0
    adj_value = 1.0

    # see here for why:
    # https://github.com/lrberge/fixest/issues/554

    # subtract one for each fixed effect, except for the first
    k_fe_adj = k_fe - (n_fe - 1) if n_fe > 1 else k_fe

    if fixef_k == "none":
        dof_k = k
    elif fixef_k == "nested":
        if n_fe == 0:
            dof_k = k
        elif k_fe_nested == 0:
            # no nested fe, so just add all fixed effects
            dof_k = k + k_fe_adj
        else:
            # subtract nested fixed effects and add one for each fully nested
            # subtracted fixed effect back
            dof_k = k + k_fe_adj - k_fe_nested + n_fe_fully_nested
    elif fixef_k == "full":
        # add all fixed effects
        dof_k = k + k_fe_adj if n_fe > 0 else k
    else:
        raise ValueError("fixef_k is neither none, nested, nor full.")

    if adj:
        adj_value = (N - 1) / (N - dof_k)

    # cluster_adj applied with G = N for hetero but not for iid
    if vcov_type in ["hetero", "CRV"] and cluster_adj:
        if cluster_df == "conventional":
            cluster_adj_value = G / (G - 1)
        elif cluster_df == "min":
            G = np.min(G)
            cluster_adj_value = G / (G - 1)
        else:
            raise ValueError("cluster_df is neither conventional nor min.")

    df_t = N - dof_k if vcov_type in ["iid", "hetero"] else G - 1

    return np.array([adj_value * cluster_adj_value * vcov_sign]), dof_k, df_t


def get_data(N=1000, seed=1234, beta_type="1", error_type="1", model="Feols"):
    """
    Create a random example data set.

    Parameters
    ----------
    N : int, optional
        Number of observations. Default is 1000.
    seed : int, optional
        Seed for the random number generator. Default is 1234.
    beta_type : str, optional
        Type of beta coefficients. Must be one of '1', '2', or '3'. Default is '1'.
    error_type : str, optional
        Type of error term. Must be one of '1', '2', or '3'. Default is '1'.
    model : str, optional
        Type of the DGP. Must be either 'Feols' or 'Fepois'. Default is 'Feols'.

    Returns
    -------
    pandas.DataFrame
        A pandas DataFrame with simulated data.

    Raises
    ------
    ValueError
        If beta_type is not '1', '2', or '3', or if error_type is not '1', '2', or '3',
        or if model is not 'Feols' or 'Fepois'.
    """
    rng = np.random.default_rng(seed)
    G = rng.choice(list(range(10, 20))).astype("int64")
    fe_dims = rng.choice(list(range(2, int(np.floor(np.sqrt(N))))), 3, True).astype(
        "int64"
    )

    # create the covariates
    X = rng.normal(0, 3, N * 5).reshape((N, 5))
    X[:, 0] = rng.choice(range(3), N, True)
    # X = pd.DataFrame(X)
    X[:, 2] = rng.choice(list(range(fe_dims[0])), N, True)
    X[:, 3] = rng.choice(list(range(fe_dims[1])), N, True)
    X[:, 4] = rng.choice(list(range(fe_dims[2])), N, True)

    X = pd.DataFrame(X)
    X.columns = ["X1", "X2", "f1", "f2", "f3"]
    # X1, X2, X3 as pd.Categorical
    X["f1"] = X["f1"].astype("category")
    X["f2"] = X["f2"].astype("category")
    X["f3"] = X["f3"].astype("category")

    mm = Formula("~ X1 + X2 + f1 + f2 + f3").get_model_matrix(data=X, output="pandas")

    k = mm.shape[1]

    # create the coefficients
    if beta_type == "1":
        beta = rng.normal(0, 1, k).reshape(k, 1)
    elif beta_type == "2":
        beta = rng.normal(0, 5, k).reshape(k, 1)
    elif beta_type == "3":
        beta = np.exp(rng.normal(0, 1, k)).reshape(k, 1)
    else:
        raise ValueError("beta_type needs to be '1', '2' or '3'.")

    # create the error term
    if error_type == "1":
        u = rng.normal(0, 1, N).reshape(N, 1)
    elif error_type == "2":
        u = rng.normal(0, 5, N).reshape(N, 1)
    elif error_type == "3":
        u = np.exp(rng.normal(0, 1, N)).reshape(N, 1)
    else:
        raise ValueError("error_type needs to be '1', '2' or '3'.")

    # create the depvar and cluster variable
    if model == "Feols":
        Y = (1 + mm.to_numpy() @ beta + u).flatten()
        Y2 = Y + rng.normal(0, 5, N)
    elif model == "Fepois":
        mu = np.exp(mm.to_numpy() @ beta).flatten()
        mu = 1 + mu / np.sum(mu)
        Y = rng.poisson(mu, N)
        Y2 = Y + rng.choice(range(10), N, True)
    else:
        raise ValueError("model needs to be 'Feols' or 'Fepois'.")

    Y, Y2 = (pd.Series(x.flatten()) for x in [Y, Y2])
    Y.name, Y2.name = "Y", "Y2"

    cluster = rng.choice(list(range(G)), N)
    cluster = pd.Series(cluster)
    cluster.name = "group_id"

    df = pd.concat([Y, Y2, X, cluster], axis=1)

    # add some NaN values
    df.loc[0, "Y"] = np.nan
    df.loc[1, "X1"] = np.nan
    df.loc[2, "f1"] = np.nan

    # compute some instruments
    df["Z1"] = df["X1"] + rng.normal(0, 1, N)
    df["Z2"] = df["X2"] + rng.normal(0, 1, N)

    # change all variables in the data frame to float
    for col in df.columns:
        df[col] = df[col].astype("float64")

    df[df == "nan"] = np.nan

    df["weights"] = rng.uniform(0, 1, N)
    # df["weights"].iloc[]

    if model == "Fepois":
        # add separation
        idx = np.array([10, 11, 12])
        df.loc[idx[0], "f1"] = np.max(df["f1"]) + 1
        df.loc[idx[1], "f2"] = np.max(df["f2"]) + 1
        df.loc[idx[2], "f3"] = np.max(df["f3"]) + 1

    return df


def simultaneous_crit_val(
    C: np.ndarray, S: int, alpha: float = 0.05, seed: Optional[int] = None
) -> float:
    """
    Simultaneous Critical Values.

    Obtain critical values for simultaneous inference on linear model parameters
    using the Multiplier bootstrap.

    Parameters
    ----------
    C: numpy.ndarray
        Covariance matrix. Symmetric, and contains as many rows/columns
        as parameters of interest.
    S: int
        Number of replications
    alpha: float
        Significance level. Defaults to 0.05
    seed: int, optional
        Seed for the random number generator. Default is None.

    Returns
    -------
    float
        Critical value, larger than 1.96
        (which is the crit-value for pointwise intervals)
    """

    def msqrt(C: np.ndarray) -> np.ndarray:
        eig_vals, eig_vecs = np.linalg.eigh(C)
        return eig_vecs @ np.diag(np.sqrt(eig_vals)) @ np.linalg.inv(eig_vecs)

    rng = _create_rng(seed)
    p = C.shape[0]
    tmaxs = np.max(np.abs(msqrt(C) @ rng.normal(size=(p, S))), axis=0)
    return np.quantile(tmaxs, 1 - alpha)


def capture_context(context: Union[int, Mapping[str, Any]]) -> Mapping[str, Any]:
    """
    Explicitly capture the context to be used by subsequent formula
    materialisations.

    Parameters
    ----------
    context: Union[int, Mapping[str, Any]]
        The context from which variables (and custom transforms/etc)
        should be inherited.

        When specified as an integer, it is interpreted as a frame offset
        from the caller's frame. Since we use this function in the context
        of a library, we need to account for the extra frames, hence
        we add 2 to the context (one for this and one for the frame the
        function is being called in).

        Otherwise, a mapping from variable name to value is expected.

        When nesting in a library, and attempting to capture user-context,
        make sure you account for the extra frames introduced by your wrappers.

    Returns
    -------
    Mapping[str, Any]
        The context that should be later passed to the Formulaic materialization
        procedure like: `.get_model_matrix(..., context=<this object>)`.
    """
    return _capture_context(context + 2) if isinstance(context, int) else context
