from typing import Optional

import numpy as np
import pandas as pd

from pyfixest.did.did import DID
from pyfixest.report.visualize import _coefplot


class ETWFE(DID):
    def __init__(
        self,
        data: pd.DataFrame,
        yname: str,
        idname: str,
        tname: str,
        gname: str,
        cluster: Optional[str] = None,
        xfml: Optional[str] = None,
        att: bool = True,
        # cgroup: Optional[str] = None,
    ):
        super().__init__(data, yname, idname, tname, gname, cluster, xfml, att)

        # if cgroup is None:
        self._cgroup = "notyet"

    def estimate(self):
        fit = _estimate_etwfe(self._fml, self._tvar, self.gvar, self.data, self._cgroup)
        self._coeftable = fit.tidy()

    def emfx(self):
        "Computes average marginal effects of interest."
        pass

    def vcov(self):
        pass

    def iplot(
        self,
        alpha: float = 0.05,
        figsize: Optional[tuple[int, int]] = None,
        yintercept: Optional[int] = None,
        xintercept: Optional[int] = None,
        rotate_xticks: int = 0,
        title: str = "LPDID Event Study Estimate",
        coord_flip: bool = False,
        plot_backend: str = "lets_plot",
    ):
        """
        Create coefficient plots.

        Parameters
        ----------
        alpha : float, optional
            Significance level for visualization options. Defaults to 0.05.
        figsize : tuple[int, int], optional
            Size of the plot (width, height) in inches. Defaults to (500, 300).
        yintercept : float, optional
            Value to set as the y-axis intercept (vertical line).
            Defaults to None.
        xintercept : float, optional
            Value to set as the x-axis intercept (horizontal line).
            Defaults to None.
        rotate_xticks : int, optional
            Rotation angle for x-axis tick labels. Defaults to 0.
        title : str, optional
            Title of the plot.
        coord_flip : bool, optional
            Whether to flip the coordinates of the plot. Defaults to False.
        plot_backend: str, optional
            The plotting backend to use between "lets_plot" (default) and "matplotlib".

        Returns
        -------
        lets-plot figure
            A lets-plot figure with coefficient estimates and confidence intervals.
        """
        df = self._coeftable
        df["fml"] = "etwfe"

        return _coefplot(
            plot_backend=plot_backend,
            df=df,
            figsize=figsize,
            alpha=alpha,
            yintercept=yintercept,
            xintercept=xintercept,
            rotate_xticks=rotate_xticks,
            title=title,
            flip_coord=coord_flip,
        )

    def tidy(self):
        pass

    def summary(self):
        pass


def _estimate_etwfe(fml, tvar, gvar, data, cgroup="notyet"):
    fml = fml.replace(" ", "")
    fml_split = fml.split("~")
    lhs = fml_split[0]

    # if fml_split["1"] != "1":
    #  raise NotImplementedError("Covariates are not yet implemented.")

    rhs = fml_split[1] if fml_split[1] != "1" else ""

    if cgroup not in ["notyet"]:
        raise NotImplementedError(f"cgroup must be 'notyet' but is {cgroup}.")

    ug = data[gvar].unique()
    ut = data[tvar].unique()

    gref = ug[ug > max(ut)]
    if gref.size == 0:
        gref = ug[ug < min(ut)]
    gref = gref[0]
    gref_min_flag = True

    # rhs = ""
    # rhs = f".Dtreat : {rhs}"
    # rhs = f"{rhs}C({gvar}, contr.treatment({gref})):C({tvar})"

    # fixef = f"{tvar} + {gvar}"
    # fml = f"{lhs} ~ {rhs} | {fixef}"

    rhs = ""
    rhs = f".Dtreat : {rhs}"
    rhs = f"{rhs}C({gvar}, contr.treatment({gref})):C({tvar})"
    lhs = "lemp"
    fixef = "first_treat + year"
    fml = f"{lhs} ~ {rhs} | {fixef}"

    data[".Dtreat"] = (data[tvar] >= data[gvar]) & (data[gvar] != gref)
    if not gref_min_flag:
        data[".Dtreat"] = np.where(data[tvar] < gref, data[".Dtreat"], np.nan)
    else:
        data[".Dtreat"] = np.where(data[tvar] > gref, data[".Dtreat"], np.nan)

    print("fml", fml)
    fit = pf.feols(fml, data=data)

    return fit
