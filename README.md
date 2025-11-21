[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=mybarnard)

🌐 Overview
Barnard’s exact test is a powerful alternative to Fisher’s exact test for 2×2 contingency tables. It computes an unconditional exact p-value by maximizing over a nuisance parameter. This implementation is fully vectorized, extremely fast, and avoids the computational burden historically associated with Barnard’s procedure.

⭐ Features
• Fully vectorized computation (no loops in the core engine)
• Wald statistic of observed table
• Exact one-tailed and two-tailed p-values
• Optimal nuisance parameter
• Optional diagnostic plots
• Reproduces the original 2009 behaviour with improved numerical stability

🛠️ Installation
Download the file mybarnard.m and place it anywhere inside your MATLAB path.

▶️ Usage
mybarnard(X)
mybarnard(X, PLTS)
mybarnard(X, PLTS, Tbx)

🔣 Inputs
X A 2×2 matrix of counts.
PLTS 0 = no plots (default), 1 = show plots.
Tbx Number of grid points for the nuisance parameter (default = 100).

📤 Outputs
If no output argument is requested, the function prints a summary table.
If an output is requested:
STATS.TX0 Wald statistic of the observed table
STATS.p_value Maximized one-tailed Barnard p-value
STATS.nuisance Optimal nuisance parameter np*

📘 Interpretation
The method evaluates the probability of observing data at least as extreme as the observed table across all possible values of the nuisance parameter np ∈ (0,1).
The reported p-value is the supremum of these conditional probabilities.

📝 Notes
Barnard’s exact test is more powerful than Fisher’s because it is unconditional.
Increasing Tbx improves accuracy but increases computational cost.
This implementation uses log-gamma expansions to avoid overflow.

📚 Citation
Cardillo G. (2009)
“MyBarnard: a very compact routine for Barnard's exact test on 2×2 matrix”.
Available from GitHub:
https://github.com/dnafinder/mybarnard

👤 Author
Giuseppe Cardillo
Email: giuseppe.cardillo.75@gmail.com

GitHub: https://github.com/dnafinder

⚖️ License
MIT License or another license of your choice.
