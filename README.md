# eds217-ipynumpy
EDS 217 group project: Plotly

Plotly website: https://plotly.com/
Plotly Github Repo: https://github.com/plotly/plotly.py


The plotly Python library is an interactive, open-source plotting library that supports over 40 unique chart types covering a wide range of statistical, financial, geographic, scientific, and 3-dimensional use-cases.

Built on top of the Plotly JavaScript library (plotly.js), plotly enables Python users to create beautiful interactive web-based visualizations that can be displayed in Jupyter notebooks, saved to standalone HTML files, or served as part of pure Python-built web applications using Dash. The plotly Python library is sometimes referred to as "plotly.py" to differentiate it from the JavaScript library.

Thanks to deep integration with our Kaleido image export utility, plotly also provides great support for non-web contexts including desktop editors (e.g. QtConsole, Spyder, PyCharm) and static document publishing (e.g. exporting notebooks to PDF with high-quality vector images).


For reproducibility, you'll need to follow these installation instructions:
1. Create the conda environnment from the ipynumpy.yml file using this command:
conda create -f ipynumpy.yml
2. Create the kernel
python -m ipykernel install --user --name ipynumpy --display-name "Python (ipynumpy)"
3. Install the following modules using pip
pip install chart-studio
pip install geopandas