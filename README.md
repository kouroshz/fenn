# Free Energy Nearest Neighbour (FENN)
FENN is an R package implementing a Von Neumann Entropy based distance metric learning model. The method identifies informative and non-informative features in the data. It scales and projects the data onto a subspace of maximal separation between classes. The corresponding optimal weights on each dimension are identified by optimizing the Von Neumann entropy. The model can be applied to single class (unsupervised) or multi class data, bridging the gap between PCA, LDA and distance metric learning. 


## Installation
To install the package directly from GitHub, you need the devtools library.
```{R}
library(devtools)
install_github("kouroshz/fenn", local = FALSE)
```
## Usage
```
fit <- fenn(y~., data)
p <- predict(fit, newdata)
```

# Citation
"A Free Energy Based Approach For Distance Metric Learning"", Inaba, Sho and Fakhry, Carl T. and Kulkarni, Rahul V. and Zarringhalam, Kourosh. The 25th ACM SIGKDD Conference on Knowledge Discovery and Data Mining ACM. August 4-8, 2019, Anchorage, AK, USA 